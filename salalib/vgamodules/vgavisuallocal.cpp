// sala - a component of the depthmapX - spatial network analysis platform
// Copyright (C) 2000-2010, University College London, Alasdair Turner
// Copyright (C) 2011-2012, Tasos Varoudis
// Copyright (C) 2017-2018, Petros Koutsolampros

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "salalib/vgamodules/vgavisuallocal.h"

#include "salalib/parallelanalysis.h"
#include "salalib/vgamodules/vgasearchstate.h"

#include "genlib/stringutils.h"

#include <algorithm>

bool VGAVisualLocal::run(Communicator *comm, PointMap &map, bool simple_version) {
    time_t atime = 0;
    if (comm) {
        qtimer(atime, 0);
        comm->CommPostMessage(Communicator::NUM_RECORDS, map.getFilledPointCount());
    }

    int cluster_col = -1, control_col = -1, controllability_col = -1;
    if (!simple_version) {
        cluster_col = map.getAttributeTable().insertOrResetColumn("Visual Clustering Coefficient");
        control_col = map.getAttributeTable().insertOrResetColumn("Visual Control");
        controllability_col = map.getAttributeTable().insertOrResetColumn("Visual Controllability");
    }

    std::vector<PixelRef> sources;
    size_t skipped = 0;
    for (size_t i = 0; i < map.getCols(); i++) {
        for (size_t j = 0; j < map.getRows(); j++) {
            PixelRef curs = PixelRef(static_cast<short>(i), static_cast<short>(j));
            if (map.getPoint(curs).filled()) {
                if ((map.getPoint(curs).contextfilled() && !curs.iseven()) || (m_gates_only)) {
                    skipped++;
                    continue;
                }
                sources.push_back(curs);
            }
        }
    }

    // Membership tests used to be linear scans through the neighbourhood lists, which made this
    // analysis roughly cubic in the number of visible points. They are now O(1) lookups in a
    // per-thread "stamp" grid: a cell is in the set when its stamp equals the current stamp.
    struct LocalState {
        std::vector<unsigned int> inHood;
        std::vector<unsigned int> inTotal;
        unsigned int stamp = 0;
        PixelRefVector neighbourhood;
        explicit LocalState(size_t cells) : inHood(cells, 0), inTotal(cells, 0) {}
    };
    const size_t cols = map.getCols();
    auto cellIndex = [cols](PixelRef p) { return size_t(p.y) * cols + size_t(p.x); };

    const unsigned int nthreads = depthmapX::threadsFor(sources.size());
    std::vector<LocalState> states;
    states.reserve(nthreads);
    for (unsigned int t = 0; t < nthreads; t++) {
        states.emplace_back(map.getRows() * map.getCols());
    }
    std::vector<depthmapX::PendingValues> results(sources.size());

    auto analyseSource = [&](PixelRef curs, LocalState &state, depthmapX::PendingValues &out) {
        state.stamp++;
        const unsigned int stamp = state.stamp;

        // This is much easier to do with a straight forward list:
        PixelRefVector &neighbourhood = state.neighbourhood;
        neighbourhood.clear();
        depthmapX::forEachVisible(map.getPoint(curs).getNode(), [&](PixelRef pix) {
            size_t idx = cellIndex(pix);
            if (state.inHood[idx] != stamp) {
                state.inHood[idx] = stamp;
                neighbourhood.push_back(pix);
            }
        });

        // only required to match previous non-stl output. Without this
        // the output differs by the last digit of the float
        std::sort(neighbourhood.begin(), neighbourhood.end());

        int cluster = 0;
        float control = 0.0f;
        size_t totalneighbourhood = 0;

        for (size_t i = 0; i < neighbourhood.size(); i++) {
            int intersect_size = 0, retro_size = 0;
            Point &retpt = map.getPoint(neighbourhood[i]);
            if (retpt.filled() && retpt.hasNode()) {
                depthmapX::forEachVisible(retpt.getNode(), [&](PixelRef pix) {
                    retro_size++;
                    size_t idx = cellIndex(pix);
                    if (state.inHood[idx] == stamp) {
                        intersect_size++;
                    }
                    if (state.inTotal[idx] != stamp) {
                        state.inTotal[idx] = stamp;
                        totalneighbourhood++;
                    }
                });
                control += 1.0f / float(retro_size);
                cluster += intersect_size;
            }
        }
#ifndef _COMPILE_dX_SIMPLE_VERSION
        if (!simple_version) {
            if (neighbourhood.size() > 1) {
                out.emplace_back(cluster_col,
                                 float(cluster / double(neighbourhood.size() * (neighbourhood.size() - 1.0))));
                out.emplace_back(control_col, float(control));
                out.emplace_back(controllability_col,
                                 float(double(neighbourhood.size()) / double(totalneighbourhood)));
            } else {
                out.emplace_back(cluster_col, -1);
                out.emplace_back(control_col, -1);
                out.emplace_back(controllability_col, -1);
            }
        }
#endif
    };

    depthmapX::runAnalysisInParallel(
        comm, sources.size(), nthreads,
        [&](size_t index, unsigned int thread) { analyseSource(sources[index], states[thread], results[index]); },
        skipped);

    for (size_t index = 0; index < sources.size(); index++) {
        AttributeRow &row = map.getAttributeTable().getRow(AttributeKey(sources[index]));
        for (const auto &value : results[index]) {
            row.setValue(value.first, value.second);
        }
    }

#ifndef _COMPILE_dX_SIMPLE_VERSION
    if (!simple_version)
        map.setDisplayedAttribute(cluster_col);
#endif

    return true;
}
