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

#include "salalib/vgamodules/vgaangular.h"

#include "salalib/parallelanalysis.h"
#include "salalib/vgamodules/vgasearchstate.h"

#include "genlib/stringutils.h"

bool VGAAngular::run(Communicator *comm, PointMap &map, bool) {
    time_t atime = 0;
    if (comm) {
        qtimer(atime, 0);
        comm->CommPostMessage(Communicator::NUM_RECORDS, map.getFilledPointCount());
    }

    std::string radius_text;
    if (m_radius != -1.0) {
        if (map.getRegion().width() > 100.0) {
            radius_text = std::string(" R") + dXstring::formatString(m_radius, "%.f");
        } else if (map.getRegion().width() < 1.0) {
            radius_text = std::string(" R") + dXstring::formatString(m_radius, "%.4f");
        } else {
            radius_text = std::string(" R") + dXstring::formatString(m_radius, "%.2f");
        }
    }

    AttributeTable &attributes = map.getAttributeTable();

    // n.b. these must be entered in alphabetical order to preserve col indexing:
    std::string mean_depth_col_text = std::string("Angular Mean Depth") + radius_text;
    int mean_depth_col = attributes.getOrInsertColumn(mean_depth_col_text.c_str());
    std::string total_detph_col_text = std::string("Angular Total Depth") + radius_text;
    int total_depth_col = attributes.getOrInsertColumn(total_detph_col_text.c_str());
    std::string count_col_text = std::string("Angular Node Count") + radius_text;
    int count_col = attributes.getOrInsertColumn(count_col_text.c_str());

    // TODO: Binary compatibility. Remove in re-examination
    total_depth_col = attributes.getOrInsertColumn(total_detph_col_text.c_str());

    std::vector<PixelRef> sources;
    size_t skipped = 0;
    for (size_t i = 0; i < map.getCols(); i++) {
        for (size_t j = 0; j < map.getRows(); j++) {
            PixelRef curs = PixelRef(static_cast<short>(i), static_cast<short>(j));
            if (map.getPoint(curs).filled()) {
                if (m_gates_only) {
                    skipped++;
                    continue;
                }
                sources.push_back(curs);
            }
        }
    }

    const unsigned int nthreads = depthmapX::threadsFor(sources.size());
    std::vector<depthmapX::VGASearchState> states;
    states.reserve(nthreads);
    for (unsigned int t = 0; t < nthreads; t++) {
        states.emplace_back(map.getRows(), map.getCols(), 0, 0.0f, -1.0f);
    }
    std::vector<depthmapX::PendingValues> results(sources.size());

    auto analyseSource = [&](PixelRef curs, depthmapX::VGASearchState &state, depthmapX::PendingValues &out) {
        state.reset();

        float total_angle = 0.0f;
        int total_nodes = 0;

        // note that misc is used in a different manner to analyseGraph / PointDepth
        // here it marks the node as used in calculation only

        depthmapX::SearchQueue<AngularTriple> search_list;
        search_list.insert(AngularTriple(0.0f, curs, NoPixel));
        state.cumangle[state.touch(curs)] = 0.0f;
        AngularTriple here;
        while (search_list.pop(here)) {
            if (m_radius != -1.0 && here.angle > m_radius) {
                break;
            }
            Point &p = map.getPoint(here.pixel);
            size_t pidx = state.index(here.pixel);
            // nb, the filled check is necessary as diagonals seem to be stored with 'gaps' left in
            if (p.filled() && state.misc[pidx] != ~0) {
                depthmapX::extractAngular(p.getNode(), search_list, map, state, here);
                state.touched.push_back(pidx);
                state.misc[pidx] = ~0;
                if (!p.getMergePixel().empty()) {
                    Point &p2 = map.getPoint(p.getMergePixel());
                    size_t p2idx = state.index(p.getMergePixel());
                    if (state.misc[p2idx] != ~0) {
                        state.touched.push_back(p2idx);
                        state.cumangle[p2idx] = state.cumangle[pidx];
                        depthmapX::extractAngular(p2.getNode(), search_list, map, state,
                                                  AngularTriple(here.angle, p.getMergePixel(), NoPixel));
                        state.misc[p2idx] = ~0;
                    }
                }
                total_angle += state.cumangle[pidx];
                total_nodes += 1;
            }
        }

        if (total_nodes > 0) {
            out.emplace_back(mean_depth_col, float(double(total_angle) / double(total_nodes)));
        }
        out.emplace_back(total_depth_col, total_angle);
        out.emplace_back(count_col, float(total_nodes));
    };

    depthmapX::runAnalysisInParallel(
        comm, sources.size(), nthreads,
        [&](size_t index, unsigned int thread) { analyseSource(sources[index], states[thread], results[index]); },
        skipped);

    for (size_t index = 0; index < sources.size(); index++) {
        AttributeRow &row = attributes.getRow(AttributeKey(sources[index]));
        for (const auto &value : results[index]) {
            row.setValue(value.first, value.second);
        }
    }

    // leave the point map's search fields as the last analysed point left them (legacy behaviour)
    if (!sources.empty()) {
        depthmapX::VGASearchState finalState(map.getRows(), map.getCols(), 0, 0.0f, -1.0f);
        depthmapX::PendingValues ignored;
        analyseSource(sources.back(), finalState, ignored);
        finalState.writeTo(map);
    }

    map.setDisplayedAttribute(-2);
    map.setDisplayedAttribute(mean_depth_col);

    return true;
}
