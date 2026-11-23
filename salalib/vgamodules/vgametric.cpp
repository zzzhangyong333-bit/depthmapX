// sala - a component of the depthmapX - spatial network analysis platform
// Copyright (C) 2000-2010, University College London, Alasdair Turner
// Copyright (C) 2011-2026, Tasos Varoudis
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

#include "salalib/vgamodules/vgametric.h"

#include "salalib/parallelanalysis.h"
#include "salalib/vgamodules/vgasearchstate.h"

#include "genlib/stringutils.h"

#include <algorithm>

// This is a slow algorithm, but should give the correct answer
// for demonstrative purposes. Each source point is searched independently, so the
// searches are spread over all CPU cores.

bool VGAMetric::run(Communicator *comm, PointMap &map, bool) {
    time_t atime = 0;
    if (comm) {
        qtimer(atime, 0);
        comm->CommPostMessage(Communicator::NUM_RECORDS, map.getFilledPointCount());
    }

    std::string radius_text;
    if (m_radius != -1.0) {
        if (m_radius > 100.0) {
            radius_text = std::string(" R") + dXstring::formatString(m_radius, "%.f");
        } else if (map.getRegion().width() < 1.0) {
            radius_text = std::string(" R") + dXstring::formatString(m_radius, "%.4f");
        } else {
            radius_text = std::string(" R") + dXstring::formatString(m_radius, "%.2f");
        }
    }
    AttributeTable &attributes = map.getAttributeTable();

    // n.b. these must be entered in alphabetical order to preserve col indexing:
    std::string penn_col_text = std::string("Metric Mean Penn Distance") + radius_text;
    int penn_col = attributes.insertOrResetColumn(penn_col_text.c_str());
    std::string mspa_col_text = std::string("Metric Mean Shortest-Path Angle") + radius_text;
    int mspa_col = attributes.insertOrResetColumn(mspa_col_text.c_str());
    std::string mspl_col_text = std::string("Metric Mean Shortest-Path Distance") + radius_text;
    int mspl_col = attributes.insertOrResetColumn(mspl_col_text.c_str());
    std::string dist_col_text = std::string("Metric Mean Straight-Line Distance") + radius_text;
    int dist_col = attributes.insertOrResetColumn(dist_col_text.c_str());
    std::string count_col_text = std::string("Metric Node Count") + radius_text;
    int count_col = attributes.insertOrResetColumn(count_col_text.c_str());

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
        states.emplace_back(map.getRows(), map.getCols(), 0, -1.0f, 0.0f);
    }
    std::vector<depthmapX::PendingValues> results(sources.size());

    auto analyseSource = [&](PixelRef curs, depthmapX::VGASearchState &state, depthmapX::PendingValues &out) {
        state.reset();

        float euclid_depth = 0.0f;
        float total_depth = 0.0f;
        float total_angle = 0.0f;
        int total_nodes = 0;

        // note that misc is used in a different manner to analyseGraph / PointDepth
        // here it marks the node as used in calculation only

        depthmapX::SearchQueue<MetricTriple> search_list;
        search_list.insert(MetricTriple(0.0f, curs, NoPixel));
        MetricTriple here;
        while (search_list.pop(here)) {
            if (m_radius != -1.0 && (here.dist * map.getSpacing()) > m_radius) {
                break;
            }
            Point &p = map.getPoint(here.pixel);
            size_t pidx = state.index(here.pixel);
            // nb, the filled check is necessary as diagonals seem to be stored with 'gaps' left in
            if (p.filled() && state.misc[pidx] != ~0) {
                depthmapX::extractMetric(p.getNode(), search_list, map, state, here);
                state.touched.push_back(pidx);
                state.misc[pidx] = ~0;
                if (!p.getMergePixel().empty()) {
                    Point &p2 = map.getPoint(p.getMergePixel());
                    size_t p2idx = state.index(p.getMergePixel());
                    if (state.misc[p2idx] != ~0) {
                        state.touched.push_back(p2idx);
                        state.cumangle[p2idx] = state.cumangle[pidx];
                        depthmapX::extractMetric(p2.getNode(), search_list, map, state,
                                                 MetricTriple(here.dist, p.getMergePixel(), NoPixel));
                        state.misc[p2idx] = ~0;
                    }
                }
                total_depth += float(here.dist * map.getSpacing());
                total_angle += state.cumangle[pidx];
                euclid_depth += float(map.getSpacing() * dist(here.pixel, curs));
                total_nodes += 1;
            }
        }

        out.emplace_back(mspa_col, float(double(total_angle) / double(total_nodes)));
        out.emplace_back(mspl_col, float(double(total_depth) / double(total_nodes)));
        out.emplace_back(dist_col, float(double(euclid_depth) / double(total_nodes)));
        // This is the legacy Mean Penn Distance formula. The root contributes zero to both
        // totals and is included in total_nodes for compatibility with historic output.
        out.emplace_back(penn_col,
                         float(std::max(0.0, double(total_depth - euclid_depth)) / double(total_nodes)));
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
        depthmapX::VGASearchState finalState(map.getRows(), map.getCols(), 0, -1.0f, 0.0f);
        depthmapX::PendingValues ignored;
        analyseSource(sources.back(), finalState, ignored);
        finalState.writeTo(map);
    }

    map.overrideDisplayedAttribute(-2);
    map.setDisplayedAttribute(mspl_col);

    return true;
}
