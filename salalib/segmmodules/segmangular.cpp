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

#include "salalib/segmmodules/segmangular.h"
#include "salalib/options.h"
#include "salalib/parallelanalysis.h"

#include <algorithm>

#include "genlib/stringutils.h"

bool SegmentAngular::run(Communicator *comm, ShapeGraph &map, bool) {

    if (map.getMapType() != ShapeMap::SEGMENTMAP) {
        return false;
    }

    AttributeTable &attributes = map.getAttributeTable();

    time_t atime = 0;
    if (comm) {
        qtimer(atime, 0);
        comm->CommPostMessage(Communicator::NUM_RECORDS, map.getConnections().size());
    }

    // note: radius must be sorted lowest to highest, but if -1 occurs ("radius n") it needs to be last...
    // ...to ensure no mess ups, we'll re-sort here:
    bool radius_n = false;
    std::vector<double> radii;
    for (double radius : m_radius_set) {
        if (radius < 0) {
            radius_n = true;
        } else {
            radii.push_back(radius);
        }
    }
    if (radius_n) {
        radii.push_back(-1.0);
    }

    std::vector<int> depth_col, count_col, total_col;
    // first enter table values
    for (int radius : radii) {
        std::string radius_text = makeRadiusText(Options::RADIUS_ANGULAR, radius);
        std::string depth_col_text = std::string("Angular Mean Depth") + radius_text;
        attributes.insertOrResetColumn(depth_col_text.c_str());
        std::string count_col_text = std::string("Angular Node Count") + radius_text;
        attributes.insertOrResetColumn(count_col_text.c_str());
        std::string total_col_text = std::string("Angular Total Depth") + radius_text;
        attributes.insertOrResetColumn(total_col_text.c_str());
    }

    for (int radius : radii) {
        std::string radius_text = makeRadiusText(Options::RADIUS_ANGULAR, radius);
        std::string depth_col_text = std::string("Angular Mean Depth") + radius_text;
        depth_col.push_back(attributes.getColumnIndex(depth_col_text.c_str()));
        std::string count_col_text = std::string("Angular Node Count") + radius_text;
        count_col.push_back(attributes.getColumnIndex(count_col_text.c_str()));
        std::string total_col_text = std::string("Angular Total Depth") + radius_text;
        total_col.push_back(attributes.getColumnIndex(total_col_text.c_str()));
    }

    // Each source segment is an independent shortest-angle search, so the searches are spread
    // over all CPU cores. The results are written back afterwards in the original order.
    std::vector<AttributeRow *> rows;
    for (auto &iter : attributes) {
        rows.push_back(&iter.getRow());
    }

    struct SearchState {
        std::vector<char> covered;
        explicit SearchState(size_t n) : covered(n, 0) {}
    };
    const unsigned int nthreads = depthmapX::threadsFor(rows.size());
    std::vector<SearchState> states;
    states.reserve(nthreads);
    for (unsigned int t = 0; t < nthreads; t++) {
        states.emplace_back(map.getShapeCount());
    }
    std::vector<depthmapX::PendingValues> results(rows.size());
    const std::vector<Connector> &connections = map.getConnections();

    auto analyseSource = [&](size_t i, SearchState &state, depthmapX::PendingValues &out) {
        std::vector<char> &covered = state.covered;
        std::fill(covered.begin(), covered.end(), 0);

        depthmapX::StableMinQueue<float, SegmentData> anglebins;
        anglebins.push(0.0f, SegmentData(0, i, SegmentRef(), 0, 0.0, 0));

        std::vector<double> total_depth(radii.size(), 0.0);
        std::vector<int> node_count(radii.size(), 0);

        auto addConnections = [&](const std::map<SegmentRef, float> &segconns, double depth_to_line,
                                  size_t coverage) {
            for (auto &segconn : segconns) {
                if (!covered[segconn.first.ref]) {
                    double angle = depth_to_line + segconn.second;
                    size_t rbin = coverage;
                    while (rbin != radii.size() && radii[rbin] != -1 && angle > radii[rbin]) {
                        rbin++;
                    }
                    if (rbin != radii.size()) {
                        anglebins.push(float(angle), SegmentData(segconn.first, SegmentRef(), 0, 0.0, rbin));
                    }
                }
            }
        };

        // node_count includes this one, but will be added in next algo:
        float key;
        SegmentData lineindex;
        while (!anglebins.empty()) {
            anglebins.pop(key, lineindex);
            if (!covered[lineindex.ref]) {
                covered[lineindex.ref] = 1;
                double depth_to_line = key;
                total_depth[lineindex.coverage] += depth_to_line;
                node_count[lineindex.coverage] += 1;
                const Connector &line = connections[lineindex.ref];
                if (lineindex.dir != -1) {
                    addConnections(line.m_forward_segconns, depth_to_line, lineindex.coverage);
                }
                if (lineindex.dir != 1) {
                    addConnections(line.m_back_segconns, depth_to_line, lineindex.coverage);
                }
            }
        }
        // set the attributes for this node:
        int curs_node_count = 0;
        double curs_total_depth = 0.0;
        for (size_t r = 0; r < radii.size(); r++) {
            curs_node_count += node_count[r];
            curs_total_depth += total_depth[r];
            out.emplace_back(count_col[r], float(curs_node_count));
            if (curs_node_count > 1) {
                // note -- node_count includes this one -- mean depth as per p.108 Social Logic of Space
                double mean_depth = curs_total_depth / double(curs_node_count - 1);
                out.emplace_back(depth_col[r], float(mean_depth));
                out.emplace_back(total_col[r], float(curs_total_depth));
            } else {
                out.emplace_back(depth_col[r], -1);
                out.emplace_back(total_col[r], -1);
            }
        }
    };

    depthmapX::runAnalysisInParallel(comm, rows.size(), nthreads, [&](size_t index, unsigned int thread) {
        analyseSource(index, states[thread], results[index]);
    });

    for (size_t index = 0; index < rows.size(); index++) {
        for (const auto &value : results[index]) {
            rows[index]->setValue(value.first, value.second);
        }
    }

    map.setDisplayedAttribute(-2); // <- override if it's already showing
    map.setDisplayedAttribute(depth_col.back());

    return true;
}
