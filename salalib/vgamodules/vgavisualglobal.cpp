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

#include "salalib/vgamodules/vgavisualglobal.h"

#include "salalib/parallelanalysis.h"

#include "genlib/stringutils.h"

bool VGAVisualGlobal::run(Communicator *comm, PointMap &map, bool simple_version) {
    time_t atime = 0;
    if (comm) {
        qtimer(atime, 0);
        comm->CommPostMessage(Communicator::NUM_RECORDS, map.getFilledPointCount());
    }
    AttributeTable &attributes = map.getAttributeTable();

    int entropy_col = -1, rel_entropy_col = -1, integ_dv_col = -1, integ_pv_col = -1, integ_tk_col = -1,
        depth_col = -1, count_col = -1;
    std::string radius_text;
    if (m_radius != -1) {
        radius_text = std::string(" R") + dXstring::formatString(int(m_radius), "%d");
    }

    // n.b. these must be entered in alphabetical order to preserve col indexing:
    // dX simple version test // TV
#ifndef _COMPILE_dX_SIMPLE_VERSION
    if (!simple_version) {
        std::string entropy_col_text = std::string("Visual Entropy") + radius_text;
        entropy_col = attributes.insertOrResetColumn(entropy_col_text.c_str());
    }
#endif

    std::string integ_dv_col_text = std::string("Visual Integration [HH]") + radius_text;
    integ_dv_col = attributes.insertOrResetColumn(integ_dv_col_text.c_str());

#ifndef _COMPILE_dX_SIMPLE_VERSION
    if (!simple_version) {
        std::string integ_pv_col_text = std::string("Visual Integration [P-value]") + radius_text;
        integ_pv_col = attributes.insertOrResetColumn(integ_pv_col_text.c_str());
        std::string integ_tk_col_text = std::string("Visual Integration [Tekl]") + radius_text;
        integ_tk_col = attributes.insertOrResetColumn(integ_tk_col_text.c_str());
        std::string depth_col_text = std::string("Visual Mean Depth") + radius_text;
        depth_col = attributes.insertOrResetColumn(depth_col_text.c_str());
        std::string count_col_text = std::string("Visual Node Count") + radius_text;
        count_col = attributes.insertOrResetColumn(count_col_text.c_str());
        std::string rel_entropy_col_text = std::string("Visual Relativised Entropy") + radius_text;
        rel_entropy_col = attributes.insertOrResetColumn(rel_entropy_col_text.c_str());
    }
#endif

    // Collect the source points up-front (in the original column-major order) so that the
    // searches can be shared between threads. Skipped points still count towards progress.
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

    const unsigned int nthreads = depthmapX::threadsFor(sources.size());
    std::vector<SearchState> states;
    states.reserve(nthreads);
    for (unsigned int t = 0; t < nthreads; t++) {
        states.emplace_back(map.getRows(), map.getCols());
    }
    std::vector<depthmapX::PendingValues> results(sources.size());

    depthmapX::runAnalysisInParallel(
        comm, sources.size(), nthreads,
        [&](size_t index, unsigned int thread) {
            analyseSource(map, sources[index], states[thread], results[index], simple_version, entropy_col,
                          rel_entropy_col, integ_dv_col, integ_pv_col, integ_tk_col, depth_col, count_col);
        },
        skipped);

    for (size_t index = 0; index < sources.size(); index++) {
        AttributeRow &row = attributes.getRow(AttributeKey(sources[index]));
        for (const auto &value : results[index]) {
            row.setValue(value.first, value.second);
        }
    }

    // Historically the search state of the last analysed point was left behind in the point map.
    // Re-create it so that anything relying on m_misc / m_extent sees exactly what it used to.
    SearchState finalState(map.getRows(), map.getCols());
    if (!sources.empty()) {
        depthmapX::PendingValues ignored;
        analyseSource(map, sources.back(), finalState, ignored, simple_version, entropy_col, rel_entropy_col,
                      integ_dv_col, integ_pv_col, integ_tk_col, depth_col, count_col);
    }
    for (size_t i = 0; i < map.getCols(); i++) {
        for (size_t j = 0; j < map.getRows(); j++) {
            PixelRef curs = PixelRef(static_cast<short>(i), static_cast<short>(j));
            map.getPoint(curs).m_misc = finalState.miscs(j, i);
            map.getPoint(curs).m_extent = finalState.extents(j, i);
        }
    }
    map.setDisplayedAttribute(integ_dv_col);

    return true;
}

VGAVisualGlobal::SearchState::SearchState(size_t rows, size_t cols) : miscs(rows, cols), extents(rows, cols) {
    for (size_t ii = 0; ii < cols; ii++) {
        for (size_t jj = 0; jj < rows; jj++) {
            miscs(jj, ii) = 0;
            extents(jj, ii) = PixelRef(static_cast<short>(ii), static_cast<short>(jj));
        }
    }
}

void VGAVisualGlobal::SearchState::reset() {
    // Only the cells written by the previous search need restoring; this replaces a full
    // clear of the whole grid for every single source point.
    for (const PixelRef &pix : touched) {
        miscs(pix.y, pix.x) = 0;
        extents(pix.y, pix.x) = pix;
    }
    touched.clear();
}

void VGAVisualGlobal::analyseSource(PointMap &map, PixelRef curs, SearchState &state,
                                    depthmapX::PendingValues &out, bool simple_version, int entropy_col,
                                    int rel_entropy_col, int integ_dv_col, int integ_pv_col, int integ_tk_col,
                                    int depth_col, int count_col) {
    state.reset();
    depthmapX::RowMatrix<int> &miscs = state.miscs;
    state.touched.push_back(curs);

    int total_depth = 0;
    int total_nodes = 0;

    std::vector<int> distribution;
    std::vector<PixelRefVector> search_tree;
    search_tree.push_back(PixelRefVector());
    search_tree.back().push_back(curs);

    int level = 0;
    while (search_tree[level].size()) {
        search_tree.push_back(PixelRefVector());
        const PixelRefVector &searchTreeAtLevel = search_tree[level];
        distribution.push_back(0);
        for (auto currLvlIter = searchTreeAtLevel.rbegin(); currLvlIter != searchTreeAtLevel.rend(); currLvlIter++) {
            int &pmisc = miscs(currLvlIter->y, currLvlIter->x);
            Point &p = map.getPoint(*currLvlIter);
            if (p.filled() && pmisc != ~0) {
                total_depth += level;
                total_nodes += 1;
                distribution.back() += 1;
                if ((int)m_radius == -1 ||
                    (level < (int)m_radius && (!p.contextfilled() || currLvlIter->iseven()))) {
                    extractUnseen(p.getNode(), search_tree[level + 1], state);
                    pmisc = ~0;
                    if (!p.getMergePixel().empty()) {
                        PixelRef mergePixel = p.getMergePixel();
                        int &p2misc = miscs(mergePixel.y, mergePixel.x);
                        Point &p2 = map.getPoint(mergePixel);
                        if (p2misc != ~0) {
                            extractUnseen(p2.getNode(), search_tree[level + 1], state); // did say p.misc
                            p2misc = ~0;
                            state.touched.push_back(mergePixel);
                        }
                    }
                } else {
                    pmisc = ~0;
                }
            }
            search_tree[level].pop_back();
        }
        level++;
    }
    auto setValue = [&out](int col, float value) { out.emplace_back(col, value); };
    // only set to single float precision after divide
    // note -- total_nodes includes this one -- mean depth as per p.108 Social Logic of Space
    if (!simple_version) {
        setValue(count_col, float(total_nodes)); // note: total nodes includes this one
    }
    // ERROR !!!!!!
    if (total_nodes > 1) {
        double mean_depth = double(total_depth) / double(total_nodes - 1);
        if (!simple_version) {
            setValue(depth_col, float(mean_depth));
        }
        // total nodes > 2 to avoid divide by 0 (was > 3)
        if (total_nodes > 2 && mean_depth > 1.0) {
            double ra = 2.0 * (mean_depth - 1.0) / double(total_nodes - 2);
            // d-value / p-values from Depthmap 4 manual, note: node_count includes this one
            double rra_d = ra / dvalue(total_nodes);
            double rra_p = ra / pvalue(total_nodes);
            double integ_tk = teklinteg(total_nodes, total_depth);
            setValue(integ_dv_col, float(1.0 / rra_d));
            if (!simple_version) {
                setValue(integ_pv_col, float(1.0 / rra_p));
            }
            if (total_depth - total_nodes + 1 > 1) {
                if (!simple_version) {
                    setValue(integ_tk_col, float(integ_tk));
                }
            } else {
                if (!simple_version) {
                    setValue(integ_tk_col, -1.0f);
                }
            }
        } else {
            setValue(integ_dv_col, (float)-1);
            if (!simple_version) {
                setValue(integ_pv_col, (float)-1);
                setValue(integ_tk_col, (float)-1);
            }
        }
        double entropy = 0.0, rel_entropy = 0.0, factorial = 1.0;
        // n.b., this distribution contains the root node itself in distribution[0]
        // -> chopped from entropy to avoid divide by zero if only one node
        for (size_t k = 1; k < distribution.size(); k++) {
            if (distribution[k] > 0) {
                double prob = double(distribution[k]) / double(total_nodes - 1);
                entropy -= prob * pafLog2(prob);
                // Formula from Turner 2001, "Depthmap"
                factorial *= double(k + 1);
                double q = (pow(mean_depth, double(k)) / double(factorial)) * exp(-mean_depth);
                rel_entropy += (float)prob * pafLog2(prob / q);
            }
        }
        if (!simple_version) {
            setValue(entropy_col, float(entropy));
            setValue(rel_entropy_col, float(rel_entropy));
        }
    } else {
        if (!simple_version) {
            setValue(depth_col, (float)-1);
            setValue(entropy_col, (float)-1);
            setValue(rel_entropy_col, (float)-1);
        }
    }
}

void VGAVisualGlobal::extractUnseen(Node &node, PixelRefVector &pixels, SearchState &state) {
    depthmapX::RowMatrix<int> &miscs = state.miscs;
    depthmapX::RowMatrix<PixelRef> &extents = state.extents;
    for (int i = 0; i < 32; i++) {
        Bin &bin = node.bin(i);
        for (auto pixVec : bin.m_pixel_vecs) {
            for (PixelRef pix = pixVec.start(); pix.col(bin.m_dir) <= pixVec.end().col(bin.m_dir);) {
                int &misc = miscs(pix.y, pix.x);
                PixelRef &extent = extents(pix.y, pix.x);
                if (misc == 0) {
                    pixels.push_back(pix);
                    misc |= (1 << i);
                    state.touched.push_back(pix);
                }
                // 10.2.02 revised --- diagonal was breaking this as it was extent in diagonal or horizontal
                if (!(bin.m_dir & PixelRef::DIAGONAL)) {
                    if (extent.col(bin.m_dir) >= pixVec.end().col(bin.m_dir))
                        break;
                    extent.col(bin.m_dir) = pixVec.end().col(bin.m_dir);
                    state.touched.push_back(pix);
                }
                pix.move(bin.m_dir);
            }
        }
    }
}
