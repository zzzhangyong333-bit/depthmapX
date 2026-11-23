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

#include "salalib/segmmodules/segmtulip.h"

#include "salalib/parallelanalysis.h"

#include "genlib/stringutils.h"

bool SegmentTulip::run(Communicator *comm, ShapeGraph &map, bool) {

    if (map.getMapType() != ShapeMap::SEGMENTMAP) {
        return false;
    }

    // TODO: Understand what these parameters do. They were never truly provided in the original function
    int weighting_col2 = m_weighted_measure_col2;
    int routeweight_col = m_routeweight_col;
    bool interactive = m_interactive;

    AttributeTable &attributes = map.getAttributeTable();

    int processed_rows = 0;

    time_t atime = 0;

    if (comm) {
        qtimer(atime, 0);
        comm->CommPostMessage(Communicator::NUM_RECORDS,
                              (m_sel_only ? map.getSelSet().size() : map.getConnections().size()));
    }

    // note: radius must be sorted lowest to highest, but if -1 occurs ("radius n") it needs to be last...
    // ...to ensure no mess ups, we'll re-sort here:
    bool radius_n = false;
    std::vector<double> radius_unconverted;
    for (int radius : m_radius_set) {
        if (radius == -1.0) {
            radius_n = true;
        } else {
            radius_unconverted.push_back(radius);
        }
    }
    if (radius_n) {
        radius_unconverted.push_back(-1.0);
    }

    // retrieve weighted col data, as this may well be overwritten in the new analysis:
    std::vector<float> weights;
    std::vector<float> routeweights; // EF
    std::string weighting_col_text;

    int tulip_bins = m_tulip_bins;

    if (m_weighted_measure_col != -1) {
        weighting_col_text = attributes.getColumnName(m_weighted_measure_col);
        for (size_t i = 0; i < map.getConnections().size(); i++) {
            weights.push_back(map.getAttributeRowFromShapeIndex(i).getValue(m_weighted_measure_col));
        }
    } else { // Normal run // TV
        for (size_t i = 0; i < map.getConnections().size(); i++) {
            weights.push_back(1.0f);
        }
    }
    // EF routeweight*
    std::string routeweight_col_text;
    if (routeweight_col != -1) {
        // we normalise the column values between 0 and 1 and reverse it so that high values can be treated as a 'low
        // cost' - similar to the angular cost
        double max_value = attributes.getColumn(routeweight_col).getStats().max;
        routeweight_col_text = attributes.getColumnName(routeweight_col);
        for (size_t i = 0; i < map.getConnections().size(); i++) {
            routeweights.push_back(1.0 - (map.getAttributeRowFromShapeIndex(i).getValue(routeweight_col) /
                                          max_value)); // scale and revert!
        }
    } else { // Normal run // TV
        for (size_t i = 0; i < map.getConnections().size(); i++) {
            routeweights.push_back(1.0f);
        }
    }
    //*EF routeweight

    // EFEF*
    // for origin-destination weighting
    std::vector<float> weights2;
    std::string weighting_col_text2;
    if (weighting_col2 != -1) {
        weighting_col_text2 = attributes.getColumnName(weighting_col2);
        for (size_t i = 0; i < map.getConnections().size(); i++) {
            weights2.push_back(map.getAttributeRowFromShapeIndex(i).getValue(weighting_col2));
        }
    } else { // Normal run // TV
        for (size_t i = 0; i < map.getConnections().size(); i++) {
            weights2.push_back(1.0f);
        }
    }
    //*EFEF

    std::string tulip_text = std::string("T") + dXstring::formatString(tulip_bins, "%d");

    // first enter the required attribute columns:
    size_t r;
    for (r = 0; r < radius_unconverted.size(); r++) {
        std::string radius_text = makeRadiusText(m_radius_type, radius_unconverted[r]);
        if (m_choice) {
            // EF routeweight *
            if (routeweight_col != -1) {
                std::string choice_col_text =
                    tulip_text + " Choice [Route weight by " + routeweight_col_text + "]" + radius_text;
                attributes.insertOrResetColumn(choice_col_text.c_str());
                if (m_weighted_measure_col != -1) {
                    std::string w_choice_col_text = tulip_text + " Choice [[Route weight by " + routeweight_col_text +
                                                    "][" + weighting_col_text + " Wgt]]" + radius_text;

                    attributes.insertOrResetColumn(w_choice_col_text.c_str());
                }
                // EFEF*
                if (weighting_col2 != -1) {
                    std::string w_choice_col_text2 = tulip_text + " Choice [[Route weight by " + routeweight_col_text +
                                                     "][" + weighting_col_text + "-" + weighting_col_text2 + " Wgt]]" +
                                                     radius_text;

                    attributes.insertOrResetColumn(w_choice_col_text2.c_str());
                }
                //*EFEF
            }
            //*EF routeweight
            else { // Normal run // TV
                std::string choice_col_text = tulip_text + " Choice" + radius_text;
                attributes.insertOrResetColumn(choice_col_text.c_str());
                if (m_weighted_measure_col != -1) {
                    std::string w_choice_col_text =
                        tulip_text + " Choice [" + weighting_col_text + " Wgt]" + radius_text;
                    attributes.insertOrResetColumn(w_choice_col_text.c_str());
                }
                // EFEF*
                if (weighting_col2 != -1) {
                    std::string w_choice_col_text2 = tulip_text + " Choice [" + weighting_col_text + "-" +
                                                     weighting_col_text2 + " Wgt]" + radius_text;
                    attributes.insertOrResetColumn(w_choice_col_text2.c_str());
                }
                //*EFEF
            }
        }

        // EF routeweight *
        if (routeweight_col != -1) {
            std::string integ_col_text = tulip_text + " Integration [Route weight by " + routeweight_col_text + "]" +
                                         radius_text; // <- note, the fact this is a tulip is unnecessary
            std::string w_integ_col_text = tulip_text + " Integration [[Route weight by " + routeweight_col_text +
                                           "][" + weighting_col_text + " Wgt]]" + radius_text;

            std::string count_col_text = tulip_text + " Node Count [Route weight by " + routeweight_col_text + "]" +
                                         radius_text; // <- note, the fact this is a tulip is unnecessary
            std::string td_col_text = tulip_text + " Total Depth [Route weight by " + routeweight_col_text + "]" +
                                      radius_text; // <- note, the fact this is a tulip is unnecessary
            // '[' comes after 'R' in ASCII, so this column will come after Mean Depth R...
            std::string w_td_text = tulip_text + " Total Depth [[Route weight by " + routeweight_col_text + "][" +
                                    weighting_col_text + " Wgt]]" + radius_text;
            std::string total_weight_text = tulip_text + " Total " + weighting_col_text + " [Route weight by " +
                                            routeweight_col_text + "]" + radius_text;

            attributes.insertOrResetColumn(integ_col_text.c_str());
            attributes.insertOrResetColumn(count_col_text.c_str());
            attributes.insertOrResetColumn(td_col_text.c_str());
            if (m_weighted_measure_col != -1) {
                attributes.insertOrResetColumn(w_integ_col_text.c_str());
                attributes.insertOrResetColumn(w_td_text.c_str());
                attributes.insertOrResetColumn(total_weight_text.c_str());
            }
        }
        //*EF routeweight
        else { // Normal run // TV
            std::string integ_col_text =
                tulip_text + " Integration" + radius_text; // <- note, the fact this is a tulip is unnecessary
            std::string w_integ_col_text = tulip_text + " Integration [" + weighting_col_text + " Wgt]" + radius_text;

            std::string count_col_text =
                tulip_text + " Node Count" + radius_text; // <- note, the fact this is a tulip is unnecessary
            std::string td_col_text =
                tulip_text + " Total Depth" + radius_text; // <- note, the fact this is a tulip is unnecessary
            // '[' comes after 'R' in ASCII, so this column will come after Mean Depth R...
            std::string w_td_text = tulip_text + " Total Depth [" + weighting_col_text + " Wgt]" + radius_text;
            std::string total_weight_text = tulip_text + " Total " + weighting_col_text + radius_text;

            attributes.insertOrResetColumn(integ_col_text.c_str());
            attributes.insertOrResetColumn(count_col_text.c_str());
            attributes.insertOrResetColumn(td_col_text.c_str());
            if (m_weighted_measure_col != -1) {
                attributes.insertOrResetColumn(w_integ_col_text.c_str());
                attributes.insertOrResetColumn(w_td_text.c_str());
                attributes.insertOrResetColumn(total_weight_text.c_str());
            }
        }
    }
    std::vector<int> choice_col, w_choice_col, w_choice_col2, count_col, integ_col, w_integ_col, td_col, w_td_col,
        total_weight_col;
    // then look them up! eek....
    for (r = 0; r < radius_unconverted.size(); r++) {
        std::string radius_text = makeRadiusText(m_radius_type, radius_unconverted[r]);
        if (m_choice) {
            // EF routeweight *
            if (routeweight_col != -1) {
                std::string choice_col_text =
                    tulip_text + " Choice [Route weight by " + routeweight_col_text + "]" + radius_text;
                choice_col.push_back(attributes.getColumnIndex(choice_col_text.c_str()));
                if (m_weighted_measure_col != -1) {
                    std::string w_choice_col_text = tulip_text + " Choice [[Route weight by " + routeweight_col_text +
                                                    "][" + weighting_col_text + " Wgt]]" + radius_text;
                    w_choice_col.push_back(attributes.getColumnIndex(w_choice_col_text.c_str()));
                }
                // EFEF*
                if (weighting_col2 != -1) {
                    std::string w_choice_col_text2 = tulip_text + " Choice [[Route weight by " + routeweight_col_text +
                                                     "][" + weighting_col_text + "-" + weighting_col_text2 + " Wgt]]" +
                                                     radius_text;
                    w_choice_col2.push_back(attributes.getColumnIndex(w_choice_col_text2.c_str()));
                }
                //*EFEF
            }
            //* EF routeweight
            else { // Normal run // TV
                std::string choice_col_text = tulip_text + " Choice" + radius_text;
                choice_col.push_back(attributes.getColumnIndex(choice_col_text.c_str()));
                if (m_weighted_measure_col != -1) {
                    std::string w_choice_col_text =
                        tulip_text + " Choice [" + weighting_col_text + " Wgt]" + radius_text;
                    w_choice_col.push_back(attributes.getColumnIndex(w_choice_col_text.c_str()));
                }
                // EFEF*
                if (weighting_col2 != -1) {
                    std::string w_choice_col_text2 = tulip_text + " Choice [" + weighting_col_text + "-" +
                                                     weighting_col_text2 + " Wgt]" + radius_text;
                    w_choice_col2.push_back(attributes.getColumnIndex(w_choice_col_text2.c_str()));
                }
                //*EFEF
            }
        }
        // EF routeweight *
        if (routeweight_col != -1) {
            std::string integ_col_text = tulip_text + " Integration [Route weight by " + routeweight_col_text + "]" +
                                         radius_text; // <- note, the fact this is a tulip is unnecessary
            std::string w_integ_col_text = tulip_text + " Integration [[Route weight by " + routeweight_col_text +
                                           "][" + weighting_col_text + " Wgt]]" + radius_text;

            std::string count_col_text = tulip_text + " Node Count [Route weight by " + routeweight_col_text + "]" +
                                         radius_text; // <- note, the fact this is a tulip is unnecessary
            std::string td_col_text = tulip_text + " Total Depth [Route weight by " + routeweight_col_text + "]" +
                                      radius_text; // <- note, the fact this is a tulip is unnecessary
            std::string w_td_text = tulip_text + " Total Depth [[Route weight by " + routeweight_col_text + "][" +
                                    weighting_col_text + " Wgt]]" + radius_text;
            std::string total_weight_col_text = tulip_text + " Total " + weighting_col_text + " [Route weight by " +
                                                routeweight_col_text + "]" + radius_text;

            integ_col.push_back(attributes.getColumnIndex(integ_col_text.c_str()));
            count_col.push_back(attributes.getColumnIndex(count_col_text.c_str()));
            td_col.push_back(attributes.getColumnIndex(td_col_text.c_str()));
            if (m_weighted_measure_col != -1) {
                // '[' comes after 'R' in ASCII, so this column will come after Mean Depth R...
                w_integ_col.push_back(attributes.getColumnIndex(w_integ_col_text.c_str()));
                w_td_col.push_back(attributes.getColumnIndex(w_td_text.c_str()));
                total_weight_col.push_back(attributes.getColumnIndex(total_weight_col_text.c_str()));
            }
        }
        //* EF routeweight
        else { // Normal run // TV
            std::string integ_col_text =
                tulip_text + " Integration" + radius_text; // <- note, the fact this is a tulip is unnecessary
            std::string w_integ_col_text = tulip_text + " Integration [" + weighting_col_text + " Wgt]" + radius_text;

            std::string count_col_text =
                tulip_text + " Node Count" + radius_text; // <- note, the fact this is a tulip is unnecessary
            std::string td_col_text =
                tulip_text + " Total Depth" + radius_text; // <- note, the fact this is a tulip is unnecessary
            std::string w_td_text = tulip_text + " Total Depth [" + weighting_col_text + " Wgt]" + radius_text;
            std::string total_weight_col_text = tulip_text + " Total " + weighting_col_text + radius_text;

            integ_col.push_back(attributes.getColumnIndex(integ_col_text.c_str()));
            count_col.push_back(attributes.getColumnIndex(count_col_text.c_str()));
            td_col.push_back(attributes.getColumnIndex(td_col_text.c_str()));
            if (m_weighted_measure_col != -1) {
                // '[' comes after 'R' in ASCII, so this column will come after Mean Depth R...
                w_integ_col.push_back(attributes.getColumnIndex(w_integ_col_text.c_str()));
                w_td_col.push_back(attributes.getColumnIndex(w_td_text.c_str()));
                total_weight_col.push_back(attributes.getColumnIndex(total_weight_col_text.c_str()));
            }
        }
    }

    tulip_bins /= 2; // <- actually use semicircle of tulip bins
    tulip_bins += 1;

    std::vector<double> radius;
    for (r = 0; r < radius_unconverted.size(); r++) {
        if (m_radius_type == Options::RADIUS_ANGULAR && radius_unconverted[r] != -1) {
            radius.push_back(floor(radius_unconverted[r] * tulip_bins * 0.5));
        } else {
            radius.push_back(radius_unconverted[r]);
        }
    }
    // entered once for each segment
    int length_col = attributes.getColumnIndex("Segment Length");
    std::vector<float> lengths;
    if (length_col != -1) {
        for (size_t i = 0; i < map.getConnections().size(); i++) {
            AttributeRow& row = map.getAttributeRowFromShapeIndex(i);
            lengths.push_back(row.getValue(length_col));
        }
    }

    int radiussize = radius.size();
    int radiusmask = 0;
    for (int i = 0; i < radiussize; i++) {
        radiusmask |= (1 << i);
    }

    const size_t nconn = map.getConnections().size();

    // Per-thread search data (formerly heap-allocated 3D arrays shared by every search).
    struct TulipState {
        int radiussize;
        std::vector<AnalysisInfo> audit;
        std::vector<unsigned int> uncovered;
        std::vector<std::vector<SegmentData>> bins;
        TulipState(size_t n, int rs, int tb) : radiussize(rs), audit(n * size_t(rs) * 2), uncovered(n * 2), bins(tb) {}
        AnalysisInfo &at(size_t j, int k, int dir) { return audit[(j * size_t(radiussize) + size_t(k)) * 2 + size_t(dir)]; }
        unsigned int &unc(size_t j, int dir) { return uncovered[j * 2 + size_t(dir)]; }
    };

    std::vector<AttributeRow *> rows(nconn);
    std::vector<size_t> sources;
    for (size_t cursor = 0; cursor < nconn; cursor++) {
        rows[cursor] = &map.getAttributeRowFromShapeIndex(cursor);
        // could use m_selection_set.searchindex(rowid) to find
        // if this row is selected as m_selection_set is ordered for axial and segment maps, etc
        // BUT, actually quicker to check the tag in the attributes that shows it's selected
        if (!m_sel_only || rows[cursor]->isSelected()) {
            sources.push_back(cursor);
        }
    }

    auto analyseSource = [&](size_t cursor, TulipState &state, depthmapX::PendingValues &out) {
        auto &bins = state.bins;
        const AttributeRow &row = *rows[cursor];
        auto setValue = [&out](int col, float value) { out.emplace_back(col, value); };

        for (int k = 0; k < tulip_bins; k++) {
            bins[k].clear();
        }
        for (size_t j = 0; j < nconn; j++) {
            for (int dir = 0; dir < 2; dir++) {
                for (int k = 0; k < radiussize; k++) {
                    state.at(j, k, dir).clearLine();
                }
                state.unc(j, dir) = radiusmask;
            }
        }

        double rootseglength = row.getValue(length_col);
        double rootweight = (m_weighted_measure_col != -1) ? weights[cursor] : 0.0;

        // setup: direction 0 (both ways), segment i, previous -1, segdepth (step depth) 0, metricdepth 0.5 *
        // rootseglength, bin 0
        SegmentData segmentData(0, cursor, SegmentRef(), 0, 0.5 * rootseglength, radiusmask);
        auto it = std::lower_bound(bins[0].begin(), bins[0].end(), segmentData);
        if (it == bins[0].end() || segmentData != *it) {
            bins[0].insert(it, segmentData);
        }
        // this version below is only designed to be used temporarily --
        // could be on an option?
        // bins[0].push_back(SegmentData(0,rowid,SegmentRef(),0,0.0,radiusmask));
        int depthlevel = 0;
        int opencount = 1;
        size_t currentbin = 0;
        while (opencount) {
            while (!bins[currentbin].size()) {
                depthlevel++;
                currentbin++;
                if (currentbin == static_cast<size_t>(tulip_bins)) {
                    currentbin = 0;
                }
            }
            SegmentData lineindex = bins[currentbin].back();
            bins[currentbin].pop_back();
            //
            opencount--;

            int ref = lineindex.ref;
            int dir = (lineindex.dir == 1) ? 0 : 1;
            int coverage = lineindex.coverage & state.unc(ref, dir);
            if (coverage != 0) {
                int rbin = 0;
                int rbinbase;
                if (lineindex.previous.ref != -1) {
                    state.unc(ref, dir) &= ~coverage;
                    while (((coverage >> rbin) & 0x1) == 0)
                        rbin++;
                    rbinbase = rbin;
                    while (rbin < radiussize) {
                        if (((coverage >> rbin) & 0x1) == 1) {
                            state.at(ref, rbin, dir).depth = depthlevel;
                            state.at(ref, rbin, dir).previous = lineindex.previous;
                            state.at(lineindex.previous.ref, rbin, (lineindex.previous.dir == 1) ? 0 : 1).leaf =
                                false;
                        }
                        rbin++;
                    }
                } else {
                    rbinbase = 0;
                    state.unc(ref, 0) &= ~coverage;
                    state.unc(ref, 1) &= ~coverage;
                }
                const Connector &line = map.getConnections()[ref];
                float seglength;
                int extradepth;
                if (lineindex.dir != -1) {
                    for (auto &segconn : line.m_forward_segconns) {
                        rbin = rbinbase;
                        SegmentRef conn = segconn.first;
                        if ((state.unc(conn.ref, (conn.dir == 1 ? 0 : 1)) & coverage) != 0) {
                            // EF routeweight*
                            if (routeweight_col != -1) { // EF here we do the weighting of the angular cost by the
                                                         // weight of the next segment
                                // note that the content of the routeweights array is scaled between 0 and 1 and is
                                // reversed
                                // such that: = 1.0-(attributes.getValue(i, routeweight_col)/max_value)
                                extradepth = (int)floor(segconn.second * tulip_bins * 0.5 * routeweights[conn.ref]);
                            }
                            //*EF routeweight
                            else {
                                extradepth = (int)floor(segconn.second * tulip_bins * 0.5);
                            }
                            seglength = lengths[conn.ref];
                            switch (m_radius_type) {
                            case Options::RADIUS_ANGULAR:
                                while (rbin != radiussize && radius[rbin] != -1 &&
                                       depthlevel + extradepth > (int)radius[rbin]) {
                                    rbin++;
                                }
                                break;
                            case Options::RADIUS_METRIC:
                                while (rbin != radiussize && radius[rbin] != -1 &&
                                       lineindex.metricdepth + seglength * 0.5 > radius[rbin]) {
                                    rbin++;
                                }
                                break;
                            case Options::RADIUS_STEPS:
                                if (rbin != radiussize && radius[rbin] != -1 &&
                                    lineindex.segdepth >= (int)radius[rbin]) {
                                    rbin++;
                                }
                                break;
                            }
                            if ((coverage >> rbin) != 0) {
                                SegmentData sd(conn, SegmentRef(1, lineindex.ref), lineindex.segdepth + 1,
                                               lineindex.metricdepth + seglength, (coverage >> rbin) << rbin);
                                size_t bin = (currentbin + tulip_bins + extradepth) % tulip_bins;
                                depthmapX::insert_sorted(bins[bin], sd);
                                opencount++;
                            }
                        }
                    }
                }
                if (lineindex.dir != 1) {
                    for (auto &segconn : line.m_back_segconns) {
                        rbin = rbinbase;
                        SegmentRef conn = segconn.first;
                        if ((state.unc(conn.ref, (conn.dir == 1 ? 0 : 1)) & coverage) != 0) {
                            // EF routeweight*
                            if (routeweight_col != -1) { // EF here we do the weighting of the angular cost by the
                                                         // weight of the next segment
                                // note that the content of the routeweights array is scaled between 0 and 1 and is
                                // reversed
                                // such that: = 1.0-(attributes.getValue(i, routeweight_col)/max_value)
                                extradepth = (int)floor(segconn.second * tulip_bins * 0.5 * routeweights[conn.ref]);
                            }
                            //*EF routeweight
                            else {
                                extradepth = (int)floor(segconn.second * tulip_bins * 0.5);
                            }
                            seglength = lengths[conn.ref];
                            switch (m_radius_type) {
                            case Options::RADIUS_ANGULAR:
                                while (rbin != radiussize && radius[rbin] != -1 &&
                                       depthlevel + extradepth > (int)radius[rbin]) {
                                    rbin++;
                                }
                                break;
                            case Options::RADIUS_METRIC:
                                while (rbin != radiussize && radius[rbin] != -1 &&
                                       lineindex.metricdepth + seglength * 0.5 > radius[rbin]) {
                                    rbin++;
                                }
                                break;
                            case Options::RADIUS_STEPS:
                                if (rbin != radiussize && radius[rbin] != -1 &&
                                    lineindex.segdepth >= (int)radius[rbin]) {
                                    rbin++;
                                }
                                break;
                            }
                            if ((coverage >> rbin) != 0) {
                                SegmentData sd(conn, SegmentRef(-1, lineindex.ref), lineindex.segdepth + 1,
                                               lineindex.metricdepth + seglength, (coverage >> rbin) << rbin);
                                size_t bin = (currentbin + tulip_bins + extradepth) % tulip_bins;
                                depthmapX::insert_sorted(bins[bin], sd);
                                opencount++;
                            }
                        }
                    }
                }
            }
        }
        // set the attributes for this node:
        for (int k = 0; k < radiussize; k++) {
            // note, curs_total_depth must use double as mantissa can get too long for int in large systems
            double curs_node_count = 0.0, curs_total_depth = 0.0;
            double curs_total_weight = 0.0, curs_total_weighted_depth = 0.0;
            size_t j;
            for (j = 0; j < nconn; j++) {
                // find dir according
                bool m0 = ((state.unc(j, 0) >> k) & 0x1) == 0;
                bool m1 = ((state.unc(j, 1) >> k) & 0x1) == 0;
                if ((m0 | m1) != 0) {
                    int dir;
                    if (m0 & m1) {
                        // dir is the one with the lowest depth:
                        if (state.at(j, k, 0).depth < state.at(j, k, 1).depth)
                            dir = 0;
                        else
                            dir = 1;
                    } else {
                        // dir is simply the one that's filled in:
                        dir = m0 ? 0 : 1;
                    }
                    curs_node_count++;
                    curs_total_depth += state.at(j, k, dir).depth;
                    curs_total_weight += weights[j];
                    curs_total_weighted_depth += state.at(j, k, dir).depth * weights[j];
                    //
                    if (m_choice && state.at(j, k, dir).leaf) {
                        // note, graph may be directed (e.g., for one way streets), so both ways must be included from
                        // now on:
                        SegmentRef here = SegmentRef(dir == 0 ? 1 : -1, j);
                        if (here.ref != static_cast<int>(cursor)) {
                            int choicecount = 0;
                            double choiceweight = 0.0;
                            // EFEF*
                            double choiceweight2 = 0.0;
                            //*EFEF
                            while (here.ref != static_cast<int>(cursor)) { // not rowid means not the current root for the path
                                int heredir = (here.dir == 1) ? 0 : 1;
                                // each node has the existing choicecount and choiceweight from previously encountered
                                // nodes added to it
                                state.at(here.ref, k, heredir).choice += choicecount;
                                // nb, weighted values calculated anyway to save time on 'if'
                                state.at(here.ref, k, heredir).weighted_choice += choiceweight;
                                // EFEF*
                                state.at(here.ref, k, heredir).weighted_choice2 += choiceweight2;
                                //*EFEF
                                // if the node hasn't been encountered before, the choicecount and choiceweight is
                                // incremented for all remaining nodes to be encountered on the backwards route from it
                                if (!state.at(here.ref, k, heredir).choicecovered) {
                                    // this node has not been encountered before: this adds the choicecount and weight
                                    // for this node, and flags it as visited
                                    choicecount++;
                                    choiceweight += weights[here.ref] * rootweight;
                                    // EFEF*
                                    choiceweight2 += weights2[here.ref] * rootweight; // rootweight!
                                    //*EFEF

                                    state.at(here.ref, k, heredir).choicecovered = true;
                                    // note, for weighted choice, the start and end points have choice added to them:
                                    if (m_weighted_measure_col != -1) {
                                        state.at(here.ref, k, heredir).weighted_choice +=
                                            (weights[here.ref] * rootweight) / 2.0;
                                        // EFEF*
                                        if (weighting_col2 != -1) {
                                            state.at(here.ref, k, heredir).weighted_choice2 +=
                                                (weights2[here.ref] * rootweight) / 2.0; // rootweight!
                                        }
                                        //*EFEF
                                    }
                                }
                                here = state.at(here.ref, k, heredir).previous;
                            }
                            // note, for weighted choice, the start and end points have choice added to them:
                            // (this is the summed weight for all starting nodes encountered in this path)
                            if (m_weighted_measure_col != -1) {
                                state.at(here.ref, k, (here.dir == 1) ? 0 : 1).weighted_choice += choiceweight / 2.0;
                                // EFEF*
                                if (weighting_col2 != -1) {
                                    state.at(here.ref, k, (here.dir == 1) ? 0 : 1).weighted_choice2 +=
                                        choiceweight2 / 2.0;
                                }
                                //*EFEF
                            }
                        }
                    }
                }
            }
            double total_depth_conv = curs_total_depth / ((tulip_bins - 1.0f) * 0.5f);
            double total_weighted_depth_conv = curs_total_weighted_depth / ((tulip_bins - 1.0f) * 0.5f);
            //
            setValue(count_col[k], float(curs_node_count));
            if (curs_node_count > 1) {
                // for dmap 8 and above, mean depth simply isn't calculated as for radius measures it is meaningless
                setValue(td_col[k], total_depth_conv);
                if (m_weighted_measure_col != -1) {
                    setValue(total_weight_col[k], float(curs_total_weight));
                    setValue(w_td_col[k], float(total_weighted_depth_conv));
                }
            } else {
                setValue(td_col[k], -1);
                if (m_weighted_measure_col != -1) {
                    setValue(total_weight_col[k], -1.0f);
                    setValue(w_td_col[k], -1.0f);
                }
            }
            // for dmap 10 an above, integration is included!
            if (total_depth_conv > 1e-9) {
                setValue(integ_col[k], (float)(curs_node_count * curs_node_count / total_depth_conv));
                if (m_weighted_measure_col != -1) {
                    setValue(w_integ_col[k],
                                 (float)(curs_total_weight * curs_total_weight / total_weighted_depth_conv));
                }
            } else {
                setValue(integ_col[k], -1);
                if (m_weighted_measure_col != -1) {
                    setValue(w_integ_col[k], -1.0f);
                }
            }
        }
    };

    // Each source is an independent search. Choice values are accumulated per thread and summed
    // at the end (integer choice counts are exact; weighted choice can differ from a sequential
    // run only in the last bits of floating point rounding).
    const unsigned int nthreads = depthmapX::threadsFor(sources.size());
    std::vector<TulipState> states;
    states.reserve(nthreads);
    for (unsigned int t = 0; t < nthreads; t++) {
        states.emplace_back(nconn, radiussize, tulip_bins);
    }
    std::vector<depthmapX::PendingValues> results(sources.size());
    std::vector<char> finished(sources.size(), 0);

    bool completed = depthmapX::parallelFor(
        sources.size(), nthreads,
        [&](size_t index, unsigned int thread) {
            analyseSource(sources[index], states[thread], results[index]);
            finished[index] = 1;
        },
        [comm](size_t done) {
            if (comm) {
                if (comm->IsCancelled()) {
                    return false;
                }
                comm->CommPostMessage(Communicator::CURRENT_RECORD, static_cast<int>(done));
            }
            return true;
        });
    if (!completed && interactive) {
        // interactive is usual Depthmap: throw an exception if cancelled
        throw Communicator::CancelledException();
    }
    // in non-interactive mode, retain what's been processed already
    for (size_t index = 0; index < sources.size(); index++) {
        if (!finished[index]) {
            continue;
        }
        for (const auto &value : results[index]) {
            rows[sources[index]]->setValue(value.first, value.second);
        }
        processed_rows++;
    }

    // merge the per-thread choice totals
    auto choiceTotal = [&](size_t cursor, size_t r, int dir, double AnalysisInfo::*field) {
        double total = 0.0;
        for (auto &state : states) {
            total += state.at(cursor, int(r), dir).*field;
        }
        return total;
    };

    if (m_choice) {
        for (size_t cursor = 0; cursor < map.getConnections().size(); cursor++) {
            AttributeRow &row =
                attributes.getRow(AttributeKey(depthmapX::getMapAtIndex(map.getAllShapes(), cursor)->first));
            for (size_t r = 0; r < radius.size(); r++) {
                // according to Eva's correction, total choice and total weighted choice
                // should already have been accumulated by radius at this stage
                double total_choice =
                    choiceTotal(cursor, r, 0, &AnalysisInfo::choice) + choiceTotal(cursor, r, 1, &AnalysisInfo::choice);
                double total_weighted_choice = choiceTotal(cursor, r, 0, &AnalysisInfo::weighted_choice) +
                                               choiceTotal(cursor, r, 1, &AnalysisInfo::weighted_choice);
                // EFEF*
                double total_weighted_choice2 = choiceTotal(cursor, r, 0, &AnalysisInfo::weighted_choice2) +
                                                choiceTotal(cursor, r, 1, &AnalysisInfo::weighted_choice2);
                //*EFEF

                // normalised choice now excluded for two reasons:
                // a) not useful measure, b) in parallel calculations, cannot be calculated at this stage
                // n.b., it is possible through the front end: the new choice takes into account bidirectional routes,
                // so it should be normalised according to (n-1)(n-2) (maximum possible through routes) not
                // (n-1)(n-2)/2 the relativised segment length weighted choice equation was
                // (total_seg_length*total_seg_length-seg_length*seg_length)/2 again, drop the divide by 2 for the new
                // implementation
                //
                //
                row.setValue(choice_col[r], float(total_choice));
                if (m_weighted_measure_col != -1) {
                    row.setValue(w_choice_col[r], float(total_weighted_choice));
                    // EFEF*
                    if (weighting_col2 != -1) {
                        row.setValue(w_choice_col2[r], float(total_weighted_choice2));
                    }
                    //*EFEF
                }
            }
        }
    }
    map.setDisplayedAttribute(-2); // <- override if it's already showing
    if (m_choice) {
        map.setDisplayedAttribute(choice_col.back());
    } else {
        map.setDisplayedAttribute(td_col.back());
    }
    return processed_rows > 0;
}
