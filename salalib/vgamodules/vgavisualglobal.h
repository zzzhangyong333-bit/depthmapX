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

#pragma once

#include "salalib/ivga.h"
#include "salalib/pixelref.h"
#include "salalib/pointdata.h"

#include "genlib/simplematrix.h"

#include <utility>
#include <vector>

class VGAVisualGlobal : IVGA {
  private:
    double m_radius;
    bool m_gates_only;

    // Scratch data for one search; one instance per worker thread.
    struct SearchState {
        depthmapX::RowMatrix<int> miscs;
        depthmapX::RowMatrix<PixelRef> extents;
        std::vector<PixelRef> touched; // cells written since the last reset
        SearchState(size_t rows, size_t cols);
        void reset();
    };
    void analyseSource(PointMap &map, PixelRef curs, SearchState &state, std::vector<std::pair<int, float>> &out,
                       bool simple_version, int entropy_col, int rel_entropy_col, int integ_dv_col,
                       int integ_pv_col, int integ_tk_col, int depth_col, int count_col);
    void extractUnseen(Node &node, PixelRefVector &pixels, SearchState &state);

  public:
    std::string getAnalysisName() const override { return "Global Visibility Analysis"; }
    bool run(Communicator *comm, PointMap &map, bool simple_version) override;
    VGAVisualGlobal(double radius, bool gates_only) : m_radius(radius), m_gates_only(gates_only) {}
};
