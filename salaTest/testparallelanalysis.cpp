// Copyright (C) 2026, depthmapX contributors

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

#include "catch.hpp"

#include "salalib/displayparams.h"
#include "salalib/pafcolor.h"
#include "salalib/parallelanalysis.h"

#include "genlib/containerutils.h"

#include <utility>
#include <vector>

TEST_CASE("StableMinQueue pops in the same order as a sorted vector drained from the front", "") {
    // the order the old code got from insert_sorted (upper_bound) + erase(begin())
    std::vector<std::pair<float, int>> inserts = {{3.0f, 0}, {1.0f, 1}, {3.0f, 2}, {2.0f, 3},
                                                  {1.0f, 4}, {3.0f, 5}, {0.5f, 6}, {2.0f, 7}};
    struct KeyOnly {
        float key;
        int id;
        bool operator<(const KeyOnly &other) const { return key < other.key; }
    };
    std::vector<KeyOnly> sorted;
    depthmapX::StableMinQueue<float, int> queue;
    for (auto &item : inserts) {
        depthmapX::insert_sorted(sorted, KeyOnly{item.first, item.second});
        queue.push(item.first, item.second);
    }
    for (auto &expected : sorted) {
        REQUIRE_FALSE(queue.empty());
        float key;
        int id;
        queue.pop(key, id);
        REQUIRE(key == expected.key);
        REQUIRE(id == expected.id);
    }
    REQUIRE(queue.empty());
}

TEST_CASE("Viridis and magma colour scales", "") {
    DisplayParams dp;
    dp.colorscale = DisplayParams::VIRIDIS;
    PafColor low, high, mid;
    low.makeColor(0.0, dp);
    high.makeColor(1.0, dp);
    REQUIRE((low.m_color & 0xFFFFFF) == 0x440154);
    REQUIRE((high.m_color & 0xFFFFFF) == 0xFDE725);
    REQUIRE((low.m_color >> 24) == 0xFF);

    // midpoint is exactly the 6th stop
    mid.makeColor(0.5, dp);
    REQUIRE((mid.m_color & 0xFFFFFF) == 0x21918C);

    dp.colorscale = DisplayParams::MAGMA;
    low.makeColor(0.0, dp);
    high.makeColor(1.0, dp);
    REQUIRE((low.m_color & 0xFFFFFF) == 0x000004);
    REQUIRE((high.m_color & 0xFFFFFF) == 0xFCFDBF);

    // values outside the range are clamped, -1 stays the "no value" grey
    PafColor clamped, none;
    clamped.makeColor(5.0, dp);
    REQUIRE((clamped.m_color & 0xFFFFFF) == 0xFCFDBF);
    none.makeColor(-1.0, dp);
    REQUIRE(none.m_color == 0x007f7f7f);
}
