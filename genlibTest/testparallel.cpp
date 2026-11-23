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
#include <genlib/parallel.h>

#include <atomic>
#include <chrono>
#include <stdexcept>
#include <thread>
#include <vector>

TEST_CASE("parallelFor visits every index exactly once", "") {
    for (unsigned int threads : {1u, 2u, 7u}) {
        const size_t count = 1000;
        std::vector<std::atomic<int>> visits(count);
        for (auto &v : visits) {
            v = 0;
        }
        std::atomic<unsigned int> maxThread(0);
        bool completed = depthmapX::parallelFor(
            count, threads,
            [&](size_t index, unsigned int thread) {
                visits[index]++;
                unsigned int seen = maxThread.load();
                while (thread > seen && !maxThread.compare_exchange_weak(seen, thread)) {
                }
            },
            nullptr);
        REQUIRE(completed);
        for (auto &v : visits) {
            REQUIRE(v.load() == 1);
        }
        REQUIRE(maxThread.load() < threads);
    }
}

TEST_CASE("parallelFor with no work", "") {
    bool called = false;
    REQUIRE(depthmapX::parallelFor(
        0, 4, [&](size_t, unsigned int) { called = true; }, nullptr));
    REQUIRE_FALSE(called);
}

TEST_CASE("parallelFor stops when progress asks to cancel", "") {
    std::atomic<size_t> done(0);
    bool completed = depthmapX::parallelFor(
        100000, 2,
        [&](size_t, unsigned int) {
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
            done++;
        },
        [](size_t) { return false; });
    REQUIRE_FALSE(completed);
    REQUIRE(done.load() < 100000);
}

TEST_CASE("parallelFor rethrows worker exceptions on the calling thread", "") {
    REQUIRE_THROWS_AS(depthmapX::parallelFor(
                          50, 3,
                          [](size_t index, unsigned int) {
                              if (index == 17) {
                                  throw std::runtime_error("boom");
                              }
                          },
                          nullptr),
                      std::runtime_error);
}

TEST_CASE("threadsFor never exceeds the amount of work", "") {
    REQUIRE(depthmapX::threadsFor(0) == 1);
    REQUIRE(depthmapX::threadsFor(1) == 1);
    REQUIRE(depthmapX::threadsFor(1000000) >= 1);
}
