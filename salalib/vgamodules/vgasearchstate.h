// sala - a component of the depthmapX - spatial network analysis platform
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

// Per-thread replacement for the search bookkeeping that used to live inside every Point
// (m_misc / m_dist / m_cumangle). Keeping it outside the map lets several searches run at
// once. The extract functions below are line-for-line copies of Bin::extractMetric and
// Bin::extractAngular operating on this state, so results are bit-identical to the old code.

#pragma once

#include "salalib/ngraph.h"
#include "salalib/pixelref.h"
#include "salalib/pointdata.h"

#include <cstdint>
#include <queue>
#include <set>
#include <vector>

namespace depthmapX {

    // Drop-in replacement for the std::set<Triple> used as a priority queue by the metric and
    // angular searches. A binary heap is much cheaper than a balanced tree here. To reproduce the
    // set exactly, entries that compare equal are popped in insertion order (the set kept the
    // first one inserted); the later duplicates are then skipped by the search's "already
    // visited" check, exactly as if they had never been inserted.
    template <typename Triple> class SearchQueue {
        struct Entry {
            Triple triple;
            uint64_t seq;
        };
        struct Later {
            bool operator()(const Entry &a, const Entry &b) const {
                if (a.triple < b.triple)
                    return false;
                if (b.triple < a.triple)
                    return true;
                return a.seq > b.seq;
            }
        };
        std::priority_queue<Entry, std::vector<Entry>, Later> m_heap;
        uint64_t m_seq = 0;
        bool m_hasLast = false;
        Triple m_last;

      public:
        void insert(const Triple &t) { m_heap.push(Entry{t, m_seq++}); }
        bool empty() const { return m_heap.empty(); }
        // Pops the smallest entry, skipping entries equal to one already popped (set semantics).
        // Returns false when the queue is exhausted.
        bool pop(Triple &out) {
            while (!m_heap.empty()) {
                Triple t = m_heap.top().triple;
                m_heap.pop();
                if (m_hasLast && !(m_last < t) && !(t < m_last)) {
                    continue; // duplicate of the entry just popped: the set would never have held it
                }
                m_last = t;
                m_hasLast = true;
                out = t;
                return true;
            }
            return false;
        }
    };

    struct VGASearchState {
        size_t cols;
        int defaultMisc;
        float defaultDist;
        float defaultCumangle;
        std::vector<int> misc;
        std::vector<float> dist;
        std::vector<float> cumangle;
        std::vector<size_t> touched;

        VGASearchState(size_t rows, size_t columns, int miscInit, float distInit, float cumangleInit)
            : cols(columns), defaultMisc(miscInit), defaultDist(distInit), defaultCumangle(cumangleInit),
              misc(rows * columns, miscInit), dist(rows * columns, distInit),
              cumangle(rows * columns, cumangleInit) {}

        size_t index(PixelRef p) const { return size_t(p.y) * cols + size_t(p.x); }

        // mark a cell as modified so reset() can restore it cheaply
        size_t touch(PixelRef p) {
            size_t i = index(p);
            touched.push_back(i);
            return i;
        }

        void reset() {
            for (size_t i : touched) {
                misc[i] = defaultMisc;
                dist[i] = defaultDist;
                cumangle[i] = defaultCumangle;
            }
            touched.clear();
        }

        // Copy this state into the point map's legacy per-point fields
        void writeTo(PointMap &map) const {
            for (size_t i = 0; i < map.getCols(); i++) {
                for (size_t j = 0; j < map.getRows(); j++) {
                    PixelRef p(static_cast<short>(i), static_cast<short>(j));
                    Point &pt = map.getPoint(p);
                    size_t idx = index(p);
                    pt.m_misc = misc[idx];
                    pt.m_dist = dist[idx];
                    pt.m_cumangle = cumangle[idx];
                }
            }
        }
    };

    inline void extractMetric(const Node &node, SearchQueue<MetricTriple> &pixels, const PointMap &map,
                              VGASearchState &state, const MetricTriple &curs) {
        if (curs.dist == 0.0f || map.getPoint(curs.pixel).blocked() || map.blockedAdjacent(curs.pixel)) {
            for (int b = 0; b < 32; b++) {
                const Bin &bin = node.bin(b);
                for (auto pixVec : bin.m_pixel_vecs) {
                    for (PixelRef pix = pixVec.start(); pix.col(bin.m_dir) <= pixVec.end().col(bin.m_dir);) {
                        size_t idx = state.index(pix);
                        if (state.misc[idx] == 0 &&
                            (state.dist[idx] == -1.0 || (curs.dist + dist(pix, curs.pixel) < state.dist[idx]))) {
                            state.touched.push_back(idx);
                            state.dist[idx] = curs.dist + (float)dist(pix, curs.pixel);
                            // n.b. dmap v4.06r now sets angle in range 0 to 4 (1 = 90 degrees)
                            state.cumangle[idx] =
                                state.cumangle[state.index(curs.pixel)] +
                                (curs.lastpixel == NoPixel
                                     ? 0.0f
                                     : (float)(angle(pix, curs.pixel, curs.lastpixel) / (M_PI * 0.5)));
                            pixels.insert(MetricTriple(state.dist[idx], pix, curs.pixel));
                        }
                        pix.move(bin.m_dir);
                    }
                }
            }
        }
    }

    inline void extractAngular(const Node &node, SearchQueue<AngularTriple> &pixels, const PointMap &map,
                               VGASearchState &state, const AngularTriple &curs) {
        if (curs.angle == 0.0f || map.getPoint(curs.pixel).blocked() || map.blockedAdjacent(curs.pixel)) {
            for (int b = 0; b < 32; b++) {
                const Bin &bin = node.bin(b);
                for (auto pixVec : bin.m_pixel_vecs) {
                    for (PixelRef pix = pixVec.start(); pix.col(bin.m_dir) <= pixVec.end().col(bin.m_dir);) {
                        size_t idx = state.index(pix);
                        if (state.misc[idx] == 0) {
                            // n.b. dmap v4.06r now sets angle in range 0 to 4 (1 = 90 degrees)
                            float ang = (curs.lastpixel == NoPixel)
                                            ? 0.0f
                                            : (float)(angle(pix, curs.pixel, curs.lastpixel) / (M_PI * 0.5));
                            if (state.cumangle[idx] == -1.0 || curs.angle + ang < state.cumangle[idx]) {
                                state.touched.push_back(idx);
                                state.cumangle[idx] = state.cumangle[state.index(curs.pixel)] + ang;
                                pixels.insert(AngularTriple(state.cumangle[idx], pix, curs.pixel));
                            }
                        }
                        pix.move(bin.m_dir);
                    }
                }
            }
        }
    }

    // Visits every pixel a node can see, in the same order (and with the same repeats) as the
    // old Node::first()/next() iterator, but without touching the node's shared cursor.
    template <typename F> void forEachVisible(const Node &node, F f) {
        for (int b = 0; b < 32; b++) {
            const Bin &bin = node.bin(b);
            for (const auto &pixVec : bin.m_pixel_vecs) {
                PixelRef pix = pixVec.start();
                do {
                    f(pix);
                    pix.move(bin.m_dir);
                } while (pix.col(bin.m_dir) <= pixVec.end().col(bin.m_dir));
            }
        }
    }
} // namespace depthmapX
