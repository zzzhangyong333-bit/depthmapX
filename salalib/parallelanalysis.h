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

#pragma once

#include "genlib/comm.h"
#include "genlib/parallel.h"

#include <cstdint>
#include <queue>
#include <utility>
#include <vector>

namespace depthmapX {

    // Values an analysis wants to write into one attribute row, recorded by a worker thread
    // and applied afterwards on the calling thread (the attribute table is not thread safe).
    using PendingValues = std::vector<std::pair<int, float>>;

    // parallelFor wired to a Communicator: reports CURRENT_RECORD progress (offset by
    // `alreadyDone`) and throws Communicator::CancelledException if the user cancels.
    inline void runAnalysisInParallel(Communicator *comm, size_t count, unsigned int nthreads,
                                      const std::function<void(size_t, unsigned int)> &work,
                                      size_t alreadyDone = 0) {
        bool completed = parallelFor(count, nthreads, work, [comm, alreadyDone](size_t done) {
            if (comm) {
                if (comm->IsCancelled()) {
                    return false;
                }
                comm->CommPostMessage(Communicator::CURRENT_RECORD, static_cast<int>(alreadyDone + done));
            }
            return true;
        });
        if (!completed) {
            throw Communicator::CancelledException();
        }
    }

    // Min-priority queue that pops equal keys in insertion order. This is exactly the order the
    // old code got from a sorted std::vector filled with insert_sorted (upper_bound) and drained
    // from the front, but each operation is O(log n) instead of O(n).
    template <typename Key, typename Value> class StableMinQueue {
        struct Entry {
            Key key;
            uint64_t seq;
            Value value;
        };
        struct Later {
            bool operator()(const Entry &a, const Entry &b) const {
                if (a.key < b.key)
                    return false;
                if (b.key < a.key)
                    return true;
                return a.seq > b.seq;
            }
        };
        std::priority_queue<Entry, std::vector<Entry>, Later> m_heap;
        uint64_t m_seq = 0;

      public:
        void push(const Key &key, const Value &value) { m_heap.push(Entry{key, m_seq++, value}); }
        bool empty() const { return m_heap.empty(); }
        void pop(Key &key, Value &value) {
            key = m_heap.top().key;
            value = m_heap.top().value;
            m_heap.pop();
        }
    };
} // namespace depthmapX
