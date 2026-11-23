// genlib - a component of the depthmapX - spatial network analysis platform
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

// Small, dependency-free parallel-for used by the analysis modules.
//
// Design notes:
//  * Work items are handed out one at a time through an atomic counter, which balances well
//    for the analyses here (each item is a full graph search and costs vary a lot by position).
//  * The calling thread never runs work itself; it only reports progress and checks for
//    cancellation, so the Communicator (progress dialog) is only ever touched from the thread
//    that called run(), exactly as before.
//  * The number of threads can be forced with the DEPTHMAPX_THREADS environment variable
//    (DEPTHMAPX_THREADS=1 gives the old single-threaded behaviour).

#pragma once

#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstdlib>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

namespace depthmapX {

    inline unsigned int analysisThreadCount() {
        if (const char *env = std::getenv("DEPTHMAPX_THREADS")) {
            int forced = std::atoi(env);
            if (forced > 0) {
                return static_cast<unsigned int>(forced);
            }
        }
        unsigned int hw = std::thread::hardware_concurrency();
        return hw == 0 ? 1 : hw;
    }

    // Number of worker threads to use for `count` work items.
    inline unsigned int threadsFor(size_t count) {
        return static_cast<unsigned int>(std::max<size_t>(1, std::min<size_t>(analysisThreadCount(), count)));
    }

    // Runs work(index, threadIndex) for every index in [0, count) on `nthreads` threads.
    // Callers keep per-thread scratch state in a vector of size nthreads, indexed by threadIndex.
    // progress(completed) is called on the calling thread roughly every 100ms; if it returns
    // false the remaining work is abandoned and parallelFor returns false.
    // Any exception thrown by a work item is re-thrown on the calling thread.
    inline bool parallelFor(size_t count, unsigned int nthreads, const std::function<void(size_t, unsigned int)> &work,
                            const std::function<bool(size_t)> &progress) {
        if (count == 0) {
            return true;
        }
        nthreads = std::max(1u, nthreads);

        std::atomic<size_t> next(0);
        std::atomic<size_t> done(0);
        std::atomic<bool> stop(false);
        std::exception_ptr error;
        std::mutex errorMutex;
        std::mutex doneMutex;
        std::condition_variable doneCv;
        unsigned int finishedThreads = 0;

        auto worker = [&](unsigned int threadIndex) {
            try {
                while (!stop.load(std::memory_order_relaxed)) {
                    size_t i = next.fetch_add(1);
                    if (i >= count) {
                        break;
                    }
                    work(i, threadIndex);
                    done.fetch_add(1, std::memory_order_relaxed);
                }
            } catch (...) {
                std::lock_guard<std::mutex> lock(errorMutex);
                if (!error) {
                    error = std::current_exception();
                }
                stop = true;
            }
            {
                std::lock_guard<std::mutex> lock(doneMutex);
                finishedThreads++;
            }
            doneCv.notify_one();
        };

        std::vector<std::thread> threads;
        threads.reserve(nthreads);
        for (unsigned int t = 0; t < nthreads; t++) {
            threads.emplace_back(worker, t);
        }

        bool cancelled = false;
        {
            std::unique_lock<std::mutex> lock(doneMutex);
            while (finishedThreads < nthreads) {
                doneCv.wait_for(lock, std::chrono::milliseconds(100));
                if (finishedThreads < nthreads && !cancelled) {
                    lock.unlock();
                    bool keepGoing = progress ? progress(done.load()) : true;
                    lock.lock();
                    if (!keepGoing) {
                        cancelled = true;
                        stop = true;
                    }
                }
            }
        }
        for (auto &t : threads) {
            t.join();
        }
        if (error) {
            std::rethrow_exception(error);
        }
        return !cancelled;
    }
} // namespace depthmapX
