// Copyright (C) 2023 by Mark Melton
//

#pragma once
#include <atomic>
#include <chrono>

namespace utils {

/// The Timer class template implements a timer designed for minimal
/// overhead, ad-hoc timing of code regions including micro-timing
/// down to single machine instructions.
template<class Duration = std::chrono::nanoseconds,
         class Clock = std::chrono::high_resolution_clock>
class Timer {
public:
    using TimePoint = typename Clock::time_point;

    /// Run the supplied `code` in a loop `n` times.
    template<class Code>
    Timer& run(size_t n, Code&& code) {
        start();
        for (auto i = 0ul; i < n; ++i) {
            code();
        }
        stop(n);
        return *this;
    }

    /// Start the timer.
    void start() {
        start_ = Clock::now();
    }

    /// Stop the timer indicating `n` operations.
    auto stop(size_t n = 1) {
        auto end = Clock::now();
        iterations_ += n;
        elapsed_ += std::chrono::duration_cast<Duration>(end - start_);
        return elapsed_;
    }

    /// Return the average number of nanoseconds per operation.
    auto elapsed_per_iteration() const {
        return iterations_ > 0 ? (double)elapsed_.count() / iterations_ : 0.0;
    }

    /// Return the elapsed duration.
    auto elapsed() const {
        return elapsed_;
    }

    auto cur_time() const {
        return Clock::now() - start_;
    }

private:
    TimePoint start_{};
    Duration elapsed_{};
    size_t iterations_{};
};

}; // core::timer
