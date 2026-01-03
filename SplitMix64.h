#pragma once
// File: SplitMix64.h

#include <cstdint>
#include <random>
#include <limits>

namespace RNG {

    // tags
    class Deterministic {};
    class NonDeterministic {};

    // -----------------------------------------------------------------
    // Very small, fast SplitMix64 (fully constexpr-friendly)
    // -----------------------------------------------------------------

    class SplitMix64 {
        using u64 = uint64_t;
        static constexpr u64 INCREMENT = 0x9e3779b97f4a7c15ULL;
        u64 state;
    public:
        using result_type = u64;

        constexpr SplitMix64(u64 seed = 0) noexcept : state(seed) {} // default = deterministic
        constexpr SplitMix64(Deterministic, u64 seed) noexcept : state(seed) {}

        // Cryptographically secure non-deterministic seeding, cannot be noexcept
        SplitMix64(NonDeterministic) noexcept(false) {
            std::random_device rd;
            state = (static_cast<uint64_t>(rd()) >> 32) | (static_cast<uint64_t>(rd()));
        }

        constexpr u64 operator()() noexcept {
            u64 z = (state += INCREMENT);
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
            z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
            return z ^ (z >> 31);
        }

        constexpr SplitMix64& discard(u64 n) noexcept {
            state += INCREMENT * n;
            return *this;
        }

        static constexpr result_type min()  noexcept { return 0; }
        static constexpr result_type max()  noexcept { return UINT64_MAX; }
    };// class SplitMix64
}// namespace RNG
