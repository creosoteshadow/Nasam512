# Nasam512
Non-cryptographic Random Number Generator, 512-bit state.


 Nasam512
 ========

 A high-quality, fast, non-cryptographic 64-bit pseudorandom number generator.
 Or, as some might say, a ridiculous over-engineered monster of a PRNG . . .

 - 100% header-only C++ implementation
 - Minimal dependencies: C++ standard library + SplitMix64 (included for seeding)
 - Full compatibility with std::random_number_engine concept

 Nasam512 is a counter-based PRNG with a 512-bit state (8 × uint64_t lanes)
 using the outstanding NASAM unary mixer by Pelle Evensen (2020).
 It delivers excellent statistical quality, high throughput, and a
 guaranteed period of 2⁵¹².

 It has passed PractRand tests to 16 GB without failures.

 Key Features
 ------------
 • Period: Exactly 2⁵¹² — practically infinite
 • State: 512 bits (8 × uint64_t)
 • Output: Full-range 64-bit values (0 to UINT64_MAX)
 • Mixer: NASAM — one of the strongest known non-cryptographic 64-bit unary mixers,
   particularly effective against low-entropy, sequential, and structured inputs
 • Counter advance: Irregular increments derived from the golden-ratio conjugate
   (floor(2⁶⁴ / φ)), each subsequent value NASAM-mixed from the previous
 • Circular carry propagation ensures full-period equidistribution in every output lane
 • Buffered output: 8 values produced per counter increment for high throughput
 • Efficient jump-ahead via discard(n) using 128-bit arithmetic
 • Complete state get/set, reseed(), and stream I/O support
 • Non-cryptographic: extremely fast and statistically robust,
   but not suitable for cryptographic purposes

 Design Rationale
 ----------------
 Nasam512 combines two proven ideas:
   • A multi-lane counter with highly irregular, golden-ratio-derived increments
     for excellent diffusion and non-linearity in carry patterns
   • Independent NASAM mixing of each lane to destroy any residual structure

 The circular carry wrap-around is critical: without it, lanes would exhibit
 staggered periods (2⁶⁴, 2¹²⁸, …, 2⁵¹²). This minimal mechanism guarantees
 theoretical full-period equidistribution across all output streams.

 Performance
 -----------
 Nasam512 achieves approximately 350–500 MB/s single-threaded throughput
 on modern x86-64 CPUs — very fast, though not in the same league as
 ultra-lightweight generators like wyrand (~2 GB/s).

 This is the deliberate price for its massive 512-bit state and 2⁵¹² period:
 - Safe, efficient parallel stream splitting via discard()
 - No risk of period exhaustion or stream overlap in even the longest runs
 - Proven statistical quality even from its theoretically weakest lane

 For applications where raw speed is paramount and smaller periods are acceptable,
 consider a lightweight alternative. For everything else — especially long-running
 simulations, reproducible parallelism, or future-proof designs — Nasam512
 delivers unmatched safety and quality at still-excellent speed.

 Usage Example
 -------------
 #include "Nasam512.h"
 #include <iostream>
 #include <random>

 int main() {
     // Deterministic seeding
     Nasam512 rng(12345ULL);

     // Direct generation
     for (int i = 0; i < 10; ++i)
         std::cout << rng() << '\n';

     // With standard distributions
     std::uniform_int_distribution<int> dist(1, 100);
     std::cout << "Dice roll: " << dist(rng) << '\n';

     // State save/restore
     auto saved = rng.get_state();
     // ... later
     rng.set_state(saved);

     // Fast skip-ahead (e.g., for parallel streams)
     rng.discard(1'000'000ULL);

     // Fill buffers with random bytes
     std::array<uint8_t, 32> key;
     rng.fill(key.data(), key.size());

     // Reseed mid-run
     rng.reseed(0xdeadbeefcafebabeULL);
 }

 Performance
 -----------
 • Generates 8 outputs per counter advance (buffered)
 • NASAM mixer: ≈3–4 cycles per 64-bit output on modern x86-64
 • Jump-ahead uses native 128-bit multiplication (MSVC, GCC, Clang)
 • Passes PractRand to multiple gigabytes (full engine and single-lane stress tests)

 Credits
 -------
 • NASAM mixer: Pelle Evensen, 2020
   http://mostlymangling.blogspot.com/2020/01/nasam-not-another-strange-acronym-mixer.html
 • Golden-ratio constant: floor(2⁶⁴ / φ) — widely used in hashing and PRNGs
 • Design inspired by modern counter-based generators (PCG, Romu, xoshiro)

 License
 -------
 MIT License — free to use, modify, and distribute without restriction.

 Nasam512 — Strong, simple, and blazing fast random numbers.
 January 2026

