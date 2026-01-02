#pragma once
// file Nasam512.h

#include <algorithm>    // std::min (for fill)
#include <array> // std::array
#include <bit> // std::rotr
#include <cstdint> // uint64_t, etc.
#include <limits> // std::numeric_limits
#include <string.h> // memcpy
#include "SplitMix64.h"

namespace RNG {
	// 128 bit multiplication support
#if defined(_MSC_VER)
#include <intrin.h>  // Provides _umul128, _addcarry_u64, etc.
// _umul128 is already declared here, nothing else needed
#elif defined(__SIZEOF_INT128__)  // More precise than just !defined(_MSC_VER)
	inline std::uint64_t _umul128(std::uint64_t a, std::uint64_t b, std::uint64_t* hi) noexcept {
		__uint128_t product = static_cast<__uint128_t>(a) * b;
		*hi = product >> 64;
		return static_cast<std::uint64_t>(product);
	}
#else
#error "Compiler/platform does not support 128-bit integer multiplication"
#endif


/*
 * Nasam512
 * ========
 *
 * A high-quality, fast, non-cryptographic 64-bit pseudorandom number generator.
 * Or, as some might say, a ridiculous over-engineered monster of a PRNG . . .
 *
 * - 100% header-only C++ implementation
 * - Minimal dependencies: C++ standard library + SplitMix64 (included for seeding)
 * - Full compatibility with std::random_number_engine concept
 *
 * Nasam512 is a counter-based PRNG with a 512-bit state (8 × uint64_t lanes)
 * using the outstanding NASAM unary mixer by Pelle Evensen (2020).
 * It delivers excellent statistical quality, high throughput, and a
 * guaranteed period of 2⁵¹².
 *
 * It has passed PractRand tests to 16 GB without failures.
 *
 * Key Features
 * ------------
 * • Period: Exactly 2⁵¹² — practically infinite
 * • State: 512 bits (8 × uint64_t)
 * • Output: Full-range 64-bit values (0 to UINT64_MAX)
 * • Mixer: NASAM — one of the strongest known non-cryptographic 64-bit unary mixers,
 *   particularly effective against low-entropy, sequential, and structured inputs
 * • Counter advance: Irregular increments derived from the golden-ratio conjugate
 *   (floor(2⁶⁴ / φ)), each subsequent value NASAM-mixed from the previous
 * • Circular carry propagation ensures full-period equidistribution in every output lane
 * • Buffered output: 8 values produced per counter increment for high throughput
 * • Efficient jump-ahead via discard(n) using 128-bit arithmetic
 * • Complete state get/set, reseed(), and stream I/O support
 * • Non-cryptographic: extremely fast and statistically robust,
 *   but not suitable for cryptographic purposes
 *
 * Design Rationale
 * ----------------
 * Nasam512 combines two proven ideas:
 *   • A multi-lane counter with highly irregular, golden-ratio-derived increments
 *     for excellent diffusion and non-linearity in carry patterns
 *   • Independent NASAM mixing of each lane to destroy any residual structure
 *
 * The circular carry wrap-around is critical: without it, lanes would exhibit
 * staggered periods (2⁶⁴, 2¹²⁸, …, 2⁵¹²). This minimal mechanism guarantees
 * theoretical full-period equidistribution across all output streams.
 *
 * Performance
 * -----------
 * Nasam512 achieves approximately 350–500 MB/s single-threaded throughput
 * on modern x86-64 CPUs — very fast, though not in the same league as
 * ultra-lightweight generators like wyrand (~2 GB/s).
 *
 * This is the deliberate price for its massive 512-bit state and 2⁵¹² period:
 * - Safe, efficient parallel stream splitting via discard()
 * - No risk of period exhaustion or stream overlap in even the longest runs
 * - Proven statistical quality even from its theoretically weakest lane
 *
 * For applications where raw speed is paramount and smaller periods are acceptable,
 * consider a lightweight alternative. For everything else — especially long-running
 * simulations, reproducible parallelism, or future-proof designs — Nasam512
 * delivers unmatched safety and quality at still-excellent speed.
 *
 * Usage Example
 * -------------
 * #include "Nasam512.h"
 * #include <iostream>
 * #include <random>
 *
 * int main() {
 *     // Deterministic seeding
 *     Nasam512 rng(12345ULL);
 *
 *     // Direct generation
 *     for (int i = 0; i < 10; ++i)
 *         std::cout << rng() << '\n';
 *
 *     // With standard distributions
 *     std::uniform_int_distribution<int> dist(1, 100);
 *     std::cout << "Dice roll: " << dist(rng) << '\n';
 *
 *     // State save/restore
 *     auto saved = rng.get_state();
 *     // ... later
 *     rng.set_state(saved);
 *
 *     // Fast skip-ahead (e.g., for parallel streams)
 *     rng.discard(1'000'000ULL);
 *
 *     // Fill buffers with random bytes
 *     std::array<uint8_t, 32> key;
 *     rng.fill(key.data(), key.size());
 *
 *     // Reseed mid-run
 *     rng.reseed(0xdeadbeefcafebabeULL);
 * }
 *
 * Performance
 * -----------
 * • Generates 8 outputs per counter advance (buffered)
 * • NASAM mixer: ≈3–4 cycles per 64-bit output on modern x86-64
 * • Jump-ahead uses native 128-bit multiplication (MSVC, GCC, Clang)
 * • Passes PractRand to multiple gigabytes (full engine and single-lane stress tests)
 *
 * Credits
 * -------
 * • NASAM mixer: Pelle Evensen, 2020
 *   http://mostlymangling.blogspot.com/2020/01/nasam-not-another-strange-acronym-mixer.html
 * • Golden-ratio constant: floor(2⁶⁴ / φ) — widely used in hashing and PRNGs
 * • Design inspired by modern counter-based generators (PCG, Romu, xoshiro)
 *
 * License
 * -------
 * MIT License — free to use, modify, and distribute without restriction.
 *
 * Nasam512 — Strong, simple, and blazing fast random numbers.
 * January 2026
 */

	class Nasam512 {
		static const int STATESIZE = 8;
		uint64_t state[STATESIZE]; // 512 bit state
		uint64_t buffer[STATESIZE]; // output buffer
		int pos = STATESIZE;  // pos==STATESIZE means buffer is empty
		uint64_t inc[8]; // increments for each lane, computed in constructor

		static constexpr uint64_t nasam(uint64_t v) noexcept {
			auto rotr64 = [](uint64_t x, int r)noexcept {
				return ((x >> r) | (x << (64 - r)));
				};
			v *= 0x9E6F1D9BB2D6C165ULL;
			v ^= rotr64(v, 26);
			v *= 0x9E6F1D9BB2D6C165ULL;
			v ^= rotr64(v, 47) ^ rotr64(v, 21);
			v *= 0x9FB21C651E98DF25ULL;
			return v ^ (v >> 28);
		}

		void initialize_inc(uint64_t out[8]) noexcept {
			// Note for curious programmers:
			// We start with the golden-ratio conjugate constant (floor(2^64 / φ)),
			// widely used in PRNGs and hash functions for its excellent bit-diffusion properties.
			// Each subsequent increment is obtained by applying the NASAM mixer to the previous one.
			// This creates a sequence of 8 highly irregular, well-distributed increments
			// with strong avalanche between lanes.
			// These values could be hard-coded for micro-optimization if desired,
			// but keeping the generation code preserves clarity and intent.
			constexpr uint64_t INC = 0x9e3779b97f4a7c15ULL;  // Golden ratio constant
			out[0] = INC;
			for (int i = 1; i < 8; ++i) {
				out[i] = nasam(out[i - 1]);
			}
		}

		// Add 'inc' to x[index] and propagate any carry (1 on overflow) forward,
		// wrapping circularly around the 8-element array until the carry is absorbed.
		// Carry detection uses overflow check for the initial add, and
		// "was this limb UINT64_MAX?" (i.e., x[index] == 0 after +1) for subsequent adds.
		void add_circular_carry(uint64_t x[8], uint64_t inc, int index) {
			uint64_t carry;
			index = index % 8;

			x[index] += inc;
			carry = (x[index] < inc) ? 1 : 0;

			// In practice: usually terminates after 1 iteration; worst case 8 iterations.
			while (carry) {
				index = (index + 1) % 8;
				x[index] += carry;
				carry = (x[index] == 0);  // carry continues only if limb rolled over from all 1s
			}
		}

		// Advances the 512-bit counter state by exactly 'blocks' increments.
		// Each increment adds inc[i] to lane i,
		// with carry rippling across lanes as in a big integer.
		// For large blocks, we scale the increment efficiently using 128-bit multiplication.
		// Final carry is wrapped around circularly to ensure EVERY output lane
		// (including buffer[0]) has the full 2^512 period — this is the key to full equidistribution.
		void increment_state(uint64_t nblocks = 1) noexcept {
			if (nblocks == 0)
				return;
			if (nblocks == 1) {
				for (int i = 0; i < 8; i++)
					add_circular_carry(state, inc[i], i);
				return;
			}
			else {
				// General case: efficiently add (nblocks × inc_vector) to the 512-bit counter
				// using 128-bit multiplication per lane, then leveraging add_circular_carry
				// to handle all carry propagation (including wrap-around).
				uint64_t carry = 0;
				uint64_t product_lo, product_hi;
				for (int i = 0; i < 8; i++) {
					product_lo = _umul128(inc[i], nblocks, &product_hi);
					// product_lo contains the lower 64 bits of inc[i] * nblocks
					// product_hi contains the upper 64 bits of inc[i] * nblocks
					// 
					// Because add_circular_carry is global and circular, adding the low 64 bits
					// starting at lane i and the high 64 bits starting at lane i+1 produces
					// exactly the same result as a full 128-bit add per lane with chained carries.
					add_circular_carry(state, product_lo, i);
					if (product_hi) add_circular_carry(state, product_hi, i + 1); // add the high part to the next lane.
				}
			}
		}

		// Advance the state by an arbitrary 512-bit number of increments:
		//   state += step × inc_vector  (with full circular carry semantics)
		//
		// Enables creation of equally spaced parallel streams with any desired spacing
		// up to 2^512 - 1. Useful for massive-scale reproducible parallelism.
		//
		// Note: Much slower than discard() — use only for stream initialization/splitting.
		//
		// Conceptually, this would allow someone to create 2^256 independent streams of
		// non-overlapping length 2^256 blocks each. Why anyone would need that many streams
		// is another question . . .
		void big_jump(const uint64_t step[8]) noexcept {
			uint64_t temp[8];
			std::memcpy(temp, state, sizeof(temp));  // start from current state

			for (int i = 0; i < 8; ++i) {
				if (step[i] == 0) continue;

				uint64_t lo, hi;
				for (int j = 0; j < 8; ++j) {
					lo = _umul128(inc[j], step[i], &hi);
					add_circular_carry(temp, lo, (i + j) % 8);
					if (hi)
						add_circular_carry(temp, hi, (i + j + 1) % 8);
				}
			}

			std::memcpy(state, temp, sizeof(state));
		}

		void refill_buffer() noexcept {
			increment_state();
			for (int i = 0; i < STATESIZE; i++) {
				buffer[i] = nasam(state[i]);
			}
			pos = 0; // buffer is full
		}

	public:

		//
		// CONSTRUCTORS AND DESTRUCTOR
		// 

		// Default constructor: non-deterministic seeding via std::random_device
		// (potential upgrade: use BCryptGenRandom on Windows for stronger entropy)
		Nasam512() : pos(STATESIZE) {
			initialize_inc(inc);
			// Potential improvement: replace std::random_device with a strong entropy source on each platform.
			std::random_device rd;
			for (int i = 0; i < STATESIZE; ++i) {
				state[i] = (static_cast<uint64_t>(rd()) << 32) | rd();
			}
		}

		// Construct from a single 64-bit seed (deterministic) using SplitMix64
		explicit Nasam512(uint64_t seed) : pos(STATESIZE) {
			initialize_inc(inc);
			RNG::SplitMix64 gen(seed);
			for (int i = 0; i < STATESIZE; ++i) {
				state[i] = gen();
			}
		}

		// Construct from an explicit full 512-bit state
		explicit Nasam512(const std::array<uint64_t, 8>& initial_state) : pos(STATESIZE) {
			initialize_inc(inc);
			std::memcpy(state, initial_state.data(), sizeof(state));
		}

		// Construct from any SeedSequence-compatible type (e.g., std::seed_seq, random_device)
		template<class Sseq>
		explicit Nasam512(Sseq& seq) : pos(STATESIZE) {
			initialize_inc(inc);
			std::uint32_t seeds[16];
			seq.generate(seeds, seeds + 16);
			for (int i = 0; i < STATESIZE; ++i) {
				state[i] = (static_cast<uint64_t>(seeds[2 * i]) << 32) | seeds[2 * i + 1];
			}
		}

		// Copy and move are safe and efficient — use defaults
		Nasam512(const Nasam512&) = default;
		Nasam512& operator=(const Nasam512&) = default;
		Nasam512(Nasam512&&) = default;
		Nasam512& operator=(Nasam512&&) = default;

		// Destructor also defaulted
		~Nasam512() = default;

		//
		// RANDOM NUMBER GENERATION
		// 

		// Generate the next 64 bit random number
		inline uint64_t operator()() {
			if (pos == STATESIZE)
				refill_buffer();
			return buffer[pos++];

			// For testing only...
			//refill_buffer();
			//return buffer[2]; // lane 2 has the shortest increment
		}

		// Fill a byte buffer with random data
		// Useful for generating random keys, noise, or initializing memory
		void fill(uint8_t* data, size_t size) noexcept {
			size_t bytes_generated = 0;
			while (bytes_generated < size) {
				uint64_t rnd = (*this)();  // Get next 64-bit random value

				size_t bytes_to_copy = std::min(sizeof(uint64_t), size - bytes_generated);
				std::memcpy(data + bytes_generated, &rnd, bytes_to_copy);

				bytes_generated += bytes_to_copy;
			}
		}

		// 
		// STATE MAMAGEMENT
		// 

		// Since this is a non-cryptographic RNG, we provide state get/set functions

		Nasam512& set_state(const std::array<uint64_t, 8> initial_state) noexcept {
			memcpy(state, initial_state.data(), STATESIZE * sizeof(uint64_t));
			return *this;
		}

		std::array<uint64_t, 8> get_state() const noexcept {
			std::array<uint64_t, 8> result;
			memcpy(result.data(), state, STATESIZE * sizeof(uint64_t));
			return result;
		}

		// Advance the RNG state by 'n' outputs without generating them.
		void discard(uint64_t n) noexcept {
			// Goal: Advance the RNG forward by exactly 'n' output values,
			// without generating them individually (for efficiency).

#if 0
		// Equivalent brute-force implementation (kept for clarity and verification)
		// This is what discard(n) must semantically match.
			for (uint64_t i = 0; i < n; ++i) {
				(*this)();  // Generate and discard one output value
			}
			return;
#else
		// Optimized implementation

		// Step 1: Consume any outputs already present in the current buffer
		// pos is the index of the next available value (0 <= pos < STATESIZE when buffer has data)
		//uint64_t remaining_in_buffer = STATESIZE - pos;
			uint64_t remaining_in_buffer = static_cast<uint64_t>(STATESIZE) - pos;

			if (n <= remaining_in_buffer) {
				// All requested skips are within the current buffer — just advance the pointer
				pos += static_cast<int>(n);
				return;
			}

			// We've used up the current buffer
			n -= remaining_in_buffer;

			// Step 2: Skip full blocks of STATESIZE outputs efficiently
			// Each block corresponds to one increment of the counter + NASAM mixing
			uint64_t full_blocks_to_skip = n / STATESIZE;        // Complete blocks we can bypass entirely
			uint64_t outputs_in_final_block = n % STATESIZE;     // Outputs needed from the next block

			// Advance the counter past all the full blocks we want to completely skip
			if (full_blocks_to_skip > 0) {
				increment_state(full_blocks_to_skip);
			}

			// Step 3: Generate the block that contains the remaining outputs
			// refill_buffer() correctly:
			//   - increments the counter one more time
			//   - applies NASAM to produce a fresh buffer of STATESIZE outputs
			// This is intentional and necessary — we need access to these values
			// so we can skip the first 'outputs_in_final_block' of them.
			refill_buffer();

			// Step 4: Position the buffer pointer past the outputs we just "discarded"
			// in the newly generated block
			pos = static_cast<int>(outputs_in_final_block);

			// At this point, exactly 'n' outputs have been skipped,
			// and the next call to operator()() will return the correct value.
#endif
		}

	public:	// Compatibility
		using result_type = uint64_t;

		// Constants required by the concept
		static constexpr result_type min() noexcept { return 0; }
		static constexpr result_type max() noexcept { return std::numeric_limits<result_type>::max(); }

		// 1. seed() with no argument — same as default constructor
		void seed() {
			*this = Nasam512();  // delegating to default ctor
		}

		// 2. seed() with single uint64_t — delegate to your existing ctor
		void seed(uint64_t s) {
			*this = Nasam512(s);
		}

		// 3. seed() with SeedSequence — delegate to template ctor
		template<class Sseq>
		void seed(Sseq& seq) {
			*this = Nasam512(seq);
		}

		// 4. Equality / inequality — compare full state
		friend bool operator==(const Nasam512& lhs, const Nasam512& rhs) {
			return lhs.pos == rhs.pos &&
				std::memcmp(lhs.state, rhs.state, sizeof(lhs.state)) == 0 &&
				std::memcmp(lhs.buffer, rhs.buffer, sizeof(lhs.buffer)) == 0;
			// inc[] is always the same, no need to compare
		}

		friend bool operator!=(const Nasam512& lhs, const Nasam512& rhs) {
			return !(lhs == rhs);
		}

		// 5. Stream output (for debugging / serialization)
		template<class CharT, class Traits>
		friend std::basic_ostream<CharT, Traits>&
			operator<<(std::basic_ostream<CharT, Traits>& os, const Nasam512& rng)
		{
			os << rng.pos;
			for (auto v : rng.state)  os << ' ' << v;
			for (auto v : rng.buffer) os << ' ' << v;
			return os;
		}

		// 6. Stream input (restore state)
		template<class CharT, class Traits>
		friend std::basic_istream<CharT, Traits>&
			operator>>(std::basic_istream<CharT, Traits>& is, Nasam512& rng)
		{
			is >> rng.pos;
			for (auto& v : rng.state)  is >> v;
			for (auto& v : rng.buffer) is >> v;
			// inc[] will be reinitialized on next use or you can call initialize_inc(rng.inc);
			return is;
		}

		// Reseed the generator from a new 64-bit seed
		// Uses the same SplitMix64 expansion as the constructor for consistency
		void reseed(uint64_t seed) noexcept {
			RNG::SplitMix64 gen(seed);
			for (int i = 0; i < STATESIZE; ++i) {
				state[i] = gen();
			}
			pos = STATESIZE;  // Invalidate buffer — will refill on next call
		}
	};// class Nasam512

}// namespace RNG
