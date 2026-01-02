#pragma once
#define _CRT_SECURE_NO_WARNINGS
#include <cstdint>
#include <iostream>
#include <string>

using u8 = std::uint8_t;
using u16 = std::uint16_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;

template <class GENERATOR>
inline void write_PractRand_file(int nGB, GENERATOR& gen) {
    //
    // Warning!! This function assumes RNG_test.exe is in the current directory.
    //

    if (nGB == 0 || nGB > 1024) {
        throw std::runtime_error("nGB must be between 1 and 1024");
    }

    constexpr size_t buffer_size = 64;
    u64 buffer[buffer_size];

    FILE* outfile = fopen("test.bin", "wb");
    if (!outfile) {
        printf("Error opening test.bin\n");
        system("pause");
        exit(EXIT_FAILURE);
    }

    u64 nwords = nGB * 1024ull * 1024ull * 1024ull / 8;
    std::cout << "Size: " << nGB << " GB\n";
    std::cout << "Words: " << nwords << "\n";

    u64 offset = 0;
    size_t remaining = nwords;
    while (remaining) {
        size_t n_used = std::min(remaining, buffer_size);
        for (size_t i = 0; i < n_used; i++)
            //buffer[i] = gen.draw64();
            buffer[i] = gen();

        size_t nb = fwrite(buffer, sizeof(u64), n_used, outfile);
        if (nb < n_used) {
            printf("Error writing to test.bin at offset=%llu\n", offset);
            system("pause");
            fclose(outfile);
            exit(EXIT_FAILURE);
        }
        remaining -= n_used;
        offset += n_used;
        if (offset % (1024ull * 1024ull * 1024ull / 8) == 0) {
            std::cout << "Wrote " << (offset * 8) / (1024ull * 1024ull * 1024ull) << " GB\n";
        }
    }

    fflush(outfile);
    fclose(outfile);

    // Construct and execute PractRand command (for Windows; for Unix use 'cat' instead of 'type')
    std::string command = "type test.bin | RNG_test.exe stdin64 -tf 2 -te 1 -tlmax " + std::to_string(nGB) + "GB -multithreaded";
    std::cout << "Executing PractRand command: " << command << "\n";
    int result = system(command.c_str());
    if (result != 0) {
        std::cerr << "PractRand command failed with return code " << result << "\n";
        system("pause");
        // Optionally exit or handle the error differently
    }
    else {
        std::cout << "PractRand command completed successfully\n";
    }

    //printf("Suggested PractRand command line (for reference):\n");
    //printf("\tWindows: %s\n", command.c_str());
    //printf("\tUnix: cat test.bin | ./RNG_test stdin64 -tf 2 -te 1 -tlmax %dGB -multithreaded\n", nGB);
}


#include "Nasam512.h"

int main() {
	RNG::Nasam512 gen(123456789ull);
	write_PractRand_file(16, gen);
	return 0;
}
