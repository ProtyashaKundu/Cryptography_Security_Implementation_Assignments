#include <stdio.h>
#include <stdint.h>
#include <immintrin.h>
#include <time.h>
#include <x86intrin.h>

//AVX: Instead of operating on one number at a time (scalar execution), AVX allows the CPU to load wide vector registers and perform identical mathematical operations on multiple data elements simultaneously.

//__m256i: array of unsinged integers (e.g., 32 X 8-bit, 16 X 16-bit, 8 X 32-bit, or 4 X 64-bit)

//_m256_* Intrinsic Functions used in this program
//_mm256_setr_epi32: load initial key, constants, counter, and nonce
//_mm256_add_epi32: add step in ARX & final state accumulation
//_mm256_xor_si256: XOR step in ARX
//_mm256_slli_epi32/srli_epi32: vector building blocks for Rotate in ARX
//_mm256_or_si256: combine shifted vectors for Rotate in ARX
//_mm256_shuffle_epi32: re-align matrix lanes between Column and Diagonal rounds
//_mm256_storeu_si256:	write computed keystream out to RAM


//macro to left-rotate 32-bit lane values within 256-bit AVX2 registers
#define ROTL32(v, n) \
    _mm256_or_si256(_mm256_slli_epi32((v), (n)), \
    _mm256_srli_epi32((v), 32 - (n)))

//macro for Intel AVX2-optimized chacha20 quarter round 
#define QUARTER_ROUND(a, b, c, d) \
    a = _mm256_add_epi32(a, b); \
    d = _mm256_xor_si256(d, a); \
    d = ROTL32(d, 16); \
    \
    c = _mm256_add_epi32(c, d); \
    b = _mm256_xor_si256(b, c); \
    b = ROTL32(b, 12); \
    \
    a = _mm256_add_epi32(a, b); \
    d = _mm256_xor_si256(d, a); \
    d = ROTL32(d, 8);  \
    \
    c = _mm256_add_epi32(c, d); \
    b = _mm256_xor_si256(b, c); \
    b = ROTL32(b, 7);

// processes two 64-byte chacha20 blocks (128 bytes total) simultaneously using 256-bit registers
void chacha20_avx2_double_block(const uint32_t key[8], const uint32_t nonce[3], uint32_t counter, uint8_t output[128]) 
{
    // row 0: constants ("expand 32-byte k")
    __m256i row0 = _mm256_setr_epi32(0x61707865, 0x3320646e, 0x7920622d, 0x6b206574, 0x61707865, 0x3320646e, 0x7920622d, 0x6b206574);

    // row 1 & 2: key duplicated across both 128-bit halves
    __m256i row1 = _mm256_setr_epi32(key[0], key[1], key[2], key[3], key[0], key[1], key[2], key[3]);
    __m256i row2 = _mm256_setr_epi32(key[4], key[5], key[6], key[7], key[4], key[5], key[6], key[7]);

    // row 3: counter and nonce (block 0 gets counter, block 1 gets counter + 1)
    __m256i row3 = _mm256_setr_epi32(counter, nonce[0], nonce[1], nonce[2],counter + 1, nonce[0], nonce[1], nonce[2]);

    // preserve initial state for final addition
    __m256i orig0 = row0, orig1 = row1, orig2 = row2, orig3 = row3;

    // 10 iterations = 20 rounds total because we have two rounds in parallel per iteration
    for (int i = 0; i < 10; i++) 
    {
	  // column round
	  QUARTER_ROUND(row0, row1, row2, row3);

	  // diagonal shuffle: shift row1, row2, row3 across 32-bit lanes to align diagonals
	  row1 = _mm256_shuffle_epi32(row1, _MM_SHUFFLE(0, 3, 2, 1));
	  row2 = _mm256_shuffle_epi32(row2, _MM_SHUFFLE(1, 0, 3, 2));
	  row3 = _mm256_shuffle_epi32(row3, _MM_SHUFFLE(2, 1, 0, 3));

	  // diagonal round
	  QUARTER_ROUND(row0, row1, row2, row3);

	  // undo diagonal shuffle
	  row1 = _mm256_shuffle_epi32(row1, _MM_SHUFFLE(2, 1, 0, 3));
	  row2 = _mm256_shuffle_epi32(row2, _MM_SHUFFLE(1, 0, 3, 2));
	  row3 = _mm256_shuffle_epi32(row3, _MM_SHUFFLE(0, 3, 2, 1));
    }

    // add initial state modulo 2^32
    row0 = _mm256_add_epi32(row0, orig0);
    row1 = _mm256_add_epi32(row1, orig1);
    row2 = _mm256_add_epi32(row2, orig2);
    row3 = _mm256_add_epi32(row3, orig3);

    // store generated keystream (128 bytes) back into output memory
    _mm256_storeu_si256((__m256i*)(output + 0), row0);
    _mm256_storeu_si256((__m256i*)(output + 32), row1);
    _mm256_storeu_si256((__m256i*)(output + 64), row2);
    _mm256_storeu_si256((__m256i*)(output + 96), row3);
}

int main() {
    uint32_t key[8] = {0x03020100, 0x07060504, 0x0b0a0908, 0x0f0e0d0c, 0x13121110, 0x17161514, 0x1b1a1918, 0x1f1e1d1c};
    uint32_t nonce[3] = {0x00000000, 0x4a000000, 0x00000000};
    uint32_t counter = 1;
    uint8_t keystream[128];

    chacha20_avx2_double_block(key, nonce, counter, keystream);//taken from claude

    const uint64_t ITERATIONS = 1000000; // 1 Million iterations (128 MB) otherwise for smaller iterations time is 0.0000 //taken from Claude

    struct timespec start_time, end_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);

    unsigned int ui;
    uint64_t start_cycles = __rdtscp(&ui);

    for (uint64_t i = 0; i < ITERATIONS; i++) 
    {
        chacha20_avx2_double_block(key, nonce, counter, keystream);
    }

    asm volatile("" : : "g"(keystream) : "memory");

    uint64_t end_cycles = __rdtscp(&ui);
    clock_gettime(CLOCK_MONOTONIC, &end_time);

    uint64_t total_time_ns = (uint64_t)(end_time.tv_sec - start_time.tv_sec) * 1000000000ULL + (uint64_t)(end_time.tv_nsec - start_time.tv_nsec);

    uint64_t total_bytes = ITERATIONS * 128;
    double ns_per_iteration = (double)total_time_ns / ITERATIONS;
    double ns_per_block  = ns_per_iteration / 2.0;
    double ns_per_byte = (double)total_time_ns / total_bytes;
    double cycles_per_byte = (double)(end_cycles - start_cycles) / total_bytes;
    double throughput_gbps = ((double)total_bytes / (1024.0 * 1024.0 * 1024.0)) / ((double)total_time_ns / 1e9);

    printf("=== ChaCha20 AVX2 Benchmark (Nanoseconds) ===\n");
    printf("Total Processed Data : %lu MB\n", total_bytes / (1024 * 1024));
    printf("Total Execution Time : %lu ns\n", total_time_ns);
    printf("Time per Iteration   : %.2f ns\n", ns_per_iteration);
    printf("Time per Block (64B) : %.2f ns\n", ns_per_block);
    printf("Time per Byte        : %.4f ns/byte\n", ns_per_byte);
    printf("Cycles per Byte      : %.2f cpb\n", cycles_per_byte);
    printf("Throughput           : %.2f GiB/s\n", throughput_gbps);

    return 0;
}
//to compile: gcc -O3 -march=native chacha20.c -o chacha20
//reference
//website: https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html#techs=AVX_ALL
