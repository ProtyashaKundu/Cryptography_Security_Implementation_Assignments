
#include <immintrin.h>   // AES-NI intrinsics (AESENC, AESENCLAST, AESKEYGENASSIST)
#include <cpuid.h>        // __get_cpuid: serializing barrier for RDTSC
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

/* ==================================================================
 * Cycle counting (same pattern used for the RSA benchmark)
 * ================================================================== */
static inline uint64_t rdtsc_start(void) 
{
    unsigned int a, b, c, d;
    __get_cpuid(0, &a, &b, &c, &d);    // serialize: drain the pipeline
    return __rdtsc();                  // start timestamp
}
static inline uint64_t rdtscp_end(unsigned int *aux) 
{
    uint64_t t = __rdtscp(aux);        // end timestamp, waits for retire
    unsigned int a, b, c, d;
    __get_cpuid(0, &a, &b, &c, &d);    // serialize again
    return t;
}

/* ==================================================================
 * Key expansion
* ================================================================== */

/* ---- AES-128: 11 round keys (Nr = 10 rounds) -------------------- */
static inline __m128i expand_step_128(__m128i key, __m128i assist) 
{
    // assist already holds {SubWord(RotWord(temp)) ^ Rcon} broadcast
    // to all 4 dwords by the caller's 0xff shuffle below.
    __m128i t = key;
    assist = _mm_shuffle_epi32(assist, 0xff);      // broadcast last dword
    t = _mm_xor_si128(t, _mm_slli_si128(t, 4));     // propagate word 0->1
    t = _mm_xor_si128(t, _mm_slli_si128(t, 4));     // ->2
    t = _mm_xor_si128(t, _mm_slli_si128(t, 4));     // ->3
    return _mm_xor_si128(t, assist);
}

static void aes128_key_expansion(const uint8_t *key, __m128i round_keys[11]) 
{
    round_keys[0] = _mm_loadu_si128((const __m128i *)key);

    // AESKEYGENASSIST's immediate is the round constant (Rcon) for
    // that key-schedule step: 0x01,0x02,0x04,...,0x36 (xtime powers).
    round_keys[1]  = expand_step_128(round_keys[0], _mm_aeskeygenassist_si128(round_keys[0], 0x01));
    round_keys[2]  = expand_step_128(round_keys[1], _mm_aeskeygenassist_si128(round_keys[1], 0x02));
    round_keys[3]  = expand_step_128(round_keys[2], _mm_aeskeygenassist_si128(round_keys[2], 0x04));
    round_keys[4]  = expand_step_128(round_keys[3], _mm_aeskeygenassist_si128(round_keys[3], 0x08));
    round_keys[5]  = expand_step_128(round_keys[4], _mm_aeskeygenassist_si128(round_keys[4], 0x10));
    round_keys[6]  = expand_step_128(round_keys[5], _mm_aeskeygenassist_si128(round_keys[5], 0x20));
    round_keys[7]  = expand_step_128(round_keys[6], _mm_aeskeygenassist_si128(round_keys[6], 0x40));
    round_keys[8]  = expand_step_128(round_keys[7], _mm_aeskeygenassist_si128(round_keys[7], 0x80));
    round_keys[9]  = expand_step_128(round_keys[8], _mm_aeskeygenassist_si128(round_keys[8], 0x1B));
    round_keys[10] = expand_step_128(round_keys[9], _mm_aeskeygenassist_si128(round_keys[9], 0x36));
}

/* ---- AES-192: 13 round keys (Nr = 12 rounds) -------------------- */

static inline void expand_step_192(__m128i *temp1, __m128i *temp2, __m128i *temp3) 
{
    __m128i t4;
    *temp2 = _mm_shuffle_epi32(*temp2, 0x55);
    t4 = _mm_slli_si128(*temp1, 4);  *temp1 = _mm_xor_si128(*temp1, t4);
    t4 = _mm_slli_si128(t4, 4);      *temp1 = _mm_xor_si128(*temp1, t4);
    t4 = _mm_slli_si128(t4, 4);      *temp1 = _mm_xor_si128(*temp1, t4);
    *temp1 = _mm_xor_si128(*temp1, *temp2);
    *temp2 = _mm_shuffle_epi32(*temp1, 0xff);
    t4 = _mm_slli_si128(*temp3, 4);  *temp3 = _mm_xor_si128(*temp3, t4);
    *temp3 = _mm_xor_si128(*temp3, *temp2);
}

static void aes192_key_expansion(const uint8_t *key, __m128i round_keys[13]) 
{
    __m128i temp1 = _mm_loadu_si128((const __m128i *)key);          // key[0..15]
    __m128i temp3 = _mm_loadu_si128((const __m128i *)(key + 16));   // low 64 bits used = key[16..23]
    __m128i temp2;

    round_keys[0] = temp1;
    round_keys[1] = temp3;

    temp2 = _mm_aeskeygenassist_si128(temp3, 0x01);
    expand_step_192(&temp1, &temp2, &temp3);
    round_keys[1] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(round_keys[1]), _mm_castsi128_pd(temp1), 0));
    round_keys[2] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(temp1), _mm_castsi128_pd(temp3), 1));

    temp2 = _mm_aeskeygenassist_si128(temp3, 0x02);
    expand_step_192(&temp1, &temp2, &temp3);
    round_keys[3] = temp1;
    round_keys[4] = temp3;

    temp2 = _mm_aeskeygenassist_si128(temp3, 0x04);
    expand_step_192(&temp1, &temp2, &temp3);
    round_keys[4] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(round_keys[4]), _mm_castsi128_pd(temp1), 0));
    round_keys[5] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(temp1), _mm_castsi128_pd(temp3), 1));

    temp2 = _mm_aeskeygenassist_si128(temp3, 0x08);
    expand_step_192(&temp1, &temp2, &temp3);
    round_keys[6] = temp1;
    round_keys[7] = temp3;

    temp2 = _mm_aeskeygenassist_si128(temp3, 0x10);
    expand_step_192(&temp1, &temp2, &temp3);
    round_keys[7] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(round_keys[7]), _mm_castsi128_pd(temp1), 0));
    round_keys[8] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(temp1), _mm_castsi128_pd(temp3), 1));

    temp2 = _mm_aeskeygenassist_si128(temp3, 0x20);
    expand_step_192(&temp1, &temp2, &temp3);
    round_keys[9]  = temp1;
    round_keys[10] = temp3;

    temp2 = _mm_aeskeygenassist_si128(temp3, 0x40);
    expand_step_192(&temp1, &temp2, &temp3);
    round_keys[10] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(round_keys[10]), _mm_castsi128_pd(temp1), 0));
    round_keys[11] = _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(temp1), _mm_castsi128_pd(temp3), 1));

    temp2 = _mm_aeskeygenassist_si128(temp3, 0x80);
    expand_step_192(&temp1, &temp2, &temp3);
    round_keys[12] = temp1;
}

/* ---- AES-256: 15 round keys (Nr = 14 rounds) -------------------- */
static inline __m128i expand_step_256_even(__m128i key, __m128i assist) 
{
    return expand_step_128(key, assist);   // identical recurrence to AES-128
}
static inline __m128i expand_step_256_odd(__m128i key_lo, __m128i key_hi) 
{
    __m128i t = _mm_aeskeygenassist_si128(key_hi, 0x00);
    __m128i sub = _mm_shuffle_epi32(t, 0xaa);        // take SubWord(word[3]), no Rcon
    __m128i r = key_lo;
    r = _mm_xor_si128(r, _mm_slli_si128(r, 4));
    r = _mm_xor_si128(r, _mm_slli_si128(r, 4));
    r = _mm_xor_si128(r, _mm_slli_si128(r, 4));
    return _mm_xor_si128(r, sub);
}

static void aes256_key_expansion(const uint8_t *key, __m128i round_keys[15]) 
{
    round_keys[0] = _mm_loadu_si128((const __m128i *)key);
    round_keys[1] = _mm_loadu_si128((const __m128i *)(key + 16));

    round_keys[2]  = expand_step_256_even(round_keys[0], _mm_aeskeygenassist_si128(round_keys[1], 0x01));
    round_keys[3]  = expand_step_256_odd (round_keys[1], round_keys[2]);
    round_keys[4]  = expand_step_256_even(round_keys[2], _mm_aeskeygenassist_si128(round_keys[3], 0x02));
    round_keys[5]  = expand_step_256_odd (round_keys[3], round_keys[4]);
    round_keys[6]  = expand_step_256_even(round_keys[4], _mm_aeskeygenassist_si128(round_keys[5], 0x04));
    round_keys[7]  = expand_step_256_odd (round_keys[5], round_keys[6]);
    round_keys[8]  = expand_step_256_even(round_keys[6], _mm_aeskeygenassist_si128(round_keys[7], 0x08));
    round_keys[9]  = expand_step_256_odd (round_keys[7], round_keys[8]);
    round_keys[10] = expand_step_256_even(round_keys[8], _mm_aeskeygenassist_si128(round_keys[9], 0x10));
    round_keys[11] = expand_step_256_odd (round_keys[9], round_keys[10]);
    round_keys[12] = expand_step_256_even(round_keys[10], _mm_aeskeygenassist_si128(round_keys[11], 0x20));
    round_keys[13] = expand_step_256_odd (round_keys[11], round_keys[12]);
    round_keys[14] = expand_step_256_even(round_keys[12], _mm_aeskeygenassist_si128(round_keys[13], 0x40));
}

/* ===============================================================================================================
 * Single-block AES encryption
 * Purpose : Encrypt one 128-bit block. This IS the CTR-mode keystream
 *           generator (CTR encrypts the COUNTER, never the plaintext,
 *           then XORs the result with plaintext/ciphertext).
 * Intrinsics:
 *   _mm_xor_si128     : initial AddRoundKey (whitening) -- plain SSE2, not AES-NI, but required before round 1.
 *   _mm_aesenc_si128  : one FULL AES round (SubBytes+ShiftRows+MixColumns+AddRoundKey) in a single instruction.
 *   _mm_aesenclast_si128 : the FINAL round, which omits MixColumns per the AES spec -- a separate instruction
 *                        because its datapath differs from AESENC.
 * ===============================================================================================================*/
static inline __m128i aes_encrypt_block(__m128i block, const __m128i *rk, int rounds)
{
    block = _mm_xor_si128(block, rk[0]);          // initial whitening
    for (int i = 1; i < rounds; i++)
        block = _mm_aesenc_si128(block, rk[i]);   // rounds 1..Nr-1
    return _mm_aesenclast_si128(block, rk[rounds]); // final round Nr
}

/* ==================================================================
 * CTR mode
 * Purpose : Turns the AES block cipher into a stream cipher.
 *           keystream_i = AES_encrypt(counter_i); ciphertext = pt ^ ks.
 *           Decryption is IDENTICAL (XOR is self-inverse), so the same
 *           function does both directions.
 * Note: counter increment is done on a plain byte array (big-endian,
 * matching the common CTR convention / NIST SP 800-38A) since it is
 * not the bottleneck -- the AES rounds above are.
 * ================================================================== */
static void increment_counter(uint8_t ctr[16]) 
{
    for (int i = 15; i >= 0; i--) 
    {
        if (++ctr[i] != 0) break;   // stop unless this byte wrapped to 0
    }
}

static void aes_ctr_crypt(uint8_t *out, const uint8_t *in, size_t len, const __m128i *round_keys, int rounds,const uint8_t nonce_ctr[16])
 {
    uint8_t ctr[16];
    memcpy(ctr, nonce_ctr, 16);

    size_t offset = 0;
    while (offset < len) {
        __m128i ctr_block = _mm_loadu_si128((const __m128i *)ctr);
        __m128i keystream = aes_encrypt_block(ctr_block, round_keys, rounds);

        size_t chunk = (len - offset < 16) ? (len - offset) : 16;
        uint8_t ks_bytes[16];
        _mm_storeu_si128((__m128i *)ks_bytes, keystream);
        for (size_t i = 0; i < chunk; i++)
            out[offset + i] = in[offset + i] ^ ks_bytes[i];

        increment_counter(ctr);
        offset += chunk;
    }
}

/* ==================================================================
 * Correctness self-test using FIPS-197 Appendix test vectors
 * (single-block ECB encryption -- validates key expansion + rounds)
 * ================================================================== */
static void print_hex(const char *label, const uint8_t *b, int n) 
{
    printf("%-12s", label);
    for (int i = 0; i < n; i++) printf("%02x", b[i]);
    printf("\n");
}

static int check_block(const uint8_t *key_bytes, int key_bits, const uint8_t *pt, const uint8_t *expected_ct) 
{
    __m128i rk[15];
    int rounds;
    if (key_bits == 128)      { aes128_key_expansion(key_bytes, rk); rounds = 10; }
    else if (key_bits == 192) { aes192_key_expansion(key_bytes, rk); rounds = 12; }
    else                      { aes256_key_expansion(key_bytes, rk); rounds = 14; }

    __m128i block = _mm_loadu_si128((const __m128i *)pt);
    __m128i ct = aes_encrypt_block(block, rk, rounds);
    uint8_t out[16];
    _mm_storeu_si128((__m128i *)out, ct);

    int ok = memcmp(out, expected_ct, 16) == 0;
    printf("AES-%d ECB test: %s\n", key_bits, ok ? "PASS" : "FAIL");
    if (!ok) print_hex("got:", out, 16);
    return ok;
}

/* ==================================================================
 * Benchmark: CTR-mode throughput / cycles-per-byte for each key size
 * ================================================================== */
static void benchmark_ctr(int key_bits, size_t buf_len) 
{
    uint8_t key[32];
    for (int i = 0; i < 32; i++) 
    	key[i] = (uint8_t)(i * 0x11 + key_bits);

    __m128i rk[15];
    int rounds;
    if (key_bits == 128)      
    { 
    	aes128_key_expansion(key, rk); rounds = 10;
     }
    else if (key_bits == 192) 
    { 
      aes192_key_expansion(key, rk); rounds = 12; 
    }
    else                      
    { 
    	aes256_key_expansion(key, rk); rounds = 14; 
    }

    uint8_t nonce_ctr[16] = {0};
    uint8_t *plaintext = malloc(buf_len);
    uint8_t *ciphertext = malloc(buf_len);
    uint8_t *decrypted = malloc(buf_len);
    for (size_t i = 0; i < buf_len; i++) 
    	plaintext[i] = (uint8_t)(rand() & 0xff);

    unsigned int aux;
    uint64_t t0 = rdtsc_start();
    aes_ctr_crypt(ciphertext, plaintext, buf_len, rk, rounds, nonce_ctr);
    uint64_t t1 = rdtscp_end(&aux);

    // Round-trip check: CTR decrypt = CTR encrypt again with same counter
    aes_ctr_crypt(decrypted, ciphertext, buf_len, rk, rounds, nonce_ctr);
    int roundtrip_ok = memcmp(plaintext, decrypted, buf_len) == 0;

    uint64_t cycles = t1 - t0;
    double cycles_per_byte = (double)cycles / (double)buf_len;

    printf("\n=== AES-%d-CTR, %zu KiB buffer ===\n", key_bits, buf_len / 1024);
    printf("round-trip decrypt : %s\n", roundtrip_ok ? "OK" : "MISMATCH");
    printf("total cycles        : %llu\n", (unsigned long long)cycles);
    printf("cycles per byte      : %.2f\n", cycles_per_byte);

    free(plaintext); free(ciphertext); free(decrypted);
}

int main(void) 
{
    /* ---- FIPS-197 known-answer tests, one per key size ---- */
    uint8_t pt[16] = {0x00,0x11,0x22,0x33,0x44,0x55,0x66,0x77,
                       0x88,0x99,0xaa,0xbb,0xcc,0xdd,0xee,0xff};

    uint8_t key128[16] = {0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
                           0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f};
    uint8_t ct128[16]  = {0x69,0xc4,0xe0,0xd8,0x6a,0x7b,0x04,0x30,
                           0xd8,0xcd,0xb7,0x80,0x70,0xb4,0xc5,0x5a};

    uint8_t key192[24] = {0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
                           0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,
                           0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17};
    uint8_t ct192[16]  = {0xdd,0xa9,0x7c,0xa4,0x86,0x4c,0xdf,0xe0,
                           0x6e,0xaf,0x70,0xa0,0xec,0x0d,0x71,0x91};

    uint8_t key256[32] = {0x00,0x01,0x02,0x03,0x04,0x05,0x06,0x07,
                           0x08,0x09,0x0a,0x0b,0x0c,0x0d,0x0e,0x0f,
                           0x10,0x11,0x12,0x13,0x14,0x15,0x16,0x17,
                           0x18,0x19,0x1a,0x1b,0x1c,0x1d,0x1e,0x1f};
    uint8_t ct256[16]  = {0x8e,0xa2,0xb7,0xca,0x51,0x67,0x45,0xbf,
                           0xea,0xfc,0x49,0x90,0x4b,0x49,0x60,0x89};

    printf("=== Correctness (FIPS-197 known-answer vectors) ===\n");
    check_block(key128, 128, pt, ct128);
    check_block(key192, 192, pt, ct192);
    check_block(key256, 256, pt, ct256);

    /* ---- Cycle-count benchmarks: 1 MiB CTR encrypt for each key size --- */
    srand(42);
    const size_t BUF_LEN = 1024 * 1024;   // 1 MiB
    benchmark_ctr(128, BUF_LEN);
    benchmark_ctr(192, BUF_LEN);
    benchmark_ctr(256, BUF_LEN);

    return 0;
}
/**
protyasha@protyashaLinux:~/Documents$ ./AES
=== Correctness (FIPS-197 known-answer vectors) ===
AES-128 ECB test: PASS
AES-192 ECB test: PASS
AES-256 ECB test: PASS

=== AES-128-CTR, 1024 KiB buffer ===
round-trip decrypt : OK
total cycles        : 5703520
cycles per byte      : 5.44

=== AES-192-CTR, 1024 KiB buffer ===
round-trip decrypt : OK
total cycles        : 5624943
cycles per byte      : 5.36

=== AES-256-CTR, 1024 KiB buffer ===
round-trip decrypt : OK
total cycles        : 5881534
cycles per byte      : 5.61

*/
/*
 * ------------------------------------------------------------------
 * Function/intrinsic summary
 * ------------------------------------------------------------------
 * _mm_aeskeygenassist_si128 : hardware key-schedule step (SubWord +
 *                              RotWord + Rcon XOR in one instruction).
 * _mm_aesenc_si128          : one full AES encryption round.
 * _mm_aesenclast_si128      : final AES round (no MixColumns).
 * _mm_xor_si128 / _mm_slli_si128 / _mm_shuffle_epi32 / _mm_shuffle_pd :
 *      plain SSE2/SSE4.1 data-movement ops used to route key-schedule
 *      words into the right lanes; not AES-specific but required
 *      scaffolding around the AES-NI instructions.
 * __rdtsc / __rdtscp / __get_cpuid : cycle-accurate timing, same
 *      serializing pattern used in the RSA benchmark file.
 *
 * aes128/192/256_key_expansion -> build the round-key schedule
 * aes_encrypt_block            -> encrypt one 16-byte block (= CTR
 *                                  keystream generator)
 * increment_counter            -> big-endian 128-bit counter increment
 * aes_ctr_crypt                -> full CTR-mode encrypt/decrypt loop
 * benchmark_ctr                -> times a 1 MiB CTR pass, reports
 *                                  total cycles and cycles/byte
 
 */
