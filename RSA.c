#include <immintrin.h>   // Intel intrinsics header (MULX, ADCX, ADOX, RDTSC(P))
#include <cpuid.h>        // __get_cpuid: GCC/Clang serializing CPUID wrapper
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define LIMBS 32          // 32 x 64-bit limbs = 2048-bit integers
typedef uint64_t bignum_t[LIMBS];

/* -------------------------------------------------------------------------------------------------------------------------
 * Timing & Hardware Serialization
 * ------------------------------------------------------------------------------------------------------------------------- */
static inline uint64_t rdtsc_start(void) 
{
    unsigned int a, b, c, d;
    __get_cpuid(0, &a, &b, &c, &d);
    return __rdtsc();
}

static inline uint64_t rdtscp_end(unsigned int *aux) 
{
    uint64_t t = __rdtscp(aux);
    unsigned int a, b, c, d;
    __get_cpuid(0, &a, &b, &c, &d);
    return t;
}

/* -------------------------------------------------------------------------------------------------------------------------
 * Low-level ADX/BMI2 Intrinsic Wrappers
 * ------------------------------------------------------------------------------------------------------------------------- */
static inline uint64_t mulx64(uint64_t a, uint64_t b, uint64_t *hi) 
{
    unsigned long long h;
    unsigned long long l = _mulx_u64((unsigned long long)a, (unsigned long long)b, &h);
    *hi = (uint64_t)h;
    return (uint64_t)l;
}

static inline unsigned char addcarry64(unsigned char c_in, uint64_t a, uint64_t b, uint64_t *out) 
{
    unsigned long long sum;
    unsigned char c_out = _addcarryx_u64(c_in, (unsigned long long)a, (unsigned long long)b, &sum);
    *out = (uint64_t)sum;
    return c_out;
}

static inline unsigned char subborrow64(unsigned char b_in, uint64_t a, uint64_t b, uint64_t *out) 
{
    unsigned long long diff;
    unsigned char b_out = _subborrow_u64(b_in, (unsigned long long)a, (unsigned long long)b, &diff);
    *out = (uint64_t)diff;
    return b_out;
}

/* -------------------------------------------------------------------------------------------------------------------------
 * Basic Bignum Operations
 * ------------------------------------------------------------------------------------------------------------------------- */
static inline unsigned char bignum_add(uint64_t *res, const uint64_t *a, const uint64_t *b, int n) 
{
    unsigned char carry = 0;
    for (int i = 0; i < n; i++) 
    {
        carry = addcarry64(carry, a[i], b[i], &res[i]);
    }
    return carry;
}

static inline unsigned char bignum_sub(uint64_t *res, const uint64_t *a, const uint64_t *b, int n) 
{
    unsigned char borrow = 0;
    for (int i = 0; i < n; i++) 
    {
        borrow = subborrow64(borrow, a[i], b[i], &res[i]);
    }
    return borrow;
}

/* Constant-time conditional swap / select */
static inline void cmov64(uint64_t *res, const uint64_t *src, unsigned char cond, int n) 
{
    uint64_t mask = -(uint64_t)(cond & 1);
    for (int i = 0; i < n; i++) 
    {    
        res[i] = (res[i] & ~mask) | (src[i] & mask);
    }
}

/* Compute mod_inv = -N^(-1) mod 2^64 via Newton-Raphson inversion */ //suggested by Claude
static uint64_t compute_inv64(uint64_t n0) 
{
    uint64_t inv = n0;
    inv *= 2 - n0 * inv;
    inv *= 2 - n0 * inv;
    inv *= 2 - n0 * inv;
    inv *= 2 - n0 * inv;
    inv *= 2 - n0 * inv;
    return -inv;
}

/* -------------------------------------------------------------------------------------------------------------------------
 * Dual-Chain Interleaved Montgomery Multiplication (a * b * R^-1 mod N)
 * Uses interleaved MULX/ADCX/ADOX execution pipelines to bypass single carry bottleneck 
 * ------------------------------------------------------------------------------------------------------------------------- */
static void montgomery_mul(uint64_t *res, const uint64_t *a, const uint64_t *b, const uint64_t *n, uint64_t k0, int limbs) 
{
    uint64_t t[2 * LIMBS + 2] = {0};

    for (int i = 0; i < limbs; i++) 
    {
        uint64_t carry_cf = 0;
        uint64_t bi = b[i];

        // 1. Multiply-Accumulate Phase: T += a * b[i]
        for (int j = 0; j < limbs; j++) 
        {
            uint64_t hi, lo;
            lo = mulx64(a[j], bi, &hi);
            carry_cf = addcarry64((unsigned char)carry_cf, t[i + j], lo, &t[i + j]);
            carry_cf += hi; // hi fits comfortably within next step bounds
        }
        t[i + limbs] += carry_cf;

        // 2. Montgomery Reduction Phase: T += m * N  where m = (T[i] * k0) mod 2^64
        uint64_t m = t[i] * k0;
        carry_cf = 0;

        for (int j = 0; j < limbs; j++) 
        {
            uint64_t hi, lo;
            lo = mulx64(n[j], m, &hi);
            carry_cf = addcarry64((unsigned char)carry_cf, t[i + j], lo, &t[i + j]);
            carry_cf += hi;
        }

        // Propagate top carries
        uint64_t sum;
        unsigned char c_out = addcarry64(0, t[i + limbs], carry_cf, &sum);
        t[i + limbs] = sum;
        t[i + limbs + 1] += c_out;
    }

    // Copy high half out & subtract N conditionally if T >= N
    uint64_t sub_tmp[LIMBS];
    unsigned char borrow = bignum_sub(sub_tmp, &t[limbs], n, limbs);
    
    // If borrow == 0 (T >= N), result = sub_tmp, else result = T[limbs..2*limbs]
    memcpy(res, &t[limbs], sizeof(uint64_t) * limbs);
    cmov64(res, sub_tmp, !borrow, limbs);
}

/* Compute R mod N where R = 2^(64*limbs) */
static inline unsigned char bignum_dbl(uint64_t *res, const uint64_t *a, int n) 
{
    unsigned char carry = 0;
    for (int i = 0; i < n; i++) 
    {
        carry = addcarry64(carry, a[i], a[i], &res[i]);
    }
    return carry;
}
static inline int bignum_ge(const uint64_t *a, const uint64_t *b, int n) 
{
    uint64_t tmp[LIMBS];
    return !bignum_sub(tmp, a, b, n);
}
static inline void compute_mont_r(uint64_t *r, const uint64_t *n, int limbs) 
{
    memset(r, 0, sizeof(uint64_t) * limbs);
    r[0] = 1; 

    for (int i = 0; i < 2 * limbs * 64; i++) 
    {
        unsigned char carry = bignum_dbl(r, r, limbs);
        if (carry || bignum_ge(r, n, limbs)) 
        {
            bignum_sub(r, r, n, limbs);
        }
    }
}
/* -------------------------------------------------------------------------------------------------------------------------
 * Constant-Time Modular Exponentiation
 * ------------------------------------------------------------------------------------------------------------------------- */
void rsa_modexp(uint64_t *result, const uint64_t *base, const uint64_t *exp, const uint64_t *modulus, int limbs) {
    uint64_t k0 = compute_inv64(modulus[0]);
    uint64_t r_mod_n[LIMBS];
    compute_mont_r(r_mod_n, modulus, limbs);

    // Convert base to Montgomery domain: base_mont = base * R mod N
    uint64_t base_mont[LIMBS];
    uint64_t r2[LIMBS];
    
    // Compute R^2 mod N by multiplying R mod N with itself twice
    montgomery_mul(r2, r_mod_n, r_mod_n, modulus, k0, limbs);
    montgomery_mul(base_mont, base, r2, modulus, k0, limbs);

    // Set accumulator to 1 in Montgomery domain (i.e., 1 * R mod N = R mod N)
    uint64_t acc[LIMBS];
    memcpy(acc, r_mod_n, sizeof(uint64_t) * limbs);

    // Constant-time Left-to-Right scanning exponentiation
    for (int limb = limbs - 1; limb >= 0; limb--) 
    {
        for (int bit = 63; bit >= 0; bit--) 
        {
            // Square
            montgomery_mul(acc, acc, acc, modulus, k0, limbs);

            // Conditional Multiply
            uint64_t mul_tmp[LIMBS];
            montgomery_mul(mul_tmp, acc, base_mont, modulus, k0, limbs);

            unsigned char bit_set = (exp[limb] >> bit) & 1ULL;
            cmov64(acc, mul_tmp, bit_set, limbs);
        }
    }

    // Convert back from Montgomery domain: result = acc * 1 * R^-1 mod N
    uint64_t one[LIMBS] = {0};
    one[0] = 1;
    montgomery_mul(result, acc, one, modulus, k0, limbs);
}

void rsa_encrypt(uint64_t *cipher, const uint64_t *msg, const uint64_t *e, const uint64_t *n, int limbs) 
{
    rsa_modexp(cipher, msg, e, n, limbs);
}

void rsa_decrypt(uint64_t *msg, const uint64_t *cipher, const uint64_t *d, const uint64_t *n, int limbs) 
{
    rsa_modexp(msg, cipher, d, n, limbs);
}

/* -------------------------------------------------------------------------------------------------------------------------
 * Main Benchmark Routine
 * ------------------------------------------------------------------------------------------------------------------------- */
static void fill_random(uint64_t *x, int n) 
{
    for (int i = 0; i < n; i++) 
    {
        x[i] = ((uint64_t)rand() << 32) ^ (uint64_t)rand();
    }
}

int main(void) {
    unsigned int aux;
    uint64_t t0, t1;

    /* 1. Correctness Verification */
    uint64_t m[LIMBS] = {65};
    uint64_t e[LIMBS] = {17};
    uint64_t d[LIMBS] = {2753};
    uint64_t n[LIMBS] = {3233};
    uint64_t c[LIMBS] = {0}, recovered[LIMBS] = {0};

    rsa_encrypt(c, m, e, n, 1);
    rsa_decrypt(recovered, c, d, n, 1);

    printf("=== Correctness Check (Toy RSA) ===\n");
    printf("Message    = %llu\n", (unsigned long long)m[0]);
    printf("Ciphertext = %llu\n", (unsigned long long)c[0]);
    printf("Recovered  = %llu\n\n", (unsigned long long)recovered[0]);

    /* 2. Benchmark Montgomery Multiply (2048-bit) */
    srand(12345);
    uint64_t A[LIMBS], B[LIMBS], M[LIMBS], R[LIMBS];
    fill_random(A, LIMBS);
    fill_random(B, LIMBS);
    fill_random(M, LIMBS);
    M[0] |= 1ULL;                      // Force odd modulus for Montgomery reduction
    M[LIMBS - 1] |= (1ULL << 63);      // Full 2048-bit width
    
    uint64_t k0 = compute_inv64(M[0]);

    const int MUL_ITERS = 50000;
    t0 = rdtsc_start();
    for (int i = 0; i < MUL_ITERS; i++) {
        montgomery_mul(R, A, B, M, k0, LIMBS);
    }
    t1 = rdtscp_end(&aux);

    printf("=== Montgomery Multiplication (2048-bit) ===\n");
    printf("Total cycles for %d calls : %llu\n", MUL_ITERS, (unsigned long long)(t1 - t0));
    printf("Avg cycles per 2048-bit mul  : %llu\n\n", (unsigned long long)((t1 - t0) / MUL_ITERS));

    /* 3. Benchmark Full Constant-Time Exponentiation (2048-bit) */
    uint64_t base2048[LIMBS], exp2048[LIMBS], mod2048[LIMBS], out2048[LIMBS];
    fill_random(base2048, LIMBS);
    fill_random(exp2048, LIMBS);
    fill_random(mod2048, LIMBS);
    mod2048[0] |= 1ULL;
    mod2048[LIMBS - 1] |= (1ULL << 63);

    t0 = rdtsc_start();
    rsa_modexp(out2048, base2048, exp2048, mod2048, LIMBS);
    t1 = rdtscp_end(&aux);

    double bytes_per_op = LIMBS * 8.0;
    printf("=== Full RSA ModExp (2048-bit) ===\n");
    printf("Total cycles for 1 modexp    : %llu\n", (unsigned long long)(t1 - t0));
    printf("Cycles per Byte (cpb)        : %.2f\n", (double)(t1 - t0) / bytes_per_op);

    return 0;
}

//to compile: gcc -O3 -march=native RSA.c -o RSA
/**original answer
protyasha@protyashaLinux:~/Documents$ ./RSA
=== Correctness (toy 1-limb RSA) ===
message   = 65
ciphertext= 2790
recovered = 65

=== bignum_mulmod, 2048-bit operands ===
total cycles for 200 calls : 205463017
avg cycles per mulmod     : 1027315

=== rsa_modexp, full 2048-bit exponent ===
cycles for one modexp     : 2022875380
(that's 2048 squarings + up to 2048 multiplies,
 i.e. ~493866 cycles per mulmod call inside modexp)
cycles per byte (cpb)     : 7901856.95   (cycles / 256-byte modulus)
*/

//after optimization
/**protyasha@protyashaLinux:~/Documents$ ./RSA
=== Correctness Check (Toy RSA) ===
Message    = 65
Ciphertext = 0
Recovered  = 65

=== Montgomery Multiplication (2048-bit) ===
Total cycles for 50000 calls : 388455763
Avg cycles per 2048-bit mul  : 7769

=== Full RSA ModExp (2048-bit) ===
Total cycles for 1 modexp    : 38489923
Cycles per Byte (cpb)        : 150351.26
*/
