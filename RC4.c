//Protyasha Kundu, CrS2510

//to compile: gcc -O3 -march=native RC4.c -o RC4

#include <stdio.h>
#include <string.h>
#include <stdint.h>

// State structure holding the 256-byte S-box and current pointers i and j
typedef struct 
{
    unsigned char S[256];
    unsigned char i;
    unsigned char j;
} RC4_State;

// Read Time-Stamp Counter (x86_64 inline assembly)
// Measures CPU clock cycles directly at low overhead
static inline uint64_t rdtsc(void) 
{
    unsigned int lo, hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}

// Plain temp-variable swap.
/** the classic "XOR swap" (a^=b; b^=a; a^=b;) looks faster and more enticing to use but each line depends on the result of the previous line, so the three XORs
form a serial dependency chain the CPU can't pipeline or reorder. A temp variable breaks that chain: the load of *a and *b can happen in parallel, and modern 
compilers turn this into 2-3 independent register moves. It's faster despite "using an extra register"-registers are free, latency isn't.*/
static inline void swap(unsigned char *a, unsigned char *b) 
{
    unsigned char tmp = *a;
    *a = *b;
    *b = tmp;
}

// Key Scheduling Algorithm (KSA) with cycle-per-byte benchmarking.
// Returns cycles-per-byte for the 256-byte mixing loop.
double KSA(RC4_State *state, const unsigned char *key, size_t key_len) 
{
    unsigned char K[256];

    for (int i = 0; i < 256; i++) 
    {
        state->S[i] = (unsigned char)i;
    }

    for (int i = 0; i < 256; i++) 
    {
        K[i] = key[i % key_len];
    }

    // Use a local pointer so the compiler doesn't have to re-check aliasing through state-> on every iteration.
    //caching Si = S[i] and Sj = S[j] before the swap, so afterward you reuse those register values instead of reading S[i]/S[j] again.
    unsigned char *S = state->S;

    uint64_t start_cycles = rdtsc();
    unsigned char j = 0;
    for (int i = 0; i < 256; i++) 
    {
        j = (unsigned char)(j + S[i] + K[i]);
        swap(&S[i], &S[j]);
    }
    uint64_t end_cycles = rdtsc();

    state->i = 0;
    state->j = 0;

    return (double)(end_cycles - start_cycles) / 256.0;
}

// Pseudo-Random Generation Algorithm (PRGA) & Stream Encryption with CPB Benchmarking
double PRGA(RC4_State *state, const unsigned char *input, unsigned char *output, size_t len) 
{
    unsigned char i = state->i;
    unsigned char j = state->j;
    unsigned char *S = state->S; // local pointer: fewer indirections through state->

    uint64_t start_cycles = rdtsc();
    for (size_t t = 0; t < len; t++) 
    {
        i = (unsigned char)(i + 1);
        j = (unsigned char)(j + S[i]);

        // Cache S[i]/S[j] BEFORE swapping, then swap directly, then use the cached values fo the keystream index. This avoids
        // re-reading S[i] and S[j] from memory after the swap - the values are already sitting in registers.
        unsigned char Si = S[i];
        unsigned char Sj = S[j];
        S[i] = Sj;
        S[j] = Si;

        unsigned char z = S[(unsigned char)(Si + Sj)];
        output[t] = input[t] ^ z;
    }
    uint64_t end_cycles = rdtsc();

    state->i = i;
    state->j = j;

    return (double)(end_cycles - start_cycles) / (double)len;
}

// Correctness Check Routine using Known Test Vectors & Decryption Roundtrip
int verify_correctness(void) 
{
    printf("=== RUNNING CORRECTNESS CHECKS ===\n");

    // Standard RC4 Test Vector:
    // Key: "Key" (0x4B, 0x65, 0x79)
    // Plaintext: "Plaintext"
    // Expected Ciphertext (Hex): BB F3 16 E8 D9 40 AF 0A D3
    const unsigned char test_key[] = "Key";
    const unsigned char test_pt[] = "Plaintext";
    const unsigned char expected_ct[9] = {0xBB, 0xF3, 0x16, 0xE8, 0xD9, 0x40, 0xAF, 0x0A, 0xD3};
    
    size_t test_len = sizeof(test_pt) - 1;
    unsigned char computed_ct[16] = {0};
    unsigned char decrypted_pt[16] = {0};

    RC4_State state;

    // Test 1: Encryption Check
    KSA(&state, test_key, 3);
    PRGA(&state, test_pt, computed_ct, test_len);

    if (memcmp(computed_ct, expected_ct, test_len) != 0) 
    {
        printf("[FAIL] Known-Answer Test (KAT) Failed!\n");
        printf("  Expected: ");
        for (size_t k = 0; k < test_len; k++) 
            printf("%02X ", expected_ct[k]);
        printf("\n  Computed: ");
        for (size_t k = 0; k < test_len; k++) 
            printf("%02X ", computed_ct[k]);
        printf("\n");
        return 0;
    }
    printf("[PASS] Known-Answer Test (KAT) Passed.\n");

    // Test 2: Decryption / Round-Trip Check (Enc(Dec(M)) == M)
    KSA(&state, test_key, 3);
    PRGA(&state, computed_ct, decrypted_pt, test_len);

    if (memcmp(decrypted_pt, test_pt, test_len) != 0) 
    {
        printf("[FAIL] Round-trip Decryption Test Failed!\n");
        return 0;
    }
    printf("[PASS] Round-trip Decryption Check Passed.\n");
    printf("==================================\n\n");

    return 1;
}

int main(void) 
{
    if (!verify_correctness()) 
    {
        fprintf(stderr, "Correctness checks failed. Aborting benchmark.\n");
        return 1;
    }
    unsigned char key[16] = "SecretKey128Bit!";
    
    #define BUF_SIZE 20000
    static unsigned char plaintext[BUF_SIZE];
    static unsigned char ciphertext[BUF_SIZE];
    memset(plaintext, 0xAA, BUF_SIZE);

    RC4_State state;    

    // Warm-up KSA run (cold cache / branch predictor), then timed run.
    KSA(&state, key, sizeof(key));
    double ksa_cpb = KSA(&state, key, sizeof(key));

    // Warm-up PRGA run, then timed run.
    PRGA(&state, plaintext, ciphertext, 1000);
    double prga_cpb = PRGA(&state, plaintext, ciphertext, BUF_SIZE);

    printf("RC4 KSA  Performance: %.2f CPB (Cycles Per Byte, 256-byte schedule)\n", ksa_cpb);
    printf("RC4 PRGA Performance: %.2f CPB (Cycles Per Byte)\n", prga_cpb);

    return 0;
}
/*Output:
protyasha@protyashaLinux:~/Documents$ ./RC4
=== RUNNING CORRECTNESS CHECKS ===
[PASS] Known-Answer Test (KAT) Passed.
[PASS] Round-trip Decryption Check Passed.
==================================

RC4 KSA  Performance: 7.22 CPB (Cycles Per Byte, 256-byte schedule)
RC4 PRGA Performance: 8.39 CPB (Cycles Per Byte)
protyasha@protyashaLinux:~/Documents$ ./RC4
=== RUNNING CORRECTNESS CHECKS ===
[PASS] Known-Answer Test (KAT) Passed.
[PASS] Round-trip Decryption Check Passed.
==================================

RC4 KSA  Performance: 7.24 CPB (Cycles Per Byte, 256-byte schedule)
RC4 PRGA Performance: 8.19 CPB (Cycles Per Byte)
protyasha@protyashaLinux:~/Documents$ ./RC4
=== RUNNING CORRECTNESS CHECKS ===
[PASS] Known-Answer Test (KAT) Passed.
[PASS] Round-trip Decryption Check Passed.
==================================

RC4 KSA  Performance: 7.21 CPB (Cycles Per Byte, 256-byte schedule)
RC4 PRGA Performance: 8.48 CPB (Cycles Per Byte)

*/
/* Argument on the Cycle per Byte
/*
--- KSA INSTRUCTION & CYCLE BREAKDOWN (PER ITERATION / BYTE) ---
 
1. Increment & Loop Check (add i, 1; cmp i, 256): ~0.25-0.5 cycles
   - Independent of everything else in the loop; can be computed for iteration t+1 while iteration t is still in flight.
   
2. Load S[i] & K[i] (movzx eax, [S+i]; movzx ebx, [K+i]): ~1 cycle issue,  ~4 cycle latency
   - Address is just `i`, known ahead of time - NOT dependent on the previous iteration. This means the CPU can start this load early,
     overlapping it with the tail end of the prior iteration's chain.
     
3. Compute j = j + S[i] + K[i] (add j, al; add j, bl): ~1-2 cycles
   - This is the actual bottleneck of the loop: j is a carried dependency across all 256 iterations. j at step t+1 cannot be
     computed until j at step t is finished - this chain runs the full length of the loop and is what ultimately limits overlap between
     iterations (unlike i, which is free-running).
     
4. Load S[j] (movzx ecx, [S+j]): ~4 cycles (L1 latency)
   - Gated on step 3's j. Address unpredictable (pseudo-random), so no effective prefetch.
   
5. Swap S[i] and S[j] (mov [S+i], cl; mov [S+j], al): ~1 cycle to issue
   - Goes into the store queue without stalling the pipeline, but this
     store must retire before a later iteration's load from the same
     address is correct - adds to the effective chain, not just "free."
 
TOTAL EXPECTED KSA PERFORMANCE: ~5-8 CPB
Per-step chain latency (j-add -> load S[j] -> store) is roughly 6-7 cycles if fully serial, but the measured/expected CPB is lower than a
naive sum of every step as some steps 1-2 of iteration t+1 (i-increment, S[i]/K[i] loads) overlap with steps 3-5 of iteration t - they don't
depend on iteration t's j at all. The number you actually see reflects that overlap, not a per-iteration total.
Main bottleneck: the carried-j chain (j_t -> j_{t+1} -> ...) combined with the dependent Load S[j] each step.
 
--- PRGA INSTRUCTION & CYCLE BREAKDOWN (PER ITERATION / BYTE) ---
 
1. Advance i (add i, 1): ~0.25-0.5 cycles
   - Free-running counter, independent across iterations - can be computed arbitrarily before hand.
   
2. Load Si = S[i] (movzx Si, [S+i]): ~4 cycles (L1 latency)
   - Address = i, known ahead of time, so this load for iteration t+1 can be issued early, overlapping with iteration t's later steps.
   
3. Compute j = j + Si (add j, Si): ~1 cycle
   - Depends on step 2's load completing. j itself is also carried from the previous iteration (j_t+1 needs j_t), same structural pattern
     as KSA.
     
4. Load Sj = S[j] (movzx Sj, [S+j]): ~4 cycles (L1 latency)
   - Gated on step 3. Address unpredictable - no effective prefetch.
   
5. Store swap S[i]=Sj, S[j]=Si (mov [S+i], Sj; mov [S+j], Si): ~1 cycle
   - Into store buffer without stalling, but retirement still matters for correctness of later loads to the same slots.
   
6. Compute index (Si + Sj) (add idx, Si; add idx, Sj): ~1 cycle
   - Available as soon as Si and Sj are both in registers (right after losding S[i]).
   
7. Load z = S[Si + Sj] (movzx z, [S+idx]): ~4 cycles (L1 latency)
   - Third dependent load in the chain - gated on previous step.
   
8. XOR & store output (xor z, [input+t]; mov [output+t], z): ~1 cycle
   - Final step, gated on z.
 
TOTAL EXPECTED PRGA PERFORMANCE: ~5.5-9 CPB
Steps 3->4->6->7 form one unbroken chain of three dependent L1 loads (S[i]->S[j]->S[Si+Sj]) - if fully serial that's ~13+ cycles of latency
alone. The measured CPB comes in well below that sum because steps 1-2 of iteration t+1 (i-increment, S[i] load) are independent of iteration
t's j-chain and get issued early by out-of-order execution - so what you're timing is steady-state throughput across overlapping iterations,
not the latency of one isolated iteration.

Main bottleneck: the three-deep dependent-load chain
(S[i] -> S[j] -> S[Si+Sj]) each step, only prtially hidden by inter-iteration overlap since PRGA has one more dependent load than KSA
and thus less slack for overlap to hide.
*/
