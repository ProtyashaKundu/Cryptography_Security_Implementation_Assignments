/*
 * Simple RSA-2048 encrypt/decrypt round trip using Intel IPP Cryptography.
 *   gcc RSA.c -o RSA -lippcp 
 */
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>
#include "ippcp.h"

/* ---- cycle counting -------------------------------------------------*/

static inline unsigned long long start_tsc(void) 
{
    _mm_lfence();//It prevents the CPU out-of-order execution engine from speculatively executing instructions across the timing boundary, ensuring precise timing measurements.
    return __rdtsc();
}
static inline unsigned long long stop_tsc(void) 
{
    unsigned int aux;
    unsigned long long tsc = __rdtscp(&aux);
    _mm_lfence();
    return tsc;
}

/* ---- helpers ------------------------------------------------------- */

static void check(IppStatus st, const char* what)//Validates every IPP function return status (IppStatus). If an API call fails, it prints the error string via ippcpGetStatusString() and exits.
{
    if (st != ippStsNoErr) 
    {
        fprintf(stderr, "IPP error in %s: %d (%s)\n", what, (int)st, ippcpGetStatusString(st));
        exit(EXIT_FAILURE);
    }
}

static Ipp32u rdrand32(void) //Wraps ippsPRNGenRDRAND to obtain hardware-generated 32-bit random integers directly from the x86 RDRAND hardware instruction.
{
    Ipp32u v;
    check(ippsPRNGenRDRAND(&v, 32, NULL), "ippsPRNGenRDRAND");
    return v;
}

/* Fill buf[0..len) with random bytes */
static void fill_nonzero_random(Ipp8u* buf, int len) 
{
    int nWords = (len + 3) / 4;
    Ipp32u* words = (Ipp32u*)malloc((size_t)nWords * sizeof(Ipp32u));
    if (!words) 
    { 
        fprintf(stderr, "malloc failed\n"); exit(EXIT_FAILURE); 
    }
    for (int i = 0; i < nWords; i++) 
        words[i] = rdrand32();
    memcpy(buf, words, (size_t)len);
    free(words);
}

int main(void) 
{
    const int rsaBits    = 2048;
    const int rsaBytes   = rsaBits / 8;   /* 256 */
    const int factorBits = rsaBits / 2;   /* 1024, size of p and q */
    const int pubExpBits = 32;            /* room for e = 65537 */
    const int bnWordsN   = (rsaBits + 31) / 32;
    const int crtWords   = (5 * factorBits + 31) / 32; /* 160 words required for 5 CRT components */

    IppStatus st;

    /* ---- public key context ---- */ 
    //Intel IPP uses user-managed memory allocation. First, ippsRSA_GetSize* is called to determine the exact byte size needed. Then standard malloc() allocates that memory, and ippsRSA_Init* formats that memory buffer into an active state object.
    int pubKeySize = 0;
    check(ippsRSA_GetSizePublicKey(rsaBits, pubExpBits, &pubKeySize),"ippsRSA_GetSizePublicKey");
    IppsRSAPublicKeyState* pPubKey = (IppsRSAPublicKeyState*)malloc((size_t)pubKeySize);
    check(ippsRSA_InitPublicKey(rsaBits, pubExpBits, pPubKey, pubKeySize),"ippsRSA_InitPublicKey");

    /* ---- private key (CRT / "type 2") context ---- */
    int privKeySize = 0;
    check(ippsRSA_GetSizePrivateKeyType2(factorBits, factorBits, &privKeySize), "ippsRSA_GetSizePrivateKeyType2");
    IppsRSAPrivateKeyState* pPrivKey = (IppsRSAPrivateKeyState*)malloc((size_t)privKeySize);
    check(ippsRSA_InitPrivateKeyType2(factorBits, factorBits, pPrivKey, privKeySize),"ippsRSA_InitPrivateKeyType2");

    /* ---- big numbers used during key generation ---- */
    int bnSize1 = 0, bnSizeN = 0, bnSizeCRT = 0;
    check(ippsBigNumGetSize(1, &bnSize1), "ippsBigNumGetSize(1)");
    check(ippsBigNumGetSize(bnWordsN, &bnSizeN), "ippsBigNumGetSize(N)");
    check(ippsBigNumGetSize(crtWords, &bnSizeCRT), "ippsBigNumGetSize(CRT)");

    IppsBigNumState* pSrcE       = (IppsBigNumState*)malloc((size_t)bnSize1);
    IppsBigNumState* pModulus   = (IppsBigNumState*)malloc((size_t)bnSizeN);//N
    IppsBigNumState* pPublicExp  = (IppsBigNumState*)malloc((size_t)bnSize1);//e
    IppsBigNumState* pPrivateExp = (IppsBigNumState*)malloc((size_t)bnSizeCRT);//d (sized for CRT components)
    
    if (!pSrcE || !pModulus || !pPublicExp || !pPrivateExp) 
    {
        fprintf(stderr, "Memory allocation failed for BigNum structures.\n");
        exit(EXIT_FAILURE);
    }

    check(ippsBigNumInit(1, pSrcE), "ippsBigNumInit(pSrcE)");
    check(ippsBigNumInit(bnWordsN, pModulus), "ippsBigNumInit(pModulus)");
    check(ippsBigNumInit(1, pPublicExp), "ippsBigNumInit(pPublicExp)");
    check(ippsBigNumInit(crtWords, pPrivateExp), "ippsBigNumInit(pPrivateExp)");

    Ipp32u eSeed = 65537;
    check(ippsSet_BN(ippBigNumPOS, 1, &eSeed, pSrcE), "ippsSet_BN(e)");

    /* ---- prime generator context (for the p, q search) ---- */
    int primeGenSize = 0;
    check(ippsPrimeGetSize(factorBits, &primeGenSize), "ippsPrimeGetSize");//Queries IPP for the buffer size (in bytes) needed to hold an internal IppsPrimeState context configured for generating primes up to factorBits long (here, 1024 bits each for p and q).
    IppsPrimeState* pPrimeGen = (IppsPrimeState*)malloc((size_t)primeGenSize);
    check(ippsPrimeInit(factorBits, pPrimeGen), "ippsPrimeInit");//Prepares the candidate search generator for 1024-bit primes p and q.

    /* ---- generate the keypair ---- */
    int genBufSize = 0;
    check(ippsRSA_GetBufferSizePrivateKey(&genBufSize, pPrivKey),"ippsRSA_GetBufferSizePrivateKey (keygen)");
    Ipp8u* pGenBuffer = (Ipp8u*)malloc((size_t)genBufSize);

    const int nTrials = 50; /* Miller-Rabin rounds for primality testing */
    // Executes prime generation, primality testing (using 50 Miller-Rabin test rounds), modulus multiplication (N = p x q), and CRT value computation using the ippsPRNGenRDRAND hardware RNG.
    check(ippsRSA_GenerateKeys(pSrcE, pModulus, pPublicExp, pPrivateExp, pPrivKey, pGenBuffer, nTrials, pPrimeGen, ippsPRNGenRDRAND, NULL), "ippsRSA_GenerateKeys");

    check(ippsRSA_SetPublicKey(pModulus, pPublicExp, pPubKey),"ippsRSA_SetPublicKey");

    /* ---- scratch buffers for encrypt/decrypt ---- */
    int encBufSize = 0, decBufSize = 0;
    check(ippsRSA_GetBufferSizePublicKey(&encBufSize, pPubKey),"ippsRSA_GetBufferSizePublicKey");
    check(ippsRSA_GetBufferSizePrivateKey(&decBufSize, pPrivKey),"ippsRSA_GetBufferSizePrivateKey (decrypt)");
    Ipp8u* pEncBuffer = (Ipp8u*)malloc((size_t)encBufSize);
    Ipp8u* pDecBuffer = (Ipp8u*)malloc((size_t)decBufSize);

    /* ---- message ---- */
    const char* msg="All my years to this moment. All my roads to this wall. All my words to this silence. All my pride to this fall.";
    int msgLen = (int)strlen(msg);
    Ipp8u pPlaintext[256];
    memset(pPlaintext, 0, sizeof(pPlaintext));
    memcpy(pPlaintext, msg, (size_t)msgLen);

    /* RSA-OAEP configurations (SHA-256) */
    const IppsHashMethod* pHash = ippsHashMethod_SHA256();
    const Ipp8u* pLabel = NULL;
    int labLen = 0;

    /* OAEP seed size matches SHA-256 digest size (32 bytes) */
	//Guarantees ciphertext integrity and prevents Bleichenbacher padding oracle attacks. Returning a generic error ensures an attacker cannot learn which specific check failed.
    Ipp8u pSeed[32];
    fill_nonzero_random(pSeed, sizeof(pSeed));

    Ipp8u pCiphertext[256];
    memset(pCiphertext, 0, sizeof(pCiphertext));
    check(ippsRSAEncrypt_OAEP_rmf(pPlaintext, msgLen, pLabel, labLen, pSeed, pCiphertext, pPubKey, pHash, pEncBuffer), "ippsRSAEncrypt_OAEP_rmf");

    Ipp8u pDecrypted[256];
    memset(pDecrypted, 0, sizeof(pDecrypted));
    int decryptedLen = 0;
    check(ippsRSADecrypt_OAEP_rmf(pCiphertext, pLabel, labLen, pDecrypted, &decryptedLen, pPrivKey, pHash, pDecBuffer), "ippsRSADecrypt_OAEP_rmf");

    /* ---- verify ---- */
    printf("Original  (%d bytes): %s\n", msgLen, msg);
    printf("Decrypted (%d bytes): %.*s\n", decryptedLen, decryptedLen, pDecrypted);
    if (decryptedLen == msgLen && memcmp(pPlaintext, pDecrypted, (size_t)msgLen) == 0) 
    {
        printf("Round trip OK.\n");
    } 
    else 
    {
        printf("Round trip FAILED.\n");
        free(pPubKey); free(pPrivKey); free(pSrcE); free(pModulus);
        free(pPublicExp); free(pPrivateExp); free(pPrimeGen);
        free(pGenBuffer); free(pEncBuffer); free(pDecBuffer);
        return EXIT_FAILURE;
    }

    /* ---- cycles-per-byte benchmark ---- */
    const int iterations = 10000;

    /* warm-up (cache, branch predictor, microcode) */
    for (int i = 0; i < 100; i++) 
    {
        fill_nonzero_random(pSeed, sizeof(pSeed));
        check(ippsRSAEncrypt_OAEP_rmf(pPlaintext, msgLen, pLabel, labLen, pSeed, pCiphertext, pPubKey, pHash, pEncBuffer), "warm-up encrypt");
    }

    unsigned long long encCycles = 0;
    for (int i = 0; i < iterations; i++) 
    {
        fill_nonzero_random(pSeed, sizeof(pSeed));
        unsigned long long t1 = start_tsc();
        check(ippsRSAEncrypt_OAEP_rmf(pPlaintext, msgLen, pLabel, labLen, pSeed, pCiphertext, pPubKey, pHash, pEncBuffer), "encrypt (timed)");
        unsigned long long t2 = stop_tsc();
        encCycles += (t2 - t1);
    }

    for (int i = 0; i < 100; i++) 
    {
        check(ippsRSADecrypt_OAEP_rmf(pCiphertext, pLabel, labLen, pDecrypted, &decryptedLen, pPrivKey, pHash, pDecBuffer), "warm-up decrypt");
    }

    unsigned long long decCycles = 0;
    for (int i = 0; i < iterations; i++) 
    {
        unsigned long long t1 = start_tsc();
        check(ippsRSADecrypt_OAEP_rmf(pCiphertext, pLabel, labLen, pDecrypted, &decryptedLen, pPrivKey, pHash, pDecBuffer), "decrypt (timed)");
        unsigned long long t2 = stop_tsc();
        decCycles += (t2 - t1);
    }

    double avgEnc = (double)encCycles / iterations;
    double avgDec = (double)decCycles / iterations;

    printf("===================================================\n");
    printf("RSA block size:                %d bytes\n", rsaBytes);
    printf("Encrypt: avg cycles            %.2f\n", avgEnc);
    printf("Encrypt: cycles/byte (block)   %.2f CPB\n", avgEnc / rsaBytes);
    printf("Decrypt: avg cycles            %.2f\n", avgDec);
    printf("Decrypt: cycles/byte (block)   %.2f CPB\n", avgDec / rsaBytes);
    printf("===================================================\n");
    /* Note: public-key ops (encrypt) touch the whole rsaBytes-sized block
     * regardless of msgLen, so CPB relative to msgLen isn't meaningful --
     * cycles/block is the metric that reflects the actual modexp cost. */

    /* ---- cleanup: only free what was malloc'd ---- */
    free(pPubKey);
    free(pPrivKey);
    free(pSrcE);
    free(pModulus);
    free(pPublicExp);
    free(pPrivateExp);
    free(pPrimeGen);
    free(pGenBuffer);
    free(pEncBuffer);
    free(pDecBuffer);

    return 0;
}
/**output
protyasha@protyashaLinux:~/Documents$ ./RSA
Original  (51 bytes): Success is 100% 0%, 0% 100%, and the Fourier Series
Decrypted (51 bytes): Success is 100% 0%, 0% 100%, and the Fourier Series
Round trip OK.
===================================================
RSA block size:                256 bytes
Encrypt: avg cycles            57274.72
Encrypt: cycles/byte (block)   223.73 CPB
Decrypt: avg cycles            1894023.49
Decrypt: cycles/byte (block)   7398.53 CPB
===================================================
protyasha@protyashaLinux:~/Documents$ gcc RSA.c -o RSA -lippcp
protyasha@protyashaLinux:~/Documents$ ./RSA
Original  (112 bytes): All my years to this moment. All my roads to this wall. All my words to this silence. All my pride to this fall.
Decrypted (112 bytes): All my years to this moment. All my roads to this wall. All my words to this silence. All my pride to this fall.
Round trip OK.
===================================================
RSA block size:                256 bytes
Encrypt: avg cycles            57428.54
Encrypt: cycles/byte (block)   224.33 CPB
Decrypt: avg cycles            1902388.30
Decrypt: cycles/byte (block)   7431.20 CPB
===================================================

**/

/**reference
https://www.intel.com/content/www/us/en/docs/ipp-crypto/developer-guide-reference/2021-9/rsa-algorithm-functions.html*/
