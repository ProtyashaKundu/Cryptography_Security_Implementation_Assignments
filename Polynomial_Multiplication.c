#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <time.h>
#include <stdint.h>

typedef double complex cplx;

// ---- Cycle-accurate timer (x86/x86-64 only) ----
static inline uint64_t rdtsc(void) 
{
#if defined(__x86_64__) || defined(__i386__)
    unsigned int lo, hi;
    __asm__ __volatile__ ("cpuid\n\t"      
                          "rdtsc" : "=a"(lo), "=d"(hi)
                          :: "%rbx", "%rcx");
    return ((uint64_t)hi << 32) | lo;
#else
    return 0; 
#endif
}

// Forward FFT: coefficient representation -> value representation
void fft(cplx *P, cplx *result, int n)
{
    if (n == 1)
    {
        result[0] = P[0];
        return;
    }
    int half = n / 2;

    cplx *Pe = malloc(half * sizeof(cplx));
    cplx *Po = malloc(half * sizeof(cplx));
    for (int i = 0; i < half; i++)
    {
        Pe[i] = P[2 * i];
        Po[i] = P[2 * i + 1];
    }

    cplx *ye = malloc(half * sizeof(cplx));
    cplx *yo = malloc(half * sizeof(cplx));
    fft(Pe, ye, half);
    fft(Po, yo, half);

    cplx w = 1;
    cplx wn = cexp(2.0 * I * M_PI / n);   // forward: e^{+2*pi*i/n}
    for (int k = 0; k < half; k++)
    {
        cplx t = w * yo[k];
        result[k] = ye[k] + t;
        result[k + half] = ye[k] - t;
        w *= wn;
    }

    free(Pe);
    free(Po);
    free(ye);
    free(yo);
}

// Inverse FFT: value representation -> coefficient representation
/**void ifft_raw(cplx *P, cplx *result, int n)
{
    if (n == 1)
    {
        result[0] = P[0];
        return;
    }
    int half = n / 2;

    cplx *Pe = malloc(half * sizeof(cplx));
    cplx *Po = malloc(half * sizeof(cplx));
    for (int i = 0; i < half; i++)
    {
        Pe[i] = P[2 * i];
        Po[i] = P[2 * i + 1];
    }

    cplx *ye = malloc(half * sizeof(cplx));
    cplx *yo = malloc(half * sizeof(cplx));
    ifft_raw(Pe, ye, half);
    ifft_raw(Po, yo, half);

    cplx w = 1;
    cplx wn = (1.0 / n) * cexp(-2.0 * I * M_PI / n); 
    for (int k = 0; k < half; k++)
    {
        cplx t = w * yo[k];
        result[k] = ye[k] + t;
        result[k + half] = ye[k] - t;
        w *= wn;
    }

    free(Pe);
    free(Po);
    free(ye);
    free(yo);
}
void ifft(cplx *P, cplx *result, int n)
{
    ifft_raw(P, result, n);
    for (int i = 0; i < n; i++)
        result[i] /= n;
}*/

void ifft(cplx *P, cplx *result, int n)
{
    cplx *tmp_in = malloc(n * sizeof(cplx));
    for (int i = 0; i < n; i++)
        tmp_in[i] = conj(P[i]);

    fft(tmp_in, result, n);   
    
    for (int i = 0; i < n; i++)
        result[i] = conj(result[i]) / n;

    free(tmp_in);
}

void multiply(int *a, int n1, int *b, int n2, long long *res) 
{
    int result_size = n1 + n2 - 1;
    int n = 1;
    while (n < result_size) 
    	n <<= 1;

    cplx *A = calloc(n, sizeof(cplx));
    cplx *B = calloc(n, sizeof(cplx));
    for (int i = 0; i < n1; i++) 
    	A[i] = a[i];
    for (int i = 0; i < n2; i++) 
    	B[i] = b[i];

    cplx *FA = malloc(n * sizeof(cplx));
    cplx *FB = malloc(n * sizeof(cplx));

    uint64_t t0, t1;

    t0 = rdtsc();
    fft(A, FA, n);
    t1 = rdtsc();
    printf("fft(A):   %llu cycles (n=%d)\n", (unsigned long long)(t1 - t0), n);

    t0 = rdtsc();
    fft(B, FB, n);
    t1 = rdtsc();
    printf("fft(B):   %llu cycles (n=%d)\n", (unsigned long long)(t1 - t0), n);

    cplx *FC = malloc(n * sizeof(cplx));

    t0 = rdtsc();
    for (int i = 0; i < n; i++)
        FC[i] = FA[i] * FB[i];
    t1 = rdtsc();
    printf("pointwise mul: %llu cycles\n", (unsigned long long)(t1 - t0));

    cplx *C = malloc(n * sizeof(cplx));

    t0 = rdtsc();
    ifft(FC, C, n);
    t1 = rdtsc();
    printf("ifft(FC): %llu cycles (n=%d)\n", (unsigned long long)(t1 - t0), n);

    for (int i = 0; i < result_size; i++)
        res[i] = (long long) round(creal(C[i]));

    free(A); 
    free(B); 
    free(FA); 
    free(FB); 
    free(FC); 
    free(C);
}

void multiply_naive(int *a, int n1, int *b, int n2, long long *res) 
{
    int result_size = n1 + n2 - 1;

    // Initialize result coefficients to 0
    for (int k = 0; k < result_size; k++)
        res[k] = 0;

    // c_k = sum over i+j=k of a[i]*b[j]
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            res[i + j] += (long long) a[i] * b[j];
        }
    }
}

int main(void) 
{
    int a[] = {1, 2, 3};
    int b[] = {4, 5};
    int n1 = 3, n2 = 2;

    long long res[10];
    uint64_t total_start_naive= rdtsc();
    clock_t c_start_naive = clock();
    multiply_naive(a, n1, b, n2, res);
    uint64_t total_end_naive = rdtsc();
    clock_t c_end_naive = clock();
    
    uint64_t total_start = rdtsc();
    clock_t c_start = clock();
    multiply(a, n1, b, n2, res);
    uint64_t total_end = rdtsc();
    clock_t c_end = clock();

    printf("\nResult: ");
    for (int i = 0; i < n1 + n2 - 1; i++)
        printf("%lld ", res[i]);
    printf("\n");
    
   
    printf("\nTotal cycles with FFT:  %llu\n",
           (unsigned long long)(total_end - total_start));
    printf("\nTotal cycles naive method:  %llu\n",
           (unsigned long long)(total_end_naive - total_start_naive));

   /** printf("Total time (clock()):  %.6f ms\n",
           1000.0 * (c_end - c_start) / CLOCKS_PER_SEC);*/

    return 0;
}
/**protyasha@protyashaLinux:~/Documents$ ./PM
fft(A):   55669 cycles (n=4)
fft(B):   4871 cycles (n=4)
pointwise mul: 910 cycles
ifft(FC): 3376 cycles (n=4)

Result: 4 13 22 15 

Total cycles with FFT:  235327

Total cycles naive method:  30125

*/
