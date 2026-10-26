// Build: g++ -O2 -std=c++17 Dilithium.cpp -o Dilithium

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <random>
#include <vector>
#include <x86intrin.h>   // __rdtsc()

// Dilithium2 parameters
const int N = 256, Q = 8380417, K = 4, L = 4;
const int ETA = 2, TAU = 39, BETA = TAU * ETA;   // BETA = 78
const int G1 = 1 << 17;                          // gamma1: range of mask y
const int G2 = (Q - 1) / 88;                     // gamma2: rounding range

typedef std::vector<int> poly;
typedef std::vector<poly> vec;

std::mt19937 rng(1);

int modq(long long a)
{
    a %= Q;
    return a < 0 ? a + Q : a;
}

int centered(int a)                 // [0,q) -> (-q/2, q/2]
{
    return a > Q / 2 ? a - Q : a;
}

// ---------------- polynomial arithmetic ----------------
// Reference: a*b in Z_q[X]/(X^256 + 1), O(n^2). Only used to check the NTT.
poly mul_schoolbook(const poly& a, const poly& b)
{
    std::vector<long long> t(N, 0);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            if (i + j < N) t[i + j]     += (long long)a[i] * b[j];
            else           t[i + j - N] -= (long long)a[i] * b[j];
        }
    poly r(N);
    for (int i = 0; i < N; i++)
        r[i] = modq(t[i]);
    return r;
}

// ---------------- NTT (Dilithium, FIPS 204) ----------------
// q = 8380417 has a 512th root of unity (1753), so X^256+1 splits all the way
// into 256 linear factors: after the NTT, multiplication is just pointwise.
int zetas[256];
int inv256;                                     // 256^-1 mod q

int power(long long b, int e)
{
    long long r = 1;
    for (b %= Q; e; e >>= 1, b = b * b % Q)
        if (e & 1) r = r * b % Q;
    return r;
}

int brv8(int x)                                 // reverse the 8 low bits
{
    int r = 0;
    for (int i = 0; i < 8; i++)
        r |= ((x >> i) & 1) << (7 - i);
    return r;
}

void init_ntt()
{
    for (int i = 0; i < 256; i++)
        zetas[i] = power(1753, brv8(i));
    inv256 = power(256, Q - 2);
}

// inputs in [0,q): one conditional correction instead of a slow %
int addq(int a, int b) 
{ 
	int r = a + b; 
	return r >= Q ? r - Q : r; 
}
int subq(int a, int b) 
{ 
	int r = a - b; 
	return r < 0 ? r + Q : r; 
}

void ntt(poly& f)                               
{
    int k = 0;
    for (int len = 128; len >= 1; len /= 2)
        for (int start = 0; start < N; start += 2 * len)
        {
            long long z = zetas[++k];
            for (int j = start; j < start + len; j++)
            {
                int t = z * f[j + len] % Q;
                f[j + len] = subq(f[j], t);
                f[j]       = addq(f[j], t);
            }
        }
}

void invntt(poly& f)                            
{
    int k = 256;
    for (int len = 1; len <= 128; len *= 2)
        for (int start = 0; start < N; start += 2 * len)
        {
            long long z = Q - zetas[--k];       // inverse uses -zeta
            for (int j = start; j < start + len; j++)
            {
                int t = f[j];
                f[j]       = addq(t, f[j + len]);
                f[j + len] = z * subq(t, f[j + len]) % Q;
            }
        }
    for (int& x : f)
        x = (long long)x * inv256 % Q;
}

// a*b via NTT: transform, multiply pointwise, transform back
poly mul(const poly& a, const poly& b)
{
    poly A = a, B = b;
    ntt(A);
    ntt(B);
    for (int i = 0; i < N; i++)
        A[i] = (long long)A[i] * B[i] % Q;
    invntt(A);
    return A;
}

poly add(const poly& a, const poly& b)
{
    poly r(N);
    for (int i = 0; i < N; i++)
        r[i] = modq(a[i] + b[i]);
    return r;
}

poly sub(const poly& a, const poly& b)
{
    poly r(N);
    for (int i = 0; i < N; i++)
        r[i] = modq(a[i] - b[i]);
    return r;
}

int norm(const poly& p)                         // infinity norm
{
    int m = 0;
    for (int x : p)
        m = std::max(m, std::abs(centered(x)));
    return m;
}

// ---------------- sampling ----------------
poly uniform()
{
    poly r(N);
    for (int i = 0; i < N; i++)
        r[i] = rng() % Q;
    return r;
}

poly small()                                    // coefficients in [-ETA, ETA]
{
    poly r(N);
    for (int i = 0; i < N; i++)
        r[i] = modq((int)(rng() % (2 * ETA + 1)) - ETA);
    return r;
}

poly mask()                                     // coefficients in (-G1, G1]
{
    poly r(N);
    for (int i = 0; i < N; i++)
        r[i] = modq((int)(rng() % (2 * G1)) - (G1 - 1));
    return r;
}

// ---------------- rounding: r = r1 * 2*G2 + r0 ----------------
void decompose(int r, int& r1, int& r0)
{
    r0 = r % (2 * G2);
    if (r0 > G2) r0 -= 2 * G2;                  // centered low part
    if (r - r0 == Q - 1) { r1 = 0; r0 -= 1; }  // wrap-around corner case
    else r1 = (r - r0) / (2 * G2);
}

poly highbits(const poly& p)
{
    poly h(N);
    int r0;
    for (int i = 0; i < N; i++)
        decompose(p[i], h[i], r0);
    return h;
}

int lowbits_norm(const poly& p)                 // max |r0|
{
    int m = 0, r1, r0;
    for (int i = 0; i < N; i++)
    {
        decompose(p[i], r1, r0);
        m = std::max(m, std::abs(r0));
    }
    return m;
}

// ---------------- hash -> challenge ----------------
// H(msg, w1) -> 64-bit seed (toy hash; real Dilithium uses SHAKE256)
uint64_t H(const uint8_t* msg, int len, const vec& w1)
{
    uint64_t h = 14695981039346656037ULL;
    auto eat = [&](uint8_t b) { h = (h ^ b) * 1099511628211ULL; };
    for (int i = 0; i < len; i++) eat(msg[i]);
    for (auto& p : w1)
        for (int x : p) { eat(x & 0xFF); eat(x >> 8); }
    return h;
}

// seed -> polynomial with exactly TAU coefficients = +-1, rest 0
poly challenge(uint64_t seed)
{
    std::mt19937_64 g(seed);
    poly c(N, 0);
    for (int i = N - TAU; i < N; i++)
    {
        int j = g() % (i + 1);
        c[i] = c[j];
        c[j] = (g() & 1) ? 1 : Q - 1;
    }
    return c;
}

// ---------------- Dilithium ----------------
std::vector<vec> A(K, vec(L));                  // public matrix
vec t(K);                                       // public key
vec s1(L), s2(K);                               // secret key
int sign_iterations = 0;

vec matvec(const vec& y)                        // A * y
{
    vec w(K, poly(N, 0));
    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
            w[i] = add(w[i], mul(A[i][j], y[j]));
    return w;
}

void keygen()
{
    for (int i = 0; i < K; i++)
        for (int j = 0; j < L; j++)
            A[i][j] = uniform();
    for (int j = 0; j < L; j++) s1[j] = small();
    for (int i = 0; i < K; i++) s2[i] = small();
    vec As1 = matvec(s1);
    for (int i = 0; i < K; i++)                 // t = A s1 + s2
        t[i] = add(As1[i], s2[i]);
}

// signature = (z, h) where h is the hash seed of the challenge c
void sign(const uint8_t* msg, int len, vec& z, uint64_t& h)
{
    while (true)
    {
        sign_iterations++;

        vec y(L);
        for (int j = 0; j < L; j++) y[j] = mask();

        vec w = matvec(y);                      // w = A y
        vec w1(K);
        for (int i = 0; i < K; i++) w1[i] = highbits(w[i]);

        h = H(msg, len, w1);
        poly c = challenge(h);

        bool ok = true;
        for (int j = 0; j < L && ok; j++)       // z = y + c s1, must not leak s1
        {
            z[j] = add(y[j], mul(c, s1[j]));
            if (norm(z[j]) >= G1 - BETA) ok = false;
        }
        for (int i = 0; i < K && ok; i++)       // w - c s2 must keep the same high bits
            if (lowbits_norm(sub(w[i], mul(c, s2[i]))) >= G2 - BETA) ok = false;

        if (ok) return;                         // otherwise reject and retry
    }
}

bool verify(const uint8_t* msg, int len, const vec& z, uint64_t h)
{
    for (int j = 0; j < L; j++)
        if (norm(z[j]) >= G1 - BETA) return false;

    poly c = challenge(h);
    vec Az = matvec(z);
    vec w1(K);
    for (int i = 0; i < K; i++)                 // A z - c t = A y - c s2
        w1[i] = highbits(sub(Az[i], mul(c, t[i])));

    return H(msg, len, w1) == h;
}

// ---------------- Benchmark ----------------
template <class F>
uint64_t cycles(F f, int runs = 101)            
{
    std::vector<uint64_t> c(runs);
    for (int i = 0; i < runs; i++)
    {
        uint64_t start = __rdtsc();
        f();
        c[i] = __rdtsc() - start;
    }
    std::sort(c.begin(), c.end());
    return c[runs / 2];
}

int main()
{
    uint8_t msg[32] = "Hello Dilithium! 32-byte msg :)";
    const int LEN = 32;
    vec z(L);
    uint64_t h;

    init_ntt();

    // check NTT multiplication against schoolbook on random inputs
    bool same = true;
    for (int i = 0; i < 100; i++)
    {
        poly a = uniform(), b = uniform();
        if (mul(a, b) != mul_schoolbook(a, b)) same = false;
    }
    std::cout << "NTT mul == schoolbook mul: " << (same ? "YES" : "NO") << "\n";

    poly a = uniform(), b = uniform();
    std::cout << "one poly mul: schoolbook " << cycles([&] { mul_schoolbook(a, b); })
              << " cycles, NTT " << cycles([&] { mul(a, b); }) << " cycles\n\n";

    keygen();
    sign(msg, LEN, z, h);
    std::cout << "Valid signature   : " << (verify(msg, LEN, z, h) ? "ACCEPT" : "REJECT") << "\n";
    msg[0] ^= 1;
    std::cout << "Tampered message  : " << (verify(msg, LEN, z, h) ? "ACCEPT" : "REJECT") << "\n\n";
    msg[0] ^= 1;

    uint64_t ck = cycles([&] { keygen(); });
    sign_iterations = 0;
    uint64_t cs = cycles([&] { sign(msg, LEN, z, h); });
    double avg_iter = sign_iterations / 101.0;
    uint64_t cv = cycles([&] { verify(msg, LEN, z, h); });

    std::cout << "            cycles      cpb (per 32-byte msg)\n";
    std::cout << "keygen : " << ck << "\t" << (double)ck / LEN << "\n";
    std::cout << "sign   : " << cs << "\t" << (double)cs / LEN << "\n";
    std::cout << "verify : " << cv << "\t" << (double)cv / LEN << "\n";
    std::cout << "\navg sign attempts (rejection sampling): " << avg_iter << "\n";
}
/**
protyasha@protyashaLinux:~/Documents$ ./Dilithium
NTT mul == schoolbook mul: YES
one poly mul: schoolbook 229562 cycles, NTT 25175 cycles

Valid signature   : ACCEPT
Tampered message  : REJECT

            cycles      cpb (per 32-byte msg)
keygen : 502820	15713.1
sign   : 1332577	41643
verify : 396094	12377.9

avg sign attempts (rejection sampling): 3.92079
protyasha@protyashaLinux:~/Documents$ ./Dilithium
NTT mul == schoolbook mul: YES
one poly mul: schoolbook 282406 cycles, NTT 30409 cycles

Valid signature   : ACCEPT
Tampered message  : REJECT

            cycles      cpb (per 32-byte msg)
keygen : 502678	15708.7
sign   : 1378034	43063.6
verify : 424011	13250.3

avg sign attempts (rejection sampling): 3.92079
*/
