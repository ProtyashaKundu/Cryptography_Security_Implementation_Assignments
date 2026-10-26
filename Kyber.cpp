// Build: g++ -O2 -std=c++17 Kyber.cpp -o Kyber

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <random>
#include <vector>
#include <x86intrin.h>   // __rdtsc()

const int N = 256, Q = 3329, K = 2, ETA = 2;

typedef std::vector<int> poly;     // N coefficients
typedef std::vector<poly> vec;     // K polynomials

std::mt19937 rng(1);

int modq(long long a)
{
    a %= Q;
    return a < 0 ? a + Q : a;
}

// Reference: a*b in Z_q[X]/(X^256 + 1), O(n^2). Only used to check the NTT.
poly mul_schoolbook(const poly& a, const poly& b)
{
    std::vector<long long> t(N, 0);
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            if (i + j < N) 
            	t[i + j]     += (long long)a[i] * b[j];
            else           
            	t[i + j - N] -= (long long)a[i] * b[j];
        }
    poly r(N);
    for (int i = 0; i < N; i++)
        r[i] = modq(t[i]);
    return r;
}


int zetas[128], gammas[128];

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

int power(long long b, int e)
{
    long long r = 1;
    for (b %= Q; e; e >>= 1, b = b * b % Q)
        if (e & 1) r = r * b % Q;
    return r;
}

int brv7(int x)                                 // reverse the 7 low bits
{
    int r = 0;
    for (int i = 0; i < 7; i++)
        r |= ((x >> i) & 1) << (6 - i);
    return r;
}

void init_ntt()
{
    for (int i = 0; i < 128; i++)
    {
        zetas[i]  = power(17, brv7(i));
        gammas[i] = power(17, 2 * brv7(i) + 1);
    }
}

void ntt(poly& f)                               
{
    int k = 1;
    for (int len = 128; len >= 2; len /= 2)
        for (int start = 0; start < N; start += 2 * len)
        {
            long long z = zetas[k++];
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
    int k = 127;
    for (int len = 2; len <= 128; len *= 2)
        for (int start = 0; start < N; start += 2 * len)
        {
            long long z = zetas[k--];
            for (int j = start; j < start + len; j++)
            {
                int t = f[j];
                f[j]       = addq(t, f[j + len]);
                f[j + len] = z * subq(f[j + len], t) % Q;
            }
        }
    for (int& x : f)
        x = (long long)x * 3303 % Q;            // 3303 = 128^-1 mod q
}

// a*b via NTT: transform, multiply pairs mod (X^2 - gamma), transform back
poly mul(const poly& a, const poly& b)
{
    poly A = a, B = b, c(N);
    ntt(A);
    ntt(B);
    for (int i = 0; i < 128; i++)
    {
        long long a0 = A[2 * i], a1 = A[2 * i + 1];
        long long b0 = B[2 * i], b1 = B[2 * i + 1];
        c[2 * i]     = (a0 * b0 + a1 * b1 % Q * gammas[i]) % Q;
        c[2 * i + 1] = (a0 * b1 + a1 * b0) % Q;
    }
    invntt(c);
    return c;
}

poly add(const poly& a, const poly& b)
{
    poly r(N);
    for (int i = 0; i < N; i++)
        r[i] = modq(a[i] + b[i]);
    return r;
}

// small noise in [-ETA, ETA] (centered binomial)
poly noise()
{
    poly r(N);
    for (int i = 0; i < N; i++)
    {
        int x = 0;
        for (int j = 0; j < ETA; j++)
            x += (rng() & 1) - (rng() & 1);
        r[i] = modq(x);
    }
    return r;
}

poly uniform()
{
    poly r(N);
    for (int i = 0; i < N; i++)
        r[i] = rng() % Q;
    return r;
}

// a . b = sum a[i]*b[i]
poly dot(const vec& a, const vec& b)
{
    poly r(N, 0);
    for (int i = 0; i < K; i++)
        r = add(r, mul(a[i], b[i]));
    return r;
}

// ---------------- Kyber ----------------
std::vector<vec> A(K, vec(K));     // public matrix
vec t(K), s(K);                    // public key t, secret key s

void keygen()
{
    for (int i = 0; i < K; i++)
        for (int j = 0; j < K; j++)
            A[i][j] = uniform();
    for (int i = 0; i < K; i++)
        s[i] = noise();
    for (int i = 0; i < K; i++)                 // t = A s + e
        t[i] = add(dot(A[i], s), noise());
}

// msg: 32 bytes = 256 bits -> one bit per coefficient
void encrypt(const uint8_t msg[32], vec& u, poly& v)
{
    vec r(K);
    for (int i = 0; i < K; i++)
        r[i] = noise();

    for (int i = 0; i < K; i++)                 // u = A^T r + e1
    {
        vec col(K);
        for (int j = 0; j < K; j++)
            col[j] = A[j][i];
        u[i] = add(dot(col, r), noise());
    }

    v = add(dot(t, r), noise());                // v = t.r + e2 + m*(q/2)
    for (int i = 0; i < N; i++)
        if ((msg[i / 8] >> (i % 8)) & 1)
            v[i] = modq(v[i] + Q / 2);
}

void decrypt(const vec& u, const poly& v, uint8_t msg[32])
{
    poly su = dot(s, u);
    for (int i = 0; i < 32; i++)
        msg[i] = 0;
    for (int i = 0; i < N; i++)
    {
        int w = modq(v[i] - su[i]);             // = m*(q/2) + small noise
        if (w > Q / 4 && w < 3 * Q / 4)         // closer to q/2 than to 0 -> bit 1
            msg[i / 8] |= 1 << (i % 8);
    }
}

// ---------------- Benchmark ----------------
// Runs f() many times, returns the median cycle count
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
    uint8_t msg[32] = "Hello Kyber! 32-byte message :)";
    uint8_t out[32];
    vec u(K);
    poly v;

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
    //std::cout << "one poly mul: schoolbook " << cycles([&] { mul_schoolbook(a, b); })
              //<< " cycles, NTT " << cycles([&] { mul(a, b); }) << " cycles\n\n";

    keygen();
    encrypt(msg, u, v);
    decrypt(u, v, out);
    std::cout << "Decrypted: " << out << "\n";
    std::cout << (std::equal(msg, msg + 32, out) ? "SUCCESS" : "FAILURE") << "\n\n";

    const int BYTES = 32;  // message size
    uint64_t ck = cycles([&] { keygen(); });
    encrypt(msg, u, v);
    uint64_t ce = cycles([&] { encrypt(msg, u, v); });
    uint64_t cd = cycles([&] { decrypt(u, v, out); });

    std::cout << "            cycles      cpb (per 32-byte msg)\n";
    std::cout << "keygen : " << ck << "\t" << (double)ck / BYTES << "\n";
    std::cout << "encrypt: " << ce << "\t" << (double)ce / BYTES << "\n";
    std::cout << "decrypt: " << cd << "\t" << (double)cd / BYTES << "\n";
}
/**
protyasha@protyashaLinux:~/Documents$ ./Kyber
NTT mul == schoolbook mul: YES
one poly mul: schoolbook 229505 cycles, NTT 23550 cycles

Decrypted: Hello Kyber! 32-byte message :)
SUCCESS

            cycles      cpb (per 32-byte msg)
keygen : 172437	5388.66
encrypt: 206617	6456.78
decrypt: 46345	1448.28
protyasha@protyashaLinux:~/Documents$ ./Kyber
NTT mul == schoolbook mul: YES
one poly mul: schoolbook 262084 cycles, NTT 28734 cycles

Decrypted: Hello Kyber! 32-byte message :)
SUCCESS

            cycles      cpb (per 32-byte msg)
keygen : 173483	5421.34
encrypt: 218447	6826.47
decrypt: 49042	1532.56
*/
