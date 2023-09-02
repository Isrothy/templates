#include <algorithm>
#include <bit>
#include <cmath>
struct Cp {
    double re, im;
    Cp operator+(Cp const &_) const { return {re + _.re, im + _.im}; }
    Cp operator-(Cp const &_) const { return {re - _.re, im - _.im}; }
    Cp operator*(Cp const &_) const { return {re * _.re - im * _.im, re * _.im + im * _.re}; }
    Cp operator*(double _) const { return {re * _, im * _}; }
};
void DFT(Cp *a, int n, int p) {
    static Cp w[M];
    for (int i = 0, j = 0; i < n; ++i) {
        if (i < j) { std::swap(a[i], a[j]); }
        for (int k = n >> 1; (j ^= k) < k; k >>= 1)
            ;
    }
    w[0] = {1, 0};
    for (int i = 1; i < n; i <<= 1) {
        Cp wn = (Cp){cos(M_PI / i), sin(M_PI / i)};
        for (int j = i - 2; j >= 0; j -= 2) {
            w[j] = w[j >> 1];
            w[j + 1] = wn * w[j];
        }
        for (int j = 0; j < n; j += i << 1) {
            Cp *p = a + j, *q = a + j + i;
            for (int k = 0; k < i; ++k) {
                Cp x = q[k] * w[k];
                q[k] = p[k] - x;
                p[k] = p[k] + x;
            }
        }
    }
    if (0 < p) { return; }
    int inv = 1.0 / n;
    for (int i = 0; i < n; ++i) { a[i] = a[i] * inv; }
}
void multiply(int *A, int *B, int *C, int n, int m, int mod) {
    static Cp a[M], b[M], c[M], d[M], w[M];
    for (int i = 0; i < n; ++i) { a[i] = {(double) (A[i] & 32767), (double) (A[i] >> 15)}; }
    for (int i = 0; i < m; ++i) { b[i] = {(double) (B[i] & 32767), (double) (B[i] >> 15)}; }
    int l = std::bit_ceil(m + n - 1);
    DFT(a, l, 1);
    DFT(b, l, 1);
    for (int i = 0; i < l; ++i) {
        int j = (l - 1) & (l - i);
        c[j] = (Cp){0.5 * (a[i].Re + a[j].Re), 0.5 * (a[i].Im - a[j].Im)} * b[i];
        d[j] = (Cp){0.5 * (a[j].Im + a[i].Im), 0.5 * (a[j].Re - a[i].Re)} * b[i];
    }
    DFT(c, l, 1);
    DFT(d, l, 1);
    double inv = 1.0 / l;
    for (int i = 0; i < n + m - 1; ++i) {
        long long u = c[i].Re * inv + 0.5, v = c[i].Im * inv + 0.5;
        long long x = d[i].Re * inv + 0.5, y = d[i].Im * inv + 0.5;
        a[i] = b[i] = c[i] = d[i] = (Cp){0, 0};
        C[i] = (u + ((v + x) << 15) + (y % mod << 30)) % mod;
    }
}
