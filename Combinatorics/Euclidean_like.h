#include <tuple>
std::tuple<int64_t, int64_t, int64_t>
Euclidean_like(int64_t n, int64_t a, int64_t b, int64_t c, int64_t mod) {
    // sum_{x=0}^n (ax + b)/c
    int64_t x = a / c % mod, y = b / c % mod;
    int64_t s0 = (n + 1) % mod;
    int64_t s1 = n * (n + 1) % mod * inv[2] % mod;
    int64_t s2 = n * (n + 1) % mod * (2 * n + 1) % mod * inv[6] % mod;
    int64_t m = (a * n + b) / c;
    int64_t _f, _g, _h, f, g, h;
    if (a == 0) {
        f = y * s0 % mod;
        g = y * s1 % mod;
        h = y * y % mod * s0 % mod;
    } else if (a >= c || b >= c) {
        std::tie(_f, _g, _h) = Euclidean_like(n, a % c, b % c, c, mod);
        f = (_f + x * s1 + y * s0) % mod;
        g = (_g + x * s2 + y * s1) % mod;
        h = (_h + 2 * y * _f % mod + 2 * x * _g % mod + x * x % mod * s2 % mod + 2 * x * y % mod * s1 % mod
             + y * y % mod * s0 % mod)
            % mod;
    } else {
        std::tie(_f, _g, _h) = Euclidean_like(m - 1, c, c - b - 1, a, mod);
        f = (m * n - _f) % mod;
        g = (m * s1 - inv[2] * _h - inv[2] * _f) % mod;
        h = (m * (m + 1) % mod * n - 2 * _g - 2 * _f - f) % mod;
    }
    return std::make_tuple(f, g, h);
}
