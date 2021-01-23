tuple<long long, long long, long long> Euclidean_like(long long n, long long a, long long b, long long c) {
    long long x = a / c % mod, y = b / c % mod, _n = n % mod;
    long long s0 = (_n + 1) % mod;
    long long s1 = _n * (_n + 1) % mod * inv[2] % mod;
    long long s2 = _n * (_n + 1) % mod * (2 * _n + 1) % mod * inv[6] % mod;
    long long m = ((__int128) a * n + b) / c, _f, _g, _h, f, g, h;
    if (a == 0) {
        f = y * s0 % mod;
        g = y * s1 % mod;
        h = y * y % mod * s0 % mod;
    } else if (a >= c || b >= c) {
        tie(_f, _g, _h) = Euclidean_like(n, a % c, b % c, c);
        f = (_f + x * s1 + y * s0) % mod;
        g = (_g + x * s2 + y * s1) % mod;
        h = (_h + 2 * y * _f % mod + 2 * x * _g % mod + x * x % mod * s2 + 2 * x * y % mod * s1 + y * y % mod * s0) % mod;
    } else {
        tie(_f, _g, _h) = Euclidean_like(m - 1, c, c - b - 1, a);
        f = (m * n - _f) % mod;
        g = (m * s1 - inv[2] * _h - inv[2] * _f) % mod;
        h = (m * (m + 1) % mod * n - 2 * _g - 2 * _f - f) % mod;
    }
    return make_tuple(f, g, h);
}
