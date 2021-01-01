long long inverse(long long a, long long mod) {
    long long x, y;
    ex_gcd(a, mod, x, y);
    return (x % mod + mod) % mod;
}
