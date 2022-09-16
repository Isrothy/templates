long long f(long long x, long long c, long long n) {
    return ((__int128) x * x + c) % n;
}

long long Pollard_Rho(long long n) {
    if (Miller_Rabin(n)) {
        return n;
    }
    for (int i = 1; i <= 10; ++i) {
        if (n % prime[i] == 0) {
            return (long long) prime[i];
        }
    }
    while (true) {
        long long c = mt_rand() % (n - 1) + 1;
        long long t = f(0, c, n), r = f(f(0, c, n), c, n);
        while (t != r) {
            long long d = gcd(llabs(t - r), n);
            if (d != 1) {
                return d;
            }
            t = f(t, c, n);
            r = f(f(r, c, n), c, n);
        }
    }
}
