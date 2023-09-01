#include <cstdlib>
#include <random>
std::mt19937_64 mt_rand(std::random_device{}());
long long f(long long x, long long c, long long n) { return ((__int128) x * x + c) % n; }
long long Pollard_Rho(long long n) {
    if (Miller_Rabin(n)) { return n; }
    int prime[11] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31};
    for (int i = 0; i <= 10; ++i) {
        if (n % prime[i] == 0) { return (long long) prime[i]; }
    }
    std::uniform_int_distribution<long long> dist(1, n - 1);
    while (true) {
        long long c = dist(mt_rand);
        long long t = f(0, c, n), r = f(f(0, c, n), c, n);
        while (t != r) {
            long long d = std::gcd(llabs(t - r), n);
            if (d != 1) { return d; }
            t = f(t, c, n);
            r = f(f(r, c, n), c, n);
        }
    }
}
