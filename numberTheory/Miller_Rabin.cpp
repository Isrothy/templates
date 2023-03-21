#include <algorithm>
bool witness(long long a, int s, long long d, long long n) {
    long long x = power(a, d, n);
    if (x == 1 || x == n - 1) {
        return false;
    }
    for (int i = 1; i < s; ++i) {
        x = (__int128) x * x % n;
        if (x == n - 1) {
            return false;
        }
    }
    return true;
}

bool Miller_Rabin(long long n) {
    if (n < M) {
        return !is_composite[n] && n != 1;
    }
    long long d = n - 1;
    int s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }
    // int p[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
    long long p[] = {2, 325, 9375, 28178, 450775, 9780504, 1795265022};
    return std::all_of(p, p + 7, [&](long long i) { return !witness(i, s, d, n); });
}
