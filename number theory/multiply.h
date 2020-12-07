long long multiply(long long x, long long k, long long mod) {
    k = (k % mod + mod) % mod;
    long long res =  0;
    while (k != 0) {
        if ((k & 1) == 1)
            res = (res + x) % mod;
        x = (x + x) % mod;
        k >>= 1;
    }
    return res;
}
