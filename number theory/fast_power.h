long long power(long long x, int k) {
    long long res = 1;
    while (k != 0) {
        if ((k & 1) == 1)
            res = res * x % mod;
        x = x * x % mod;
    }
    return res;
}
