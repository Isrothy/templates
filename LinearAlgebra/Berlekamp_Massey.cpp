int Berlekamp_Massey(int *A, int *f, int n) {
    static int g[M], tmp[M];
    int k = 0, last_k = 0, last_delta, last = -1;
    for (int i = 0; i <= n; ++i) { tmp[i] = f[i] = 0; }
    for (int i = 0; i < n; ++i) {
        long long delta = -A[i];
        for (int j = 1; j <= k; ++j) { delta = (delta + (long long) f[j] * A[i - j]) % mod; }
        if (delta == 0) { continue; }
        if (last == -1) {
            k = i + 1;
        } else {
            long long t = delta * power(last_delta, mod - 2) % mod;
            tmp[i - last] = (tmp[i - last] + t) % mod;
            for (int j = 1; j <= last_k; ++j) { tmp[i - last + j] = (tmp[i - last + j] - t * g[j]) % mod; }
            int p = last_k;
            last_k = k;
            k = max(k, i - last + p);
            for (int j = 1; j <= last_k; ++j) { g[j] = f[j]; }
            for (int j = 1; j <= k; ++j) { f[j] = tmp[j]; }
        }
        last_delta = delta;
        last = i;
    }
    return k;
}
