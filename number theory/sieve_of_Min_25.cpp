namespace Min_25 {
    long long w[M], h0[M], h1[M], h[M];
    int s, m;

    int index(long long x, long long n) {
        return x <= s ? x : m - n / x + 1;
    }

    int dfs_mu(long long x, int k, long long n) {
        if (x <= prime[k]) {
            return 0;
        }
        int res = (h[index(x, n)] - h[prime[k]]) % mod;
        for (int i = k + 1; (long long) prime[i] * prime[i] <= x; ++i) {
            res = (res - dfs_mu(x / prime[i], i, n)) % mod;
        }
        return res;
    }

    int dfs_phi(long long x, int k, long long n) {
        if (x <= prime[k]) {
            return 0;
        }
        int res = (h[index(x, n)] - h[prime[k]]) % mod;
        for (int i = k + 1; (long long) prime[i] * prime[i] <= x; ++i) {
            long long p = prime[i], d = prime[i], g = p - 1;
            while (d * p <= x) {
                res = (res + g * (dfs_phi(x / d, i, n) + p)) % mod;
                g = (g * p) % mod;
                d *= p;
            }
        }
        return res;
    }

    int sum_of_mu(long long n) {
        m = 0;
        s = (int) sqrtl(n);
        for (long long x = 1; x <= n; x = w[m] + 1) {
            w[++m] = n / (n / x);
            h0[m] = w[m];
        }
        for (int i = 1; prime[i] <= s; ++i) {
            long long p = prime[i];
            for (int j = m; p * p <= w[j]; --j) {
                int k = index(w[j] / p, n);
                h0[j] -= h0[k] - h0[p - 1];
            }
        }
        for (int i = 2; i <= m; ++i) {
            h[i] = -h0[i];
        }
        return 1 + dfs_mu(n, 0, n);
    }

    int sum_of_phi(long long n) {
        m = 0;
        s = (int) sqrtl(n);
        for (long long x = 1; x <= n; x = w[m] + 1) {
            w[++m] = n / (n / x);
            long long tmp = w[m] % mod;
            h0[m] = tmp;
            h1[m] = ((tmp * (tmp + 1)) >> 1) % mod;
        }
        for (int i = 1; prime[i] <= s; ++i) {
            long long p = prime[i];
            for (int j = m; p * p <= w[j]; --j) {
                int k = index(w[j] / p, n);
                h0[j] -= h0[k] - h0[p - 1];
                h1[j] -= (h1[k] - h1[p - 1]) * p % mod;
            }
        }
        for (int i = 2; i <= m; ++i) {
            h[i] = (h1[i] - h0[i]) % mod;
        }
        return (1 + dfs_phi(n, 0, n)) % mod;
    }
}// namespace Min_25
