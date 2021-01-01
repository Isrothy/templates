namespace ex_Lucas {
    long long S[K][M], P[K], Pk[K], Mt[K];
    int t;
    
    void initialize(long long mod) {
        t = 0;
        long long m = mod;
        for (int i = 2; i * i <= m; ++i) {
            if (mod % i == 0) {
                P[t] = i;
                Pk[t] = 1;
                while (m % i == 0) {
                    m /= i;
                    Pk[t] *= i;
                }
                S[t][0] = 1;
                for (int j = 1; j <= Pk[t]; ++j) {
                    S[t][j] = j % i == 0 ? S[t][j - 1] : S[t][j - 1] * j % Pk[t];
                }
                Mt[t] = (mod / Pk[t]) * inverse(mod / Pk[t], Pk[t]) % mod;
                ++t;
            }
        }
        if (m != 1) {
            P[t] = Pk[t] = m;
            S[t][0] = 1;
            for (int j = 1; j <= Pk[t]; ++j) {
                S[t][j] = j % m == 0 ? S[t][j - 1] : S[t][j - 1] * j % Pk[t];
            }
            Mt[t] = (mod / m) * inverse(mod / m, m) % mod;
            ++t;
        }
    }
    
    long long f(long long n, long long p, long long mod, long long *S) {
        long long res = 1;
        while (n != 0) {
            res = res * power(S[mod], n / mod, mod) % mod * S[n % mod] % mod;
            n /= p;
        }
        return res;
    }
    
    long long g(long long n, long long p) {
        long long res = 0;
        while (n != 0) {
            n = n / p;
            res += n;
        }
        return res;
    }
    
    long long combination(long long n, long long m) {
        long long res = 0;
        for (int i = 0; i < t; ++i) {
            long long x = f(n, P[i], Pk[i], S[i]);
            long long y = inverse(f(n - m, P[i], Pk[i], S[i]), Pk[i]);
            long long z = inverse(f(m, P[i], Pk[i], S[i]), Pk[i]);
            long long r = g(n, P[i]) - g(m, P[i]) - g(n - m, P[i]);
            long long b = x * y % Pk[i] * z % Pk[i] * power(P[i], r, Pk[i]) % mod;
            res = (res + b * Mt[i]) % mod;
        }
        return res;
    }
}
