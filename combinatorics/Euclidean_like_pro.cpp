long long power_sum(int m, long long n) {
    n %= mod;
    long long sum = 0, x = n;
    for (int i = m; i >= 0; --i) {
        sum = (sum + P[m][i] * x) % mod;
        x = x * n % mod;
    }
    return sum;
}

void initialize() {
    inv[1] = 1;
    for (int i = 2; i < M; ++i) {
        inv[i] = -mod / i * inv[mod % i] % mod;
    }
    for (int i = 0; i < M; ++i) {
        C[i][0] = 1;
        for (int j = 1; j <= i; ++j) {
            C[i][j] = (C[i - 1][j] + C[i - 1][j - 1]) % mod;
        }
    }
    B[0] = 1;
    for (int i = 1; i < M; ++i) {
        B[i] = 1;
        for (int j = 0; j < i; ++j) {
            B[i] = (B[i] - C[i][j] * B[j] % mod * inv[i - j + 1]) % mod;
        }
    }
    for (int i = 0; i < M - 1; ++i) {
        for (int j = 0; j <= i + 1; ++j) {
            P[i][j] = inv[i + 1] * C[i + 1][j] % mod * B[j] % mod;
        }
    }
}

vector<vector<long long>> Euclidean_like(long long n, long long a, long long b, long long c, int K) {
    //sum_{x=0}^n x^k1((ax + b)/c)^k2
    vector<vector<long long>> res;
    long long s[K + 1], power_x[K + 1], power_y[K + 1], power_m[K + 1];
    long long x = a / c % mod, y = b / c % mod;
    long long m = ((__int128) a * n + b) / c;
    for (int i = 0; i <= K; ++i) {
        s[i] = power_sum(i, n);
    }
    power_x[0] = power_y[0] = power_m[0] = 1;
    for (int i = 1; i <= K; ++i) {
        power_x[i] = power_x[i - 1] * x % mod;
        power_y[i] = power_y[i - 1] * y % mod;
        power_m[i] = power_m[i - 1] * (m % mod) % mod;
    }
    res.resize(K + 1);
    for (int i = 0; i <= K; ++i) {
        res[i].resize(K - i + 1);
    }
    if (a == 0) {
        for (int k2 = 0; k2 <= K; ++k2) {
            for (int k1 = 0; k1 <= K - k2; ++k1) {
                res[k1][k2] = power_y[k2] * (k1 == 0 ? (n + 1) : s[k1]) % mod;
            }
        }
    } else if (a >= c) {
        auto tmp = Euclidean_like(n, a % c, b, c, K);
        for (int k2 = 0; k2 <= K; ++k2) {
            for (int j = 0; j <= k2; ++j) {
                long long u = C[k2][j] * power_x[j] % mod;
                for (int k1 = 0; k1 <= K - k2; ++k1) {
                    res[k1][k2] = (res[k1][k2] + u * tmp[k1 + j][k2 - j]) % mod;
                }
            }
        }
    } else if (b >= c) {
        auto tmp = Euclidean_like(n, a, b % c, c, K);
        for (int k2 = 0; k2 <= K; ++k2) {
            for (int j = 0; j <= k2; ++j) {
                long long u = C[k2][j] * power_y[j] % mod;
                for (int k1 = 0; k1 <= K - k2; ++k1) {
                    res[k1][k2] = (res[k1][k2] + u * tmp[k1][k2 - j]) % mod;
                }
            }
        }
    } else {
        auto tmp = Euclidean_like(m - 1, c, c - b - 1, a, K), D = res;
        for (int k1 = 0; k1 <= K; ++k1) {
            res[k1][0] = (k1 == 0 ? (n + 1) : s[k1]) % mod;
        }
        for (int i = 0; i <= K; ++i) {
            for (int j = 0; j <= K - i; ++j) {
                for (int k = 0; k <= i; ++k) {
                    D[i][j] = (D[i][j] + C[i + 1][k] * tmp[k][j]) % mod;
                }
            }
        }
        for (int k2 = 1; k2 <= K; ++k2) {
            for (int k1 = 0; k1 <= K - k2; ++k1) {
                res[k1][k2] = power_m[k2] * s[k1] % mod;
                for (int j = 0; j <= k1; ++j) {
                    res[k1][k2] = (res[k1][k2] - P[k1][j] * D[k2 - 1][k1 + 1 - j]) % mod;
                }
            }
        }
    }
    return res;
}
