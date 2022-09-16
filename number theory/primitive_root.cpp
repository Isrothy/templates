long long power(long long x, int k, int mod) {
    long long res = 1;
    while (k != 0) {
        if ((k & 1) == 1) {
            res = res * x % mod;
        }
        x = x * x % mod;
        k >>= 1;
    }
    return res;
}

vector<int> prime_factors(int n) {
    vector<int> res;
    for (int i = 1; prime[i] * prime[i] <= n; ++i) {
        if (n % prime[i] == 0) {
            while (n % prime[i] == 0) {
                n /= prime[i];
            }
            res.push_back(prime[i]);
        }
    }
    if (n != 1)
        res.push_back(n);
    sort(res.begin(), res.end());
    return res;
}

int phi(int n) {
    vector<int> A = prime_factors(n);
    for (auto p: A) {
        n = n / p * (p - 1);
    }
    return n;
}

int primitive_root(int n) {
    if (n == 2 || n == 4) {
        return n - 1;
    }
    if ((n & 3) == 0) {
        return 0;
    }
    vector<int> A = prime_factors(n);
    if (2 < A.size() || (A.size() == 2 && A[0] != 2)) {
        return false;
    }
    int m = A.size() == 2 ? n / 2 / A[1] * (A[1] - 1) : n / A[0] * (A[0] - 1);
    vector<int> B = prime_factors(m);
    for (int g = 2; g <= n; ++g) {
        if (power(g, m, n) == 1) {
            bool flag = true;
            for (auto p: B) {
                if (power(g, m / p, n) == 1) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                return g;
            }
        }
    }
    return 0;
}
