#include <optional>
#include <vector>

long long ex_gcd(long long a, long long b, long long &x, long long &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    long long d = ex_gcd(b, a % b, y, x);
    y -= a / b * x;
    return d;
}

long long inv(long long a, long long mod) {
    long long x, y;
    a = (a + mod) % mod;
    ex_gcd(a, mod, x, y);
    return (x % mod + mod) % mod;
}

struct Matrix : public std::vector<std::vector<long long>> {
    static const int MOD = 1e9 + 7;
    size_t n{}, m{};
    Matrix() = default;
    explicit Matrix(std::vector<std::vector<long long>> v)
        : std::vector<std::vector<long long>>(v) {
        n = v.size();
        m = v[0].size();
    }
    Matrix(size_t n, size_t m) : n(n), m(m) {
        resize(n);
        for (int i = 0; i < n; ++i) {
            (*this)[i].resize(m);
        }
    }

    Matrix augment(std::vector<long long> v) const {
        assert(n == v.size());
        Matrix res(n, m + 1);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                res[i][j] = (*this)[i][j];
            }
            res[i][m] = v[i];
        }
        return res;
    }

    Matrix augment(Matrix B) const {
        assert(n == B.n);
        Matrix res(n, m + B.m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                res[i][j] = (*this)[i][j];
            }
            for (int j = 0; j < B.m; ++j) {
                res[i][m + j] = B[i][j];
            }
        }
        return res;
    }

    Matrix transpose() const {
        Matrix res(m, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                res[j][i] = (*this)[i][j];
            }
        }
        return res;
    }

    Matrix remove_column(size_t k) const {
        assert(m != 1 && k < m);
        Matrix res(n, m - 1);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m - 1; ++j) {
                res[i][j] = (*this)[i][j < k ? j : j + 1];
            }
        }
        return res;
    }
    Matrix remove_row(size_t k) const {
        assert(n != 1 && k < n);
        Matrix res(n - 1, m);
        for (int i = 0; i < n - 1; ++i) {
            for (int j = 0; j < m; ++j) {
                res[i][j] = (*this)[i < k ? i : i + 1][j];
            }
        }
        return res;
    }

    long long determinate() const {
        assert(n == m);
        Matrix tmp = *this;
        long long res = 1;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                while (tmp[j][i] != 0) {
                    long long t = tmp[i][i] / tmp[j][i];
                    for (int k = i; k < n; ++k) {
                        tmp[i][k] = (tmp[i][k] - t * tmp[j][k]) % MOD;
                    }
                    std::swap(tmp[i], tmp[j]);
                    res = -res;
                }
            }
            res = (res * tmp[i][i]) % MOD;
        }
        return res;
    }

    std::vector<std::vector<long long>> Gaussian_elimination(const std::vector<long long> &v
    ) const {
        assert(n == v.size());
        std::vector<long long> v0(m);
        std::vector<int> p(n, -1), f;
        Matrix tmp = this->augment(v);
        for (int i = 0, pivot = 0; i < n; ++i) {
            while (pivot < m && tmp[i][pivot] == 0) {
                for (int j = i + 1; j < n; ++j) {
                    if (tmp[j][pivot] != 0) {
                        std::swap(tmp[i], tmp[j]);
                        break;
                    }
                }
                if (tmp[i][pivot] == 0) {
                    f.push_back(pivot);
                    ++pivot;
                }
            }
            if (pivot == m) {
                break;
            }
            long long t = inv(tmp[i][pivot], MOD);
            for (int j = pivot; j <= m; ++j) {
                tmp[i][j] = tmp[i][j] * t % MOD;
            }
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    long long s = tmp[j][pivot];
                    for (int k = pivot; k <= m; ++k) {
                        tmp[j][k] = (tmp[j][k] - tmp[i][k] * s) % MOD;
                    }
                }
            }
            p[i] = pivot++;
        }
        for (int i = 0; i < n; ++i) {
            if (p[i] == -1) {
                if (tmp[i][m] != 0) {
                    return {};
                }
            } else {
                v0[p[i]] = tmp[i][m];
            }
        }
        std::vector<std::vector<long long>> res;
        res.push_back(v0);
        for (auto i: f) {
            std::vector<long long> v(m, 0);
            v[i] = 1;
            for (int j = 0; j < n; ++j) {
                if (i != j && p[j] != -1) {
                    v[p[j]] = -tmp[j][i];
                }
            }
            res.push_back(v);
        }
        return res;
    }
    std::optional<Matrix> inverse() const {
        assert(n == m);
        auto tmp = this->augment(Matrix::identity(n));
        for (int i = 0; i < n; ++i) {
            if (tmp[i][i] == 0) {
                for (int j = i + 1; j < n; ++j) {
                    if (tmp[j][i] != 0) {
                        std::swap(tmp[i], tmp[j]);
                        break;
                    }
                }
                if (tmp[i][i] == 0) {
                    return std::nullopt;
                }
            }
            long long t = inv(tmp[i][i], MOD);
            for (int j = i; j < 2 * n; ++j) {
                tmp[i][j] = tmp[i][j] * t % MOD;
            }
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    long long s = tmp[j][i];
                    for (int k = i; k < 2 * n; ++k) {
                        tmp[j][k] = (tmp[j][k] - tmp[i][k] * s) % MOD;
                    }
                }
            }
        }
        Matrix res(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                res[i][j] = tmp[i][j + n];
            }
        }
        return res;
    }
    static Matrix identity(size_t n) {
        Matrix res(n, n);
        for (int i = 0; i < n; ++i) {
            res[i][i] = 1;
        }
        return res;
    }
    static Matrix zero(size_t n) {
        return {n, n};
    }
};
