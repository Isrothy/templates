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

struct Matrix {
    static const int MOD = 998244353;
    std::vector<std::vector<long long>> mat;
    size_t n, m;
    Matrix(size_t n, size_t m) : n(n), m(m) {
        mat.resize(n);
        for (auto &v: mat) {
            v.resize(m);
        }
    }
    Matrix(std::vector<std::vector<long long>> mat) : mat(mat) {
        n = mat.size();
        m = mat[0].size();
    }
    long long &operator()(size_t i, size_t j) {
        return mat[i][j];
    }
    std::vector<long long> &operator()(size_t i) {
        return mat[i];
    }

    Matrix argument(std::vector<long long> v) const {
        assert(n == v.size());
        Matrix res(n, m + 1);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                res(i, j) = mat[i][j];
            }
            res(i, m) = v[i];
        }
        return res;
    }

    Matrix argument(Matrix B) const {
        assert(n == B.n);
        Matrix res(n, m + B.m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                res(i, j) = mat[i][j];
            }
            for (int j = 0; j < B.m; ++j) {
                res(i, m + j) = B(i, j);
            }
        }
        return res;
    }

    Matrix transpose() const {
        Matrix res(m, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) {
                res(j, i) = mat[i][j];
            }
        }
        return res;
    }

    Matrix remove_column(size_t k) const {
        assert(m != 1 && k < m);
        Matrix res(n, m - 1);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m - 1; ++j) {
                res(i, j) = mat[i][j < k ? j : j + 1];
            }
        }
        return res;
    }
    Matrix remove_row(size_t k) const {
        assert(n != 1 && k < n);
        Matrix res(n - 1, m);
        for (int i = 0; i < n - 1; ++i) {
            for (int j = 0; j < m; ++j) {
                res(i, j) = mat[i < k ? i : i + 1][j];
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
                while (tmp(j, i) != 0) {
                    long long t = tmp(i, i) / tmp(j, i);
                    for (int k = i; k < n; ++k) {
                        tmp(i, k) = (tmp(i, k) - t * tmp(j, k)) % MOD;
                    }
                    std::swap(tmp(i), tmp(j));
                    res = -res;
                }
            }
            res = (res * tmp(i, i)) % MOD;
        }
        return res;
    }

    std::vector<std::vector<long long>> Gaussian_elimination(const std::vector<long long> &v) const {
        assert(n == v.size());
        std::vector<long long> v0(m);
        Matrix tmp = this->argument(v);
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                while (tmp(j, i) != 0) {
                    long long t = tmp(i, i) / tmp(j, i);
                    for (int k = i; k <= m; ++k) {
                        tmp(i, k) = (tmp(i, k) - t * tmp(j, k)) % MOD;
                    }
                    std::swap(tmp(i), tmp(j));
                }
            }
        }

        for (int i = n - 1; i >= 0; --i) {
            if (tmp(i, i) == 0) {
                continue;
            }
            long long t = inv(tmp(i, i), MOD);
            for (int j = i; j <= m; ++j) {
                tmp(i, j) = tmp(i, j) * t % MOD;
            }
            tmp(i, i) = 1;
            for (int j = 0; j < i; ++j) {
                auto s = tmp(j, i);
                for (int k = i; k <= m; ++k) {
                    tmp(j, k) = (tmp(j, k) - tmp(i, k) * s) % MOD;
                }
            }
            v0[i] = tmp(i, m);
        }
        std::vector<std::vector<long long>> res;
        res.push_back(v0);
        for (int i = 0; i < n; ++i) {
            if (tmp(i, i) == 0) {
                std::vector<long long> v(m, 0);
                v[i] = 1;
                for (int j = 0; j < n; ++j) {
                    if (i != j) {
                        v[j] = -tmp(j, i);
                    }
                }
                res.push_back(v);
            }
        }
        return res;
    }

    std::optional<Matrix> inverse() const {
        assert(n == m);
        auto tmp = this->argument(Matrix::identity(n));
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                while (tmp(j, i) != 0) {
                    long long t = tmp(i, i) / tmp(j, i);
                    for (int k = i; k < 2 * n; ++k) {
                        tmp(i, k) = (tmp(i, k) - t * tmp(j, k)) % MOD;
                    }
                    std::swap(tmp(i), tmp(j));
                }
            }
            if (tmp(i, i) == 0) {
                return {};
            }
            long long t = inv(tmp(i, i), MOD);
            for (int j = i; j < 2 * n; ++j) {
                tmp(i, j) = tmp(i, j) * t % MOD;
            }
            tmp(i, i) = 1;
            for (int j = 0; j < i; ++j) {
                auto s = tmp(j, i);
                for (int k = i; k < 2 * n; ++k) {
                    tmp(j, k) = (tmp(j, k) - tmp(i, k) * s) % MOD;
                }
            }
        }
        Matrix res(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                res(i, j) = tmp(i, j + n);
            }
        }
        return res;
    }
    static Matrix identity(size_t n) {
        Matrix res(n, n);
        for (int i = 0; i < n; ++i) {
            res(i, i) = 1;
        }
        return res;
    }
    static Matrix zero(size_t n) {
        return {n, n};
    }
};
