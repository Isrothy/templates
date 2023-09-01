#include <optional>
#include <vector>
int64_t ex_gcd(int64_t a, int64_t b, int64_t &x, int64_t &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    int64_t d = ex_gcd(b, a % b, y, x);
    y -= a / b * x;
    return d;
}
int64_t inv(int64_t a, int64_t mod) {
    int64_t x, y;
    a = (a + mod) % mod;
    ex_gcd(a, mod, x, y);
    return (x % mod + mod) % mod;
}
template<int64_t mod>
struct Matrix : private std::vector<std::vector<int64_t>> {
    size_t n{}, m{};
    Matrix() = default;
    using std::vector<std::vector<int64_t>>::vector;
    int64_t operator()(size_t i, size_t j) const { return (*this)[i][j]; }
    int64_t &operator()(size_t i, size_t j) { return (*this)[i][j]; }
    explicit Matrix(std::vector<std::vector<int64_t>> v) : std::vector<std::vector<int64_t>>(v) {
        n = v.size();
        m = v[0].size();
    }
    Matrix(size_t n, size_t m) : n(n), m(m) {
        resize(n);
        for (int i = 0; i < n; ++i) { (*this)[i].resize(m); }
    }
    Matrix augment(std::vector<int64_t> v) const {
        assert(n == v.size());
        Matrix res(n, m + 1);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) { res[i][j] = (*this)[i][j]; }
            res[i][m] = v[i];
        }
        return res;
    }
    Matrix augment(Matrix B) const {
        assert(n == B.n);
        Matrix res(n, m + B.m);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) { res[i][j] = (*this)[i][j]; }
            for (int j = 0; j < B.m; ++j) { res[i][m + j] = B[i][j]; }
        }
        return res;
    }
    Matrix transpose() const {
        Matrix res(m, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m; ++j) { res[j][i] = (*this)[i][j]; }
        }
        return res;
    }
    Matrix remove_column(size_t k) const {
        assert(m != 1 && k < m);
        Matrix res(n, m - 1);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < m - 1; ++j) { res[i][j] = (*this)[i][j < k ? j : j + 1]; }
        }
        return res;
    }
    Matrix remove_row(size_t k) const {
        assert(n != 1 && k < n);
        Matrix res(n - 1, m);
        for (int i = 0; i < n - 1; ++i) {
            for (int j = 0; j < m; ++j) { res[i][j] = (*this)[i < k ? i : i + 1][j]; }
        }
        return res;
    }
    int64_t determinate() const {
        assert(n == m);
        Matrix tmp = *this;
        int64_t res = 1;
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                while (tmp[j][i] != 0) {
                    int64_t t = tmp[i][i] / tmp[j][i];
                    for (int k = i; k < n; ++k) { tmp[i][k] = (tmp[i][k] - t * tmp[j][k]) % mod; }
                    std::swap(tmp[i], tmp[j]);
                    res = -res;
                }
            }
            res = (res * tmp[i][i]) % mod;
        }
        return res;
    }
    std::vector<std::vector<int64_t>> Gaussian_elimination(const std::vector<int64_t> &v) const {
        assert(n == v.size());
        std::vector<int64_t> v0(m);
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
            if (pivot == m) { break; }
            int64_t t = inv(tmp[i][pivot], mod);
            for (int j = pivot; j <= m; ++j) { tmp[i][j] = tmp[i][j] * t % mod; }
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    int64_t s = tmp[j][pivot];
                    for (int k = pivot; k <= m; ++k) { tmp[j][k] = (tmp[j][k] - tmp[i][k] * s) % mod; }
                }
            }
            p[i] = pivot++;
        }
        for (int i = 0; i < n; ++i) {
            if (p[i] == -1) {
                if (tmp[i][m] != 0) { return {}; }
            } else {
                v0[p[i]] = tmp[i][m];
            }
        }
        std::vector<std::vector<int64_t>> res;
        res.push_back(v0);
        for (auto i: f) {
            std::vector<int64_t> v(m, 0);
            v[i] = 1;
            for (int j = 0; j < n; ++j) {
                if (i != j && p[j] != -1) { v[p[j]] = -tmp[j][i]; }
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
                if (tmp[i][i] == 0) { return std::nullopt; }
            }
            int64_t t = inv(tmp[i][i], mod);
            for (int j = i; j < 2 * n; ++j) { tmp[i][j] = tmp[i][j] * t % mod; }
            for (int j = 0; j < n; ++j) {
                if (i != j) {
                    int64_t s = tmp[j][i];
                    for (int k = i; k < 2 * n; ++k) { tmp[j][k] = (tmp[j][k] - tmp[i][k] * s) % mod; }
                }
            }
        }
        Matrix res(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) { res[i][j] = tmp[i][j + n]; }
        }
        return res;
    }
    static Matrix identity(size_t n) {
        Matrix res(n, n);
        for (int i = 0; i < n; ++i) { res[i][i] = 1; }
        return res;
    }
    static Matrix zero(size_t n) { return {n, n}; }
};
