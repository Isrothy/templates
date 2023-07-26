#include <cmath>
#include <iostream>
#include <optional>
#include <unordered_map>
#include <vector>

namespace polynomial {
    auto pow(int64_t x, int64_t k, int64_t mod) {
        int64_t res = 1;
        while (k != 0) {
            if (k & 1) {
                res = res * x % mod;
            }
            x = x * x % mod;
            k >>= 1;
        }
        return res;
    }
    auto gcd(int64_t a, int64_t b) -> int64_t {
        return b ? gcd(b, a % b) : a;
    }
    auto ex_gcd(int64_t a, int64_t b) -> std::tuple<int64_t, int64_t, int64_t> {
        if (b == 0) {
            return {a, 1, 0};
        }
        auto [d, y, x] = ex_gcd(b, a % b);
        return {d, x, y - a / b * x};
    }
    int64_t inv(int64_t a, int64_t p) {
        a = (a % p + p) % p;
        auto [d, x, y] = ex_gcd(a, p);
        return x;
    }
    auto BSGS(int64_t a, int64_t b, int64_t mod) -> std::optional<int64_t> {
        auto P = static_cast<int64_t>(std::sqrt(mod) + 1);
        std::unordered_map<int64_t, int64_t> map;
        b = (b % mod + mod) % mod;
        for (int k = 1; k <= P; ++k) {
            b = b * a % mod;
            map[b] = k;
        }
        auto w = pow(a, P, mod);
        int64_t x = 1;
        for (int k = 1;; ++k) {
            x = x * w % mod;
            if (map.contains(x)) {
                return k * P - map[x];
            }
        }
        return std::nullopt;
    }
    auto congruence_equation(int64_t a, int64_t b, int64_t p) {
        b = (b % p + p) % p;
        int d = static_cast<int>(gcd(a, p));
        a /= d;
        p /= d;
        b /= d;
        return (b * inv(a, p) % p + p) % p;
    }
    auto quadratic_residue(int64_t x, int64_t mod, int64_t g) -> int64_t {
        if (x == 0) {
            return 0;
        }
        x = BSGS(g, x, mod).value();
        x = pow(g, congruence_equation(2, x, mod - 1), mod);
        return std::min(x, mod - x);
    }

    auto adequate_length(size_t n) {
        size_t m = 1;
        while (m < n) {
            m <<= 1;
        }
        return m;
    }
    template<int64_t MOD, int64_t G> struct Polynmial : protected std::vector<int> {
        using std::vector<int>::vector;
        using std::vector<int>::operator[];
        using std::vector<int>::begin;
        using std::vector<int>::end;
        using std::vector<int>::size;
        using std::vector<int>::resize;

        auto &operator*=(int64_t k) {
            for (auto &x: *this) {
                x = x * k % MOD;
            }
            return *this;
        }
        auto &operator+=(const Polynmial<MOD, G> &rhs) {
            if (rhs.size() > size()) {
                resize(rhs.size());
            }
            for (int i = 0; i < rhs.size(); ++i) {
                (*this)[i] = ((*this)[i] + rhs[i]) % MOD;
            }
            return *this;
        }
        auto &operator-=(const Polynmial<MOD, G> &rhs) {
            if (rhs.size() > size()) {
                resize(rhs.size());
            }
            for (int i = 0; i < rhs.size(); ++i) {
                (*this)[i] = ((*this)[i] - rhs[i]) % MOD;
            }
            return *this;
        }
        auto &operator+=(int64_t x) {
            (*this)[0] = ((*this)[0] + x) % MOD;
            return *this;
        }
        auto &operator-=(int64_t x) {
            (*this)[0] = ((*this)[0] - x) % MOD;
            return *this;
        }
    };

    template<int64_t MOD, int64_t G> auto dft(Polynmial<MOD, G> a, int f) {
        auto n = a.size();
        std::vector<int> w(n);
        for (int i = 0, j = 0; i < n; ++i) {
            if (i < j) {
                std::swap(a[i], a[j]);
            }
            for (int l = n >> 1; (j ^= l) < l; l >>= 1)
                ;
        }
        w[0] = 1;
        for (int i = 1; i < n; i <<= 1) {
            int64_t wn = pow(G, MOD - 1 + f * (MOD - 1) / (i << 1), MOD);
            for (auto j = i - 2; j >= 0; j -= 2) {
                w[j] = w[j >> 1];
                w[j + 1] = w[j] * wn % MOD;
            }
            for (int j = 0; j < n; j += i << 1) {
                auto *p = &a[j], *q = &a[j + i], *r = &w[0];
                for (int k = 0; k < i; ++k) {
                    int64_t t = (int64_t) (*q) * (*r);
                    *q = (*p - t) % MOD;
                    *p = (*p + t) % MOD;
                    ++p;
                    ++q;
                    ++r;
                }
            }
        }
        if (f == -1) {
            int64_t in = inv(n, MOD);
            for (auto &x: a) {
                x = x * in % MOD;
            }
        }
        return a;
    }
    template<int64_t MOD, int64_t G> auto modXN(Polynmial<MOD, G> &&p, size_t n) {
        p.resize(n);
        return p;
    }
    template<int64_t MOD, int64_t G> auto modXN(const Polynmial<MOD, G> &p, size_t n) {
        Polynmial<MOD, G> res(n);
        std::copy(p.begin(), p.begin() + std::min(n, p.size()), res.begin());
        return res;
    }
    template<int64_t MOD, int64_t G> auto divXN(Polynmial<MOD, G> &&p, size_t n) {
        std::copy(p.begin() + n, p.end(), p.begin());
        p.resize(p.size() - n);
        return p;
    }
    template<int64_t MOD, int64_t G> auto divXN(const Polynmial<MOD, G> &p, size_t n) {
        Polynmial res(p.size() - n);
        std::copy(p.begin() + n, p.end(), res.begin());
        return res;
    }
    template<int64_t MOD, int64_t G> auto reverse(Polynmial<MOD, G> p) {
        std::reverse(p.begin(), p.end());
        return p;
    }
    template<int64_t MOD, int64_t G>
    auto operator+(Polynmial<MOD, G> lhs, const Polynmial<MOD, G> &rhs) {
        lhs += rhs;
        return lhs;
    }
    template<int64_t MOD, int64_t G>
    auto operator-(Polynmial<MOD, G> lhs, const Polynmial<MOD, G> &rhs) {
        lhs -= rhs;
        return lhs;
    }

    template<int64_t MOD, int64_t G> auto operator+(int64_t x, Polynmial<MOD, G> p) {
        p += x;
        return p;
    }
    template<int64_t MOD, int64_t G> auto operator+(Polynmial<MOD, G> p, int64_t x) {
        p += x;
        return p;
    }
    template<int64_t MOD, int64_t G> auto operator-(Polynmial<MOD, G> p, int64_t x) {
        p -= x;
        return p;
    }
    template<int64_t MOD, int64_t G> auto operator-(int64_t x, Polynmial<MOD, G> p) {
        p -= x;
        return p;
    }
    template<int64_t MOD, int64_t G> auto operator*(int64_t x, Polynmial<MOD, G> p) {
        p *= x;
        return p;
    }
    template<int64_t MOD, int64_t G> auto operator*(Polynmial<MOD, G> p, int64_t x) {
        p *= x;
        return p;
    }
    template<int64_t MOD, int64_t G> auto operator*(Polynmial<MOD, G> lhs, Polynmial<MOD, G> rhs) {
        auto m = lhs.size() + rhs.size() - 1;
        auto n = adequate_length(m);
        lhs = dft(modXN(std::move(lhs), n), 1);
        rhs = dft(modXN(std::move(rhs), n), 1);
        for (size_t i = 0; i < n; ++i) {
            lhs[i] = (int64_t) lhs[i] * rhs[i] % MOD;
        }
        return modXN(dft(std::move(lhs), -1), m);
    }
    template<int64_t MOD, int64_t G>
    auto &operator*=(Polynmial<MOD, G> &lhs, const Polynmial<MOD, G> &rhs) {
        lhs = lhs * rhs;
        return &lhs;
    }
    template<int64_t MOD, int64_t G> auto inverse(const Polynmial<MOD, G> &p) {
        Polynmial<MOD, G> res = {static_cast<int>(inv(p[0], MOD))};
        auto n = adequate_length(p.size());
        for (size_t i = 2; i <= n; i <<= 1) {
            auto a = dft(modXN(modXN(p, i), i << 1), 1);
            auto b = dft(modXN(std::move(res), i << 1), 1);
            for (size_t j = 0; j < i << 1; ++j) {
                b[j] = b[j] * (2 - (int64_t) a[j] * b[j] % MOD) % MOD;
            }
            res = modXN(dft(std::move(b), -1), i);
        }
        return modXN(std::move(res), p.size());
    }
    template<int64_t MOD, int64_t G> auto derivative(Polynmial<MOD, G> p) {
        for (size_t i = 1; i < p.size(); ++i) {
            p[i - 1] = (int64_t) i * p[i] % MOD;
        }
        p.resize(p.size() - 1);
        return p;
    }
    template<int64_t MOD, int64_t G> auto integral(Polynmial<MOD, G> p) {
        p.resize(p.size() + 1);
        for (size_t i = p.size(); i != 0; --i) {
            p[i] = inv(static_cast<int64_t>(i), MOD) * p[i - 1] % MOD;
        }
        p[0] = 0;
        return p;
    }
    template<int64_t MOD, int64_t G> auto log(const Polynmial<MOD, G> &p) {
        return modXN(integral(derivative(p) * inverse(p)), p.size());
    }
    template<int64_t MOD, int64_t G> auto exp(const Polynmial<MOD, G> &p) {
        Polynmial<MOD, G> res = {1};
        auto n = adequate_length(p.size());
        for (size_t i = 2; i <= n; i <<= 1) {
            auto a = dft(modXN(modXN(p, i), i << 1), 1);
            auto b = dft(modXN(res, i << 1), 1);
            auto c = dft(modXN(log(modXN(std::move(res), i)), i << 1), 1);
            for (size_t j = 0; j < i << 1; ++j) {
                b[j] = (int64_t) b[j] * (1 + a[j] - c[j]) % MOD;
            }
            res = modXN(dft(std::move(b), -1), i);
        }
        return modXN(std::move(res), p.size());
    }
    template<int64_t MOD, int64_t G> auto pow(const Polynmial<MOD, G> &p, int64_t k) {
        return exp(log(p) * k);
    }
    template<int64_t MOD, int64_t G> auto sqrt(const Polynmial<MOD, G> &p) {
        Polynmial<MOD, G> res = {static_cast<int>(quadratic_residue(p[0], MOD, G))};
        int64_t inv2 = inv(2, MOD);
        auto n = adequate_length(p.size());
        for (size_t i = 2; i <= n; i <<= 1) {
            auto a = dft(modXN(modXN(p, i), i << 1), 1);
            auto b = dft(modXN(res, i << 1), 1);
            auto c = dft(modXN(inverse(modXN(std::move(res), i)), i << 1), 1);
            for (size_t j = 0; j < i << 1; ++j) {
                b[j] = (b[j] + (int64_t) a[j] * c[j]) % MOD * inv2 % MOD;
            }
            res = modXN(dft(std::move(b), -1), i);
        }
        return modXN(std::move(res), p.size());
    }
    template<int64_t MOD, int64_t G>
    auto operator/(const Polynmial<MOD, G> &lhs, const Polynmial<MOD, G> &rhs) {
        auto n = lhs.size();
        auto m = rhs.size();
        if (n < m) {
            return Polynmial<MOD, G>{0};
        }
        auto a = modXN(reverse(lhs), n - m + 1);
        auto b = modXN(reverse(rhs), n - m + 1);
        return reverse(modXN(a * inverse(b), n - m + 1));
    }
    template<int64_t MOD, int64_t G>
    auto operator%(const Polynmial<MOD, G> &lhs, const Polynmial<MOD, G> &rhs) {
        auto m = rhs.size();
        return modXN(lhs - lhs / rhs * rhs, m - 1);
    }
    template<int64_t MOD, int64_t G>
    auto operator/=(Polynmial<MOD, G> &lhs, const Polynmial<MOD, G> &rhs) {
        lhs = lhs / rhs;
        return lhs;
    }
    template<int64_t MOD, int64_t G>
    auto operator%=(Polynmial<MOD, G> &lhs, const Polynmial<MOD, G> &rhs) {
        lhs = lhs % rhs;
        return lhs;
    }

    template<int64_t MOD, int64_t G>
    auto eva_build(
        size_t p, size_t l, size_t r, const std::vector<int> &x, std::vector<Polynmial<MOD, G>> &a
    ) {
        if (l == r) {
            a[p] = {1, l < x.size() ? -x[l] : 0};
            return;
        }
        auto mid = (l + r) >> 1;
        eva_build(p << 1, l, mid, x, a);
        eva_build(p << 1 | 1, mid + 1, r, x, a);
        a[p] = a[p << 1] * a[p << 1 | 1];
    }
    template<int64_t MOD, int64_t G>
    auto eva_work(
        size_t p,
        size_t l,
        size_t r,
        const Polynmial<MOD, G> &f,
        std::vector<Polynmial<MOD, G>> &a,
        std::vector<int> &res
    ) {
        if (l == r) {
            if (l < res.size()) {
                res[l] = f[0];
            }
            return;
        }
        size_t mid = (l + r) >> 1;
        auto helper = [](const Polynmial<MOD, G> &f, const Polynmial<MOD, G> &g) {
            size_t n = adequate_length(f.size());
            auto a = dft(modXN(f, n), 1);
            auto b = dft(modXN(g, n), 1);
            for (int i = 0; i < n; ++i) {
                a[i] = (int64_t) a[i] * b[i] % MOD;
            }
            return divXN(modXN(dft(std::move(a), -1), f.size()), g.size() - 1);
        };
        auto lf = helper(f, a[p << 1 | 1]);
        auto rf = helper(f, a[p << 1]);
        eva_work(p << 1, l, mid, lf, a, res);
        eva_work(p << 1 | 1, mid + 1, r, rf, a, res);
    }

    template<int64_t MOD, int64_t G>
    auto evaluation(const Polynmial<MOD, G> &p, const std::vector<int> &x) {
        size_t m = std::max(x.size(), p.size() - 1);
        std::vector<Polynmial<MOD, G>> a(m << 2);
        std::vector<int> res(x.size());
        eva_build(1, 0, m - 1, x, a);
        auto f = modXN(reverse(modXN(p, m + 1)) * inverse(a[1]), m + 1);
        eva_work(1, 0, m - 1, f, a, res);
        for (size_t i = 0; i < x.size(); ++i) {
            res[i] = (p[0] + (int64_t) res[i] * x[i]) % MOD;
        }
        return res;
    }
    template<int64_t MOD, int64_t G>
    Polynmial<MOD, G> interpolation_work(
        size_t p,
        size_t l,
        size_t r,
        const std::vector<int> &y,
        std::vector<Polynmial<MOD, G>> &a,
        const std::vector<int> &b
    ) {
        if (l == r) {
            return {(int) (y[l] * inv(b[l], MOD) % MOD)};
        }
        auto mid = (l + r) >> 1;
        auto lf = interpolation_work(p << 1, l, mid, y, a, b);
        auto rf = interpolation_work(p << 1 | 1, mid + 1, r, y, a, b);
        return lf * reverse(a[p << 1 | 1]) + rf * reverse(a[p << 1]);
    }

    template<int64_t MOD, int64_t G>
    auto interpolation(const std::vector<int> &x, const std::vector<int> &y) {
        auto n = x.size();
        std::vector<Polynmial<MOD, G>> a(n << 2);
        std::vector<int> b(n);
        eva_build(1, 0, n - 1, x, a);
        auto f = derivative(reverse(a[1]));
        auto g = modXN(reverse(modXN(f, n + 1)) * inverse(a[1]), n + 1);
        eva_work(1, 0, n - 1, g, a, b);
        for (int i = 0; i < n; ++i) {
            b[i] = (f[0] + (long long) b[i] * x[i]) % MOD;
        }
        return interpolation_work(1, 0, n - 1, y, a, b);
    }
}// namespace polynomial
