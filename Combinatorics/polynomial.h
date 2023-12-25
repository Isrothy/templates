#include "../NumberTheory/common.h"
#include <bit>
#include <cmath>
#include <optional>
#include <vector>
namespace polynomial {
    auto congruence_equation(int64_t a, int64_t b, int64_t p) {
        b = (b % p + p) % p;
        auto d = std::gcd(a, p);
        std::tie(a, b, p) = std::make_tuple(a / d, b / d, p / d);
        return (b * inverse(a, p) % p + p) % p;
    }
    auto quadratic_residue(int64_t x, int64_t mod, int64_t g) -> int64_t {
        if (x == 0) { return 0; }
        x = bsgs(g, x, mod).value();
        x = power(g, congruence_equation(2, x, mod - 1), mod);
        return std::min(x, mod - x);
    }
    template<int64_t Mod, int64_t G> struct polynomial : private std::vector<int> {
        using std::vector<int>::vector;
        using std::vector<int>::operator[];
        using std::vector<int>::begin;
        using std::vector<int>::end;
        using std::vector<int>::size;
        using std::vector<int>::resize;
        auto &operator*=(int64_t k) {
            for (auto &x: *this) { x = x * k % Mod; }
            return *this;
        }
        auto &operator+=(const polynomial<Mod, G> &rhs) {
            if (rhs.size() > size()) { resize(rhs.size()); }
            for (int i = 0; i < rhs.size(); ++i) { (*this)[i] = ((*this)[i] + rhs[i]) % Mod; }
            return *this;
        }
        auto &operator-=(const polynomial<Mod, G> &rhs) {
            if (rhs.size() > size()) { resize(rhs.size()); }
            for (int i = 0; i < rhs.size(); ++i) { (*this)[i] = ((*this)[i] - rhs[i]) % Mod; }
            return *this;
        }
        auto &operator+=(int64_t x) {
            (*this)[0] = ((*this)[0] + x) % Mod;
            return *this;
        }
        auto &operator-=(int64_t x) {
            (*this)[0] = ((*this)[0] - x) % Mod;
            return *this;
        }
        auto &operator*=(polynomial<Mod, G> rhs) {
            auto m = size() + rhs.size() - 1;
            auto n = std::bit_ceil(m);
            *this = dft(modXN(std::move(*this), n), 1);
            rhs = dft(modXN(std::move(rhs), n), 1);
            for (size_t i = 0; i < n; ++i) { (*this)[i] = (int64_t) (*this)[i] * rhs[i] % Mod; }
            *this = modXN(dft(std::move(*this), -1), m);
            return *this;
        }
    };
    template<int64_t Mod, int64_t G> auto dft(polynomial<Mod, G> a, int f) {
        static constexpr auto wn{[]() constexpr {
            constexpr auto len = std::countr_zero(static_cast<uint64_t>(Mod) - 1);
            std::array<std::array<int, len>, 2> wn{};
            for (size_t i = 0; i < len; ++i) {
                wn[0][i] = power(G, (Mod - 1) >> (i + 1), Mod);
                wn[1][i] = inv(wn[0][i], Mod);
            }
            return wn;
        }()};
        auto n = a.size();
        std::vector<int> w(n);
        for (size_t i = 0, j = 0; i < n; ++i) {
            if (i < j) { std::swap(a[i], a[j]); }
            for (size_t l = n >> 1; (j ^= l) < l; l >>= 1)
                ;
        }
        w[0] = 1;
        for (size_t i = 0; 1 << i < n; ++i) {
            for (size_t j = (1 << (i + 1)) - 1; j; --j) { w[j] = j & 1 ? (int64_t) w[j >> 1] * wn[(1 - f) / 2][i] % Mod : w[j >> 1]; }
            for (size_t j = 0; j < n; j += 1 << (i + 1)) {
                auto *p = &a[j], *q = &a[j | 1 << i], *r = &w[0];
                for (size_t k = 0; k < 1 << i; ++k) {
                    auto t = (int64_t) q[k] * r[k];
                    q[k] = (p[k] - t) % Mod;
                    p[k] = (p[k] + t) % Mod;
                }
            }
        }
        if (f == -1) {
            int64_t in = inv(n, Mod);
            for (auto &x: a) { x = x * in % Mod; }
        }
        return a;
    }
    template<int64_t Mod, int64_t G> auto modXN(polynomial<Mod, G> &&p, size_t n) {
        p.resize(n);
        return p;
    }
    template<int64_t Mod, int64_t G> auto modXN(const polynomial<Mod, G> &p, size_t n) {
        polynomial<Mod, G> res(n);
        std::copy(p.begin(), p.begin() + std::min(n, p.size()), res.begin());
        return res;
    }
    template<int64_t Mod, int64_t G> auto divXN(polynomial<Mod, G> &&p, size_t n) {
        std::copy(p.begin() + n, p.end(), p.begin());
        p.resize(p.size() - n);
        return p;
    }
    template<int64_t Mod, int64_t G> auto divXN(const polynomial<Mod, G> &p, size_t n) {
        polynomial res(p.size() - n);
        std::copy(p.begin() + n, p.end(), res.begin());
        return res;
    }
    template<int64_t Mod, int64_t G> auto reverse(polynomial<Mod, G> p) {
        std::reverse(p.begin(), p.end());
        return p;
    }
    template<int64_t Mod, int64_t G> auto operator+(polynomial<Mod, G> lhs, const polynomial<Mod, G> &rhs) { return lhs += rhs; }
    template<int64_t Mod, int64_t G> auto operator-(polynomial<Mod, G> lhs, const polynomial<Mod, G> &rhs) { return lhs -= rhs; }
    template<int64_t Mod, int64_t G> auto operator+(int64_t x, polynomial<Mod, G> p) { return p += x; }
    template<int64_t Mod, int64_t G> auto operator+(polynomial<Mod, G> p, int64_t x) { return p += x; }
    template<int64_t Mod, int64_t G> auto operator-(polynomial<Mod, G> p, int64_t x) { return p -= x; }
    template<int64_t Mod, int64_t G> auto operator-(int64_t x, polynomial<Mod, G> p) { return p -= x; }
    template<int64_t Mod, int64_t G> auto operator*(int64_t x, polynomial<Mod, G> p) { return p *= x; }
    template<int64_t Mod, int64_t G> auto operator*(polynomial<Mod, G> p, int64_t x) { return p *= x; }
    template<int64_t Mod, int64_t G> auto operator*(polynomial<Mod, G> lhs, const polynomial<Mod, G> &rhs) { return lhs *= rhs; }
    template<int64_t Mod, int64_t G> auto inverse(const polynomial<Mod, G> &p) {
        polynomial<Mod, G> res = {static_cast<int>(inv(p[0], Mod))};
        auto n = std::bit_ceil(p.size());
        for (size_t i = 2; i <= n; i <<= 1) {
            auto a = dft(modXN(modXN(p, i), i << 1), 1);
            auto b = dft(modXN(std::move(res), i << 1), 1);
            for (size_t j = 0; j < i << 1; ++j) { b[j] = b[j] * (2 - (int64_t) a[j] * b[j] % Mod) % Mod; }
            res = modXN(dft(std::move(b), -1), i);
        }
        return modXN(std::move(res), p.size());
    }
    template<int64_t Mod, int64_t G> auto derivative(polynomial<Mod, G> p) {
        for (size_t i = 1; i < p.size(); ++i) { p[i - 1] = (int64_t) i * p[i] % Mod; }
        p.resize(p.size() - 1);
        return p;
    }
    template<int64_t Mod, int64_t G> auto integral(polynomial<Mod, G> p) {
        p.resize(p.size() + 1);
        for (size_t i = p.size(); i-- > 0;) { p[i] = ::inverse(static_cast<int64_t>(i), Mod) * p[i - 1] % Mod; }
        p[0] = 0;
        return p;
    }
    template<int64_t Mod, int64_t G> auto log(const polynomial<Mod, G> &p) { return modXN(integral(derivative(p) * inverse(p)), p.size()); }
    template<int64_t Mod, int64_t G> auto exp(const polynomial<Mod, G> &p) {
        polynomial<Mod, G> res = {1};
        auto n = std::bit_ceil(p.size());
        for (size_t i = 2; i <= n; i <<= 1) {
            auto a = dft(modXN(modXN(p, i), i << 1), 1);
            auto b = dft(modXN(res, i << 1), 1);
            auto c = dft(modXN(log(modXN(std::move(res), i)), i << 1), 1);
            for (size_t j = 0; j < i << 1; ++j) { b[j] = (int64_t) b[j] * (1 + a[j] - c[j]) % Mod; }
            res = modXN(dft(std::move(b), -1), i);
        }
        return modXN(std::move(res), p.size());
    }
    template<int64_t Mod, int64_t G> auto pow(const polynomial<Mod, G> &p, int64_t k) { return exp(log(p) * k); }
    template<int64_t Mod, int64_t G> auto sqrt(const polynomial<Mod, G> &p) {
        polynomial<Mod, G> res = {static_cast<int>(quadratic_residue(p[0], Mod, G))};
        constexpr auto inv2 = ::inverse(static_cast<int64_t>(2), Mod);
        auto n = std::bit_ceil(p.size());
        for (size_t i = 2; i <= n; i <<= 1) {
            auto a = dft(modXN(modXN(p, i), i << 1), 1);
            auto b = dft(modXN(res, i << 1), 1);
            auto c = dft(modXN(inverse(modXN(std::move(res), i)), i << 1), 1);
            for (size_t j = 0; j < i << 1; ++j) { b[j] = (b[j] + (int64_t) a[j] * c[j]) % Mod * inv2 % Mod; }
            res = modXN(dft(std::move(b), -1), i);
        }
        return modXN(std::move(res), p.size());
    }
    template<int64_t Mod, int64_t G> auto operator/(const polynomial<Mod, G> &lhs, const polynomial<Mod, G> &rhs) {
        auto n = lhs.size();
        auto m = rhs.size();
        if (n < m) { return polynomial<Mod, G>{0}; }
        auto a = modXN(reverse(lhs), n - m + 1);
        auto b = modXN(reverse(rhs), n - m + 1);
        return reverse(modXN(a * inverse(b), n - m + 1));
    }
    template<int64_t Mod, int64_t G> auto operator%(const polynomial<Mod, G> &lhs, const polynomial<Mod, G> &rhs) { return modXN(lhs - lhs / rhs * rhs, rhs.size() - 1); }
    template<int64_t Mod, int64_t G> auto operator/=(polynomial<Mod, G> &lhs, const polynomial<Mod, G> &rhs) { return lhs = lhs / rhs; }
    template<int64_t Mod, int64_t G> auto operator%=(polynomial<Mod, G> &lhs, const polynomial<Mod, G> &rhs) { return lhs = lhs % rhs; }
    template<int64_t Mod, int64_t G> auto eva_build(size_t p, size_t l, size_t r, const std::vector<int> &x, std::vector<polynomial<Mod, G>> &a) {
        if (l == r) {
            a[p] = {1, l < x.size() ? -x[l] : 0};
            return;
        }
        auto mid = (l + r) >> 1;
        eva_build(p << 1, l, mid, x, a);
        eva_build(p << 1 | 1, mid + 1, r, x, a);
        a[p] = a[p << 1] * a[p << 1 | 1];
    }
    template<int64_t Mod, int64_t G> auto eva_work(size_t p, size_t l, size_t r, const polynomial<Mod, G> &f, std::vector<polynomial<Mod, G>> &a, std::vector<int> &res) {
        if (l == r) {
            if (l < res.size()) { res[l] = f[0]; }
            return;
        }
        size_t mid = (l + r) >> 1;
        auto fsize = f.size();
        auto n = std::bit_ceil(fsize);
        auto x = dft(modXN(f, n), 1);
        auto helper = [n, fsize](polynomial<Mod, G> x, const polynomial<Mod, G> &g) {
            auto b = dft(modXN(g, n), 1);
            for (size_t i = 0; i < n; ++i) { x[i] = (int64_t) x[i] * b[i] % Mod; }
            return divXN(modXN(dft(std::move(x), -1), fsize), g.size() - 1);
        };
        auto lf = helper(x, a[p << 1 | 1]);
        auto rf = helper(x, a[p << 1]);
        eva_work(p << 1, l, mid, lf, a, res);
        eva_work(p << 1 | 1, mid + 1, r, rf, a, res);
    }
    template<int64_t Mod, int64_t G> auto evaluation(const polynomial<Mod, G> &p, const std::vector<int> &x) {
        size_t m = std::max(x.size(), p.size() - 1);
        std::vector<polynomial<Mod, G>> a(m << 2);
        std::vector<int> res(x.size());
        eva_build(1, 0, m - 1, x, a);
        auto f = modXN(reverse(modXN(p, m + 1)) * inverse(a[1]), m + 1);
        eva_work(1, 0, m - 1, f, a, res);
        for (size_t i = 0; i < x.size(); ++i) { res[i] = (p[0] + (int64_t) res[i] * x[i]) % Mod; }
        return res;
    }
    template<int64_t Mod, int64_t G> polynomial<Mod, G> interpolation_work(size_t p, size_t l, size_t r, const std::vector<int> &y, std::vector<polynomial<Mod, G>> &a, const std::vector<int> &b) {
        if (l == r) { return {(int) (y[l] * ::inverse<int64_t>(b[l], Mod) % Mod)}; }
        auto mid = (l + r) >> 1;
        auto lf = interpolation_work(p << 1, l, mid, y, a, b);
        auto rf = interpolation_work(p << 1 | 1, mid + 1, r, y, a, b);
        return lf * reverse(a[p << 1 | 1]) + rf * reverse(a[p << 1]);
    }
    template<int64_t Mod, int64_t G> auto interpolation(const std::vector<int> &x, const std::vector<int> &y) {
        auto n = x.size();
        std::vector<polynomial<Mod, G>> a(n << 2);
        std::vector<int> b(n);
        eva_build(1, 0, n - 1, x, a);
        auto f = derivative(reverse(a[1]));
        auto g = modXN(reverse(modXN(f, n + 1)) * inverse(a[1]), n + 1);
        eva_work(1, 0, n - 1, g, a, b);
        for (int i = 0; i < n; ++i) { b[i] = (f[0] + (int64_t) b[i] * x[i]) % Mod; }
        return interpolation_work(1, 0, n - 1, y, a, b);
    }
}// namespace polynomial
