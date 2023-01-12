#include <unordered_map>
#include <vector>

struct Polynomial : public std::vector<int> {
    static const long long MOD = 998244353;
    static const long long g = 3;

    using std::vector<int>::vector;

    static long long gcd(long long a, long long b) {
        return b ? gcd(b, a % b) : a;
    }
    static long long ex_gcd(long long a, long long b, long long &x, long long &y) {
        if (b == 0) {
            x = 1;
            y = 0;
            return a;
        }
        auto res = ex_gcd(b, a % b, y, x);
        y -= a / b * x;
        return res;
    }

    static int BSGS(long long a, long long b) {
        static int P = 30000;
        std::unordered_map<long long, int> Mp;
        b = (b % MOD + MOD) % MOD;
        for (int k = 1; k <= P; ++k) {
            b = b * a % MOD;
            Mp[b] = k;
        }
        long long w = power(a, P), x = 1;
        for (int k = 1;; ++k) {
            x = x * w % MOD;
            if (Mp.count(x)) {
                return k * P - Mp[x];
            }
        }
    }
    static long long inv(long long a, long long p) {
        long long x, y;
        a = (a % p + p) % p;
        ex_gcd(a, p, x, y);
        return x;
    }
    static int congruence_equation(int a, int b, int p) {
        b = (b % p + p) % p;
        int d = gcd(a, p);
        a /= d;
        p /= d;
        b /= d;
        return (b * inv(a, p) % p + p) % p;
    }
    static int quadratic_residue(int x, int mod) {
        if (x == 0) {
            return 0;
        }
        x = BSGS(3, x);
        x = power(3, congruence_equation(2, x, mod - 1));
        return std::min(x, mod - x);
    }
    static size_t adequate_length(size_t n) {
        size_t m = 1;
        while (m < n) {
            m <<= 1;
        }
        return m;
    }
    static long long power(long long x, long long k) {
        long long res = 1;
        while (k != 0) {
            if (k & 1) {
                res = res * x % MOD;
            }
            x = x * x % MOD;
            k >>= 1;
        }
        return res;
    }
    Polynomial DFT(size_t n, int f) const {
        std::vector<int> w(n);
        Polynomial a(n);
        for (int i = 0; i < std::min(n, size()); ++i) {
            a[i] = at(i);
        }
        for (int i = 0, j = 0; i < n; ++i) {
            if (i < j) {
                std::swap(a[i], a[j]);
            }
            for (int l = n >> 1; (j ^= l) < l; l >>= 1)
                ;
        }
        w[0] = 1;
        for (int i = 1; i < n; i <<= 1) {
            long long wn = power(g, MOD - 1 + f * (MOD - 1) / (i << 1));
            for (int j = i - 2; j >= 0; j -= 2) {
                w[j] = w[j >> 1];
                w[j + 1] = w[j] * wn % MOD;
            }
            for (int j = 0; j < n; j += i << 1) {
                auto *p = &a[j], *q = &a[j + i], *r = &w[0];
                for (int k = 0; k < i; ++k) {
                    long long t = (long long) (*q) * (*r);
                    *q = (*p - t) % MOD;
                    *p = (*p + t) % MOD;
                    ++p;
                    ++q;
                    ++r;
                }
            }
        }
        if (f == -1) {
            long long in = inv(n, MOD);
            for (int i = 0; i < n; ++i) {
                a[i] = a[i] * in % MOD;
            }
        }
        return a;
    }
    Polynomial modXN(size_t n) const {
        Polynomial res = *this;
        res.resize(n);
        return res;
    }
    Polynomial divXN(size_t n) const {
        Polynomial res(size() - n);
        for (size_t i = 0; i < res.size(); ++i) {
            res[i] = at(i + n);
        }
        return res;
    }
    Polynomial reverse() const {
        Polynomial res = *this;
        std::reverse(res.begin(), res.end());
        return res;
    }
    Polynomial operator+(const Polynomial &rhs) const {
        Polynomial res = *this;
        res.resize(std::max(res.size(), rhs.size()));
        for (int i = 0; i < rhs.size(); ++i) {
            res[i] = (res[i] + rhs[i]) % MOD;
        }
        return res;
    }
    Polynomial operator-(const Polynomial &rhs) const {
        return *this + (-rhs);
    }
    Polynomial operator-() const {
        return (-1) * (*this);
    }
    Polynomial operator+(const long long x) const {
        Polynomial res = *this;
        res[0] = (res[0] + x) % MOD;
        return res;
    }
    Polynomial operator-(const long long x) const {
        return *this + -x;
    }
    friend Polynomial operator+(const long long x, const Polynomial &rhs) {
        return rhs + x;
    }
    friend Polynomial operator-(const long long x, const Polynomial &rhs) {
        return (rhs * -1) + x;
    }
    Polynomial operator*(const Polynomial &rhs) const {
        size_t m = size() + rhs.size() - 1;
        size_t n = adequate_length(m);
        auto a = DFT(n, 1);
        auto b = rhs.DFT(n, 1);
        for (int i = 0; i < n; ++i) {
            a[i] = (long long) a[i] * b[i] % MOD;
        }
        auto res = a.DFT(n, -1);
        res.resize(m);
        return res;
    }
    Polynomial operator*(const long long k) const {
        Polynomial res = *this;
        for (auto &i: res) {
            i = i * k % MOD;
        }
        return res;
    }
    friend Polynomial operator*(const long long k, const Polynomial &rhs) {
        return rhs * k;
    }
    Polynomial inverse() const {
        Polynomial res = {(int) inv(at(0), MOD)};
        auto n = adequate_length(size());
        for (size_t i = 2; i <= n; i <<= 1) {
            Polynomial a = modXN(i).DFT(i << 1, 1);
            Polynomial b = res.DFT(i << 1, 1);
            for (int j = 0; j < i << 1; ++j) {
                b[j] = b[j] * (2 - (long long) a[j] * b[j] % MOD) % MOD;
            }
            res = b.DFT(i << 1, -1).modXN(i);
        }
        return res.modXN(size());
    }
    Polynomial derivative() const {
        Polynomial res(size() - 1);
        for (int i = 1; i < size(); ++i) {
            res[i - 1] = (long long) at(i) * i % MOD;
        }
        return res;
    }
    Polynomial integral() const {
        Polynomial res(size() + 1);
        for (int i = 0; i < size(); ++i) {
            res[i + 1] = at(i) * inv(i + 1, MOD) % MOD;
        }
        return res;
    }
    Polynomial logarithm() const {
        return (derivative() * inverse()).modXN(size() - 1).integral();
    }
    Polynomial exponent() const {
        Polynomial res = {1};
        auto n = adequate_length(size());
        for (size_t i = 2; i <= n; i <<= 1) {
            auto a = modXN(i).DFT(i << 1, 1);
            auto b = res.DFT(i << 1, 1);
            auto c = res.modXN(i).logarithm().DFT(i << 1, 1);
            for (int j = 0; j < i << 1; ++j) {
                b[j] = (long long) b[j] * (1 + a[j] - c[j]) % MOD;
            }
            res = b.DFT(i << 1, -1).modXN(i);
        }
        return res.modXN(size());
    }
    Polynomial power(long long k) const {
        return (logarithm() * k).exponent();
    }
    Polynomial sqrt() const {
        Polynomial res = {quadratic_residue(at(0), MOD)};
        long long inv2 = inv(2, MOD);
        auto n = adequate_length(size());
        for (int i = 2; i <= n; i <<= 1) {
            auto a = modXN(i).DFT(i << 1, 1);
            auto b = res.DFT(i << 1, 1);
            auto c = res.modXN(i).inverse().DFT(i << 1, 1);
            for (int j = 0; j < i << 1; ++j) {
                b[j] = (b[j] + (long long) a[j] * c[j]) % MOD * inv2 % MOD;
            }
            res = b.DFT(i << 1, -1).modXN(i);
        }
        return res.modXN(size());
    }
    Polynomial operator/(const Polynomial &rhs) const {
        size_t n = size();
        size_t m = rhs.size();
        if (n < m) {
            return {1, 0};
        }
        auto a = reverse().modXN(n - m + 1);
        auto b = rhs.reverse().modXN(n - m + 1).inverse();
        return (a * b).modXN(n - m + 1).reverse();
    }
    Polynomial operator%(const Polynomial &rhs) const {
        size_t m = rhs.size();
        return (*this - *this / rhs * rhs).modXN(m - 1);
    }
    static void eva_build(int p, int l, int r, const std::vector<int> &x, std::vector<Polynomial> &a) {
        if (l == r) {
            a[p] = {1, l < x.size() ? -x[l] : 0};
            return;
        }
        int mid = (l + r) >> 1;
        eva_build(p << 1, l, mid, x, a);
        eva_build(p << 1 | 1, mid + 1, r, x, a);
        a[p] = a[p << 1] * a[p << 1 | 1];
    }
    static void eva_work(int p, int l, int r, const Polynomial &f, std::vector<Polynomial> &a, std::vector<int> &res) {
        if (l == r) {
            if (l < res.size()) {
                res[l] = f[0];
            }
            return;
        }
        int mid = (l + r) >> 1;
        auto helper = [](const Polynomial &f, const Polynomial &g) {
            size_t n = adequate_length(f.size());
            auto a = f.DFT(n, 1);
            auto b = g.DFT(n, 1);
            for (int i = 0; i < n; ++i) {
                a[i] = (long long) a[i] * b[i] % MOD;
            }
            return a.DFT(n, -1).modXN(f.size()).divXN(g.size() - 1);
        };
        auto lf = helper(f, a[p << 1 | 1]);
        auto rf = helper(f, a[p << 1]);
        eva_work(p << 1, l, mid, lf, a, res);
        eva_work(p << 1 | 1, mid + 1, r, rf, a, res);
    }
    std::vector<int> evaluation(const std::vector<int> &x) const {
        size_t m = std::max(x.size(), size() - 1);
        std::vector<Polynomial> a(m << 2);
        std::vector<int> res(x.size());
        eva_build(1, 0, (int) m - 1, x, a);
        auto f = (modXN(m + 1).reverse() * a[1].inverse()).modXN(m + 1);
        eva_work(1, 0, (int) m - 1, f, a, res);
        for (int i = 0; i < x.size(); ++i) {
            res[i] = (at(0) + (long long) res[i] * x[i]) % MOD;
        }
        return res;
    }

    static Polynomial interpolation_work(
        int p, int l, int r, const std::vector<int> &y, std::vector<Polynomial> &a, const std::vector<int> &b
    ) {
        if (l == r) {
            return {(int) (y[l] * inv(b[l], MOD) % MOD)};
        }
        int mid = (l + r) >> 1;
        auto lf = interpolation_work(p << 1, l, mid, y, a, b);
        auto rf = interpolation_work(p << 1 | 1, mid + 1, r, y, a, b);
        return lf * a[p << 1 | 1].reverse() + rf * a[p << 1].reverse();
    }

    static Polynomial interpolation(const std::vector<int> &x, const std::vector<int> &y) {
        size_t n = x.size();
        std::vector<Polynomial> a(n << 2);
        std::vector<int> b(n);
        eva_build(1, 0, (int) n - 1, x, a);
        auto f = a[1].reverse().derivative();
        auto g = (f.modXN(n + 1).reverse() * a[1].inverse()).modXN(n + 1);
        eva_work(1, 0, (int) n - 1, g, a, b);
        for (int i = 0; i < n; ++i) {
            b[i] = (f[0] + (long long) b[i] * x[i]) % MOD;
        }
        return interpolation_work(1, 0, n - 1, y, a, b);
    }
};
