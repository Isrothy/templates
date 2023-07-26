#include <functional>
#include <random>

std::random_device rd;
std::mt19937_64 gen(rd());


template<class T> T power(T x, long long k, std::function<T(T, T)> multiply) {
    auto ret = T(1);
    while (k != 0) {
        if (k & 1) {
            ret = multiply(ret, x);
        }
        x = multiply(x, x);
        k >>= 1;
    }
    return ret;
}

long long multiply(long long a, long long b, long long p) {
    b = (b % p + p) % p;
    long long ret = 0;
    while (b != 0) {
        if (b & 1) {
            ret = (ret + a) % p;
        }
        a = (a + a) % p;
        b >>= 1;
    }
    return ret;
}
int Legendre(long long a, long long p) {
    a = (a % p + p) % p;
    if (a == 0) {
        return 0;
    }
    auto t = power<long long>(a, (p - 1) >> 1, [p](long long a, long long b) {
        return multiply(a, b, p);
    });
    return t == 1 ? 1 : -1;
}

std::vector<long long> quadratic_residue(long long a, long long p) {
    auto tmp = Legendre(a, p);
    if (tmp == -1) {
        return {};
    }
    if (tmp == 0) {
        return {0};
    }
    std::uniform_int_distribution<long long> distrib(1LL, p - 1);
    auto mulP = [=](long long x, long long y) {
        return multiply(x, y, p);
    };
    while (true) {
        long long r = distrib(gen);
        long long I2 = (mulP(r, r) - a) % p;
        if (Legendre(I2, p) != 1) {
            struct Complex {
                long long re, im;
                explicit Complex(long long re = 0, long long im = 0) : re(re), im(im) {}
            };
            auto complex_multiply = [=](const Complex &z, const Complex &w) {
                return Complex(
                    (mulP(z.re, w.re) + mulP(mulP(z.im, w.im), I2)) % p,
                    (mulP(z.re, w.im) + mulP(z.im, w.re)) % p
                );
            };
            long long x1 = (power<Complex>(Complex(r, 1), (p + 1) >> 1, complex_multiply)).re;
            long long x2 = (p - x1) % p;
            return {x1, x2};
        }
    }
}
