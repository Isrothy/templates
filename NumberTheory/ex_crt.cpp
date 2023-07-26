#include <optional>
#include <tuple>
std::optional<std::pair<long long, long long>>
ex_crt(long long a1, long long m1, long long a2, long long m2) {
    long long x, y;
    long long d = ex_gcd(m1, m2, x, y);
    if ((a2 - a1) % d != 0) {
        return std::nullopt;
    }
    long long m = m1 / d * m2;
    long long t = multiply(multiply(x, (a2 - a1) / d, m2 / d), m1, m);
    long long a = (a1 + t) % m;
    if (a < 0) {
        a += m;
    }
    return std::make_pair(a, m);
}

