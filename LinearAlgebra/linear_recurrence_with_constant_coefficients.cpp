#include "../Combinatorics/polynomial.h"
#include <vector>
template<int64_t Mod, int64_t G>
auto linear_recurrence_with_constant_coefficients(const std::vector<int> &A, const std::vector<int> &f, int n) {
    size_t k = A.size();
    polynomial::polynomial<Mod, G> a(k + 1), b(k), c(k);
    for (int i = 0; i < k; ++i) { a[i] = -f[k - i]; }
    a[k] = b[1] = c[0] = 1;
    while (n != 0) {
        if ((n & 1) == 1) { c = c * b % a; }
        b = b * b % a;
        n >>= 1;
    }
    long long res = 0;
    for (int i = 0; i < k; ++i) { res = (res + (long long) c[i] * A[i]) % Mod; }
    return res;
}
