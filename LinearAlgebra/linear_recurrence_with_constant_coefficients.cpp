#include <vector>
long long
linear_recurrence_with_constant_coefficients(std::vector<int> A, std::vector<int> f, int n) {
    size_t k = A.size();
    Polynomial a(k + 1), b(k), c(k);
    for (int i = 0; i < k; ++i) {
        a[i] = -f[k - i];
    }
    a[k] = b[1] = c[0] = 1;
    while (n != 0) {
        if ((n & 1) == 1) {
            c = c * b % a;
        }
        b = b * b % a;
        n >>= 1;
    }
    long long res = 0;
    for (int i = 0; i < k; ++i) {
        res = (res + (long long) c[i] * A[i]) % MOD;
    }
    return res;
}
