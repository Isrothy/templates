int homogeneous_linear_recursion_with_constant_coefficients(int *A, int *f, int k, int n) {
    static int a[M], b[M], c[M], d[M];
    for (int i = 0; i < k; ++i) {
        a[i] = -f[k - i];
        b[i] = c[i] = 0;
    }
    a[k] = b[1] = c[0] = 1;
    while (n != 0) {
        if ((n & 1) == 1) {
            multiply(b, c, c, k, k);
            modular(c, a, d, c, 2 * k - 1, k + 1);
        }
        multiply(b, b, b, k, k);
        modular(b, a, d, b, 2 * k - 1, k + 1);
        n >>= 1;
    }
    long long res = 0;
    for (int i = 0; i < k; ++i) {
        res = (res + (long long) c[i] * A[i]) % mod;
    }
    return res;
}
