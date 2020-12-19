void division(int *A, int *B, int *C, int n, int m) {
    static int a[M], b[M];
    int l = get_length(n << 1);
    for (int i = 0; i < n; ++i)
        a[i] = A[n - i - 1];
    for (int i = 0; i < m; ++i)
        b[i] = B[m - i - 1];
    inverse(b, b, n - m + 1);
    DFT(a, l, 1);
    DFT(b, l, 1);
    for (int i = 0; i < l; ++i)
        a[i] = (long long) a[i] * b[i] % mod;
    DFT(a, l, -1);
    for (int i = 0; i <= n - m; ++i)
        C[i] = a[n - m - i];
    for (int i = 0; i < l; ++i)
        a[i] = b[i] = 0;
}

void modular(int *A, int *B, int *C, int *D, int n, int m) {
    static int a[M];
    division(A, B, C, n, m);
    multiply(B, C, a, n, n - m + 1);
    for (int i = 0; i < m - 1; ++i)
        D[i] = (A[i] - a[i]) % mod;
    int l = get_length(n + 1);
    for (int i = 0; i < l; ++i)
        a[i] = 0;
}
