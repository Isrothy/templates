void division(int *A, int *B, int *C, int n, int m) {
    static int a[M], b[M];
    if (n < m) {
        C[0] = 0;
        return;
    }
    reverse_copy(A, A + n, a);
    reverse_copy(B, B + m, b);
    inverse(b, b, n - m + 1);
    multiply(a, b, a, n - m + 1, n - m + 1);
    reverse_copy(a, a + n - m + 1, C);
}

void modular(int *A, int *B, int *C, int *D, int n, int m) {
    static int a[M];
    if (n < m) {
        C[0] = 0;
        copy(A, A + n, D);
        return;
    }
    division(A, B, C, n, m);
    multiply(B, C, a, m - 1, min(m - 1, n - m + 1));
    for (int i = 0; i < m - 1; ++i) {
        D[i] = (A[i] - a[i]) % mod;
    }
    int l = get_length(n + 1);
}

void modular(int *A, int *B, int *D, int n, int m) {
    static int a[M], c[M];
    if (n < m) {
        copy(A, A + n, D);
        return;
    }
    division(A, B, c, n, m);
    multiply(B, c, a, m - 1, min(m - 1, n - m + 1));
    for (int i = 0; i < m - 1; ++i) {
        D[i] = (A[i] - a[i]) % mod;
    }
    int l = get_length(n + 1);
}
