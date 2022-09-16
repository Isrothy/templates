void interpolation_work(int p, int l, int r, int *v, int **f, int **g) {
    static int a[M], b[M];
    g[p] = new int[r - l + 1];
    if (l == r) {
        g[p][0] = v[l];
        return;
    }
    int mid = (l + r) >> 1;
    int len = r - l + 1, lenl = mid - l + 2, lenr = r - mid + 1;
    interpolation_work(p << 1, l, mid, v, f, g);
    interpolation_work(p << 1 | 1, mid + 1, r, v, f, g);
    multiply(f[p << 1], g[p << 1 | 1], a, lenl, lenr - 1);
    multiply(f[p << 1 | 1], g[p << 1], b, lenr, lenl - 1);
    for (int i = 0; i < len; ++i) {
        g[p][i] = (a[i] + b[i]) % mod;
    }
}

void interpolation(int *X, int *Y, int *A, int n) {
    static int a[M], b[M], c[M], d[M], v[M], *f[2 * M], *g[2 * M];
    evaluation_build(1, 0, n - 1, X, f);
    derivative(f[1], c, n + 1);
    reverse_copy(c, c + n, a);
    reverse_copy(f[1], f[1] + n + 1, b);
    inverse(b, b, n + 1);
    multiply(a, b, a, n, n + 1);
    a[n - 1] = 0;
    reverse(a, a + n - 1);
    evaluation_work(1, 0, n - 1, a, d, f);
    for (int i = 0; i < n; ++i) {
        d[i] = ((long long) d[i] * X[i] + c[0]) % mod;
        v[i] = Y[i] * power(d[i], mod - 2) % mod;
    }
    interpolation_work(1, 0, n - 1, v, f, g);
    copy(g[1], g[1] + n, A);
    free_segment_tree(1, 0, n - 1, f);
    free_segment_tree(1, 0, n - 1, g);
}
