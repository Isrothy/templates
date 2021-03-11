void free_segment_tree(int p, int l, int r, int **f) {
    delete f[p];
    if (l == r) {
        return;
    }
    int mid = (l + r) >> 1;
    free_segment_tree(p << 1, l, mid, f);
    free_segment_tree(p << 1 | 1, mid + 1, r, f);
}

void evaluation_build(int p, int l, int r, int *B, int **f) {
    f[p] = new int[r - l + 2];
    if (l == r) {
        f[p][0] = -B[l];
        f[p][1] = 1;
        return;
    }
    int mid = (l + r) >> 1;
    evaluation_build(p << 1, l, mid, B, f);
    evaluation_build(p << 1 | 1, mid + 1, r, B, f);
    multiply(f[p << 1], f[p << 1 | 1], f[p], mid - l + 2, r - mid + 1);
}

void evaluation_work(int p, int l, int r, int *a, int *C, int **f) {
    static int b[M], c[M];
    if (l == r) {
        C[l] = a[0];
        return;
    }
    int mid = (l + r) >> 1;
    int len = r - l + 1, lenl = mid - l + 2, lenr = r - mid + 1;
    int *al = new int[lenl];
    int *ar = new int[lenr];
    multiply(a, f[p << 1 | 1], b, len, lenr);
    multiply(a, f[p << 1], c, len, lenl);
    copy(b + lenr - 1, b + len, al);
    copy(c + lenl - 1, c + len, ar);
    evaluation_work(p << 1, l, mid, al, C, f);
    evaluation_work(p << 1 | 1, mid + 1, r, ar, C, f);
    delete[] al;
    delete[] ar;
}

void evaluation(int *A, int *B, int *C, int n, int m) {
    static int a[M], b[M], *f[2 * M];
    int l = max(m, n - 1);
    evaluation_build(1, 0, l - 1, B, f);
    reverse_copy(A, A + n, a);
    reverse_copy(f[1], f[1] + l + 1, b);
    inverse(b, b, l + 1);
    multiply(a, b, a, n, l + 1);
    for (int i = n - 1; i < l; ++i) {
        a[i] = 0;
    }
    reverse(a, a + n - 1);
    evaluation_work(1, 0, l - 1, a, C, f);
    for (int i = 0; i < m; ++i) {
        C[i] = ((long long) C[i] * B[i] + A[0]) % mod;
    }
    free_segment_tree(1, 0, l - 1, f);
}
