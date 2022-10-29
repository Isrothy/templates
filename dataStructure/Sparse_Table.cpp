struct Sparse_Table {
    int ST[K][M];
    int Log2[M];
    void init(int *a, int n) {
        for (int i = 2; i <= n; ++i) { Log2[i] = Log2[i >> 1] + 1; }
        for (int i = 1; i <= n; ++i) { ST[0][i] = a[i]; }
        for (int k = 1; k < K; ++k) {
            for (int i = 1; i <= n - (1 << k) + 1; ++i) {
                ST[k][i] = max(ST[k - 1][i], ST[k - 1][i + (1 << (k - 1))]);
            }
        }
    }
    int query(int l, int r) {
        int k = Log2[r - l + 1];
        return max(ST[k][l], ST[k][r - (1 << k) + 1]);
    }
};
