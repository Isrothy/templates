struct binary_indexed_tree {
    long long b0[M], b1[M];
    int n;
    void update(int l, int r, int x) {
        for (int i = l; i <= n; i += i & -i) {
            b0[i] -= (long long) (l - 1) * x;
            b1[i] += x;
        }
        for (int i = r + 1; i <= n; i += i & -i) {
            b0[i] += (long long) r * x;
            b1[i] -= x;
        }
    }
    long long query(int i) {
        long long x = 0, y = 0;
        for (int j = i; j != 0; j -= j & -j) {
            x += b0[j];
            y += b1[j];
        }
        return x + y * i;
    }
    long long query(int l, int r) {
        return query(r) - query(l - 1);
    }
    void build(int *A, int n) {
        this->n = n;
        for (int i = 1; i <= n; ++i) {
            b0[i] = b1[i] = 0;
        }
        for (int i = 1; i <= n; ++i) {
            b0[i] += A[i];
            if (i + (i & -i) < M){
                b0[i + (i & -i)] += b0[i];
            }
        }
    }
};
