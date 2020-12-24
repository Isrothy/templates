int determinant(int A[M][M], int n, int mod) {
    long long res = 1;
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            while (A[j][i] != 0) {
                long long tmp = A[i][i] / A[j][i] % mod;
                for (int k = i; k < n; ++k) {
                    A[i][k] = (A[i][k] + A[j][k] * (mod - tmp)) % mod;
                    swap(A[i][k], A[j][k]);
                }
                res = -res;
            }
        }
        res = res * A[i][i] % mod;
    }
    return res;
}
