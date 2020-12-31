bool inverse_matrix(int A[M][M], int B[M][M], int n) {
    static int tmp[M][M];
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            tmp[i][j] = A[i][j];
            B[i][j] = (int) (i == j);
        }
    }
    for (int i = 0; i < n; ++i) {
        int p = -1;
        for (int j = i; j < n; ++j) {
            if (tmp[j][i] != 0) {
                p = j;
                break;
            }
        }
        if (p == -1)
            return false;
        for (int j = i; j < n; ++j) {
            swap(tmp[i][j], tmp[p][j]);
            swap(B[i][j], B[p][j]);
        }
        long long inv = power(tmp[i][i], mod - 2);
        for (int j = 0; j < n; ++j) {
            tmp[i][j] = tmp[i][j] * inv % mod;
            B[i][j] = B[i][j] * inv % mod;
        }
        for (int k = i + 1; k < n; ++k) {
            long long t = tmp[k][i];
            for (int j = 0; j < n; ++j) {
                tmp[k][j] = (tmp[k][j] - t * tmp[i][j]) % mod;
                B[k][j] = (B[k][j] - t * B[i][j]) % mod;
            }
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        for (int k = i - 1; k >= 0; --k) {
            long long t = tmp[k][i];
            for (int j = 0; j < n; ++j) {
                tmp[k][j] = (tmp[k][j] - t * tmp[i][j]) % mod;
                B[k][j] = (B[k][j] - t * B[i][j]) % mod;
            }
        }
    }
    return true;
}
