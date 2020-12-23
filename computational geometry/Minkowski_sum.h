int Minkowski_sum(Point *A, Point *B, Point *C, int n, int m) {
    int i = 0, j = 0, k = 0;
    C[k++] = A[0] + B[0];
    while (i < n && j < m) {
        Vector a = A[(i + 1) % n] - A[i], b = B[(j + 1) % m] - B[j];
        if (0 < det(a, b)) {
            C[k] = C[k - 1] + a;
            ++i;
        } else {
            C[k] = C[k - 1] + b;
            ++j;
        }
        ++k;
    }
    while (i < n) {
        C[k] = C[k - 1] + A[(i + 1) % n] - A[i];
        ++i;
        ++k;
    }
    while (j < m) {
        C[k] = C[k - 1] + B[(j + 1) % m] - B[j];
        ++j;
        ++k;
    }
    return k - 1;
}
