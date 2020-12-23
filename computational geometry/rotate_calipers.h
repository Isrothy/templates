double rotate_calipers(Point *A, int n) {
    double ans = 0;
    for (int i = 0, j = 1; i < n; ++i) {
        while (det(A[i] - A[j], A[(i + 1) % n] - A[j]) < det(A[i] - A[(j + 1) % n], A[(i + 1) % n] - A[(j + 1) % n]))
            j = (j + 1) % n;
        ans = max(ans, (A[i] - A[j]).len());
    }
    return ans;
}
