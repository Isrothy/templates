void initialize(long long *A, int n) {
    m = 1;
    for (int i = 0; i < n; ++i) {
        m *= A[i];
    }
    for (int i = 0; i < n; ++i) {
        long long mi = m / A[i];
        mt[i] = mi * inverse(mi, A[i]) % m;
    }
}

long long query(long long *B, int n) {
    long long res = 0;
    for (int i = 0; i < n; ++i) {
        res = (res + B[i] * mt[i]) % m;
    }
    return res;
}
