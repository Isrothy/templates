void initialize(const long long *A, int n, long long *mt, long long &m) {
    m = 1;
    for (int i = 0; i < n; ++i) {
        m *= A[i];
    }
    for (int i = 0; i < n; ++i) {
        long long mi = m / A[i];
        mt[i] = mi * inverse(mi, A[i]) % m;
    }
}

long long query(const long long *B, int n, long long *mt, long long m) {
    long long res = 0;
    for (int i = 0; i < n; ++i) {
        res = (res + B[i] * mt[i]) % m;
    }
    return res;
}
