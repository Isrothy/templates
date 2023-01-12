void FWT_or(int *a, int n) {
    for (int i = 1; i < 1 << n; i <<= 1)
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                q[k] = (q[k] + p[k]) % mod;
            }
        }
}
void IFWT_or(int *a, int n) {
    for (int i = 1; i < 1 << n; i <<= 1)
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                q[k] = (q[k] - p[k]) % mod;
            }
        }
}
void FWT_and(int *a, int n) {
    for (int i = 1; i < 1 << n; i <<= 1)
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                p[k] = (p[k] + q[k]) % mod;
            }
        }
}
void IFWT_and(int *a, int n) {
    for (int i = 1; i < 1 << n; i <<= 1)
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                p[k] = (p[k] - q[k]) % mod;
            }
        }
}
void FWT_xor(int *a, int n) {
    for (int i = 1; i < 1 << n; i <<= 1)
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                int x = p[k], y = q[k];
                p[k] = (x + y) % mod;
                q[k] = (x - y) % mod;
            }
        }
}
void IFWT_xor(int *a, int n) {
    int w = (mod + 1) / 2;
    for (int i = 1; i < 1 << n; i <<= 1)
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                int x = p[k], y = q[k];
                p[k] = (long long) (x + y) * w % mod;
                q[k] = (long long) (x - y) * w % mod;
            }
        }
}
