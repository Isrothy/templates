int get_length(int n) {
    int m = 1;
    while (m < n) {
        m <<= 1;
    }
    return m;
}

void DFT(int *a, int n, int p) {
    static int w[M];
    for (int i = 0, j = 0; i < n; ++i) {
        if (i < j) {
            swap(a[i], a[j]);
        }
        for (int k = n >> 1; (j ^= k) < k; k >>= 1);
    }
    w[0] = 1;
    for (int i = 1; i < n; i <<= 1) {
        long long wn = power(g, mod - 1 + p * (mod - 1) / (i << 1));
        for (int j = i - 2; 0 <= j; j -= 2) {
            w[j] = w[j >> 1];
            w[j + 1] = w[j] * wn % mod;
        }
        for (int j = 0; j < n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                long long x = (long long) w[k] * q[k];
                q[k] = (p[k] - x) % mod;
                p[k] = (p[k] + x) % mod;
            }
        }
    }
    if (0 < p) {
        return;
    }
    long long inv = power(n, mod - 2);
    for (int i = 0; i < n; ++i) {
        a[i] = (long long) a[i] * inv % mod;
    }
}

void multiply(int *A, int *B, int *C, int n, int m) {
    static int a[M], b[M];
    copy(A, A + n, a);
    copy(B, B + m, b);
    m += n - 1;
    n = get_length(m);
    DFT(a, n, 1);
    DFT(b, n, 1);
    for (int i = 0; i < n; ++i) {
        a[i] = (long long) a[i] * b[i] % mod;
    }
    DFT(a, n, -1);
    copy(a, a + m, C);
    fill(a, a + n, 0);
    fill(b, b + n, 0);
}

void derivative(int *A, int *B, int n) {
    for (int i = 1; i < n; ++i) {
        B[i - 1] = (long long) A[i] * i % mod;
    }
}

void integral(int *A, int *B, int n) {
    for (int i = n; i != 0; --i) {
        B[i] = (long long) A[i - 1] * Inv[i] % mod;
    }
    B[0] = 0;
}

void inverse(int *A, int *B, int m) {
    static int a[M], b[M];
    int n = get_length(m);
    b[0] = power(A[0], mod - 2);
    for (int i = 2; i <= n; i <<= 1) {
        for (int j = 0; j < i; ++j) {
            a[j] = j < m ? A[j] : 0;
        }
        DFT(a, i << 1, 1);
        DFT(b, i << 1, 1);
        for (int j = 0; j < i << 1; ++j) {
            b[j] = b[j] * (2 - (long long) b[j] * a[j] % mod) % mod;
        }
        DFT(b, i << 1, -1);
        fill(b + i, b + i + i, 0);
    }
    copy(b, b + m, B);
    fill(a, a + (n << 1), 0);
    fill(b, b + (n << 1), 0);
}

void logarithm(int *A, int *B, int m) {
    static int a[M], b[M];
    int n = get_length(m * 2);
    derivative(A, a, n);
    inverse(A, b, m);
    DFT(a, n, 1);
    DFT(b, n, 1);
    for (int i = 0; i < n; ++i) {
        a[i] = (long long) a[i] * b[i] % mod;
    }
    DFT(a, n, -1);
    B[0] = 0;
    integral(a, B, m - 1);
    fill(a, a + n, 0);
    fill(b, b + n, 0);
}

void exponential(int *A, int *B, int m) {
    static int a[M], b[M], c[M];
    int n = get_length(m);
    b[0] = 1;
    for (int i = 2; i <= n; i <<= 1) {
        logarithm(b, c, i);
        for (int j = 0; j < i; ++j) {
            a[j] = j < m ? A[j] : 0;
        }
        DFT(a, i << 1, 1);
        DFT(b, i << 1, 1);
        DFT(c, i << 1, 1);
        for (int j = 0; j < i << 1; ++j) {
            b[j] = (long long) b[j] * (1 + a[j] - c[j]) % mod;
        }
        DFT(b, i << 1, -1);
        fill(b + i, b + i + i, 0);
    }
    copy(b, b + m, B);
    fill(a, a + (n << 1), 0);
    fill(b, b + (n << 1), 0);
    fill(c, c + (n << 1), 0);
}

void square_root(int *A, int *B, int m) {
    static int a[M], b[M], c[M];
    int n = get_length(m);
    b[0] = quadratic_residue(A[0], mod);
    for (int i = 2; i <= n; i <<= 1) {
        for (int j = 0; j < i; ++j) {
            a[j] = j < m ? A[j] : 0;
        }
        inverse(b, c, i);
        DFT(a, i << 1, 1);
        DFT(b, i << 1, 1);
        DFT(c, i << 1, 1);
        for (int j = 0; j < i << 1; ++j) {
            b[j] = Inv[2] * (b[j] + (long long) a[j] * c[j] % mod) % mod;
        }
        DFT(b, i << 1, -1);
        fill(b + i, b + i + i, 0);
    }
    copy(b, b + m, B);
    fill(a, a + (n << 1), 0);
    fill(b, b + (n << 1), 0);
}

void power(int *A, int *B, int m, int k) {
    static int a[M];
    logarithm(A, a, m);
    for (int i = 0; i < m; ++i) {
        a[i] = (long long) a[i] * k % mod;
    }
    exponential(a, a, m);
    copy(a, a + m, B);
}
