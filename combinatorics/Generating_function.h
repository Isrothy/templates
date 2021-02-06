long long power(long long x, int k) {
    long long res = 1;
    while (k != 0) {
        if ((k & 1) == 1) {
            res = res * x % mod;
        }
        x = x * x % mod;
        k >>= 1;
    }
    return res;
}

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
    if (0 < p)
        return;
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
    for (int i = 0; i < n; ++i) {
        a[i] = b[i] = 0;
    }
}

void inverse(int *A, int *B, int m) {
    static int a[M], b[M];
    int n = get_length(m);
    b[0] = power(A[0], mod - 2);
    for (int i = 2; i <= n; i <<= 1) {
        for (int j = 0; j < i; ++j)
            a[j] = j < m ? A[j] : 0;
        DFT(a, i << 1, 1);
        DFT(b, i << 1, 1);
        for (int j = 0; j < i << 1; ++j){
            b[j] = b[j] * (2 - (long long) b[j] * a[j] % mod) % mod;
        }
        DFT(b, i << 1, -1);
        for (int j = 0; j < i; ++j) {
            b[j + i] = 0;
        }
    }
    copy(b, b + m, B);
    for (int i = 0; i < n << 1; ++i) {
        a[i] = b[i] = 0;
    }
}

void logarithm(int *A, int *B, int m) {
    static int a[M], b[M];
    int n = get_length(m * 2);
    for (int i = 1; i < n; ++i) {
        a[i - 1] = (long long) i * A[i] % mod;
    }
    inverse(A, b, m);
    DFT(a, n, 1);
    DFT(b, n, 1);
    for (int i = 0; i < n; ++i) {
        a[i] = (long long) a[i] * b[i] % mod;
    }
    DFT(a, n, -1);
    B[0] = 0;
    for (int i = 1; i < m; ++i) {
        B[i] = (long long) Inv[i] * a[i - 1] % mod;
    }
    for (int i = 0; i < n; ++i) {
        a[i] = b[i] = 0;
    }
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
        for (int j = 0; j < i; ++j) {
            b[i + j] = 0;
        }
    }
    copy(b, b + m, B);
    for (int i = 0; i < n << 1; ++i) {
        a[i] = b[i] = c[i] = 0;
    }
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
        for (int j = 0; j < i; ++j) {
            b[i + j] = 0;
        }
    }
    copy(b, b + m, B);
    for (int i = 0; i < n << 1; ++i) {
        a[i] = b[i] = 0;
    }
}

void power(int *A, int *B, int m, int k) {
    static int a[M], b[M];
    logarithm(A, a, m);
    for (int i = 0; i < m; ++i) {
        a[i] = (long long) a[i] * k % mod;
    }
    exponential(a, a, m);
    copy(a, a + m, B);
}

void division(int *A, int *B, int *C, int n, int m) {
    static int a[M], b[M];
    if (n < m) {
        C[0] = 0;
        return;
    }
    int l = get_length(n << 1);
    copy(A, A + n, a);
    reverse(a, a + n);
    copy(B, B + m, b);
    reverse(b, b + m);
    inverse(b, b, n - m + 1);
    DFT(a, l, 1);
    DFT(b, l, 1);
    for (int i = 0; i < l; ++i) {
        a[i] = (long long) a[i] * b[i] % mod;
    }
    DFT(a, l, -1);
    for (int i = 0; i <= n - m; ++i) {
        C[i] = a[n - m - i];
    }
    for (int i = 0; i < l; ++i) {
        a[i] = b[i] = 0;
    }
}

void modular(int *A, int *B, int *C, int *D, int n, int m) {
    static int a[M];
    if (n < m) {
        C[0] = 0;
        copy(A, A + n, D);
        return;
    }
    division(A, B, C, n, m);
    multiply(B, C, a, m, n - m + 1);
    for (int i = 0; i < m - 1; ++i) {
        D[i] = (A[i] - a[i]) % mod;
    }
    int l = get_length(n + 1);
    for (int i = 0; i < l; ++i) {
        a[i] = 0;
    }
}

void modular(int *A, int *B, int *D, int n, int m) {
    static int a[M], c[M];
    if (n < m) {
        copy(A, A + n, D);
        return;
    }
    division(A, B, c, n, m);
    multiply(B, c, a, m, n - m + 1);
    for (int i = 0; i < m - 1; ++i) {
        D[i] = (A[i] - a[i]) % mod;
    }
    int l = get_length(n + 1);
    for (int i = 0; i < l; ++i) {
        a[i] = 0;
    }
}
