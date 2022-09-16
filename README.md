# Templates for ICPC

## Concrete Mathematics

### polynomials

```cpp
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
```

### MTT

```cpp
struct Cp {
    double Re, Im;

    inline Cp operator+(Cp const &_) const {
        return {Re + _.Re, Im + _.Im};
    }

    inline Cp operator-(Cp const &_) const {
        return {Re - _.Re, Im - _.Im};
    }

    inline Cp operator*(Cp const &_) const {
        return {Re * _.Re - Im * _.Im, Re * _.Im + Im * _.Re};
    }

    inline Cp operator*(double _) const {
        return {Re * _, Im * _};
    }
};


void DFT(Cp *a, int n, int p) {
    static Cp w[M];
    for (int i = 0, j = 0; i < n; ++i) {
        if (i < j) {
            swap(a[i], a[j]);
        }
        for (int k = n >> 1; (j ^= k) < k; k >>= 1);
    }
    w[0] = {1, 0};
    for (int i = 1; i < n; i <<= 1) {
        Cp wn = (Cp) {cos(pi / i), sin(pi / i)};
        for (int j = i - 2; j >= 0; j -= 2) {
            w[j] = w[j >> 1];
            w[j + 1] = wn * w[j];
        }
        for (int j = 0; j < n; j += i << 1) {
            Cp *p = a + j, *q = a + j + i;
            for (int k = 0; k < i; ++k) {
                Cp x = q[k] * w[k];
                q[k] = p[k] - x;
                p[k] = p[k] + x;
            }
        }
    }
    if (0 < p) {
        return;
    }
    int inv = 1.0 / n;
    for (int i = 0; i < n; ++i){
        a[i] = a[i] * inv;
    }
}

void multiply(int *A, int *B, int *C, int n, int m, int mod) {
    static Cp a[M], b[M], c[M], d[M], w[M];
    for (int i = 0; i < n; ++i){
        a[i] = {(double) (A[i] & 32767), (double) (A[i] >> 15)};
    }
    for (int i = 0; i < m; ++i){
        b[i] = {(double) (B[i] & 32767), (double) (B[i] >> 15)};
    }
    int l = get_length(m + n - 1);
    DFT(a, l, 1);
    DFT(b, l, 1);
    for (int i = 0; i < l; ++i) {
        int j = (l - 1) & (l - i);
        c[j] = (Cp) {0.5 * (a[i].Re + a[j].Re), 0.5 * (a[i].Im - a[j].Im)} * b[i];
        d[j] = (Cp) {0.5 * (a[j].Im + a[i].Im), 0.5 * (a[j].Re - a[i].Re)} * b[i];
    }
    DFT(c, l, 1);
    DFT(d, l, 1);
    double inv = 1.0 / l;
    for (int i = 0; i < n + m - 1; ++i) {
        long long u = c[i].Re * inv + 0.5, v = c[i].Im * inv + 0.5;
        long long x = d[i].Re * inv + 0.5, y = d[i].Im * inv + 0.5;
        a[i] = b[i] = c[i] = d[i] = (Cp) {0, 0};
        C[i] = (u + ((v + x) << 15) + (y % mod << 30)) % mod;
    }
}
```

### Modular

```cpp
void division(int *A, int *B, int *C, int n, int m) {
    static int a[M], b[M];
    if (n < m) {
        C[0] = 0;
        return;
    }
    reverse_copy(A, A + n, a);
    reverse_copy(B, B + m, b);
    inverse(b, b, n - m + 1);
    multiply(a, b, a, n - m + 1, n - m + 1);
    reverse_copy(a, a + n - m + 1, C);
}

void modular(int *A, int *B, int *C, int *D, int n, int m) {
    static int a[M];
    if (n < m) {
        C[0] = 0;
        copy(A, A + n, D);
        return;
    }
    division(A, B, C, n, m);
    multiply(B, C, a, m - 1, min(m - 1, n - m + 1));
    for (int i = 0; i < m - 1; ++i) {
        D[i] = (A[i] - a[i]) % mod;
    }
    int l = get_length(n + 1);
}

void modular(int *A, int *B, int *D, int n, int m) {
    static int a[M], c[M];
    if (n < m) {
        copy(A, A + n, D);
        return;
    }
    division(A, B, c, n, m);
    multiply(B, c, a, m - 1, min(m - 1, n - m + 1));
    for (int i = 0; i < m - 1; ++i) {
        D[i] = (A[i] - a[i]) % mod;
    }
    int l = get_length(n + 1);
}
```

### FWT

```cpp
void FWT_or(int *a, int n) {
    for (int i = 1; i < 1 << n; i <<= 1) {
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                q[k] = (q[k] + p[k]) % mod;
            }
        }
    }
}

void IFWT_or(int *a, int n) {
    for (int i = 1; i < 1 << n; i <<= 1) {
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                q[k] = (q[k] - p[k]) % mod;
            }
        }
    }
}

void FWT_and(int *a, int n) {
    for (int i = 1; i < 1 << n; i <<= 1) {
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                p[k] = (p[k] + q[k]) % mod;
            }
        }
    }
}

void IFWT_and(int *a, int n) {
    for (int i = 1; i < 1 << n; i <<= 1) {
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                p[k] = (p[k] - q[k]) % mod;
            }
        }
    }
}

void FWT_xor(int *a, int n) {
    for (int i = 1; i < 1 << n; i <<= 1) {
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                int x = p[k], y = q[k];
                p[k] = (x + y) % mod;
                q[k] = (x - y) % mod;
            }
        }
    }
}

void IFWT_xor(int *a, int n) {
    int w = (mod + 1) / 2;
    for (int i = 1; i < 1 << n; i <<= 1) {
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a + j, *q = a + i + j;
            for (int k = 0; k < i; ++k) {
                int x = p[k], y = q[k];
                p[k] = (long long) (x + y) * w % mod;
                q[k] = (long long) (x - y) * w % mod;
            }
        }
    }
}
```

### evaluation

```cpp
void free_segment_tree(int p, int l, int r, int **f) {
    delete f[p];
    if (l == r) {
        return;
    }
    int mid = (l + r) >> 1;
    free_segment_tree(p << 1, l, mid, f);
    free_segment_tree(p << 1 | 1, mid + 1, r, f);
}

void evaluation_build(int p, int l, int r, int *B, int **f) {
    f[p] = new int[r - l + 2];
    if (l == r) {
        f[p][0] = -B[l];
        f[p][1] = 1;
        return;
    }
    int mid = (l + r) >> 1;
    evaluation_build(p << 1, l, mid, B, f);
    evaluation_build(p << 1 | 1, mid + 1, r, B, f);
    multiply(f[p << 1], f[p << 1 | 1], f[p], mid - l + 2, r - mid + 1);
}

void evaluation_work(int p, int l, int r, int *a, int *C, int **f) {
    static int b[M], c[M];
    if (l == r) {
        C[l] = a[0];
        return;
    }
    int mid = (l + r) >> 1;
    int len = r - l + 1, lenl = mid - l + 2, lenr = r - mid + 1;
    int *al = new int[lenl];
    int *ar = new int[lenr];
    multiply(a, f[p << 1 | 1], b, len, lenr);
    multiply(a, f[p << 1], c, len, lenl);
    copy(b + lenr - 1, b + len, al);
    copy(c + lenl - 1, c + len, ar);
    evaluation_work(p << 1, l, mid, al, C, f);
    evaluation_work(p << 1 | 1, mid + 1, r, ar, C, f);
    delete[] al;
    delete[] ar;
}

void evaluation(int *A, int *B, int *C, int n, int m) {
    static int a[M], b[M], *f[2 * M];
    int l = max(m, n - 1);
    evaluation_build(1, 0, l - 1, B, f);
    reverse_copy(A, A + n, a);
    reverse_copy(f[1], f[1] + l + 1, b);
    inverse(b, b, l + 1);
    multiply(a, b, a, n, l + 1);
    for (int i = n - 1; i < l; ++i) {
        a[i] = 0;
    }
    reverse(a, a + n - 1);
    evaluation_work(1, 0, l - 1, a, C, f);
    for (int i = 0; i < m; ++i) {
        C[i] = ((long long) C[i] * B[i] + A[0]) % mod;
    }
    free_segment_tree(1, 0, l - 1, f);
}
```

### interpolation

```cpp
void interpolation_work(int p, int l, int r, int *v, int **f, int **g) {
    static int a[M], b[M];
    g[p] = new int[r - l + 1];
    if (l == r) {
        g[p][0] = v[l];
        return;
    }
    int mid = (l + r) >> 1;
    int len = r - l + 1, lenl = mid - l + 2, lenr = r - mid + 1;
    interpolation_work(p << 1, l, mid, v, f, g);
    interpolation_work(p << 1 | 1, mid + 1, r, v, f, g);
    multiply(f[p << 1], g[p << 1 | 1], a, lenl, lenr - 1);
    multiply(f[p << 1 | 1], g[p << 1], b, lenr, lenl - 1);
    for (int i = 0; i < len; ++i) {
        g[p][i] = (a[i] + b[i]) % mod;
    }
}

void interpolation(int *X, int *Y, int *A, int n) {
    static int a[M], b[M], c[M], d[M], v[M], *f[2 * M], *g[2 * M];
    evaluation_build(1, 0, n - 1, X, f);
    derivative(f[1], c, n + 1);
    reverse_copy(c, c + n, a);
    reverse_copy(f[1], f[1] + n + 1, b);
    inverse(b, b, n + 1);
    multiply(a, b, a, n, n + 1);
    a[n - 1] = 0;
    reverse(a, a + n - 1);
    evaluation_work(1, 0, n - 1, a, d, f);
    for (int i = 0; i < n; ++i) {
        d[i] = ((long long) d[i] * X[i] + c[0]) % mod;
        v[i] = Y[i] * power(d[i], mod - 2) % mod;
    }
    interpolation_work(1, 0, n - 1, v, f, g);
    copy(g[1], g[1] + n, A);
    free_segment_tree(1, 0, n - 1, f);
    free_segment_tree(1, 0, n - 1, g);
}
```

### Euclidean like

```cpp
tuple<long long, long long, long long>
Euclidean_like(long long n, long long a, long long b, long long c) {
    long long x = a / c % mod, y = b / c % mod;
    long long s0 = (n + 1) % mod;
    long long s1 = n * (n + 1) % mod * inv[2] % mod;
    long long s2 = n * (n + 1) % mod * (2 * n + 1) % mod * inv[6] % mod;
    long long m = (a * n + b) / c;
    long long _f, _g, _h, f, g, h;
    if (a == 0) {
        f = y * s0 % mod;
        g = y * s1 % mod;
        h = y * y % mod * s0 % mod;
    } else if (a >= c || b >= c) {
        tie(_f, _g, _h) = Euclidean_like(n, a % c, b % c, c);
        f = (_f + x * s1 + y * s0) % mod;
        g = (_g + x * s2 + y * s1) % mod;
        h = (_h + 2 * y * _f % mod + 2 * x * _g % mod + x * x % mod * s2 % mod
               + 2 * x * y % mod * s1 % mod + y * y % mod * s0 % mod) % mod;
    } else {
        tie(_f, _g, _h) = Euclidean_like(m - 1, c, c - b - 1, a);
        f = (m * n - _f) % mod;
        g = (m * s1 - inv[2] * _h - inv[2] * _f) % mod;
        h = (m * (m + 1) % mod * n - 2 * _g - 2 * _f - f) % mod;
    }
    return make_tuple(f, g, h);
}
```

```cpp
long long power_sum(int m, long long n) {
    n %= mod;
    long long sum = 0, x = n;
    for (int i = m; i >= 0; --i) {
        sum = (sum + P[m][i] * x) % mod;
        x = x * n % mod;
    }
    return sum;
}

void initialize() {
    inv[1] = 1;
    for (int i = 2; i < M; ++i) {
        inv[i] = -mod / i * inv[mod % i] % mod;
    }
    for (int i = 0; i < M; ++i) {
        C[i][0] = 1;
        for (int j = 1; j <= i; ++j) {
            C[i][j] = (C[i - 1][j] + C[i - 1][j - 1]) % mod;
        }
    }
    B[0] = 1;
    for (int i = 1; i < M; ++i) {
        B[i] = 1;
        for (int j = 0; j < i; ++j) {
            B[i] = (B[i] - C[i][j] * B[j] % mod * inv[i - j + 1]) % mod;
        }
    }
    for (int i = 0; i < M - 1; ++i) {
        for (int j = 0; j <= i + 1; ++j) {
            P[i][j] = inv[i + 1] * C[i + 1][j] % mod * B[j] % mod;
        }
    }
}

vector<vector<long long>>
Euclidean_like(long long n, long long a, long long b, long long c, int K) {
    vector<vector<long long>> res;
    long long s[K + 1], power_x[K + 1], power_y[K + 1], power_m[K + 1];
    long long x = a / c % mod, y = b / c % mod;
    long long m = ((__int128) a * n + b) / c;
    for (int i = 0; i <= K; ++i) {
        s[i] = power_sum(i, n);
    }
    power_x[0] = power_y[0] = power_m[0] = 1;
    for (int i = 1; i <= K; ++i) {
        power_x[i] = power_x[i - 1] * x % mod;
        power_y[i] = power_y[i - 1] * y % mod;
        power_m[i] = power_m[i - 1] * (m % mod) % mod;
    }
    res.resize(K + 1);
    for (int i = 0; i <= K; ++i) {
        res[i].resize(K - i + 1);
    }
    if (a == 0) {
        for (int k2 = 0; k2 <= K; ++k2) {
            for (int k1 = 0; k1 <= K - k2; ++k1) {
                res[k1][k2] = power_y[k2] * (k1 == 0 ? (n + 1) : s[k1]) % mod;
            }
        }
    } else if (a >= c) {
        auto tmp = Euclidean_like(n, a % c, b, c, K);
        for (int k2 = 0; k2 <= K; ++k2) {
            for (int j = 0; j <= k2; ++j) {
                long long u = C[k2][j] * power_x[j] % mod;
                for (int k1 = 0; k1 <= K - k2; ++k1) {
                    res[k1][k2] = (res[k1][k2] + u * tmp[k1 + j][k2 - j]) % mod;
                }
            }
        }
    } else if (b >= c) {
        auto tmp = Euclidean_like(n, a, b % c, c, K);
        for (int k2 = 0; k2 <= K; ++k2) {
            for (int j = 0; j <= k2; ++j) {
                long long u = C[k2][j] * power_y[j] % mod;
                for (int k1 = 0; k1 <= K - k2; ++k1) {
                    res[k1][k2] = (res[k1][k2] + u * tmp[k1][k2 - j]) % mod;
                }
            }
        }
    } else {
        auto tmp = Euclidean_like(m - 1, c, c - b - 1, a, K), D = res;
        for (int k1 = 0; k1 <= K; ++k1) {
            res[k1][0] = (k1 == 0 ? (n + 1) : s[k1]) % mod;
        }
        for (int i = 0; i <= K; ++i) {
            for (int j = 0; j <= K - i; ++j) {
                for (int k = 0; k <= i; ++k) {
                    D[i][j] = (D[i][j] + C[i + 1][k] * tmp[k][j]) % mod;
                }
            }
        }
        for (int k2 = 1; k2 <= K; ++k2) {
            for (int k1 = 0; k1 <= K - k2; ++k1) {
                res[k1][k2] = power_m[k2] * s[k1] % mod;
                for (int j = 0; j <= k1; ++j) {
                    res[k1][k2] = (res[k1][k2] - P[k1][j] * D[k2 - 1][k1 + 1 - j]) % mod;
                }
            }
        }
    }
    return res;
}
```

## computational geometry

### 2D basic

```cpp
int dcmp(double x) {
    if (x < -EPS) {
        return -1;
    }
    if (EPS < x) {
        return 1;
    }
    return 0;
}

int dcmp(double x, double y) {
    if (x - y < -EPS) {
        return -1;
    }
    if (EPS < x - y) {
        return 1;
    }
    return 0;
}

struct Point {
    double x, y;

    double len2() const {
        return x * x + y * y;
    }

    double len() const {
        return sqrt(len2());
    }

    Point operator+(Point const &_) const {
        return (Point){x + _.x, y + _.y};
    }

    Point operator-(Point const &_) const {
        return (Point){x - _.x, y - _.y};
    }

    Point operator*(double p) const {
        return (Point){x * p, y * p};
    }

    Point operator/(double p) const {
        return (Point){x / p, y / p};
    }

    bool operator==(Point const &_) const {
        return dcmp(x, _.x) == 0 && dcmp(y, _.y) == 0;
    }

    Point unit() const {
        return *this / len();
    }

    double angle() const {
        return atan2(y, x);
    }

    Point normal() const {
        return (Point){-y, x};
    }

    void read() {
        scanf("%lf%lf", &x, &y);
    }

    void write() const {
        printf("%lf %lf ", x, y);
    }

    void writeln() const {
        printf("%lf %lf\n", x, y);
    }
};

typedef Point Vector;

Point operator*(double p, Point const &_) {
    return (Point){p * _.x, p * _.y};
}

double dot(Vector const &A, Vector const &B) {
    return A.x * B.x + A.y * B.y;
}

double det(Vector const &A, Vector const &B) {
    return A.x * B.y - A.y * B.x;
}

Point middle(Point const &A, Point const &B) {
    return 0.5 * (A + B);
}

double point_line_distance(Point const &P, Point const &A, Point const &B) {
    Vector v1 = B - A, v2 = P - A;
    return fabs(det(v1, v2) / v1.len());
}

int point_on_line_segment(Point const &P, Point const &A, Point const &B) {
    return (int) (dcmp(det(A - P, B - P)) == 0 && dcmp(dot(A - P, B - P)) <= 0);
}

double point_line_segment_distance(Point const &P, Point const &A, Point const &B) {
    if (A == B) {
        return (P - A).len();
    }
    Vector v1 = B - A, v2 = P - A, v3 = P - B;
    if (dcmp(dot(v1, v2)) < 0) {
        return v2.len();
    }
    if (dcmp(dot(v1, v3)) > 0) {
        return v3.len();
    }
    return det(v1, v2) / v1.len();
}

Point projection(Point const &P, Point const &A, Point const &B) {
    Vector v = B - A;
    return A + v * (dot(v, P - A) / v.len2());
}

Point symmetry(Point const &P, Point const &A, Point const &B) {
    return 2 * projection(P, A, B) - P;
}

int intersection(Point const &A, Point const &B, Point const &C, Point const &D, Point &O) {
    if (dcmp(det(B - A, D - C)) == 0) {
        return 0;
    }
    double s1 = det(D - A, C - A);
    double s2 = det(C - B, D - B);
    O = A + (B - A) * (s1 / (s1 + s2));
    return 1;
}

int line_segment_line_segment_intersection(Point const &A, Point const &B, Point const &C,
                                           Point const &D, Point &O) {
    if (!intersection(A, B, C, D, O))
        return 0;
    return (int) (point_on_line_segment(O, A, B) == 1 && point_on_line_segment(O, C, D) == 1);
}

int point_in_triangle(Point const &P, Point const &A, Point const &B, Point const &C) {
    double s1 = det(A - P, B - P) + det(B - P, C - P) + det(C - P, A - P);
    double s2 = det(A - B, C - B);
    return (int) (dcmp(s1, s2) == 0);
}

int circle_line_intersection(Point const &O, double r, Point const &A, Point const &B, Point &P1,
                             Point &P2) {
    Point H = projection(O, A, B);
    double tmp = r * r - (H - O).len2();
    int tmp1 = dcmp(tmp);
    if (tmp1 == -1) {
        return 0;
    } else if (tmp1 == 0) {
        P1 = H;
        return 1;
    } else {
        Vector v = (A - B) / (A - B).len() * sqrt(tmp);
        P1 = H + v;
        P2 = H - v;
        return 2;
    }
}

int circle_line_segment_intersection(Point const &O, double r, Point const &A, Point const &B,
                                     Point &P1, Point &P2) {
    int tmp = circle_line_intersection(O, r, A, B, P1, P2);
    bool o1 = tmp >= 1 && dcmp(dot(P1 - A, P1 - B)) <= 0;
    bool o2 = tmp == 2 && dcmp(dot(P2 - A, P2 - B)) <= 0;
    if (o1 && o2) {
        return 2;
    }
    if (o1) {
        return 1;
    }
    if (o2) {
        P1 = P2;
        return 1;
    }
    return 0;
}

int circle_circle_intersection(Point const &O1, double r1, Point const &O2, double r2, Point &P1,
                               Point &P2) {
    if (O2 == O1) {
        return 0;
    }
    double tmp1 = ((O2 - O1).len2() + r1 * r1 - r2 * r2) / (2 * (O2 - O1).len());
    double tmp2 = r1 * r1 - tmp1 * tmp1;
    int tmp3 = dcmp(tmp2);
    if (tmp3 == -1) {
        return 0;
    } else if (tmp3 == 0) {
        P1 = O1 + (O2 - O1).unit() * tmp1;
        return 1;
    } else {
        Point H = O1 + (O2 - O1).unit() * tmp1;
        Vector v = (O2 - O1).unit().normal() * sqrt(tmp2);
        P1 = H + v;
        P2 = H - v;
        return 2;
    }
}

double circle_point_tangent(Point const &O, double r, Point const &A, Point &P1, Point &P2) {
    double tmp = (O - A).len2();
    if (dcmp(tmp, r * r) == -1){
        return -1;
    }
    Point H = O + (A - O) * (r * r / tmp);
    tmp = r * r - (H - O).len2();
    int tmp1 = dcmp(tmp);
    if (tmp1 == -1) {
        return -1;
    } else if (tmp1 == 0) {
        P1 = H;
        return 0;
    } else {
        Vector v = (A - O).unit().normal() * sqrt(tmp);
        P1 = H + v;
        P2 = H - v;
        return (A - P1).len();
    }
}

double external_co_tangent(Point const &O1, double r1, Point const &O2, double r2, Point &P1,
                           Point &P2, Point &P3, Point &P4) {
    if (r1 < r2){
        return external_co_tangent(O2, r2, O1, r1, P3, P4, P1, P2);
    }
    double res = circle_point_tangent(O1, r1 - r2, O2, P1, P2);
    if (res <= 0){
        return res;
    }
    Vector v1 = (P1 - O1).unit() * r2;
    Vector v2 = (P2 - O1).unit() * r2;
    P1 = P1 + v1;
    P2 = P2 + v2;
    P3 = O2 + v2;
    P4 = O2 + v1;
    return res;
}
```

### Andrew

```cpp
int Andrew(Point *a, int n) {
    static Point b[M];
    sort(a, a + n, [](Point const &A, Point const &B) {
        return A.x == B.x ? A.y < B.y : A.x < B.x;
    });
    int top = 0;
    b[top++] = a[0];
    for (int i = 1; i < n; ++i) {
        while (2 <= top && det(b[top - 1] - b[top - 2], a[i] - b[top - 1]) <= 0) {
            --top;
        }
        b[top++] = a[i];
    }
    int tmp = top;
    for (int i = n - 2; 0 <= i; --i) {
        while (tmp < top && det(b[top - 1] - b[top - 2], a[i] - b[top - 1]) <= 0) {
            --top;
        }
        b[top++] = a[i];
    }
    for (int i = 0; i < top; ++i) {
        a[i] = b[i];
    }
    return top - 1;
}
```

### closest pair of points

```cpp
bool cmp_x(Point const &A, Point const &B) {
    return A.x == B.x ? A.y < B.y : A.x < B.x;
}

bool cmp_y(Point const &A, Point const &B) {
    return A.y == B.y ? A.x < B.x : A.y < B.y;
}

void divide_and_conquer(Point *A, int l, int r, double &res) {
    if (r - l < 3) {
        for (int i = l; i < r; ++i) {
            for (int j = i + 1; j < r; ++j) {
                res = min(res, (A[i] - A[j]).len2());
            }
        }
        sort(A + l, A + r, cmp_y);
        return;
    }
    int mid = (l + r) >> 1;
    double xmid = A[mid].x;
    divide_and_conquer(A, l, mid, res);
    divide_and_conquer(A, mid, r, res);
    inplace_merge(A + l, A + mid, A + r, cmp_y);
    static Point B[M];
    int k = 0;
    for (int i = l; i < r; ++i) {
        if ((xmid - A[i].x) * (xmid - A[i].x) >= res) {
            continue;
        }
        for (int j = k - 1; 0 <= j && (A[i].y - B[j].y) * (A[i].y - B[j].y) < res; --j) {
            res = min(res, (A[i] - B[j]).len2());
        }
        B[k++] = A[i];
    }
}

double closest_pair_of_points(Point *A, int n) {
    sort(A, A + n, cmp_x);
    double res = 1e30;
    divide_and_conquer(A, 0, n, res);
    return sqrt(res);
}
```

### dynamic convex hull

```cpp
struct dynamic_convex_hull {

    struct cmp1 {
        bool operator()(Point A, Point B) {
            return A.x == B.x ? A.y < B.y : A.x < B.x;
        }
    };

    struct cmp2 {
        bool operator()(Point const &A, Point const &B) {
            return A.x == B.x ? A.y > B.y : A.x > B.x;
        }
    };

    set<Point, cmp1> L;
    set<Point, cmp2> U;

    template<typename T>
    bool contain(T &S, Point const &P) {
        if (S.size() < 2) {
            return false;
        }
        typename T::iterator i = S.lower_bound(P);
        return *i == P || (i != S.end() && i != S.begin() && dcmp(det(*prev(i) - P, *i - P)) >= 0);
    }

    template<typename T>
    void insert(T &S, Point const &P) {
        if (contain(S, P)) {
            return;
        }
        S.insert(P);
        auto p = S.lower_bound(P), l_bound = S.begin(), r_bound = prev(S.end());
        if (p != l_bound) {
            for (typename T::iterator i = prev(p), j; i != l_bound; i = j) {
                j = prev(i);
                if (dcmp(det(P - *j, *i - *j)) < 0) {
                    break;
                }
                S.erase(i);
                i = j;
            }
        }
        if (p != r_bound) {
            for (typename T::iterator i = next(p), j; i != r_bound; i = j) {
                j = next(i);
                if (dcmp(det(P - *j, *i - *j)) > 0) {
                    break;
                }
                S.erase(i);
                i = j;
            }
        }
    }

    bool contain(Point const &P) {
        return contain(L, P) && contain(U, P);
    }

    void insert(Point const &P) {
        insert(L, P);
        insert(U, P);
    }
};
```

### half planes intersection

```cpp
int half_planes_intersection(pair<Point, Point> *A, Point *B, int n) {
    static Point t[M];
    static pair<Point, Point> Q[M];
    sort(A, A + n, [](pair<Point, Point> l1, pair<Point, Point> l2) {
        int d = dcmp((l1.second - l1.first).angle(), (l2.second - l2.first).angle());
        return d == 0 ? dcmp(det(l2.first - l1.first, l2.second - l1.first)) > 0 : d < 0;
    });
    int l = 0, r = 0;
    Q[r++] = A[0];
    for (int i = 1; i < n; ++i) {
        if (dcmp(det(Q[r - 1].first - Q[r - 1].second, A[i].first - A[i].second)) == 0){
            continue;
        }
        while (r - l > 1 && dcmp(det(A[i].first - t[r - 1], A[i].second - t[r - 1])) <= 0) {
            --r;
        }
        while (r - l > 1 && dcmp(det(A[i].first - t[l + 1], A[i].second - t[l + 1])) <= 0) {
            ++l;
        }
        intersection(Q[r - 1].first, Q[r - 1].second, A[i].first, A[i].second, t[r]);
        Q[r++] = A[i];
    }
    while (r - l > 1 && dcmp(det(Q[l].first - t[r - 1], Q[l].second - t[r - 1])) <= 0) {
        --r;
    }
    intersection(Q[l].first, Q[l].second, Q[r - 1].first, Q[r - 1].second, t[l]);
    for (int i = l; i < r; ++i) {
        A[i - l] = Q[i];
        B[i - l] = t[i];
    }
    return r - l;
}
```

### Minkowski sum

```cpp
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
```

### rotate calipers

```cpp
double rotate_calipers(Point *A, int n) {
    double ans = 0;
    for (int i = 0, j = 1; i < n; ++i) {
        while (det(A[i] - A[j], A[(i + 1) % n] - A[j])
            < det(A[i] - A[(j + 1) % n], A[(i + 1) % n] - A[(j + 1) % n]))
            j = (j + 1) % n;
        ans = max(ans, (A[i] - A[j]).len());
    }
    return ans;
}
```

### smallest_circle

```cpp
double smallest_circle(Point *a, int n, Point &O) {
    random_shuffle(a, a + n);
    O = a[0];
    double r = 0;
    for (int i = 0; i < n; ++i) {
        if (dcmp(r * r, (O - a[i]).len2()) < 0) {
            O = a[i];
            r = 0;
            for (int j = 0; j < i; ++j) {
                if (dcmp(r * r, (O - a[j]).len2()) < 0) {
                    O = middle(a[i], a[j]);
                    r = (O - a[j]).len();
                    for (int k = 0; k < j; ++k) {
                        if (dcmp(r * r, (O - a[k]).len2()) < 0){
                            r = circumscribed_circle(a[i], a[j], a[k], O);
                        }
                    }
                }
            }
        }
    }
    return r;
}
```

### Adaptive Simpson integral

```cpp
double Simpson(double l, double r) {
    double mid = (l + r) / 2;
    return (r - l) * (F(l) + 4 * F(mid) + F(r)) / 6;
}

double ASR(double l, double r, double tmp) {
    double mid = (l + r) * 0.5;
    double sl = Simpson(l, mid), sr = Simpson(mid, r);
    if (fabs(sl + sr - tmp) < EPS) {
        return sl + sr + (sl + sr - tmp);
    }
    return ASR(l, mid, sl) + ASR(mid, r, sr);
}
```

### 3D basic

```cpp
struct Point {
    double x, y, z;

    Point operator+(const Point &_) const {
        return (Point){x + _.x, y + _.y, z + _.z};
    }

    Point operator-(const Point &_) const {
        return (Point){x - _.x, y - _.y, z - _.z};
    }

    Point operator*(const double &p) const {
        return (Point){x * p, y * p, z * p};
    }

    Point operator/(const double &p) const {
        return (Point){x / p, y / p, z / p};
    }

    bool operator==(const Point &_) const {
        return dcmp(x, _.x) == 0 && dcmp(y, _.y) == 0 && dcmp(z, _.z) == 0;
    }

    double operator/(const Point &_) const {
        if (unit() == _.unit())
            return len() / _.len();
        return -len() / _.len();
    }

    double len2() const {
        return x * x + y * y + z * z;
    }

    double len() const {
        return sqrtl(len2());
    }

    Point unit() const {
        return *this / len();
    }

    void read() {
        scanf("%lf%lf%lf", &x, &y, &z);
    }

    void write() const {
        printf("%lf %lf %lf\n", x, y, z);
    }
};

typedef Point Vector;

Point operator*(const double &p, const Point &A) {
    return (Point){A.x * p, A.y * p, A.z * p};
}

double dot(const Vector &A, const Vector &B) {
    return A.x * B.x + A.y * B.y + A.z * B.z;
}

Vector det(const Vector &A, const Vector &B) {
    return (Vector){A.y * B.z - A.z * B.y, A.z * B.x - A.x * B.z, A.x * B.y - A.y * B.x};
}

double point_line_distance(const Point &P, const Point &A, const Point &B) {
    Vector v1 = B - A, v2 = P - A;
    return det(v1, v2).len() / v1.len();
}

int point_on_segment(const Point &P, const Point &A, const Point &B) {
    return (int) (dcmp(det(A - P, B - P).len2()) == 0 && dcmp(dot(A - P, B - P)) <= 0);
}

double point_segment_distance(const Point &P, const Point &A, const Point &B) {
    if (A == B) {
        return (P - A).len();
    }
    Vector v1 = B - A, v2 = P - A, v3 = P - B;
    if (dcmp(dot(v1, v2)) < 0) {
        return v2.len();
    }
    if (dcmp(dot(v1, v3)) > 0) {
        return v3.len();
    }
    return det(v1, v2).len() / v1.len();
}

double point_plane_diatance(const Point &P, const Point &P0, const Vector &n) {
    return fabs(dot(P - P0, n));
}

Point point_plane_projection(const Point &P, const Point &P0, const Vector &n) {
    return P - n * dot(P - P0, n);
}

int coplaner(const Point &A, const Point &B, const Point &C, const Point &D) {
    return (int) (dcmp(det(det(C - A, D - A), det(C - B, D - B)).len2()) == 0);
}

int line_line_intersection(const Point &A, const Point &B, const Point &C, const Point &D,
                           Point &O) {
    if (!coplaner(A, B, C, D) || dcmp(det(B - A, D - C).len2()) == 0) {
        return 0;
    }
    Vector s1 = det(D - A, C - A);
    Vector s2 = det(C - B, D - B);
    O = A + (B - A) * (s1 / (s1 + s2));
    return 1;
}

int line_segment_line_segment_intersection(const Point &A, const Point &B, const Point &C,
                                           const Point &D, Point &O) {
    if (!line_line_intersection(A, B, C, D, O)) {
        return 0;
    }
    return (int) (point_on_segment(O, A, B) == 1 && point_on_segment(O, C, D) == 1);
}

double line_segment_line_segment_distance(const Point &P1, const Point &P2, const Point &Q1,
                                          const Point &Q2) {
    double A1 = (P1 - P2).len2(), B1 = -dot(P2 - P1, Q2 - Q1), C1 = dot(P1 - P2, P1 - Q1);
    double A2 = -dot(P2 - P1, Q2 - Q1), B2 = (Q1 - Q2).len2(), C2 = dot(P1 - Q1, Q2 - Q1);
    double x1 = B2 * C1 - B1 * C2, y1 = A1 * B2 - A2 * B1;
    double x2 = A2 * C1 - A1 * C2, y2 = A2 * B1 - A1 * B2;

    double s = dcmp(x1) == 0 && dcmp(y1) == 0 ? 0 : x1 / y1;
    double t = dcmp(x2) == 0 && dcmp(y2) == 0 ? 0 : x2 / y2;

    if (0 <= dcmp(s) && dcmp(s, 1) <= 0 && 0 <= dcmp(t) && dcmp(t, 1) <= 0) {
        return ((P1 + s * (P2 - P1)) - (Q1 + t * (Q2 - Q1))).len();
    }

    double a = point_segment_distance(P1, Q1, Q2);
    double b = point_segment_distance(P2, Q1, Q2);
    double c = point_segment_distance(Q1, P1, P2);
    double d = point_segment_distance(Q2, P1, P2);
    return min(min(a, b), min(c, d));
}

int line_plane_intersection(const Point &P1, const Point &P2, const Point P0, const Vector &n,
                            Point &P) {
    double s = dot(n, P2 - P1);
    if (dcmp(s) == 0) {
        return 0;
    }
    Vector v = P2 - P1;
    double t = dot(n, P0 - P1) / s;
    P = P1 + v * t;
    return 1;
}

int point_in_triangle(const Point &P, const Point &A, const Point &B, const Point &C) {
    double S1 = det(A - P, B - P).len() + det(B - P, C - P).len() + det(C - P, A - P).len();
    double S2 = det(A - B, C - B).len();
    return (int) (dcmp(S1, S2) == 0);
}

double point_triangle_distance(const Point &P, const Point &A, const Point &B, const Point &C) {
    Point P0 = point_plane_projection(P, A, det(B - A, C - A).unit());
    if (point_in_triangle(P0, A, B, C)) {
        return (P - P0).len();
    }
    return min(min(point_segment_distance(P, A, B), point_segment_distance(P, A, C)),
               point_segment_distance(P, B, C));
}

bool line_segment_triangle_intersection(const Point &A, const Point &B, const Point &P0,
                                        const Point &P1, const Point &P2, Point &P) {
    if (line_segment_line_segment_intersection(A, B, P0, P1, P)) {
        return true;
    }
    if (line_segment_line_segment_intersection(A, B, P1, P2, P)) {
        return true;
    }
    if (line_segment_line_segment_intersection(A, B, P2, P0, P)) {
        return true;
    }
    Vector n = det(P1 - P0, P2 - P0);
    if (dcmp(dot(n, B - A)) == 0) {
        return false;
    }
    double t = dot(n, P0 - A) / dot(n, B - A);
    if (dcmp(t) < 0 || 0 < dcmp(t, 1)) {
        return false;
    }
    P = A + (B - A) * t;
    return point_in_triangle(P, P0, P1, P2);
}

bool triangle_triangle_intersection(Point *T1, Point *T2) {
    Point P;
    for (int i = 0; i < 3; ++i) {
        if (line_segment_triangle_intersection(T1[i], T1[(i + 1) % 3], T2[0], T2[1], T2[2], P)) {
            return true;
        }
        if (line_segment_triangle_intersection(T2[i], T2[(i + 1) % 3], T1[0], T1[1], T1[2], P)) {
            return true;
        }
    }
    return false;
}

double triangle_triangle_distrance(Point *T1, Point *T2) {
    if (triangle_triangle_intersection(T1, T2)) {
        return 0;
    }
    double res = inf;
    for (int i = 0; i < 3; ++i) {
        res = min(res, point_triangle_distance(T1[i], T2[0], T2[1], T2[2]));
        res = min(res, point_triangle_distance(T2[i], T1[0], T1[1], T1[2]));
        for (int j = 0; j < 3; ++j) {
            res = min(res, line_segment_line_segment_distance(T1[i], T1[(i + 1) % 3], T2[j],
                                                              T2[(j + 1) % 3]));
        }
    }
    return res;
}
```

## data structure

### Treap

```cpp
struct treap {
    int val, size;
    unsigned long long pri;
    treap *ch[2];

    void push_up() {
        size = 1;
        if (ch[0] != nullptr) {
            size += ch[0]->size;
        }
        if (ch[1] != nullptr) {
            size += ch[1]->size;
        }
    }
    void push_down() {}
};

typedef pair<treap *, treap *> ptt;

treap pool[M], *allc = pool;

mt19937_64 mt_rand(time(NULL));

int Size(treap *p) {
    return p == nullptr ? 0 : p->size;
}

treap *rotate(treap *p, bool f) {
    treap *q = p->ch[f];
    p->ch[f] = q->ch[!f];
    q->ch[!f] = p;
    p->push_up();
    return q;
}

treap *insert(treap *p, int x) {
    if (p == nullptr) {
        *allc = (treap) {x, 1, mt_rand(), {nullptr, nullptr}};
        return allcpp;
    }
    bool f = p->val < x;
    p->ch[f] = insert(p->ch[f], x);
    if (p->ch[f]->pri < p->pri) {
        p = rotate(p, f);
    }
    p->push_up();
    return p;
}

treap *erase(treap *p, int x) {
    if (x == p->val) {
        if (p->ch[0] == nullptr || p->ch[1] == nullptr) {
            return p->ch[p->ch[0] == nullptr];
        }
        bool f = p->ch[1]->pri < p->ch[0]->pri;
        p = rotate(p, f);
        p->ch[!f] = erase(p->ch[!f], x);
    } else {
        bool f = p->val < x;
        p->ch[f] = erase(p->ch[f], x);
    }
    p->push_up();
    return p;
}

int index(treap *p, int x) {
    int res = 1;
    while (p != nullptr) {
        if (p->val < x) {
            res += Size(p->ch[0]) + 1;
            p = p->ch[1];
        } else {
            p = p->ch[0];
        }
    }
    return res;
}

int kth(treap *p, int k) {
    for (;;) {
        int s = Size(p->ch[0]);
        if (s + 1 == k) {
            return p->val;
        }
        if (k <= s) {
            p = p->ch[0];
        }
        else {
            k -= s + 1;
            p = p->ch[1];
        }
    }
}

int pre(treap *p, int x) {
    int res = -1;
    while (p != nullptr) {
        if (p->val < x) {
            res = p->val;
            p = p->ch[1];
        } else {
            p = p->ch[0];
        }
    }
    return res;
}

int nxt(treap *p, int x) {
    int res = -1;
    while (p != nullptr) {
        if (x < p->val) {
            res = p->val;
            p = p->ch[0];
        } else {
            p = p->ch[1];
        }
    }
    return res;
}

treap *merge(treap *p, treap *q) {
    if (p == nullptr) {
        return q;
    }
    if (q == nullptr) {
        return p;
    }
    p->push_down();
    q->push_down();
    if (p->pri < q->pri) {
        p->ch[1] = merge(p->ch[1], q);
        p->push_up();
        return p;
    } else {
        q->ch[0] = merge(p, q->ch[0]);
        q->push_up();
        return q;
    }
}

ptt split(treap *p, int k) {
    if (p == nullptr) {
        return make_pair(nullptr, nullptr);
    }
    p->push_down();
    if (k <= Size(p->ch[0])) {
        ptt o = split(p->ch[0], k);
        p->ch[0] = o.second;
        p->push_up();
        return make_pair(o.first, p);
    } else {
        ptt o = split(p->ch[1], k - Size(p->ch[0]) - 1);
        p->ch[1] = o.first;
        p->push_up();
        return make_pair(p, o.second);
    }
}

ptt split_by_value(treap *p, int v) {
    if (p == nullptr) {
        return make_pair(nullptr, nullptr);
    }
    if (v < p->val) {
        ptt o = split_by_value(p->ch[0], v);
        p->ch[0] = o.second;
        p->push_up();
        return make_pair(o.first, p);
    } else {
        ptt o = split_by_value(p->ch[1], v);
        p->ch[1] = o.first;
        p->push_up();
        return make_pair(p, o.second);
    }
}

treap *heuristic_merge(treap *p, treap *q) {
    if (p == nullptr) {
        return q;
    }
    if (q == nullptr) {
        return p;
    }
    if (p->pri < q->pri) {
        swap(p, q);
    }
    ptt o = split_by_value(p, q->val);
    q->ch[0] = heuristic_merge(q->ch[0], o.first);
    q->ch[1] = heuristic_merge(q->ch[1], o.second);
    q->push_up();
    return q;
}
```

### AVL

```cpp
struct avl {
    int val, size, height;
    avl *ch[2];

    void push_up() {
        size = height = 1;
        if (ch[0] != nullptr) {
            size += ch[0]->size;
            height = max(height, ch[0]->height + 1);
        }
        if (ch[1] != nullptr) {
            size += ch[1]->size;
            height = max(height, ch[1]->height + 1);
        }
    }
};

avl pool[M], *allc = pool;

int Size(avl *p) {
    return p == nullptr ? 0 : p->size;
}

int Height(avl *p) {
    return p == nullptr ? 0 : p->height;
}

avl *rotate(avl *p, bool f) {
    avl *q = p->ch[f];
    p->ch[f] = q->ch[!f];
    q->ch[!f] = p;
    p->push_up();
    return q;
}

avl *maintain(avl *p) {
    if (Height(p->ch[0]) - Height(p->ch[1]) == 2) {
        if (Height(p->ch[0]->ch[0]) >= Height(p->ch[0]->ch[1])) {
            p = rotate(p, 0);
        } else {
            p->ch[0] = rotate(p->ch[0], 1);
            p = rotate(p, 0);
        }
    } else if (Height(p->ch[1]) - Height(p->ch[0]) == 2) {
        if (Height(p->ch[1]->ch[1]) >= Height(p->ch[1]->ch[0])) {
            p = rotate(p, 1);
        } else {
            p->ch[1] = rotate(p->ch[1], 0);
            p = rotate(p, 1);
        }
    }
    p->push_up();
    return p;
}

avl *insert(avl *p, int x) {
    if (p == nullptr) {
        *allc = (avl) {x, 1, 1, {nullptr, nullptr}};
        return allcpp;
    }
    bool f = p->val < x;
    p->ch[f] = insert(p->ch[f], x);
    p->push_up();
    return maintain(p);
}

avl *erase(avl *p, int x) {
    if (p->val == x) {
        if (p->ch[0] == nullptr || p->ch[1] == nullptr) {
            return p->ch[p->ch[0] == nullptr];
        }
        avl *q = p->ch[1];
        while (q->ch[0] != nullptr) {
            q = q->ch[0];
        }
        p->val = q->val;
        p->ch[1] = erase(p->ch[1], q->val);
    } else {
        bool f = p->val < x;
        p->ch[f] = erase(p->ch[f], x);
    }
    p->push_up();
    return maintain(p);
}

```

### BIT

```cpp
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
```

### KD Tree

```cpp
const double alpha = 0.7;

struct kd_tree {
    int x, y, val, sum, sz;
    int x_max, x_min, y_max, y_min;
    kd_tree *ch[2];

    void push_up() {
        sz = 1;
        sum = val;
        x_min = x_max = x;
        y_min = y_max = y;
        for (int k = 0; k < 2; ++k) {
            if (ch[k] != nullptr) {
                sz += ch[k]->sz;
                sum += ch[k]->sum;
                x_min = min(x_min, ch[k]->x_min);
                x_max = max(x_max, ch[k]->x_max);
                y_min = min(y_min, ch[k]->y_min);
                y_max = max(y_max, ch[k]->y_max);
            }
        }
    }
};

kd_tree pool[M], *allc = pool;

kd_tree *build(kd_tree **p, int l, int r, int D = 1) {
    if (r < l)
        return nullptr;
    int mid = (l + r) >> 1;
    nth_element(p + l, p + mid, p + r + 1, [=](kd_tree *p, kd_tree *q) {
        if (D == 1) {
            return p->x == q->x ? p->y < q->y : p->x < q->x;
        } else {
            return p->y == q->y ? p->x < q->x : p->y < q->y;
        }
    });
    kd_tree *root = p[mid];
    root->ch[0] = build(p, l, mid - 1, D ^ 1);
    root->ch[1] = build(p, mid + 1, r, D ^ 1);
    root->push_up();
    return root;
}

int query(kd_tree *p, int x1, int y1, int x2, int y2) {
    if (p == nullptr) {
        return 0;
    }
    if (x2 < p->x_min || x1 > p->x_max || y2 < p->y_min || y1 > p->y_max) {
        return 0;
    }
    if (x1 <= p->x_min && p->x_max <= x2 && y1 <= p->y_min && p->y_max <= y2) {
        return p->sum;
    }
    int s = 0;
    if (x1 <= p->x && p->x <= x2 && y1 <= p->y && p->y <= y2) {
        s += p->val;
    }
    return s + query(p->ch[0], x1, y1, x2, y2) + query(p->ch[1], x1, y1, x2, y2);
}

kd_tree *rebuild(kd_tree *p, int D) {
    static kd_tree *Q[M];
    int head = 0, tail = 0;
    Q[tail++] = p;
    while (head < tail) {
        kd_tree *q = Q[head++];
        if (q->ch[0] != nullptr) {
            Q[tail++] = q->ch[0];
        }
        if (q->ch[1] != nullptr) {
            Q[tail++] = q->ch[1];
        }
    }
    return build(Q, 0, tail - 1, D);
}

void insert(kd_tree *&p, int x, int y, int a, int D = 1) {
    if (p == nullptr) {
        p = allcpp;
        *p = (kd_tree){x, y, a, a, 1, x, x, y, y, {nullptr, nullptr}};
        return;
    }
    if (p->x == x && p->y == y) {
        p->val += a;
        p->sum += a;
        return;
    }
    bool f;
    if (D == 1) {
        f = x == p->x ? p->y < y : p->x < x;
    } else {
        f = y == p->y ? p->x < x : p->y < y;
    }
    insert(p->ch[f], x, y, a, D ^ 1);
    p->push_up();
    if (p->sz * alpha < p->ch[f]->sz) {
        p = rebuild(p, D);
    }
}
```

### Li Chao segment tree

```cpp
struct segment {
    long long x1, y1, x2, y2;
    int id;
};

segment S[4 * M];

bool cmp(segment const &A, segment const &B, long long x0) {
    long long h1
        = (B.x2 - B.x1) * (x0 - A.x1) * (A.y2 - A.y1) + A.y1 * (A.x2 - A.x1) * (B.x2 - B.x1);
    long long h2
        = (A.x2 - A.x1) * (x0 - B.x1) * (B.y2 - B.y1) + B.y1 * (B.x2 - B.x1) * (A.x2 - A.x1);
    return h1 == h2 ? A.id > B.id : h1 < h2;
}

void update(int p, int l, int r, int a, int b, segment L) {
    int mid = (l + r) >> 1;
    if (l == a && r == b) {
        if (S[p].id == 0) {
            S[p] = L;
        } else {
            if (cmp(S[p], L, mid)) {
                swap(S[p], L);
            }
            if (l != r) {
                if (cmp(S[p], L, l)) {
                    update(p << 1, l, mid, a, mid, L);
                }
                if (cmp(S[p], L, r)) {
                    update(p << 1 | 1, mid + 1, r, mid + 1, b, L);
                }
            }
        }
        return;
    }
    if (b <= mid) {
        update(p << 1, l, mid, a, b, L);
    } else if (mid < a) {
        update(p << 1 | 1, mid + 1, r, a, b, L);
    } else {
        update(p << 1, l, mid, a, mid, L);
        update(p << 1 | 1, mid + 1, r, mid + 1, b, L);
    }
}

segment query(int p, int l, int r, long long x) {
    if (l == r) {
        return S[p];
    }
    int mid = (l + r) >> 1;
    segment res;
    if (x <= mid) {
        res = query(p << 1, l, mid, x);
    } else {
        res = query(p << 1 | 1, mid + 1, r, x);
    }
    if (S[p].id != 0 && (res.id == 0 || cmp(res, S[p], x))) {
        res = S[p];
    }
    return res;
}
```

### Link Cut Tree

```cpp
struct LCT {
    static const int M = 300005;

    int par[M], ch[M][2];
    bool rev[M];

    bool is_root(int p) {
        return ch[par[p]][0] != p && ch[par[p]][1] != p;
    }

    bool dir(int p) {
        return ch[par[p]][1] == p;
    }

    void update_rev(int p) {
        rev[p] ^= true;
        swap(ch[p][0], ch[p][1]);
    }

    void push_up(int p) {}

    void push_down(int p) {
        if (rev[p]) {
            if (ch[p][0] != 0) {
                update_rev(ch[p][0]);
            }
            if (ch[p][1] != 0) {
                update_rev(ch[p][1]);
            }
            rev[p] = false;
        }
    }

    void rotate(int p) {
        bool f = dir(p);
        int q = par[p], r = ch[p][!f];
        if (!is_root(q)) {
            ch[par[q]][dir(q)] = p;
        }
        par[p] = par[q];
        if (r != 0){
            par[r] = q;
        }
        ch[q][f] = r;
        par[q] = p;
        ch[p][!f] = q;
        push_up(q);
    }

    void splay(int p) {
        static int stk[M];
        int top = 0;
        for (int q = p; !is_root(q); q = par[q]){
            stk[top++] = par[q];
        }
        while (top != 0){
            push_down(stk[--top]);
        }
        push_down(p);
        while (!is_root(p)) {
            int q = par[p];
            if (!is_root(q)){
                rotate(dir(p) == dir(q) ? q : p);
            }
            rotate(p);
        }
        push_up(p);
    }

    void access(int p) {
        int q = 0;
        while (p != 0) {
            splay(p);
            ch[p][1] = q;
            push_up(p);
            q = p;
            p = par[p];
        }
    }

    void make_root(int p) {
        access(p);
        splay(p);
        update_rev(p);
    }

    void link(int p, int q) {
        make_root(p);
        par[p] = q;
    }

    void cut(int p, int q) {
        make_root(p);
        access(q);
        splay(q);
        par[ch[q][0]] = 0;
        ch[q][0] = 0;
        push_up(q);
    }

    int find_root(int p) {
        access(p);
        splay(p);
        while (ch[p][0] != 0) {
            push_down(p);
            p = ch[p][0];
        }
        splay(p);
        return p;
    }

    bool cotree(int p, int q) {
        return find_root(p) == find_root(q);
    }
};

```

### segment tree beats

```cpp
int mx[4 * M], mi[4 * M], second_max[4 * M], second_min[4 * M], cnt_max[4 * M], cnt_min[4 * M],
    len[4 * M];
int lazy_min[4 * M], lazy_max[4 * M], lazy_add[4 * M];

void push_up(int p) {
    sum[p] = sum[p << 1] + sum[p << 1 | 1];
    mx[p] = max(mx[p << 1], mx[p << 1 | 1]);
    mi[p] = min(mi[p << 1], mi[p << 1 | 1]);
    second_max[p] = max(second_max[p << 1], second_max[p << 1 | 1]);
    second_min[p] = min(second_min[p << 1], second_min[p << 1 | 1]);
    cnt_max[p] = cnt_min[p] = 0;
    if (mx[p] == mx[p << 1]) {
        cnt_max[p] += cnt_max[p << 1];
    } else {
        second_max[p] = max(second_max[p], mx[p << 1]);
    }
    if (mx[p] == mx[p << 1 | 1]) {
        cnt_max[p] += cnt_max[p << 1 | 1];
    } else {
        second_max[p] = max(second_max[p], mx[p << 1 | 1]);
    }
    if (mi[p] == mi[p << 1]) {
        cnt_min[p] += cnt_min[p << 1];
    } else {
        second_min[p] = min(second_min[p], mi[p << 1]);
    }
    if (mi[p] == mi[p << 1 | 1]) {
        cnt_min[p] += cnt_min[p << 1 | 1];
    } else {
        second_min[p] = min(second_min[p], mi[p << 1 | 1]);
    }
}

void add(int p, int x) {
    mx[p] += x;
    mi[p] += x;
    sum[p] += (long long) x * len[p];
    if (second_max[p] != -INF) {
        second_max[p] += x;
    }
    if (second_min[p] != INF) {
        second_min[p] += x;
    }
    lazy_add[p] += x;
    if (lazy_min[p] != INF) {
        lazy_min[p] += x;
    }
    if (lazy_max[p] != -INF) {
        lazy_max[p] += x;
    }
}

void check_max(int p, int x) {
    if (x <= mi[p]) {
        return;
    }
    sum[p] += (long long) cnt_min[p] * (x - mi[p]);
    if (mx[p] == mi[p]) {
        mx[p] = lazy_min[p] = x;
    } else if (second_max[p] == mi[p]) {
        second_max[p] = x;
    }
    mi[p] = lazy_max[p] = x;
}

void check_min(int p, int x) {
    if (mx[p] <= x) {
        return;
    }
    sum[p] -= (long long) cnt_max[p] * (mx[p] - x);
    if (mx[p] == mi[p]) {
        mi[p] = lazy_max[p] = x;
    } else if (second_min[p] == mx[p]) {
        second_min[p] = x;
    }
    mx[p] = lazy_min[p] = x;
}

void push_down(int p) {
    if (lazy_add[p] != 0) {
        add(p << 1, lazy_add[p]);
        add(p << 1 | 1, lazy_add[p]);
        lazy_add[p] = 0;
    }
    if (lazy_max[p] != -INF) {
        check_max(p << 1, lazy_max[p]);
        check_max(p << 1 | 1, lazy_max[p]);
        lazy_max[p] = -INF;
    }
    if (lazy_min[p] != INF) {
        check_min(p << 1, lazy_min[p]);
        check_min(p << 1 | 1, lazy_min[p]);
        lazy_min[p] = INF;
    }
}

void build(int p, int l, int r, int *A) {
    len[p] = r - l + 1;
    lazy_add[p] = 0;
    lazy_max[p] = -INF;
    lazy_min[p] = INF;
    if (l == r) {
        sum[p] = mx[p] = mi[p] = A[l];
        cnt_max[p] = cnt_min[p] = 1;
        lazy_add[p] = 0;
        second_max[p] = -INF;
        second_min[p] = INF;
        return;
    }
    int mid = (l + r) >> 1;
    build(p << 1, l, mid, A);
    build(p << 1 | 1, mid + 1, r, A);
    push_up(p);
}

void update_add(int p, int l, int r, int a, int b, int x) {
    if (l == a && r == b) {
        add(p, x);
        return;
    }
    int mid = (l + r) >> 1;
    push_down(p);
    if (b <= mid) {
        update_add(p << 1, l, mid, a, b, x);
    } else if (mid < a) {
        update_add(p << 1 | 1, mid + 1, r, a, b, x);
    } else {
        update_add(p << 1, l, mid, a, mid, x);
        update_add(p << 1 | 1, mid + 1, r, mid + 1, b, x);
    }
    push_up(p);
}

void check_max(int p, int l, int r, int a, int b, int x) {
    if (x <= mi[p]) {
        return;
    }
    if (l == a && r == b && x < second_min[p]) {
        check_max(p, x);
        return;
    }
    int mid = (l + r) >> 1;
    push_down(p);
    if (b <= mid) {
        check_max(p << 1, l, mid, a, b, x);
    } else if (mid < a) {
        check_max(p << 1 | 1, mid + 1, r, a, b, x);
    } else {
        check_max(p << 1, l, mid, a, mid, x);
        check_max(p << 1 | 1, mid + 1, r, mid + 1, b, x);
    }
    push_up(p);
}

void check_min(int p, int l, int r, int a, int b, int x) {
    if (mx[p] <= x) {
        return;
    }
    if (l == a && r == b && second_max[p] < x) {
        check_min(p, x);
        return;
    }
    int mid = (l + r) >> 1;
    push_down(p);
    if (b <= mid) {
        check_min(p << 1, l, mid, a, b, x);
    } else if (mid < a) {
        check_min(p << 1 | 1, mid + 1, r, a, b, x);
    } else {
        check_min(p << 1, l, mid, a, mid, x);
        check_min(p << 1 | 1, mid + 1, r, mid + 1, b, x);
    }
    push_up(p);
}

long long query_sum(int p, int l, int r, int a, int b) {
    if (l == a && r == b) {
        return sum[p];
    }
    int mid = (l + r) >> 1;
    push_down(p);
    if (b <= mid) {
        return query_sum(p << 1, l, mid, a, b);
    }
    if (mid < a) {
        return query_sum(p << 1 | 1, mid + 1, r, a, b);
    }
    return query_sum(p << 1, l, mid, a, mid) + query_sum(p << 1 | 1, mid + 1, r, mid + 1, b);
}

int query_max(int p, int l, int r, int a, int b) {
    if (l == a && r == b) {
        return mx[p];
    }
    int mid = (l + r) >> 1;
    push_down(p);
    if (b <= mid) {
        return query_max(p << 1, l, mid, a, b);
    }
    if (mid < a) {
        return query_max(p << 1 | 1, mid + 1, r, a, b);
    }
    return max(query_max(p << 1, l, mid, a, mid), query_max(p << 1 | 1, mid + 1, r, mid + 1, b));
}

int query_min(int p, int l, int r, int a, int b) {
    if (l == a && r == b) {
        return mi[p];
    }
    int mid = (l + r) >> 1;
    push_down(p);
    if (b <= mid)
        return query_min(p << 1, l, mid, a, b);
    if (mid < a)
        return query_min(p << 1 | 1, mid + 1, r, a, b);
    return min(query_min(p << 1, l, mid, a, mid), query_min(p << 1 | 1, mid + 1, r, mid + 1, b));
}
```

### Sparse Table

```cpp
struct Sparse_Table {
    int ST[K][M];
    int Log2[M];

    void init(int *a, int n) {
        for (int i = 2; i <= n; ++i) {
            Log2[i] = Log2[i >> 1] + 1;
        }
        for (int i = 1; i <= n; ++i) {
            ST[0][i] = a[i];
        }
        for (int k = 1; k < K; ++k) {
            for (int i = 1; i <= n - (1 << k) + 1; ++i) {
                ST[k][i] = max(ST[k - 1][i], ST[k - 1][i + (1 << (k - 1))]);
            }
        }
    }

    int query(int l, int r) {
        int k = Log2[r - l + 1];
        return max(ST[k][l], ST[k][r - (1 << k) + 1]);
    }
};
```

## graph theory

### Euler path

#### directed

```cpp
pair<int, int> stk[N];
int cur[M];
bool vis[N];

bool Euler(int S, int* ans, int m) {
    memset(cur, 0, sizeof cur);
    memset(vis, 0, sizeof vis);
    int sz = 0, top = 0;
    stk[top++] = make_pair(S, 0);
    while (top != 0) {
        pair<int,int> p = stk[top - 1];
        int u = p.first, &i = cur[u];
        while (i < (int) E[u].size() && vis[abs(E[u][i].second)]) {
            ++i;
        }
        if (i < (int) E[u].size()) {
            stk[top++] = E[u][i];
            vis[abs(E[u][i].second)] = true;
            ++i;
        } else {
            ans[++sz] = p.second;
            --top;
        }
    }
    return sz - 1 == m;
}
```

#### undirected

```cpp
bool Euler(int S, int* ans, int m) {
    memset(cur, 0, sizeof cur);
    memset(vis, 0, sizeof vis);
    int sz = 0, top = 0;
    stk[top++] = make_pair(S, 0);
    while (top != 0) {
        pair<int,int> p = stk[top - 1];
        int u = p.first, &i = cur[u];
        while (i < (int) E[u].size() && vis[E[u][i].second]) {
            ++i;
        }
        if (i < (int) E[u].size()) {
            stk[top++] = E[u][i];
            vis[E[u][i].second] = true;
            ++i;
        } else {
            ans[++sz] = p.second;
            --top;
        }
    }
    return sz - 1 == m;
}

```

### 2-SAt

```cpp
namespace Two_SAT {
    vector<int> E[2 * M];
    int stk[M];
    bool mark[2 * M];
    int top;

    void add_clause(int u, bool f1, int v, bool f2) {
        u = u << 1 | f1;
        v = v << 1 | f2;
        E[u ^ 1].push_back(v);
        E[v ^ 1].push_back(u);
    }

    bool dfs(int x) {
        if (mark[x ^ 1]) {
            return false;
        }
        if (mark[x]) {
            return true;
        }
        mark[x] = true;
        stk[top++] = x;
        for (auto y: E[x]) {
            if (!dfs(y))
                return false;
        }
        return true;
    }

    bool check(int n) {
        for (int i = 0; i < 2 * n; i += 2) {
            if (!mark[i] && !mark[i ^ 1]) {
                top = 0;
                if (!dfs(i)) {
                    while (top != 0) {
                        mark[stk[--top]] = false;
                    }
                    if (!dfs(i ^ 1)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
}// namespace Two_SAT
```

```cpp
namespace Two_SAT {
    vector<int> E[2 * M], S[2 * M];
    int dfn[2 * M], low[2 * M], sccno[2 * M], stk[2 * M];
    int deg[2 * M], topono[2 * M], Q[2 * M];
    int dfs_clock, scc_cnt, top;

    void add_clause(int u, bool f1, int v, bool f2) {
        u = u << 1 | f1;
        v = v << 1 | f2;
        E[u ^ 1].push_back(v);
        E[v ^ 1].push_back(u);
    }

    void Tarjan(int u) {
        dfn[u] = low[u] = ++dfs_clock;
        stk[top++] = u;
        for (auto v : E[u]) {
            if (dfn[v] == 0) {
                Tarjan(v);
                low[u] = min(low[u], low[v]);
            } else if (sccno[v] == 0) {
                low[u] = min(low[u], dfn[v]);
            }
        }
        if (dfn[u] == low[u]) {
            ++scc_cnt;
            for (;;) {
                int v = stk[--top];
                sccno[v] = scc_cnt;
                if (v == u) {
                    break;
                }
            }
        }
    }

    bool query(int u) {
        return topono[sccno[u << 1]] < topono[sccno[u << 1 | 1]];
    }

    bool check(int n) {
        for (int u = 0; u < 2 * n; ++u) {
            if (dfn[u] == 0) {
                Tarjan(u);
            }
        }
        for (int u = 0; u < n; ++u) {
            if (sccno[u << 1] == sccno[u << 1 | 1]) {
                return false;
            }
        }
        for (int u = 0; u < 2 * n; ++u) {
            for (auto v : E[u]) {
                if (sccno[u] != sccno[v]) {
                    S[sccno[u]].push_back(sccno[v]);
                    ++deg[sccno[v]];
                }
            }
        }
        int head = 0, tail = 0;
        for (int u = 1; u <= scc_cnt; ++u) {
            if (deg[u] == 0) {
                Q[tail++] = u;
                topono[u] = tail;
            }
        }
        while (head < tail) {
            int u = Q[head++];
            for (auto v : S[u]) {
                if (--deg[v] == 0) {
                    Q[tail++] = v;
                    topono[v] = tail;
                }
            }
        }
        return true;
    }
}
```

### Bipartite graph matching

```cpp
namespace HK {
    vector<int> E[M];
    int dx[M], dy[M], S[M], T[M], Q[M];

    void add_edge(int u, int v) {
        E[u].push_back(v);
    }

    bool BFS(int nl, int nr) {
        memset(dx, 0, sizeof(dx));
        memset(dy, 0, sizeof(dy));
        int head = 0, tail = 0;
        bool res = false;
        for (int u = 1; u <= nl; ++u) {
            if (S[u] == 0)
                Q[tail++] = u;
        }
        while (head < tail) {
            int u = Q[head++];
            for (auto v : E[u]) {
                if (dy[v] == 0) {
                    dy[v] = dx[u] + 1;
                    if (T[v] != 0) {
                        dx[T[v]] = dy[v];
                        Q[tail++] = T[v];
                    } else
                        res = true;
                }
            }
        }
        return res != 0;
    }

    bool DFS(int u) {
        for (auto v : E[u]) {
            if (dy[v] == dx[u] + 1) {
                dy[v] = 0;
                if (T[v] == 0 || DFS(T[v])) {
                    S[u] = v;
                    T[v] = u;
                    return true;
                }
            }
        }
        return false;
    }

    int maximum_matching(int nl, int nr) {
        int res = 0;
        while (BFS(nl, nr)) {
            for (int u = 1; u <= nl; ++u) {
                if (S[u] == 0 && DFS(u))
                    ++res;
            }
        }
        return res;
    }
}

```

```cpp
namespace KM {
    int d[M][M], Q[M], S[M], T[M], lx[M], ly[M], slack[M], pre[M];
    bool visx[M], visy[M];

    void add_edge(int u, int v, int c) {
        d[u][v] = max(d[u][v], c);
        lx[u] = max(lx[u], c);
    }

    int BFS(int u, int n) {
        memset(slack, 0x3f, sizeof slack);
        memset(pre, 0, sizeof pre);
        memset(visx, 0, sizeof visx);
        memset(visy, 0, sizeof visy);
        int head = 0, tail = 0;
        Q[tail++] = u;
        visx[u] = true;
        for (;;) {
            while (head < tail) {
                int x = Q[head++];
                for (int v = 1; v <= n; ++v) {
                    int gap = lx[x] + ly[v] - d[x][v];
                    if (slack[v] < gap or visy[v]) {
                        continue;
                    }
                    pre[v] = x;
                    if (gap == 0) {
                        if (T[v] == 0) {
                            return v;
                        }
                        visy[v] = visx[T[v]] = true;
                        Q[tail++] = T[v];
                    } else {
                        slack[v] = gap;
                    }
                }
            }
            int gap = INF;
            for (int v = 1; v <= n; ++v) {
                if (visy[v] == 0) {
                    gap = min(gap, slack[v]);
                }
            }
            for (int v = 1; v <= n; ++v) {
                if (visx[v]) {
                    lx[v] -= gap;
                }
                if (visy[v]) {
                    ly[v] += gap;
                } else {
                    slack[v] -= gap;
                }
            }
            for (int v = 1; v <= n; ++v) {
                if (visy[v] != 0 || slack[v] != 0) {
                    continue;
                }
                if (T[v] == 0) {
                    return v;
                }
                visy[v] = visx[T[v]] = true;
                Q[tail++] = T[v];
            }
        }
    }

    long long maximum_weight_matching(int n) {
        for (int i = 1; i <= n; ++i) {
            if (S[i] != 0) {
                continue;
            }
            int u = BFS(i, n);
            while (u != 0) {
                T[u] = pre[u];
                swap(u, S[T[u]]);
            }
        }
        long long res = 0;
        for (int i = 1; i <= n; ++i) {
            res += d[i][S[i]];
        }
        return res;
    }
}// namespace KM

```

### Blossom

```cpp
namespace Blossom{
    vector<int> E[M];
    int mate[M], link[M], label[M], fa[M], Q[M];
    int head, tail;

    int find(int x) {
        if (x == fa[x]) {
            return x;
        }
        fa[x] = find(fa[x]);
        return fa[x];
    }
    void add_edge(int u, int v) {
        E[u].push_back(v);
        E[v].push_back(u);
    }
    int LCA(int u, int v) {
        static int last_time[M], time;
        ++time;
        while (last_time[u] != time) {
            if (u != 0) {
                last_time[u] = time;
                u = find(link[mate[u]]);
            }
            swap(u, v);
        }
        return u;
    }

    void blossom(int u, int v, int w) {
        while (find(u) != w) {
            link[u] = v;
            v = mate[u];
            fa[v] = fa[u] = w;
            if (label[v] == 1) {
                label[v] = 2;
                Q[tail++] = v;
            }
            u = link[v];
        }
    }

    bool BFS(int S, int n) {
        head = tail = 0;
        for (int i = 1; i <= n; ++i) {
            fa[i] = i;
            label[i] = 0;
        }
        Q[tail++] = S;
        label[S] = 2;
        while (head < tail) {
            int u = Q[head++];
            for (auto v : E[u]) {
                if (label[v] == 0) {
                    label[v] = 1;
                    link[v] = u;
                    if (mate[v] == 0) {
                        while (u != 0) {
                            u = mate[link[v]];
                            mate[v] = link[v];
                            mate[mate[v]] = v;
                            v = u;
                        }
                        return true;
                    } else {
                        Q[tail++] = mate[v];
                        label[mate[v]] = 2;
                    }
                } else if (label[v] == 2 && find(v) != find(u)) {
                    int w = LCA(u, v);
                    blossom(u, v, w);
                    blossom(v, u, w);
                }
            }
        }
        return false;
    }

    int maximum_matching(int n) {
        int res = 0;
        for (int u = 1; u <= n; ++u) {
            if (mate[u] == 0 && BFS(u, n)) {
                ++res;
            }
        }
        return res;
    }
}
```

### Dinic

```cpp
namespace Dinic {
    struct edge {
        int from, to, cap, flow;
    };
    edge edges[2 * M];
    vector<int> E[M];
    int dis[N], Q[N], cur[M];
    bool vis[N];
    int edge_cnt;

    void add_edge(int u, int v, int c) {
        edges[edge_cnt] = (edge){u, v, c, 0};
        E[u].push_back(edge_cnt++);
        edges[edge_cnt] = (edge){v, u, c, 0};
        E[v].push_back(edge_cnt++);
    }

    void clear() {
        for (int i = 0; i < edge_cnt; ++i) {
            edges[i].flow = 0;
        }
    }

    bool BFS(int S, int T, int n) {
        int head = 0, tail = 0;
        for (int i = 1; i <= n; ++i) {
            vis[i] = false;
        }
        dis[T] = 0;
        vis[T] = true;
        Q[tail++] = T;
        while (head < tail) {
            int u = Q[head++];
            for (int i = 0; i < (int) E[u].size(); ++i) {
                edge e = edges[E[u][i] ^ 1];
                if (e.flow < e.cap && !vis[e.from]) {
                    vis[e.from] = true;
                    dis[e.from] = dis[u] + 1;
                    Q[tail++] = e.from;
                }
            }
        }
        return vis[S];
    }

    int DFS(int u, int T, int a) {
        if (u == T) {
            return a;
        }
        int m = a;
        for (int &i = cur[u]; i < (int) E[u].size(); ++i) {
            edge &e = edges[E[u][i]];
            if (e.flow < e.cap && vis[e.to] && dis[e.to] == dis[u] - 1) {
                int f = DFS(e.to, T, min(a, e.cap - e.flow));
                e.flow += f;
                edges[E[u][i] ^ 1].flow -= f;
                a -= f;
                if (a == 0) {
                    break;
                }
            }
        }
        return m - a;
    }

    long long max_flow(int S, int T, int n) {
        long long flow = 0;
        while (BFS(S, T, n)) {
            for (int i = 1; i <= n; ++i) {
                cur[i] = 0;
            }
            flow += DFS(S, T, INF);
        }
        return flow;
    }
}// namespace Dinic
```

### SSP

```cpp
namespace SSP {
    struct edge {
        int from, to, cap, cost, flow;
    };

    edge edges[2 * M];
    vector<int> E[N];
    int Q[N], dis[N], cur[N];
    bool in_queue[N], vis[N];
    int edge_cnt;

    void add_edge(int u, int v, int cap, int cost) {
        edges[edge_cnt] = (edge){u, v, cap, cost, 0};
        E[u].push_back(edge_cnt++);
        edges[edge_cnt] = (edge){v, u, 0, -cost, 0};
        E[v].push_back(edge_cnt++);
    }

    bool SPFA(int S, int T) {
        int head = 0, tail = 0;
        memset(in_queue, 0, sizeof in_queue);
        memset(dis, 0x3f, sizeof dis);
        dis[T] = 0;
        in_queue[T] = true;
        Q[tail++] = T;
        while (head != tail) {
            int u = Q[head++ % N];
            in_queue[u] = false;
            for (int i = 0; i < (int) E[u].size(); ++i) {
                edge e = edges[E[u][i] ^ 1];
                if (e.flow != e.cap && dis[u] + e.cost < dis[e.from]) {
                    dis[e.from] = dis[u] + e.cost;
                    if (!in_queue[e.from]) {
                        in_queue[e.from] = true;
                        Q[tail++ % N] = e.from;
                    }
                }
            }
        }
        return dis[S] != INF;
    }

    int DFS(int u, int T, int a) {
        if (u == T) {
            return a;
        }
        vis[u] = true;
        int m = a;
        for (int i = 0; i < (int) E[u].size(); ++i) {
            edge &e = edges[E[u][i]];
            if (e.flow < e.cap && !vis[e.to] && dis[e.to] == dis[u] - e.cost) {
                int f = DFS(e.to, T, min(a, e.cap - e.flow));
                e.flow += f;
                edges[E[u][i] ^ 1].flow -= f;
                a -= f;
                if (a == 0) {
                    return m;
                }
            }
        }
        return m - a;
    }

    pair<int, int> minimum_cost_flow(int S, int T) {
        int flow = 0, cost = 0;
        while (SPFA(S, T)) {
            memset(vis, 0, sizeof vis);
            int f = DFS(S, T, INF);
            flow += f;
            cost += f * dis[S];
        }
        return make_pair(flow, cost);
    }
}// namespace SSP

```

### Johnson

```cpp
bool Johnson(int n, vector<pair<int, int>> *E, long long dis[M][M]) {
    queue<int> Q;
    static int h[M], cnt[M];
    static bool in_queue[M];
    for (int u = 1; u <= n; ++u) {
        h[u] = cnt[u] = 0;
        in_queue[u] = true;
        Q.push(u);
    }
    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();
        in_queue[u] = false;
        for (auto e: E[u]) {
            int v = e.first;
            long long d = h[u] + e.second;
            if (d < h[v]) {
                h[v] = d;
                if (!in_queue[v]) {
                    in_queue[v] = true;
                    Q.push(v);
                    if (++cnt[v] == n) {
                        return false;
                    }
                }
            }
        }
    }
    for (int i = 1; i <= n; ++i) {
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> heap;
        for (int u = 1; u <= n; ++u) {
            dis[i][u] = INF;
        }
        dis[i][i] = 0;
        heap.push(make_pair(0, i));
        while (!heap.empty()) {
            pair<int, int> p = heap.top();
            heap.pop();
            int u = p.second;
            if (dis[i][u] < p.first) {
                continue;
            }
            for (auto e: E[u]) {
                int v = e.first;
                long long d = dis[i][u] + h[u] - h[v] + e.second;
                if (d < dis[i][v]) {
                    dis[i][v] = d;
                    heap.push(make_pair(d, e.first));
                }
            }
        }
        for (int u = 1; u <= n; ++u) {
            if (dis[i][u] != INF) {
                dis[i][u] += h[u] - h[i];
            }
        }
    }
    return true;
}

```

### kth shortest path

```cpp
namespace Kth_SSP {
    struct heap {
        int val, to, dis;
        heap *ch[2];
    };

    struct node {
        heap *h;
        int val;

        bool operator>(node const &_) const {
            return val > _.val;
        }
    };

    struct edge {
        int to, id, cost;
    };

    heap pool[17 * M], *allc, *H[N], *tmp[N];
    vector<edge> E[N], R[N];
    int A[N], dis[N], pre[N], deg[N], tree_edge[N];
    int edge_cnt;

    void add_edge(int u, int v, int c) {
        ++edge_cnt;
        E[u].push_back((edge){v, edge_cnt, c});
        R[v].push_back((edge){u, edge_cnt, c});
    }

    heap *merge(heap *p, heap *q) {
        if (p == nullptr) {
            return q;
        }
        if (q == nullptr) {
            return p;
        }
        if (q->val < p->val) {
            swap(p, q);
        }
        heap *r = allcpp;
        *r = *p;
        r->ch[1] = merge(r->ch[1], q);
        if (r->ch[0] == nullptr || r->ch[0]->dis < r->ch[1]->dis) {
            swap(r->ch[0], r->ch[1]);
        }
        r->dis = r->ch[0]->dis + 1;
        return r;
    }

    int Dijkstra(int S, int T, int n) {
        for (int u = 1; u <= n; ++u) {
            dis[u] = INF;
            pre[u] = 0;
        }
        dis[T] = 0;
        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> Q;
        Q.push(make_pair(0, T));
        int sz = 0;
        while (!Q.empty()) {
            pair<double, int> p = Q.top();
            Q.pop();
            int u = p.second;
            if (dis[u] != p.first)
                continue;
            A[sz++] = u;
            for (auto e: R[u]) {
                int v = e.to, d = p.first + e.cost;
                if (d < dis[v]) {
                    dis[v] = d;
                    pre[v] = u;
                    tree_edge[v] = e.id;
                    Q.push(make_pair(d, v));
                }
            }
        }
        return sz;
    }

    int calc(int S, int T, int n, int k) {
        int m = Dijkstra(S, T, n);
        if (dis[S] == INF) {
            return -1;
        }
        if (S != T && --k == 0) {
            return dis[S];
        }
        allc = pool;
        for (int i = 0; i < m; ++i) {
            int u = A[i], sz = 0;
            for (auto e: E[u]) {
                if (dis[e.to] < INF && e.id != tree_edge[u]) {
                    tmp[++sz] = allcpp;
                    *tmp[sz] = (heap){dis[e.to] + e.cost - dis[u], e.to, 1, {nullptr, nullptr}};
                }
            }
            make_heap(tmp + 1, tmp + sz + 1, [](heap *p, heap *q) {
                return p->val > q->val;
            });
            for (int j = sz; 0 < j; --j) {
                if ((j << 1) <= sz) {
                    tmp[j]->ch[0] = tmp[j << 1];
                }
                if ((j << 1 | 1) <= sz) {
                    tmp[j]->ch[1] = tmp[j << 1 | 1];
                }
                tmp[j]->dis = tmp[j]->ch[0] == nullptr ? 1 : tmp[j]->ch[0]->dis + 1;
            }
            H[u] = sz == 0 ? nullptr : tmp[1];
            if (pre[u] != 0) {
                H[u] = merge(H[u], H[pre[u]]);
            }
        }
        priority_queue<node, vector<node>, greater<node>> Q;
        if (H[S] != nullptr) {
            Q.push((node){H[S], dis[S] + H[S]->val});
        }
        while (!Q.empty()) {
            node p = Q.top();
            Q.pop();
            if (--k == 0) {
                return p.val;
            }
            if (p.h->ch[0] != nullptr) {
                Q.push((node){p.h->ch[0], p.val - p.h->val + p.h->ch[0]->val});
            }
            if (p.h->ch[1] != nullptr) {
                Q.push((node){p.h->ch[1], p.val - p.h->val + p.h->ch[1]->val});
            }
            int v = p.h->to;
            if (H[v] != nullptr) {
                Q.push((node){H[v], p.val + H[v]->val});
            }
        }
        return -1;
    }
}// namespace Kth_SSP

```

### DMST

```cpp
namespace Edmonds {
    struct heap {
        int from, c, lazy;
        heap *ch[2];

        void push_down() {
            if (lazy != 0) {
                if (ch[0] != nullptr) {
                    ch[0]->lazy += lazy;
                    ch[0]->c += lazy;
                }
                if (ch[1] != nullptr) {
                    ch[1]->lazy += lazy;
                    ch[1]->c += lazy;
                }
            }
            lazy = 0;
        }
    };

    heap pool[M], *allc = pool, *H[M];
    int fa[M], mark[M];

    heap *merge(heap *p, heap *q) {
        if (p == nullptr) {
            return q;
        }
        if (q == nullptr) {
            return p;
        }
        if (q->c < p->c) {
            swap(p, q);
        }
        p->push_down();
        p->ch[1] = merge(p->ch[1], q);
        swap(p->ch[0], p->ch[1]);
        return p;
    }

    void pop(heap *&p) {
        p->push_down();
        p = merge(p->ch[0], p->ch[1]);
    }

    void add_edge(int u, int v, int c) {
        if (u == v) {
            return;
        }
        *allc = (heap) {u, c, 0, {nullptr, nullptr}};
        H[v] = merge(H[v], allc);
        ++allc;
    }

    int find(int u) {
        if (fa[u] == u) {
            return u;
        }
        fa[u] = find(fa[u]);
        return fa[u];
    }

    int DMST(int n, int root) {
        for (int i = 1; i <= n; ++i) {
            fa[i] = i;
        }
        mark[root] = root;
        int res = 0;
        for (int i = 1; i <= n; ++i) {
            if (mark[i] != 0 || i == root)
                continue;
            int u = i;
            mark[u] = i;
            for (;;) {
                while (H[u] != nullptr && find(H[u]->from) == find(u)) {
                    pop(H[u]);
                }
                if (H[u] == nullptr) {
                    return -1;
                }
                int v = find(H[u]->from);
                res += H[u]->c;
                if (mark[v] == i) {
                    H[u]->lazy -= H[u]->c;
                    pop(H[u]);
                    while (v != u) {
                        int w = find(H[v]->from);
                        H[v]->lazy -= H[v]->c;
                        pop(H[v]);
                        fa[find(v)] = u;
                        H[u] = merge(H[u], H[v]);
                        v = w;
                    }
                } else if (mark[v] == 0) {
                    u = v;
                    mark[u] = i;
                } else {
                    break;
                }
            }
        }
        return res;
    }
}
```

### Domination tree

```cpp
vector<int> E[M], R[M], Q[M];
int par[M], dfn[M], vertices[M], Fa[M], idom[M], sdom[M], val[M];
int dfs_clock;

bool cmp(int u, int v) {
    return dfn[u] < dfn[v];
}

int Find(int u) {
    if (u == Fa[u]) {
        return u;
    }
    int &v = Fa[u], w = Find(Fa[u]);
    if (v != w) {
        if (cmp(sdom[val[v]], sdom[val[u]])) {
            val[u] = val[v];
        }
        v = w;
    }
    return Fa[u];
}

void dfs(int u) {
    dfn[u] = ++dfs_clock;
    vertices[dfs_clock] = u;
    for (auto v : E[u]) {
        if (dfn[v] == 0) {
            par[v] = u;
            dfs(v);
        }
    }
}

void Lengauer_Tarjan(int root) {
    dfs(root);
    for (int u = 1; u <= n; ++u) {
        Fa[u] = sdom[u] = val[u] = u;
    }
    for (int i = dfs_clock; i != 0; --i) {
        int u = vertices[i];
        for (auto v : R[u]) {
            if (dfn[v] == 0) {
                continue;
            }
            Find(v);
            sdom[u] = min(sdom[u], sdom[val[v]], cmp);
        }
        for (auto v : Q[u]) {
            Find(v);
            idom[v] = sdom[v] == sdom[val[v]] ? sdom[v] : val[v];
        }
        if (i != 1) {
            Fa[u] = par[u];
            Q[sdom[u]].push_back(u);
        }
    }
    for (int i = 2; i <= dfs_clock; ++i) {
        int u = vertices[i];
        if (idom[u] != sdom[u]) {
            idom[u] = idom[idom[u]];
        }
    }
}
```

### Circle-square tree

```cpp
vector<int> E[M], R[2 * M];
int dfn[M], low[M], stk[M], bccno[M], A[M];
int n, m, bcc_cnt, dfs_clock, top;

void Tarjan(int u, int fa) {
    dfn[u] = low[u] = ++dfs_clock;
    stk[top++] = u;
    for (auto v: E[u]) {
        if (v == fa) {
            continue;
        }
        if (dfn[v] == 0) {
            Tarjan(v, u);
            low[u] = min(low[u], low[v]);
            if (dfn[u] <= low[v]) {
                ++bcc_cnt;
                int x, w = n + bcc_cnt;
                do {
                    x = stk[--top];
                    R[w].push_back(x);
                } while (x != v);
                R[u].push_back(w);
            }
        } else {
            low[u] = min(low[u], dfn[v]);
        }
    }
}

```

### 3-membered rings

```cpp
int three_membered_rings(vector<int> *E, int n) {
    static long long Rank[M];
    static int vis_time[M];
    static vector<int> F[M];
    for (int u = 1; u <= n; ++u) {
        Rank[u] = (long long) E[u].size() * (n + 1) + u;
        vis_time[u] = 0;
        F[u].clear();
    }
    for (int u = 1; u <= n; ++u) {
        for (auto v: E[u]) {
            if (Rank[u] < Rank[v]) {
                F[u].push_back(v);
            }
        }
    }
    int res = 0;
    for (int u = 1; u <= n; ++u) {
        for (auto v: F[u]) {
            vis_time[v] = u;
        }
        for (auto v: F[u]) {
            for (auto w: F[v]) {
                if (vis_time[w] == u) {
                    ++res;
                }
            }
        }
    }
    return res;
}

```

### 4-membered ring

```cpp
int four_membered_rings(vector<int> *E, int n) {
    static long long Rank[M];
    static int vis_time[M], cnt[M];
    static vector<int> F[M];
    for (int u = 1; u <= n; ++u) {
        Rank[u] = (long long) E[u].size() * (n + 1) + u;
        vis_time[u] = 0;
        F[u].clear();
    }
    for (int u = 1; u <= n; ++u) {
        for (auto v: E[u]) {
            if (Rank[u] < Rank[v]) {
                F[u].push_back(v);
            }
        }
    }
    int res = 0;
    for (int u = 1; u <= n; ++u) {
        for (auto v: E[u]) {
            for (auto w: F[v]) {
                if (Rank[u] < Rank[w]) {
                    if (vis_time[w] < u) {
                        vis_time[w] = u;
                        cnt[w] = 0;
                    }
                    res = (res + cnt[w]) % mod;
                    ++cnt[w];
                }
            }
        }
    }
    return res;
}
```

### Gomory Hu Tree

```cpp
void Gomory_Hu_Tree(int *p, int l, int r, int n) {
    if (l == r) {
        return;
    }
    Dinic::clear();
    int u = p[l], v = p[r];
    add_edge(u, v, Dinic::max_flow(u, v, n));
    Dinic::BFS(u, v, n);
    int i = l + 1, j = r - 1;
    while (i <= j) {
        while (i <= j && !Dinic::vis[p[i]]) {
            ++i;
        }
        while (i <= j && Dinic::vis[p[j]]) {
            --j;
        }
        if (i <= j) {
            swap(p[i], p[j]);
        }
    }
    Gomory_Hu_Tree(p, l, i - 1, n);
    Gomory_Hu_Tree(p, j + 1, r, n);
}
```

### Stoer Wagner

```cpp
int Stoer_Wagner(int d[M][M], int n) {
    static int w[M];
    static bool vis[M], del[M];
    int res = INF;
    for (int u = 1; u <= n; ++u) {
        del[u] = false;
    }
    for (int i = 1; i < n; ++i) {
        for (int u = 1; u <= n; ++u) {
            w[u] = 0;
            vis[u] = false;
        }
        int s = -1, t = -1;
        for (int j = 1; j <= n - i + 1; ++j) {
            int v = -1;
            for (int u = 1; u <= n; ++u) {
                if (!del[u] && !vis[u] && (v == -1 || w[v] < w[u])) {
                    v = u;
                }
            }
            vis[v] = true;
            for (int u = 1; u <= n; ++u) {
                if (!del[u] && !vis[u]) {
                    w[u] += d[u][v];
                }
            }
            s = t;
            t = v;
        }
        res = min(res, w[t]);
        del[t] = true;
        for (int u = 1; u <= n; ++u) {
            d[u][s] += d[u][t];
            d[s][u] += d[t][u];
        }
    }
    return res;
}
```

## Tarjan

### BCC

```cpp
vector<int> E[M];
int dfn[M], low[M], stk[M], bccno[M];
int n, m, bcc_cnt, dfs_clock, top;

void Tarjan(int u, int fa) {
    dfn[u] = low[u] = ++dfs_clock;
    stk[top++] = u;
    for (auto v: E[u]) {
        if (v == fa) {
            continue;
        }
        if (dfn[v] == 0) {
            Tarjan(v, u);
            low[u] = min(low[u], low[v]);
            if (dfn[u] <= low[v]) {
                ++bcc_cnt;
                int x;
                do {
                    x = stk[--top];
                    bccno[x] = bcc_cnt;
                } while (x != v);
            }
        } else {
            low[u] = min(low[u], dfn[v]);
        }
    }
}
```

### EBC

```cpp
vector<pair<int, int>> E[M];
int dfn[M], low[M], ebcno[M];
bool is_bridge[M];
int n, m, dfs_clock, ebc_cnt;

void Tarjan(int u, int pre) {
    dfn[u] = low[u] = ++dfs_clock;
    for (int i = 0; i < (int) E[u].size(); ++i) {
        int v = E[u][i].first, e = E[u][i].second;
        if (e == pre) {
            continue;
        }
        if (dfn[v] == 0) {
            Tarjan(v, e);
            low[u] = min(low[u], low[v]);
            is_bridge[e] = dfn[u] < low[v];
        } else {
            low[u] = min(low[u], dfn[v]);
        }
    }
}

void dfs(int u) {
    ebcno[u] = ebc_cnt;
    for (int i = 0; i < (int) E[u].size(); ++i) {
        int v = E[u][i].first;
        if (!is_bridge[E[u][i].second] && ebcno[v] == 0) {
            dfs(v);
        }
    }
}
```

### SCC

```cpp
vector<int> E[M];
int dfn[M], low[M], sccno[M], stk[M];
int dfs_clock, scc_cnt, top;

void Tarjan(int u) {
    dfn[u] = low[u] = ++dfs_clock;
    stk[top++] = u;
    for (auto v: E[u]) {
        if (dfn[v] == 0) {
            Tarjan(v);
            low[u] = min(low[u], low[v]);
        } else if (sccno[v] == 0) {
            low[u] = min(low[u], dfn[v]);
        }
    }
    if (dfn[u] == low[u]) {
        ++scc_cnt;
        for (;;) {
            int v = stk[--top];
            sccno[v] = scc_cnt;
            if (v == u) {
                break;
            }
        }
    }
}
```

### Steiner Minimum Tree

```cpp
int Steiner_Minimum_Tree(int *s, int k) {
    for (int S = 0; S < 1 << k; ++S) {
        for (int u = 1; u <= n; ++u) {
            dp[S][u] = INF;
        }
    }
    for (int i = 0; i < k; ++i) {
        dp[1 << i][s[i]] = 0;
    }
    for (int S = 1; S < 1 << k; ++S) {
        for (int T = (S - 1) & S; T != 0; --T &= S) {
            for (int u = 1; u <= n; ++u) {
                dp[S][u] = min(dp[S][u], dp[T][u] + dp[S ^ T][u]);
            }
        }
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> Q;
        for (int u = 1; u <= n; ++u) {
            Q.push(make_pair(dp[S][u], u));
        }
        while (!Q.empty()) {
            pair<int, int> p = Q.top();
            Q.pop();
            int u = p.second;
            if (p.first != dp[S][u]) {
                continue;
            }
            for (auto e: E[u]) {
                int v = e.first;
                if (p.first + e.second < dp[S][v]) {
                    dp[S][v] = p.first + e.second;
                    Q.push(make_pair(dp[S][v], v));
                }
            }
        }
    }
    int res = INF;
    for (int u = 1; u <= n; ++u) {
        res = min(res, dp[(1 << k) - 1][u]);
    }
    return res;
}

```

## Linear algebra

### Linear bases

```cpp
struct Bases {
    unsigned long long A[K];

    void insert(unsigned long long x) {
        for (int k = K - 1; k >= 0; --k) {
            if (((x >> k) & 1) == 1) {
                if (A[k] == 0) {
                    A[k] = x;
                    break;
                } else {
                    x ^= A[k];
                }
            }
        }
    }

    unsigned long long maximum_xor_sum(unsigned long long res = 0) {
        for (int k = K - 1; k >= 0; --k) {
            res = max(res, res ^ A[k]);
        }
        return res;
    }
};

```

### determinant

```cpp
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
```

### matrix inverse

```cpp
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
        if (p == -1) {
            return false;
        }
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

```

### homogeneous linear recursion with constant coefficients

```cpp
int homogeneous_linear_recursion_with_constant_coefficients(int *A, int *f, int k, int n) {
    static int a[M], b[M], c[M];
    for (int i = 0; i < k; ++i) {
        a[i] = -f[k - i];
        b[i] = c[i] = 0;
    }
    a[k] = b[1] = c[0] = 1;
    while (n != 0) {
        if ((n & 1) == 1) {
            multiply(b, c, c, k, k);
            modular(c, a, c, 2 * k - 1, k + 1);
        }
        multiply(b, b, b, k, k);
        modular(b, a, b, 2 * k - 1, k + 1);
        n >>= 1;
    }
    long long res = 0;
    for (int i = 0; i < k; ++i) {
        res = (res + (long long) c[i] * A[i]) % mod;
    }
    return res;
}

```

### Berlekamp Massey

```cpp
int Berlekamp_Massey(int *A, int *f, int n) {
    static int g[M], tmp[M];
    int k = 0, last_k = 0, last_delta, last = -1;
    for (int i = 0; i <= n; ++i) {
        tmp[i] = f[i] = 0;
    }
    for (int i = 0; i < n; ++i) {
        long long delta = -A[i];
        for (int j = 1; j <= k; ++j) {
            delta = (delta + (long long) f[j] * A[i - j]) % mod;
        }
        if (delta == 0) {
            continue;
        }
        if (last == -1) {
            k = i + 1;
        } else {
            long long t = delta * power(last_delta, mod - 2) % mod;
            tmp[i - last] = (tmp[i - last] + t) % mod;
            for (int j = 1; j <= last_k; ++j) {
                tmp[i - last + j] = (tmp[i - last + j] - t * g[j]) % mod;
            }
            int p = last_k;
            last_k = k;
            k = max(k, i - last + p);
            for (int j = 1; j <= last_k; ++j) {
                g[j] = f[j];
            }
            for (int j = 1; j <= k; ++j) {
                f[j] = tmp[j];
            }
        }
        last_delta = delta;
        last = i;
    }
    return k;
}
```

## number theory

### ex_gcd

```cpp
long long ex_gcd(long long a, long long b, long long &x, long long &y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    }
    long long res = ex_gcd(b, a % b, y, x);
    y -= a / b * x;
    return res;
}
```

### CRT

```cpp
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

```

### Ex-CTR

```cpp
void ex_crt(long long a1, long long b1, long long a2, long long b2, long long &a, long long &b) {
    long long x, y;
    long long d = ex_gcd(a2, a1, x, y);
    long long _a = a1 / d;
    long long _b = multiply(x, (b1 - b2) / d, _a);
    a = a2 * _a;
    b = (multiply(a2, _b, a) + b2) % a;
}

```

### BSGS

```cpp
int BSGS(long long a, long long b, int mod) {
    if (mod == 1) {
        return 0;
    }
    unordered_map<long long, int> Mp;
    long long w = 1, x = 1;
    int S = sqrt(mod) + 1;
    for (int k = 1; k <= S; ++k) {
        b = b * a % mod;
        w = w * a % mod;
        Mp[b] = k;
    }
    for (int k = 1; k <= S; ++k) {
        x = x * w % mod;
        if (Mp.count(x) != 0) {
            return (k * S - Mp[x]) % (mod - 1);
        }
    }
    return -1;
}

```

### Ex-BSGS

```cpp
int ex_BSGS(long long A, long long B, long long mod) {
    A %= mod;
    B %= mod;
    if (B == 1 || mod == 1)
        return 0;
    int G = gcd(A, mod);
    if (B % G != 0) {
        return -1;
    }
    if (G == 1) {
        return BSGS(A, B, mod);
    }
    int g = inv(A / G, mod / G);
    int x = ex_BSGS(A, B / G * g, mod / G);
    return x == -1 ? -1 : x + 1;
}

```

### ex-Lucas

```cpp
namespace ex_Lucas {
    long long S[K][M], P[K], Pk[K], Mt[K];
    int t;

    void initialize(long long mod) {
        t = 0;
        long long m = mod;
        for (int i = 2; i * i <= m; ++i) {
            if (mod % i == 0) {
                P[t] = i;
                Pk[t] = 1;
                while (m % i == 0) {
                    m /= i;
                    Pk[t] *= i;
                }
                S[t][0] = 1;
                for (int j = 1; j <= Pk[t]; ++j) {
                    S[t][j] = j % i == 0 ? S[t][j - 1] : S[t][j - 1] * j % Pk[t];
                }
                Mt[t] = (mod / Pk[t]) * inverse(mod / Pk[t], Pk[t]) % mod;
                ++t;
            }
        }
        if (m != 1) {
            P[t] = Pk[t] = m;
            S[t][0] = 1;
            for (int j = 1; j <= Pk[t]; ++j) {
                S[t][j] = j % m == 0 ? S[t][j - 1] : S[t][j - 1] * j % Pk[t];
            }
            Mt[t] = (mod / m) * inverse(mod / m, m) % mod;
            ++t;
        }
    }

    long long f(long long n, long long p, long long mod, long long *S) {
        long long res = 1;
        while (n != 0) {
            res = res * power(S[mod], n / mod, mod) % mod * S[n % mod] % mod;
            n /= p;
        }
        return res;
    }

    long long g(long long n, long long p) {
        long long res = 0;
        while (n != 0) {
            n = n / p;
            res += n;
        }
        return res;
    }

    long long combination(long long n, long long m) {
        long long res = 0;
        for (int i = 0; i < t; ++i) {
            long long x = f(n, P[i], Pk[i], S[i]);
            long long y = inverse(f(n - m, P[i], Pk[i], S[i]), Pk[i]);
            long long z = inverse(f(m, P[i], Pk[i], S[i]), Pk[i]);
            long long r = g(n, P[i]) - g(m, P[i]) - g(n - m, P[i]);
            long long b = x * y % Pk[i] * z % Pk[i] * power(P[i], r, Pk[i]) % mod;
            res = (res + b * Mt[i]) % mod;
        }
        return res;
    }
}

```

### Miller-Rabin

```cpp
bool witness(long long a, int s, long long d, long long n) {
    long long x = power(a, d, n);
    if (x == 1 || x == n - 1) {
        return false;
    }
    for (int i = 1; i < s; ++i) {
        x = (__int128) x * x % n;
        if (x == n - 1) {
            return false;
        }
    }
    return true;
}

bool Miller_Rabin(long long n) {
    if (n < M) {
        return !is_composite[n] && n != 1;
    }
    long long d = n - 1;
    int s = 0;
    while ((d & 1) == 0) {
        d >>= 1;
        ++s;
    }
    int p[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29};
    for (int i = 0; i < 10; ++i) {
        if (witness(p[i], s, d, n)) {
            return false;
        }
    }
    return true;
}
```

### Pollard-rho

```cpp
long long f(long long x, long long c, long long n) {
    return ((__int128) x * x + c) % n;
}

long long Pollard_Rho(long long n) {
    if (Miller_Rabin(n)) {
        return n;
    }
    for (int i = 1; i <= 10; ++i) {
        if (n % prime[i] == 0) {
            return (long long) prime[i];
        }
    }
    while (true) {
        long long c = mt_rand() % (n - 1) + 1;
        long long t = f(0, c, n), r = f(f(0, c, n), c, n);
        while (t != r) {
            long long d = gcd(llabs(t - r), n);
            if (d != 1) {
                return d;
            }
            t = f(t, c, n);
            r = f(f(r, c, n), c, n);
        }
    }
}

```

### primitive root

```cpp
vector<int> prime_factors(int n) {
    vector<int> res;
    for (int i = 1; prime[i] * prime[i] <= n; ++i) {
        if (n % prime[i] == 0) {
            while (n % prime[i] == 0) {
                n /= prime[i];
            }
            res.push_back(prime[i]);
        }
    }
    if (n != 1) {
        res.push_back(n);
    }
    sort(res.begin(), res.end());
    return res;
}

int phi(int n) {
    vector<int> A = prime_factors(n);
    for (auto p: A) {
        n = n / p * (p - 1);
    }
    return n;
}

int primitive_root(int n) {
    if (n == 2 || n == 4) {
        return n - 1;
    }
    if ((n & 3) == 0) {
        return 0;
    }
    vector<int> A = prime_factors(n);
    if (2 < A.size() || (A.size() == 2 && A[0] != 2)) {
        return false;
    }
    int m = A.size() == 2 ? n / 2 / A[1] * (A[1] - 1) : n / A[0] * (A[0] - 1);
    vector<int> B = prime_factors(m);
    for (int g = 2; g <= n; ++g) {
        if (power(g, m, n) == 1) {
            bool flag = true;
            for (auto p: B) {
                if (power(g, m / p, n) == 1) {
                    flag = false;
                    break;
                }
            }
            if (flag) {
                return g;
            }
        }
    }
    return 0;
}

```

### sieve of Euler

```cpp
void sieve_of_Euler(int n) {
    for (int i = 2; i <= n; ++i) {
        if (!is_composite[i]) {
            prime[++prime_cnt] = i;
        }
        for (int j = 1; j <= prime_cnt && i * prime[j] <= n; ++j) {
            is_composite[i * prime[j]] = true;
            if (i % prime[j] == 0) {
                break;
            }
        }
    }
}

```

### sieve of Min_25

```cpp
namespace Min_25 {
    long long w[M], h0[M], h1[M], h[M];
    int s, m;

    int index(long long x, long long n) {
        return x <= s ? x : m - n / x + 1;
    }

    int dfs_mu(long long x, int k, long long n) {
        if (x <= prime[k]) {
            return 0;
        }
        int res = (h[index(x, n)] - h[prime[k]]) % mod;
        for (int i = k + 1; (long long) prime[i] * prime[i] <= x; ++i) {
            res = (res - dfs_mu(x / prime[i], i, n)) % mod;
        }
        return res;
    }

    int dfs_phi(long long x, int k, long long n) {
        if (x <= prime[k]) {
            return 0;
        }
        int res = (h[index(x, n)] - h[prime[k]]) % mod;
        for (int i = k + 1; (long long) prime[i] * prime[i] <= x; ++i) {
            long long p = prime[i], d = prime[i], g = p - 1;
            while (d * p <= x) {
                res = (res + g * (dfs_phi(x / d, i, n) + p)) % mod;
                g = (g * p) % mod;
                d *= p;
            }
        }
        return res;
    }

    int sum_of_mu(long long n) {
        m = 0;
        s = (int) sqrtl(n);
        for (long long x = 1; x <= n; x = w[m] + 1) {
            w[++m] = n / (n / x);
            h0[m] = w[m];
        }
        for (int i = 1; prime[i] <= s; ++i) {
            long long p = prime[i];
            for (int j = m; p * p <= w[j]; --j) {
                int k = index(w[j] / p, n);
                h0[j] -= h0[k] - h0[p - 1];
            }
        }
        for (int i = 2; i <= m; ++i) {
            h[i] = -h0[i];
        }
        return 1 + dfs_mu(n, 0, n);
    }

    int sum_of_phi(long long n) {
        m = 0;
        s = (int) sqrtl(n);
        for (long long x = 1; x <= n; x = w[m] + 1) {
            w[++m] = n / (n / x);
            long long tmp = w[m] % mod;
            h0[m] = tmp;
            h1[m] = ((tmp * (tmp + 1)) >> 1) % mod;
        }
        for (int i = 1; prime[i] <= s; ++i) {
            long long p = prime[i];
            for (int j = m; p * p <= w[j]; --j) {
                int k = index(w[j] / p, n);
                h0[j] -= h0[k] - h0[p - 1];
                h1[j] -= (h1[k] - h1[p - 1]) * p % mod;
            }
        }
        for (int i = 2; i <= m; ++i) {
            h[i] = (h1[i] - h0[i]) % mod;
        }
        return (1 + dfs_phi(n, 0, n)) % mod;
    }
}// namespace Min_25

```

## String

### AC automaton

```cpp
struct AC_ automaton {

    static const int M = 100005;
    static const int C = 26;

    int Next[M][C], fail[M], Q[M], dfn[M];
    int tot;

    void clear() {
        tot = 0;
        memset(Next[0], 0, sizeof Next[0]);
    }

    int insert(char *S, int l) {
        int p = 0;
        for (int i = 1; i <= l; ++i) {
            int &q = Next[p][S[i] - 'a'];
            if (q == 0) {
                q = ++tot;
                fail[q] = 0;
                memset(Next[q], 0, sizeof Next[q]);
            }
            p = q;
        }
        return p;
    }

    void bfs() {
        int l = 0, r = 0;
        for (int i = 0; i < C; ++i) {
            if (Next[0][i] != 0) {
                Q[r++] = Next[0][i];
            }
        }
        while (l < r) {
            int p = Q[l++];
            for (int i = 0; i < C; ++i) {
                int &q = Next[p][i];
                if (q == 0) {
                    q = Next[fail[p]][i];
                } else {
                    Q[r++] = q;
                    fail[q] = Next[fail[p]][i];
                }
            }
        }
    }
};
```

### minimal representation

```cpp
int minimal_representation(char *S, int n) {
    int i = 0, j = 1, k = 0;
    while (i < n && j < n && k < n) {
        if (S[(i + k) % n] == S[(j + k) % n]) {
            ++k;
        } else {
            if (S[(i + k) % n] < S[(j + k) % n]) {
                j += k + 1;
            } else {
                i += k + 1;
            }
            if (i == j){
                ++i;
            }
            k = 0;
        }
    }
    return min(i, j);
}
```

### Manacher

```cpp
void Manacher(char *S, int *p, int n) {
    static char T[2 * M];
    T[1] = '#';
    for (int i = 1; i <= n; ++i) {
        T[2 * i] = S[i];
        T[2 * i + 1] = '#';
    }
    for (int i = 1, l = 1, r = 0; i <= 2 * n + 1; ++i) {
        int k = r < i ? 0 : min(p[l + r - i], r - i);
        while (0 < i - k - 1 && i + k + 1 <= 2 * n + 1 && T[i - k - 1] == T[i + k + 1]) {
            ++k;
        }
        p[i] = k;
        if (r < i + k) {
            l = i - k;
            r = i + k;
        }
    }
}
```

### palindromic automation

```cpp
struct palindromic_automaton {
    int Next[M][C], fail[M], len[M], cnt[M];
    int tot;

    int get_fail(char *S, int i, int p) {
        while (i - len[p] - 1 == 0 || S[i - len[p] - 1] != S[i]) {
            p = fail[p];
        }
        return p;
    }

    void build(char *S, int n) {
        len[1] = -1;
        fail[0] = 1;
        tot = 1;
        int p = 0;
        for (int i = 1; i <= n; ++i) {
            int q = get_fail(S, i, p);
            if (Next[q][S[i] - 'a'] == 0) {
                int r = ++tot;
                len[r] = len[q] + 2;
                fail[r] = Next[get_fail(S, i, fail[q])][S[i] - 'a'];
                Next[q][S[i] - 'a'] = r;
            }
            p = Next[q][S[i] - 'a'];
        }
    }
};

```

### suffix sort

```cpp
int sa[M], cnt[M], tmp[M], Rank[M], height[M];

void suffix_sort(char *S, int n) {
    int m = 128;
    for (int i = 1; i <= n; ++i) {
        Rank[i] = S[i];
        ++cnt[Rank[i]];
    }
    for (int i = 1; i <= m; ++i) {
        cnt[i] += cnt[i - 1];
    }
    for (int i = n; 1 <= i; --i) {
        sa[cnt[Rank[i]]--] = i;
    }
    for (int k = 1; k <= n; k <<= 1) {
        int p = 0;
        for (int i = n - k + 1; i <= n; ++i) {
            tmp[++p] = i;
        }
        for (int i = 1; i <= n; ++i) {
            if (k < sa[i]) {
                tmp[++p] = sa[i] - k;
            }
        }
        for (int i = 1; i <= m; ++i) {
            cnt[i] = 0;
        }
        for (int i = 1; i <= n; ++i) {
            ++cnt[Rank[i]];
        }
        for (int i = 1; i <= m; ++i) {
            cnt[i] += cnt[i - 1];
        }
        for (int i = n; 1 <= i; --i) {
            sa[cnt[Rank[tmp[i]]]--] = tmp[i];
        }
        m = tmp[sa[1]] = 1;
        for (int i = 2; i <= n; ++i) {
            tmp[sa[i]] = Rank[sa[i]] != Rank[sa[i - 1]] || Rank[sa[i] + k] != Rank[sa[i - 1] + k]
                             ? ++m
                             : m;
        }
        for (int i = 1; i <= n; ++i) {
            Rank[i] = tmp[i];
        }
        if (n <= m) {
            break;
        }
    }
    for (int i = 1, h = 0; i <= n; ++i) {
        if (h != 0) {
            --h;
        }
        int j = sa[Rank[i] - 1];
        while (S[i + h] == S[j + h]) {
            ++h;
        }
        height[Rank[i]] = h;
    }
}

```

### suffix balanced tree

```cpp
mt19937_64 mt_rand(0);

struct treap {
    double tag;
    int pos, size;
    unsigned long long pri;
    treap *ch[2];
    void push_up() {
        size = 1;
        if (ch[0] != nullptr) {
            size += ch[0]->size;
        }
        if (ch[1] != nullptr) {
            size += ch[1]->size;
        }
    }
};

treap pool[M], *allc = pool, *root;
double tag[M];
char S[M];

void rotate(treap *&p, bool f) {
    treap *q = p->ch[f];
    p->ch[f] = q->ch[!f];
    q->ch[!f] = p;
    p->push_up();
    p = q;
}

int Size(treap *p) {
    return p == nullptr ? 0 : p->size;
}

bool comp(int i, int j) {
    return S[i] == S[j] ? tag[i - 1] < tag[j - 1] : S[i] < S[j];
}

bool comp(char *S, char *Q, int l1, int l2) {
    for (int i = 0; i < l1 && i < l2; ++i) {
        if (S[l1 - i] != Q[l2 - i]) {
            return S[l1 - i] < Q[l2 - i];
        }
    }
    return l1 < l2;
}

void re_tag(treap *p, double l, double r) {
    if (p == nullptr) {
        return;
    }
    double mid = (l + r) * 0.5;
    tag[p->pos] = p->tag = mid;
    re_tag(p->ch[0], l, mid);
    re_tag(p->ch[1], mid, r);
}

void insert(treap *&p, int i, double l, double r) {
    double mid = (l + r) * 0.5;
    if (p == nullptr) {
        p = allcpp;
        tag[i] = mid;
        *p = (treap){mid, i, 1, mt_rand(), {nullptr, nullptr}};
        return;
    }
    bool f = comp(p->pos, i);
    if (f) {
        insert(p->ch[1], i, mid, r);
    } else {
        insert(p->ch[0], i, l, mid);
    }
    if (p->ch[f]->pri < p->pri) {
        rotate(p, f);
        re_tag(p, l, r);
    }
    p->push_up();
}

void remove(treap *&p, int i, double l, double r) {
    double mid = (l + r) * 0.5;
    if (p->pos == i) {
        if (p->ch[0] == nullptr || p->ch[1] == nullptr) {
            p = p->ch[p->ch[0] == nullptr];
            if (p != nullptr) {
                re_tag(p, l, r);
            }
            return;
        }
        bool f = p->ch[1]->pri < p->ch[0]->pri;
        rotate(p, f);
        if (f) {
            remove(p->ch[!f], i, l, mid);
            re_tag(p->ch[f], mid, r);
        } else {
            remove(p->ch[!f], i, mid, r);
            re_tag(p->ch[f], l, mid);
        }
        tag[p->pos] = p->tag = mid;
    } else {
        bool f = comp(p->pos, i);
        if (f) {
            remove(p->ch[f], i, mid, r);
        } else {
            remove(p->ch[f], i, l, mid);
        }
    }
    p->push_up();
}
```

### suffix automaton

```cpp
struct suffix_automaton {
    static const int M = 500005;
    static const int C = 26;
    int trans[2 * M][C], mxlen[2 * M], slink[2 * M], deg[2 * M], Q[2 * M];
    int tot;

    void clear() {
        tot = 1;
        memset(trans[1], 0, sizeof trans[1]);
    }

    int extend(int p, int c) {
        int q = ++tot;
        mxlen[q] = mxlen[p] + 1;
        memset(trans[q], 0, sizeof trans[q]);
        while (p != 0 && trans[p][c] == 0) {
            trans[p][c] = q;
            p = slink[p];
        }
        if (p == 0) {
            slink[q] = 1;
        } else {
            int r = trans[p][c];
            if (mxlen[r] == mxlen[p] + 1) {
                slink[q] = r;
            } else {
                int o = ++tot;
                slink[o] = slink[r];
                mxlen[o] = mxlen[p] + 1;
                memcpy(trans[o], trans[r], sizeof trans[o]);
                while (trans[p][c] == r) {
                    trans[p][c] = o;
                    p = slink[p];
                }
                slink[q] = slink[r] = o;
            }
        }
        return q;
    }

    void build(char *S, int n) {
        int p = 1, l = 0, r = 0;
        tot = 1;
        for (int i = 1; i <= n; ++i) {
            p = extend(p, S[i] - 'a');
        }
    }
};

```

### general suffix automaton

```cpp
struct general_suffix_automaton {
    static const int M = 1000005;
    static const int C = 26;
    int trans[2 * M][C], mxlen[2 * M], slink[2 * M], deg[2 * M], Q[2 * M];
    int tot;

    void initialize() {
        tot = 1;
        memset(trans[1], 0, sizeof trans[1]);
    }

    general_suffix_automaton() {
        initialize();
    }

    int extend(int p, int c) {
        if (trans[p][c] != 0) {
            int r = trans[p][c];
            if (mxlen[r] == mxlen[p] + 1) {
                return r;
            }
            int o = ++tot;
            slink[o] = slink[r];
            mxlen[o] = mxlen[p] + 1;
            memcpy(trans[o], trans[r], sizeof trans[o]);
            while (trans[p][c] == r) {
                trans[p][c] = o;
                p = slink[p];
            }
            slink[r] = o;
            return o;
        }
        int q = ++tot;
        mxlen[q] = mxlen[p] + 1;
        memset(trans[q], 0, sizeof trans[q]);
        while (p != 0 && trans[p][c] == 0) {
            trans[p][c] = q;
            p = slink[p];
        }
        if (p == 0) {
            slink[q] = 1;
        } else {
            int r = trans[p][c];
            if (mxlen[r] == mxlen[p] + 1) {
                slink[q] = r;
            } else {
                int o = ++tot;
                slink[o] = slink[r];
                mxlen[o] = mxlen[p] + 1;
                memcpy(trans[o], trans[r], sizeof trans[o]);
                while (trans[p][c] == r) {
                    trans[p][c] = o;
                    p = slink[p];
                }
                slink[q] = slink[r] = o;
            }
        }
        return q;
    }

    void insert(char *S, int n) {
        int p = 1;
        for (int i = 1; i <= n; ++i) {
            p = extend(p, S[i] - 'a');
        }
    }
};

```

### Lyndon decomposition

```cpp
vector<string> duval(string const& s) {
    int n = s.size(), i = 0;
    vector<string> decomposition;
    while (i < n) {
        int j = i + 1, k = i;
        while (j < n && s[k] <= s[j]) {
            if (s[k] < s[j]) {
                k = i;
            } else {
                ++k;
            }
            j++;
        }
        while (i <= k) {
            decomposition.push_back(s.substr(i, j - k));
            i += j - k;
        }
    }
  return decomposition;
}
```
