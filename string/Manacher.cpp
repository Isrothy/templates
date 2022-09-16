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
