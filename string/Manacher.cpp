#include <algorithm>
#include <cstring>

void Manacher(char *S, int *p) {
    static char T[2 * M];
    T[0] = '#';
    int n = strlen(S);
    for (int i = 0; i < n; ++i) {
        T[2 * i + 1] = S[i];
        T[2 * i + 2] = '#';
    }
    for (int i = 0, l = 0, r = -1; i <= 2 * n; ++i) {
        int k = r < i ? 0 : std::min(p[l + r - i], r - i);
        while (0 <= i - k - 1 && i + k + 1 <= 2 * n && T[i - k - 1] == T[i + k + 1]) {
            ++k;
        }
        p[i] = k;
        if (r < i + k) {
            l = i - k;
            r = i + k;
        }
    }
}
