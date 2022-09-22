#include <numeric>
#include <string>

void suffix_sort(char *S, int *sa, int *rank, int *height) {
    static int tmp[M], cnt[M];
    int n = strlen(S) + 1, m = 128;
    for (int i = 0; i < m; ++i) {
        cnt[i] = 0;
    }
    for (int i = 0; i < n; ++i) {
        rank[i] = S[i];
        ++cnt[rank[i]];
    }
    std::partial_sum(cnt, cnt + m, cnt);
    for (int i = n - 1; i >= 0; --i) {
        sa[--cnt[rank[i]]] = i;
    }
    for (int k = 1; k <= n; k <<= 1) {
        int p = 0;
        for (int i = n - k; i < n; ++i) {
            tmp[p++] = i;
        }
        for (int i = 0; i < n; ++i) {
            if (k <= sa[i]) {
                tmp[p++] = sa[i] - k;
            }
        }
        for (int i = 0; i < m; ++i) {
            cnt[i] = 0;
        }
        for (int i = 0; i < n; ++i) {
            ++cnt[rank[i]];
        }
        std::partial_sum(cnt, cnt + m, cnt);
        for (int i = n - 1; i >= 0; --i) {
            sa[--cnt[rank[tmp[i]]]] = tmp[i];
        }
        m = 1;
        tmp[sa[0]] = 0;
        for (int i = 1; i < n; ++i) {
            if (rank[sa[i]] != rank[sa[i - 1]] || rank[sa[i] + k] != rank[sa[i - 1] + k]) {
                ++m;
            }
            tmp[sa[i]] = m - 1;
        }
        std::copy(tmp, tmp + n, rank);
        if (n == m) {
            break;
        }
    }
    --n;
    int h = 0;
    for (int i = 0; i < n; ++i) {
        if (h != 0) {
            --h;
        }
        int j = sa[rank[i] - 1];
        while (h <= n && S[i + h] == S[j + h]) {
            ++h;
        }
        height[rank[i] - 1] = h;
    }
    height[n] = 0;
}
