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
    for (int k = 1; k < n; k <<= 1) {
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
        while (S[i + h] == S[sa[Rank[i] - 1] + h]) {
            ++h;
        }
        height[Rank[i]] = h;
    }
}
