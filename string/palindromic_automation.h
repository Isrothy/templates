struct Palindromic_Automation {
    int Next[M][C], len[M], fail[M];
    int tot;

    void build(char *S, int n) {
        fail[0] = tot = 1;
        len[1] = -1;
        int p = 1;
        for (int i = 1; i <= n; ++i) {
            while (S[i - len[p] - 1] != S[i])
                p = fail[p];
            int &q = Next[p][S[i] - 'a'];
            if (!q) {
                int r = ++tot;
                len[r] = len[p] + 2;
                p = fail[p];
                while (S[i - len[p] - 1] != S[i])
                    p = fail[p];
                fail[r] = Next[p][S[i] - 'a'];
                q = r;
            }
            p = q;
        }
    }
};
