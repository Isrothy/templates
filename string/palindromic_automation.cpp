struct palindromic_automaton {
    int Next[M][C], fail[M], len[M], cnt[M];
    int tot;

    int get_fail(char *S, int i, int p) {
        while (i - len[p] - 1 < 0 || S[i - len[p] - 1] != S[i]) {
            p = fail[p];
        }
        return p;
    }

    void build(char *S) {
        int n = strlen(S);
        len[1] = -1;
        fail[0] = 1;
        tot = 1;
        int p = 0;
        for (int i = 0; i < n; ++i) {
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
