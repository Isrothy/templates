struct AC_automaton {

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
