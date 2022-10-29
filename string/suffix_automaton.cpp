struct suffix_automaton {
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
    void build(char *S) {
        int p = 1;
        int n = strlen(S);
        tot = 1;
        for (int i = 0; i < n; ++i) { p = extend(p, S[i] - 'a'); }
    }
};
