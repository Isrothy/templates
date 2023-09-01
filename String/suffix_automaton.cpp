#include <cstring>
template<size_t SIGMA, size_t M>
struct SuffixAutomaton {
    int trans[2 * M][SIGMA], mxlen[2 * M], slink[2 * M];
    int tot;
    SuffixAutomaton() : tot(1) {
        slink[0] = -1;
        memset(trans[0], -1, sizeof(trans[0]));
    }
    int extend(int p, int c) {
        int q = tot++;
        mxlen[q] = mxlen[p] + 1;
        memset(trans[q], -1, sizeof trans[q]);
        while (p != -1 && trans[p][c] == -1) {
            trans[p][c] = q;
            p = slink[p];
        }
        if (p == -1) {
            slink[q] = 0;
        } else {
            int r = trans[p][c];
            if (mxlen[r] == mxlen[p] + 1) {
                slink[q] = r;
            } else {
                int o = tot++;
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
    void insert(char *S) {
        int p = 0, n = (int) strlen(S);
        for (int i = 0; i < n; ++i) { p = extend(p, S[i] - 'a'); }
    }
};
