#include <cstring>
#include <vector>

template<size_t SIGMA, size_t M> struct ACAutomaton {
    int next[M][SIGMA], fail[M], Q[M], cnt[M];
    std::vector<int> position;
    int tot;
    ACAutomaton() {
        tot = 1;
        memset(next[0], -1, sizeof(next[0]));
    }
    void insert(char *S) {
        int p = 0;
        while (*S != '\0') {
            int c = *S - 'a';
            if (next[p][c] == -1) {
                memset(next[tot], -1, sizeof(next[tot]));
                next[p][c] = tot++;
            }
            p = next[p][c];
            ++S;
        }
        position.push_back(p);
    }
    void build() {
        int head = 0, tail = 0;
        fail[0] = 0;
        for (int i = 0; i < SIGMA; ++i) {
            if (next[0][i] != -1) {
                fail[next[0][i]] = 0;
                Q[tail++] = next[0][i];
            } else {
                next[0][i] = 0;
            }
        }
        while (head < tail) {
            int p = Q[head++];
            for (int i = 0; i < SIGMA; ++i) {
                int q = next[p][i];
                if (q != -1) {
                    fail[q] = next[fail[p]][i];
                    Q[tail++] = q;
                } else {
                    next[p][i] = next[fail[p]][i];
                }
            }
        }
    }

    std::vector<int> query(char *s) {
        memset(cnt, 0, sizeof(cnt));
        int p = 0;
        while (*s != '\0') {
            p = next[p][*s - 'a'];
            ++cnt[p];
            ++s;
        }
        for (int i = tot - 1; i >= 0; i--) {
            if (Q[i] != 0) {
                cnt[fail[Q[i]]] += cnt[Q[i]];
            }
        }
        std::vector<int> ret(position.size());
        for (int i = 0; i < position.size(); i++) {
            ret[i] = cnt[position[i]];
        }
        return ret;
    }
};
