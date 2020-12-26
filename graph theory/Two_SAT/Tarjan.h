namespace Two_SAT {
    vector<int> E[2 * M], S[2 * M];
    int dfn[2 * M], low[2 * M], sccno[2 * M], stk[2 * M];
    int deg[2 * M], topono[2 * M], Q[2 * M];
    int dfs_clock, scc_cnt, top;
    
    void add_clause(int u, bool f1, int v, bool f2) {
        u = u << 1 | f1;
        v = v << 1 | f2;
        E[u ^ 1].push_back(v);
        E[v ^ 1].push_back(u);
    }
    
    void Tarjan(int u) {
        dfn[u] = low[u] = ++dfs_clock;
        stk[top++] = u;
        for (auto v : E[u]) {
            if (dfn[v] == 0) {
                Tarjan(v);
                low[u] = min(low[u], low[v]);
            } else if (sccno[v] == 0) {
                low[u] = min(low[u], dfn[v]);
            }
        }
        if (dfn[u] == low[u]) {
            ++scc_cnt;
            for (;;) {
                int v = stk[--top];
                sccno[v] = scc_cnt;
                if (v == u)
                    break;
            }
        }W
    }
    
    bool query(int u) {
        return topono[sccno[u << 1]] < topono[sccno[u << 1 | 1]];
    }
    
    bool check(int n) {
        for (int u = 0; u < 2 * n; ++u) {
            if (dfn[u] == 0) {
                Tarjan(u);
            }
        }
        for (int u = 0; u < n; ++u) {
            if (sccno[u << 1] == sccno[u << 1 | 1])
                return false;
        }
        for (int u = 0; u < 2 * n; ++u) {
            for (auto v : E[u]) {
                if (sccno[u] != sccno[v]) {
                    S[sccno[u]].push_back(sccno[v]);
                    ++deg[sccno[v]];
                }
            }
        }
        int head = 0, tail = 0;
        for (int u = 1; u <= scc_cnt; ++u) {
            if (deg[u] == 0) {
                Q[tail++] = u;
                topono[u] = tail;
            }
        }
        while (head < tail) {
            int u = Q[head++];
            for (auto v : S[u]) {
                if (--deg[v] == 0) {
                    Q[tail++] = v;
                    topono[v] = tail;
                }
            }
        }
        return true;
    }
}
