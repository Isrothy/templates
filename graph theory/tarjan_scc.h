vector <int> E[M];
int dfn[M], low[M], sccno[M], stk[M];
int dfs_clock, scc_cnt, top;

void Tarjan(int u) {
    dfn[u] = low[u] = ++dfs_clock;
    stk[top++] = u;
    for (int i = 0; i < (int) E[u].size(); ++i) {
        int v = E[u][i];
        if (dfn[v] == 0) {
            Tarjan(v);
            low[u] = min(low[u], low[v]);
        } else if (sccno[v] == 0)
            low[u] = min(low[u], dfn[v]);
    }
    if (dfn[u] == low[u]) {
        ++scc_cnt;
        for (;;) {
            int v = stk[--top];
            sccno[v] = scc_cnt;
            if (v == u)
                break;
        }
    }
}
