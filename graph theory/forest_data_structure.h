vector<int>E[M], R[2 * M];
int dfn[M], low[M], stk[M], bccno[M], A[M];
int n, m, bcc_cnt, dfs_clock, top;

void Tarjan(int u, int fa) {
    dfn[u] = low[u] = ++dfs_clock;
    stk[top++] = u;
    for (int i = 0; i < (int) E[u].size(); ++i) {
        int v = E[u][i];
        if (v == fa)
            continue;
        if (dfn[v] == 0) {
            Tarjan(v, u);
            low[u] = min(low[u], low[v]);
            if (dfn[u] <= low[v]) {
                ++bcc_cnt;
                int x, w = n + bcc_cnt;
                do {
                    x = stk[--top];
                    R[w].push_back(x);
                } while (x != v);
                R[u].push_back(w);
            }
        } else
            low[u] = min(low[u], dfn[v]);
    }
}
