vector<int> E[M];
int dfn[M], low[M], stk[M], bccno[M];
int n, m, bcc_cnt, dfs_clock, top;
void Tarjan(int u, int fa) {
    dfn[u] = low[u] = ++dfs_clock;
    stk[top++] = u;
    for (auto v: E[u]) {
        if (v == fa) { continue; }
        if (dfn[v] == 0) {
            Tarjan(v, u);
            low[u] = min(low[u], low[v]);
            if (dfn[u] <= low[v]) {
                ++bcc_cnt;
                int x;
                do {
                    x = stk[--top];
                    bccno[x] = bcc_cnt;
                } while (x != v);
            }
        } else {
            low[u] = min(low[u], dfn[v]);
        }
    }
}
