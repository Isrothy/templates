vector<pair<int, int>> E[M];
int dfn[M], low[M], ebcno[M];
bool is_bridge[M];
int n, m, dfs_clock, ebc_cnt;
void Tarjan(int u, int pre) {
    dfn[u] = low[u] = ++dfs_clock;
    for (int i = 0; i < (int) E[u].size(); ++i) {
        int v = E[u][i].first, e = E[u][i].second;
        if (e == pre) { continue; }
        if (dfn[v] == 0) {
            Tarjan(v, e);
            low[u] = min(low[u], low[v]);
            is_bridge[e] = dfn[u] < low[v];
        } else {
            low[u] = min(low[u], dfn[v]);
        }
    }
}
void dfs(int u) {
    ebcno[u] = ebc_cnt;
    for (int i = 0; i < (int) E[u].size(); ++i) {
        int v = E[u][i].first;
        if (!is_bridge[E[u][i].second] && ebcno[v] == 0) {
            dfs(v);
        }
    }
}
