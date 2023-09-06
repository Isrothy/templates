vector<int> E[M], R[M], Q[M];
int par[M], dfn[M], vertices[M], Fa[M], idom[M], sdom[M], val[M];
int dfs_clock;
bool cmp(int u, int v) { return dfn[u] < dfn[v]; }
int Find(int u) {
    if (u == Fa[u]) { return u; }
    int &v = Fa[u], w = Find(Fa[u]);
    if (v != w) {
        if (cmp(sdom[val[v]], sdom[val[u]])) { val[u] = val[v]; }
        v = w;
    }
    return Fa[u];
}
void dfs(int u) {
    dfn[u] = ++dfs_clock;
    vertices[dfs_clock] = u;
    for (auto v: E[u]) {
        if (dfn[v] == 0) {
            par[v] = u;
            dfs(v);
        }
    }
}
void Lengauer_Tarjan(int root) {
    dfs(root);
    for (int u = 1; u <= n; ++u) { Fa[u] = sdom[u] = val[u] = u; }
    for (int i = dfs_clock; i; --i) {
        int u = vertices[i];
        for (auto v: R[u]) {
            if (dfn[v] == 0) { continue; }
            Find(v);
            sdom[u] = min(sdom[u], sdom[val[v]], cmp);
        }
        for (auto v: Q[u]) {
            Find(v);
            idom[v] = sdom[v] == sdom[val[v]] ? sdom[v] : val[v];
        }
        if (i != 1) {
            Fa[u] = par[u];
            Q[sdom[u]].push_back(u);
        }
    }
    for (int i = 2; i <= dfs_clock; ++i) {
        int u = vertices[i];
        if (idom[u] != sdom[u]) { idom[u] = idom[idom[u]]; }
    }
}
