void dfs(int u) {
    len[u] = 1;
    for (auto v : E[u]) {
        if (v == par[u])
            continue;
        anc[v][0] = par[v] = u;
        dep[v] = dep[u] + 1;
        for (int k = 1; k < K; ++k) {
            anc[v][k] = anc[anc[v][k - 1]][k - 1];
        }
        dfs(v);
        if (len[u] < len[v] + 1) {
            len[u] = len[v] + 1;
            son[u] = v;
        }
    }
}

void re_dfs(int u) {
    ladder[++ladder_sz] = u;
    id[u] = ladder_sz;
    if (son[u] != 0) {
        re_dfs(son[u]);
    }
    for (auto v : E[u]) {
        if (v == par[u] || v == son[u])
            continue;
        int tmp = ladder_sz;
        for (int j = 1, w = u; j < len[v] && w != 0; ++j) {
            ladder[++ladder_sz] = w;
            w = par[w];
        }
        reverse(ladder + tmp + 1, ladder + ladder_sz + 1);
        re_dfs(v);
    }
}

int Kth_anc(int u, int k) {
    if (k == 0)
        return u;
    int v = anc[u][Log2[k]];
    return ladder[id[v] - k + (1 << Log2[k])];
}
