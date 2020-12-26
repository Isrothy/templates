int four_membered_rings(vector<int> *E, int n) {
    static long long Rank[M];
    static int vis_time[M], cnt[M];
    static vector<int> F[M];
    for (int u = 1; u <= n; ++u) {
        Rank[u] = (long long) E[u].size() * (n + 1) + u;
        vis_time[u] = 0;
        F[u].clear();
    }
    for (int u = 1; u <= n; ++u) {
        for (auto v : E[u]) {
            if (Rank[u] < Rank[v])
                F[u].push_back(v);
        }
    }
    int res = 0;
    for (int u = 1; u <= n; ++u) {
            for (auto v : E[u]) {
                for (auto w : F[v]) {
                    if (Rank[u] < Rank[w]) {
                        if (vis_time[w] < u) {
                            vis_time[w] = u;
                            cnt[w] = 0;
                        }
                        res = (res + cnt[w]) % mod;
                        ++cnt[w];
                    }
                }
            }
        }
    return res;
}
