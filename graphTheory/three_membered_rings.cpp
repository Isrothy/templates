int three_membered_rings(vector<int> *E, int n) {
    static long long Rank[M];
    static int vis_time[M];
    static vector<int> F[M];
    for (int u = 1; u <= n; ++u) {
        Rank[u] = (long long) E[u].size() * (n + 1) + u;
        vis_time[u] = 0;
        F[u].clear();
    }
    for (int u = 1; u <= n; ++u) {
        for (auto v: E[u]) {
            if (Rank[u] < Rank[v]) {
                F[u].push_back(v);
            }
        }
    }
    int res = 0;
    for (int u = 1; u <= n; ++u) {
        for (auto v: F[u]) {
            vis_time[v] = u;
        }
        for (auto v: F[u]) {
            for (auto w: F[v]) {
                if (vis_time[w] == u) {
                    ++res;
                }
            }
        }
    }
    return res;
}
