#include <span>
auto four_membered_rings(std::span<std::vector<size_t>> adj) {
    auto n = adj.size() - 1;
    std::vector<size_t> rank(n + 1), vis_time(n + 1), cnt(n + 1);
    std::vector<std::vector<size_t>> f(n + 1);
    for (size_t u = 1; u <= n; ++u) { rank[u] = adj[u].size() * (n + 1) + u; }
    for (size_t u = 1; u <= n; ++u) {
        for (auto v: adj[u]) {
            if (rank[u] < rank[v]) { f[u].push_back(v); }
        }
    }
    size_t res = 0;
    for (size_t u = 1; u <= n; ++u) {
        for (auto v: adj[u]) {
            for (auto w: f[v]) {
                if (rank[u] < rank[w]) {
                    if (vis_time[w] < u) {
                        vis_time[w] = u;
                        cnt[w] = 0;
                    }
                    res += cnt[w];
                    ++cnt[w];
                }
            }
        }
    }
    return res;
}
