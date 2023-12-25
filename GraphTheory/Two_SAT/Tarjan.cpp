#include <optional>
#include <stack>
#include <vector>
struct TwoSAT {
    std::stack<size_t> scc{};
    std::vector<size_t> sccno, dfn, low;
    std::vector<std::vector<size_t>> adj;
    std::vector<bool> mark;
    size_t n, scc_cnt, idx;
    explicit TwoSAT(size_t n) : dfn(2 * n), low(2 * n), adj(2 * n), mark(2 * n), n(n), scc_cnt(0), idx(0) {}
    void add_clause(size_t x, size_t xval, size_t y, size_t yval) {
        x = x * 2 + xval;
        y = y * 2 + yval;
        adj[x ^ 1].push_back(y);
        adj[y ^ 1].push_back(x);
    }
    void Tarjan(size_t u) {
        dfn[u] = low[u] = ++idx;
        scc.push(u);
        mark[u] = true;
        for (auto v: adj[u]) {
            if (!dfn[v]) {
                Tarjan(v);
                low[u] = std::min(low[u], low[v]);
            } else if (mark[v]) {
                low[u] = std::min(low[u], dfn[v]);
            }
        }
        if (dfn[u] == low[u]) {
            ++scc_cnt;
            while (true) {
                auto x = scc.top();
                scc.pop();
                mark[x] = false;
                sccno[x] = scc_cnt;
                if (x == u) { break; }
            }
        }
    }
    auto solve() -> std::optional<std::vector<bool>> {
        sccno.resize(2 * n);
        for (size_t i = 0; i < 2 * n; ++i) {
            if (!dfn[i]) { Tarjan(i); }
        }
        std::vector<bool> res(n);
        for (size_t i = 0; i < n; ++i) {
            if (sccno[i * 2] == sccno[i * 2 + 1]) { return std::nullopt; }
            res[i] = sccno[i * 2] > sccno[i * 2 + 1];
        }
        return res;
    }
};
