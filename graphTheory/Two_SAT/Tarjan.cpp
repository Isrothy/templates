struct TwoSAT {
    std::stack<int> scc;
    std::vector<int> sccno, dfn, low;
    std::vector<std::vector<int>> adj;
    std::vector<bool> mark;
    int n, scc_cnt, idx;
    explicit TwoSAT(int n)
        : dfn(2 * n), low(2 * n), adj(2 * n), mark(2 * n), n(n), scc_cnt(0), idx(0) {}
    void add_clause(int x, int xval, int y, int yval) {
        x = x * 2 + xval;
        y = y * 2 + yval;
        adj[x ^ 1].push_back(y);
        adj[y ^ 1].push_back(x);
    }
    void Tarjan(int u) {
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
            int x;
            ++scc_cnt;
            do {
                x = scc.top();
                scc.pop();
                mark[x] = false;
                sccno[x] = scc_cnt;
            } while (x != u);
        }
    }
    std::vector<bool> solve() {
        sccno.resize(2 * n);
        for (int i = 0; i < 2 * n; ++i) {
            if (!dfn[i]) {
                Tarjan(i);
            }
        }
        std::vector<bool> res(n);
        for (int i = 0; i < n; ++i) {
            if (sccno[i * 2] == sccno[i * 2 + 1]) {
                return {};
            }
            res[i] = sccno[i * 2] > sccno[i * 2 + 1];
        }
        return res;
    }
};
