#include <vector>
template<size_t M> struct Two_SAT {
    std::vector<int> E[2 * M];
    int stk[M]{};
    bool mark[2 * M]{};
    int top{};
    void add_clause(int u, bool f1, int v, bool f2) {
        u = u << 1 | f1;
        v = v << 1 | f2;
        E[u ^ 1].push_back(v);
        E[v ^ 1].push_back(u);
    }
    bool dfs(int x) {
        if (mark[x ^ 1]) {
            return false;
        }
        if (mark[x]) {
            return true;
        }
        mark[x] = true;
        stk[top++] = x;
        for (auto y: E[x]) {
            if (!dfs(y))
                return false;
        }
        return true;
    }
    bool check(int n) {
        for (int i = 0; i < 2 * n; i += 2) {
            if (!mark[i] && !mark[i ^ 1]) {
                top = 0;
                if (!dfs(i)) {
                    while (top != 0) {
                        mark[stk[--top]] = false;
                    }
                    if (!dfs(i ^ 1)) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
}
