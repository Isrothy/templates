#include <queue>
#include <vector>

struct BipartiteGraph {
    std::vector<std::vector<int>> E;
    std::vector<int> dx, dy, S, T;
    size_t nl, nr;
    BipartiteGraph(size_t nl, size_t nr)
        : E(nl), dx(nl), dy(nr), S(nl, -1), T(nr, -1), nl(nl), nr(nr) {}
    void add_edge(int u, int v) {
        assert(0 <= u && u < nl && 0 <= v && v < nr);
        E[u].push_back(v);
    }
    bool BFS() {
        fill(dx.begin(), dx.end(), -1);
        fill(dy.begin(), dy.end(), -1);
        int ret = -1;
        std::queue<int> Q;
        for (int u = 0; u < nl; ++u) {
            if (S[u] == -1) {
                dx[u] = 0;
                Q.push(u);
            }
        }
        while (!Q.empty()) {
            int u = Q.front();
            Q.pop();
            if (ret != -1 && dx[u] > ret) {
                break;
            }
            for (auto v: E[u]) {
                if (dy[v] == -1) {
                    dy[v] = dx[u];
                    if (T[v] != -1) {
                        dx[T[v]] = dy[v] + 1;
                        Q.push(T[v]);
                    } else {
                        ret = dx[u];
                    }
                }
            }
        }
        return ret != -1;
    }

    bool DFS(int u) {
        for (auto v: E[u]) {
            if (dy[v] == dx[u]) {
                dy[v] = -1;
                if (T[v] == -1 || DFS(T[v])) {
                    S[u] = v;
                    T[v] = u;
                    return true;
                }
            }
        }
        return false;
    }
    int maximum_matching() {
        int res = 0;
        while (BFS()) {
            for (int u = 0; u < nl; ++u) {
                if (S[u] == -1 && DFS(u)) {
                    ++res;
                }
            }
        }
        return res;
    }
};
