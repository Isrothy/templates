#include <cstdio>
#include <queue>
#include <vector>
struct Network {
    static constexpr int INF = 0x3f3f3f3f;
    struct Edge {
        size_t from, to;
        int cap, flow;
        Edge(size_t u, size_t v, int c) : from(u), to(v), cap(c), flow(0) {}
    };
    std::vector<Edge> edges;
    std::vector<std::vector<size_t>> adj;
    std::vector<int> dis;
    std::vector<size_t> cur;
    std::vector<bool> vis;
    size_t n;
    explicit Network(size_t n) : adj(n), dis(n), cur(n), vis(n), n(n) {}
    void add_edge(size_t u, size_t v, int c) {
        adj[u].push_back(edges.size());
        edges.emplace_back(u, v, c);
        adj[v].push_back(edges.size());
        edges.emplace_back(v, u, 0);
    }
    bool BFS(size_t s, size_t t) {
        std::queue<size_t> q;
        fill(vis.begin(), vis.end(), false);
        dis[t] = 0;
        vis[t] = true;
        q.push(t);
        while (!q.empty()) {
            auto u = q.front();
            q.pop();
            for (auto i: adj[u]) {
                Edge e = edges[i ^ 1];
                if (e.flow < e.cap && !vis[e.from]) {
                    vis[e.from] = true;
                    dis[e.from] = dis[u] + 1;
                    q.push(e.from);
                }
            }
        }
        return vis[s];
    }
    auto DFS(size_t u, size_t t, int a) {
        if (u == t) { return a; }
        int m = a;
        for (auto &i = cur[u]; i < adj[u].size(); ++i) {
            Edge &e = edges[adj[u][i]];
            if (e.flow < e.cap && vis[e.to] && dis[e.to] == dis[u] - 1) {
                int f = DFS(e.to, t, std::min(m, e.cap - e.flow));
                e.flow += f;
                edges[adj[u][i] ^ 1].flow -= f;
                m -= f;
                if (a == 0) { break; }
            }
        }
        return a - m;
    }
    auto max_flow(size_t S, size_t T) {
        int flow = 0;
        while (BFS(S, T)) {
            fill(cur.begin(), cur.end(), 0);
            flow += DFS(S, T, INF);
        }
        return flow;
    }
};
