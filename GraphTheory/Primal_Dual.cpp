#include <queue>
#include <vector>
struct Network {
    static constexpr int INF = 0x3f3f3f3f;
    struct Edge {
        size_t from, to;
        int cap, cost, flow;
        Edge(size_t from, size_t to, int cap, int cost) : from(from), to(to), cap(cap), cost(cost), flow(0) {}
    };
    std::vector<Edge> edges;
    std::vector<std::vector<size_t>> adj;
    std::vector<int> dis, h;
    std::vector<bool> vis, in_queue;
    size_t n;
    explicit Network(size_t n) : adj(n), dis(n), h(n), vis(n), in_queue(n), n(n) {}
    void add_edge(size_t u, size_t v, int cap, int cost) {
        adj[u].push_back(edges.size());
        edges.emplace_back(u, v, cap, cost);
        adj[v].push_back(edges.size());
        edges.emplace_back(v, u, 0, -cost);
    }
    bool spfa(size_t s, size_t t) {
        std::queue<size_t> q;
        std::vector<bool> in_queue(n + 1, false);
        std::fill(dis.begin(), dis.end(), INF);
        dis[t] = 0;
        in_queue[t] = true;
        q.push(t);
        while (!q.empty()) {
            auto u = q.front();
            q.pop();
            in_queue[u] = false;
            for (auto i: adj[u]) {
                const auto &e = edges[i ^ 1];
                if (e.flow != e.cap && dis[u] + e.cost < dis[e.from]) {
                    dis[e.from] = dis[u] + e.cost;
                    if (!in_queue[e.from]) {
                        in_queue[e.from] = true;
                        q.push(e.from);
                    }
                }
            }
        }
        return dis[s] != INF;
    }
    bool dijkstra(size_t s, size_t t) {
        std::priority_queue<std::pair<int, int>> q;
        std::fill(dis.begin(), dis.end(), INF);
        std::fill(vis.begin(), vis.end(), false);
        dis[t] = 0;
        q.emplace(0, t);
        while (!q.empty()) {
            auto [d, u] = q.top();
            q.pop();
            if (dis[u] != -d) { continue; }
            for (auto i: adj[u]) {
                const auto &e = edges[i ^ 1];
                int c = dis[u] + e.cost + h[u] - h[e.from];
                if (e.flow < e.cap && c < dis[e.from]) {
                    dis[e.from] = c;
                    q.emplace(-c, e.from);
                }
            }
        }
        return dis[s] != INF;
    }
    auto dfs(size_t u, size_t t, int a) {
        if (u == t) { return a; }
        vis[u] = true;
        int m = a;
        for (auto i: adj[u]) {
            auto &e = edges[i];
            if (e.flow < e.cap && !vis[e.to] && h[e.to] == h[u] - e.cost) {
                int f = dfs(e.to, t, std::min(m, e.cap - e.flow));
                e.flow += f;
                edges[i ^ 1].flow -= f;
                m -= f;
                if (m == 0) { break; }
            }
        }
        return a - m;
    }
    auto minimum_cost_flow(size_t s, size_t t) {
        int flow = 0;
        int cost = 0;
        bool first = true;
        while (first ? spfa(s, t) : dijkstra(s, t)) {
            first = false;
            for (size_t i = 0; i < n; ++i) { h[i] += dis[i]; }
            while (true) {
                std::fill(vis.begin(), vis.end(), false);
                int f = dfs(s, t, INF);
                if (f == 0) { break; }
                flow += f;
                cost += f * h[s];
            }
        }
        return std::make_pair(flow, cost);
    }
};
