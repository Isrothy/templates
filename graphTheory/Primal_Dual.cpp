#include <vector>
#include <queue>
struct Network {
    static const int INF = 0x3f3f3f3f;
    struct edge {
        int from, to, cap, cost, flow;
        edge(int from, int to, int cap, int cost)
            : from(from), to(to), cap(cap), cost(cost), flow(0) {}
    };
    std::vector<edge> edges;
    std::vector<std::vector<size_t>> adj;
    std::vector<int> dis, h;
    std::vector<bool> vis, in_queue;
    int n;
    explicit Network(int n) : adj(n), dis(n), h(n), vis(n), in_queue(n), n(n) {}
    void add_edge(int u, int v, int cap, int cost) {
        adj[u].push_back(edges.size());
        edges.emplace_back(u, v, cap, cost);
        adj[v].push_back(edges.size());
        edges.emplace_back(v, u, 0, -cost);
    }
    bool SPFA(int S, int T) {
        std::queue<int> q;
        std::vector<bool> in_queue(n + 1, false);
        fill(dis.begin(), dis.end(), INF);
        dis[T] = 0;
        in_queue[T] = true;
        q.push(T);
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            in_queue[u] = false;
            for (auto i: adj[u]) {
                edge e = edges[i ^ 1];
                if (e.flow != e.cap && dis[u] + e.cost < dis[e.from]) {
                    dis[e.from] = dis[u] + e.cost;
                    if (!in_queue[e.from]) {
                        in_queue[e.from] = true;
                        q.push(e.from);
                    }
                }
            }
        }
        return dis[S] != INF;
    }
    bool Dijkstra(int S, int T) {
        std::priority_queue<std::pair<int, int>> q;
        fill(dis.begin(), dis.end(), INF);
        fill(vis.begin(), vis.end(), false);
        dis[T] = 0;
        q.emplace(0, T);
        while (!q.empty()) {
            auto [d, u] = q.top();
            q.pop();
            if (dis[u] != -d) {
                continue;
            }
            for (auto i: adj[u]) {
                edge e = edges[i ^ 1];
                int c = dis[u] + e.cost + h[u] - h[e.from];
                if (e.flow < e.cap && c < dis[e.from]) {
                    dis[e.from] = c;
                    q.emplace(-c, e.from);
                }
            }
        }
        return dis[S] != INF;
    }
    int DFS(int u, int T, int a) {
        if (u == T) {
            return a;
        }
        vis[u] = true;
        int m = a;
        for (auto i: adj[u]) {
            edge &e = edges[i];
            if (e.flow < e.cap && !vis[e.to] && h[e.to] == h[u] - e.cost) {
                int f = DFS(e.to, T, std::min(a, e.cap - e.flow));
                e.flow += f;
                edges[i ^ 1].flow -= f;
                a -= f;
                if (a == 0) {
                    return m;
                }
            }
        }
        return m - a;
    }
    std::pair<int, int> minimum_cost_flow(int S, int T) {
        int flow = 0, cost = 0;
        bool first = true;
        while (first ? SPFA(S, T) : Dijkstra(S, T)) {
            first = false;
            for (int i = 0; i < n; ++i) {
                h[i] += dis[i];
            }
            while (true) {
                fill(vis.begin(), vis.end(), false);
                int f = DFS(S, T, INF);
                if (f == 0) {
                    break;
                }
                flow += f;
                cost += f * h[S];
            }
        }
        return std::make_pair(flow, cost);
    }
}
