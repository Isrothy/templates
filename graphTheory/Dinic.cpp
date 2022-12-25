#include <cstdio>
#include <queue>
#include <vector>

struct Network {
    static const int INF = 0x3f3f3f3f;
    struct edge {
        int from, to, cap, flow;
        edge(int u, int v, int c, int f) : from(u), to(v), cap(c), flow(f) {}
    };
    std::vector<edge> edges;
    std::vector<std::vector<size_t>> adj;
    std::vector<int> dis;
    std::vector<size_t> cur;
    std::vector<bool> vis;
    int n;
    explicit Network(int n) : adj(n), dis(n), cur(n), vis(n), n(n) {}
    void add_edge(int u, int v, int c) {
        adj[u].push_back(edges.size());
        edges.emplace_back(u, v, c, 0);
        adj[v].push_back(edges.size());
        edges.emplace_back(v, u, 0, 0);
    }
    void clear() {
        for (auto &e: edges) {
            e.flow = 0;
        }
    }
    bool BFS(int S, int T) {
        std::queue<int> q;
        fill(vis.begin(), vis.end(), false);
        dis[T] = 0;
        vis[T] = true;
        q.push(T);
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (auto i: adj[u]) {
                edge e = edges[i ^ 1];
                if (e.flow < e.cap && !vis[e.from]) {
                    vis[e.from] = true;
                    dis[e.from] = dis[u] + 1;
                    q.push(e.from);
                }
            }
        }
        return vis[S];
    }
    int DFS(int u, int T, int a) {
        if (u == T) {
            return a;
        }
        int m = a;
        for (auto &i = cur[u]; i < adj[u].size(); ++i) {
            edge &e = edges[adj[u][i]];
            if (e.flow < e.cap && vis[e.to] && dis[e.to] == dis[u] - 1) {
                int f = DFS(e.to, T, std::min(a, e.cap - e.flow));
                e.flow += f;
                edges[adj[u][i] ^ 1].flow -= f;
                a -= f;
                if (a == 0) {
                    break;
                }
            }
        }
        return m - a;
    }
    long long max_flow(int S, int T) {
        long long flow = 0;
        while (BFS(S, T)) {
            fill(cur.begin(), cur.end(), 0);
            flow += DFS(S, T, INF);
        }
        return flow;
    }
};

int main() {
    int n, m, S, T;
    scanf("%d%d%d%d", &n, &m, &S, &T);
    Network network(n + 1);
    for (int i = 0; i < m; ++i) {
        int u, v, c;
        scanf("%d%d%d", &u, &v, &c);
        network.add_edge(u, v, c);
    }
    printf("%lld\n", network.max_flow(S, T));
    return 0;
};
