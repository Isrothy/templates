namespace Dinic {
    const int INF = INT_MAX;
    struct edge {
        int from, to, cap, flow;
    };
    edge edges[M * 2];
    vector<int> E[N];
    int dis[N], Q[N], cur[N];
    bool vis[N];
    int edge_cnt;
    
    void add_edge(int u, int v, int c) {
        edges[edge_cnt] = (edge) {u, v, c, 0};
        E[u].push_back(edge_cnt++);
        edges[edge_cnt] = (edge) {v, u, 0, 0};
        E[v].push_back(edge_cnt++);
    }

    bool BFS(int S, int T) {
        int head = 0, tail = 0;
        memset(vis, 0, sizeof vis);
        dis[T] = 0;
        vis[T] = true;
        Q[tail++] = T;
        while (head < tail) {
            int u = Q[head++];
            for (int i = 0; i < (int) E[u].size(); ++i) {
                edge e = edges[E[u][i] ^ 1];
                if (e.flow < e.cap && !vis[e.from]) {
                    vis[e.from] = true;
                    dis[e.from] = dis[u] + 1;
                    Q[tail++] = e.from;
                }
            }
        }
        return vis[S];
    }

    int DFS(int u, int T, int a) {
        if (u == T)
            return a;
        int m = a;
        for (int &i = cur[u]; i < (int) E[u].size(); ++i) {
            edge &e = edges[E[u][i]];
            if (e.flow < e.cap && vis[e.to] && dis[e.to] == dis[u] - 1) {
                int f = DFS(e.to, T, min(a, e.cap - e.flow));
                e.flow += f;
                edges[E[u][i] ^ 1].flow -= f;
                a -= f;
                if (a == 0)
                    return m;
            }
        }
        return m - a;
    }

    long long max_flow(int S, int T) {
        long long flow = 0;
        while (BFS(S, T)) {
            memset(cur, 0, sizeof cur);
            flow += DFS(S, T, INF);
        }
        return flow;
    }
}
