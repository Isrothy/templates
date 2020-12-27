namespace Dinic {
    struct edge {
        int from, to, cap, flow;
    };
    edge edges[2 * M];
    vector<int> E[M];
    int dis[N], Q[N], cur[M];
    bool vis[N];
    int edge_cnt;
    
    void add_edge(int u, int v, int c) {
        edges[edge_cnt] = (edge) {u, v, c, 0};
        E[u].push_back(edge_cnt++);
        edges[edge_cnt] = (edge) {v, u, c, 0};
        E[v].push_back(edge_cnt++);
    }
    
    void clear() {
        for (int i = 0; i < edge_cnt; ++i) {
            edges[i].flow = 0;
        }
    }
    
    bool BFS(int S, int T, int n) {
        int head = 0, tail = 0;
        for (int i = 1; i <= n; ++i) {
            vis[i] = false;
        }
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
                    break;
            }
        }
        return m - a;
    }
    
    long long max_flow(int S, int T, int n) {
        long long flow = 0;
        while (BFS(S, T, n)) {
            for (int i = 1; i <= n; ++i) {
                cur[i] = 0;
            }
            flow += DFS(S, T, INF);
        }
        return flow;
    }
}
