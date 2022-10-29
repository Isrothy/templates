namespace SSP {
    struct edge {
        int from, to, cap, cost, flow;
    };
    edge edges[2 * M];
    vector<int> E[N];
    int Q[N], dis[N], cur[N];
    bool in_queue[N], vis[N];
    int edge_cnt;
    void add_edge(int u, int v, int cap, int cost) {
        edges[edge_cnt] = (edge){u, v, cap, cost, 0};
        E[u].push_back(edge_cnt++);
        edges[edge_cnt] = (edge){v, u, 0, -cost, 0};
        E[v].push_back(edge_cnt++);
    }
    bool SPFA(int S, int T) {
        int head = 0, tail = 0;
        memset(in_queue, 0, sizeof in_queue);
        memset(dis, 0x3f, sizeof dis);
        dis[T] = 0;
        in_queue[T] = true;
        Q[tail++] = T;
        while (head != tail) {
            int u = Q[head++ % N];
            in_queue[u] = false;
            for (int i = 0; i < (int) E[u].size(); ++i) {
                edge e = edges[E[u][i] ^ 1];
                if (e.flow != e.cap && dis[u] + e.cost < dis[e.from]) {
                    dis[e.from] = dis[u] + e.cost;
                    if (!in_queue[e.from]) {
                        in_queue[e.from] = true;
                        Q[tail++ % N] = e.from;
                    }
                }
            }
        }
        return dis[S] != INF;
    }
    int DFS(int u, int T, int a) {
        if (u == T) { return a; }
        vis[u] = true;
        int m = a;
        for (int i = 0; i < (int) E[u].size(); ++i) {
            edge &e = edges[E[u][i]];
            if (e.flow < e.cap && !vis[e.to] && dis[e.to] == dis[u] - e.cost) {
                int f = DFS(e.to, T, min(a, e.cap - e.flow));
                e.flow += f;
                edges[E[u][i] ^ 1].flow -= f;
                a -= f;
                if (a == 0) { return m; }
            }
        }
        return m - a;
    }
    pair<int, int> minimum_cost_flow(int S, int T) {
        int flow = 0, cost = 0;
        while (SPFA(S, T)) {
            memset(vis, 0, sizeof vis);
            int f = DFS(S, T, INF);
            flow += f;
            cost += f * dis[S];
        }
        return make_pair(flow, cost);
    }
}// namespace SSP
