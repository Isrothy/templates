namespace Blossom {
    vector<int> E[M];
    int mate[M], link[M], label[M], fa[M], Q[M];
    int head, tail;
    int find(int x) {
        if (x == fa[x]) { return x; }
        fa[x] = find(fa[x]);
        return fa[x];
    }
    void add_edge(int u, int v) {
        E[u].push_back(v);
        E[v].push_back(u);
    }
    int LCA(int u, int v) {
        static int last_time[M], time;
        ++time;
        while (last_time[u] != time) {
            if (u != 0) {
                last_time[u] = time;
                u = find(link[mate[u]]);
            }
            swap(u, v);
        }
        return u;
    }
    void blossom(int u, int v, int w) {
        while (find(u) != w) {
            link[u] = v;
            v = mate[u];
            fa[v] = fa[u] = w;
            if (label[v] == 1) {
                label[v] = 2;
                Q[tail++] = v;
            }
            u = link[v];
        }
    }
    bool BFS(int S, int n) {
        head = tail = 0;
        for (int i = 1; i <= n; ++i) {
            fa[i] = i;
            label[i] = 0;
        }
        Q[tail++] = S;
        label[S] = 2;
        while (head < tail) {
            int u = Q[head++];
            for (auto v: E[u]) {
                if (label[v] == 0) {
                    label[v] = 1;
                    link[v] = u;
                    if (mate[v] == 0) {
                        while (u != 0) {
                            u = mate[link[v]];
                            mate[v] = link[v];
                            mate[mate[v]] = v;
                            v = u;
                        }
                        return true;
                    } else {
                        Q[tail++] = mate[v];
                        label[mate[v]] = 2;
                    }
                } else if (label[v] == 2 && find(v) != find(u)) {
                    int w = LCA(u, v);
                    blossom(u, v, w);
                    blossom(v, u, w);
                }
            }
        }
        return false;
    }
    int maximum_matching(int n) {
        int res = 0;
        for (int u = 1; u <= n; ++u) {
            if (mate[u] == 0 && BFS(u, n)) { ++res; }
        }
        return res;
    }
}// namespace Blossom
