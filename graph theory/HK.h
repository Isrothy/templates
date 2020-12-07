namespace HK {
    vector<int> E[M];
    int dx[M], dy[M], S[M], T[M], Q[M];

    void add_edge(int u, int v) {
        E[u].push_back(v);
    }

    bool BFS(int nl, int nr) {
        memset(dx, 0, sizeof(dx));
        memset(dy, 0, sizeof(dy));
        int head = 0, tail = 0;
        bool res = false;
        for (int u = 1; u <= nl; ++u) {
            if (S[u] == 0)
                Q[tail++] = u;
        }
        while (head < tail) {
            int u = Q[head++];
            for (int i = 0; i < (int) E[u].size(); ++i) {
                int v = E[u][i];
                if (dy[v] == 0) {
                    dy[v] = dx[u] + 1;
                    if (T[v] != 0) {
                        dx[T[v]] = dy[v];
                        Q[tail++] = T[v];
                    } else
                        res = true;
                }
            }
        }
        return res != 0;
    }

    bool DFS(int u) {
        for (int i = 0; i < (int) E[u].size(); ++i) {
            int v = E[u][i];
            if (dy[v] == dx[u] + 1) {
                dy[v] = 0;
                if (T[v] == 0 || DFS(T[v])) {
                    S[u] = v;
                    T[v] = u;
                    return true;
                }
            }
        }
        return false;
    }

    int maximum_matching(int nl, int nr) {
        int res = 0;
        while (BFS(nl, nr)) {
            for (int u = 1; u <= nl; ++u) {
                if (S[u] == 0 && DFS(u))
                    ++res;
            }
        }
        return res;
    }
}
