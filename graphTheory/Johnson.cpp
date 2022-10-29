bool Johnson(int n, vector<pair<int, int>> *E, long long dis[M][M]) {
    queue<int> Q;
    static int h[M], cnt[M];
    static bool in_queue[M];
    for (int u = 1; u <= n; ++u) {
        h[u] = cnt[u] = 0;
        in_queue[u] = true;
        Q.push(u);
    }
    while (!Q.empty()) {
        int u = Q.front();
        Q.pop();
        in_queue[u] = false;
        for (auto e: E[u]) {
            int v = e.first;
            long long d = h[u] + e.second;
            if (d < h[v]) {
                h[v] = d;
                if (!in_queue[v]) {
                    in_queue[v] = true;
                    Q.push(v);
                    if (++cnt[v] == n) { return false; }
                }
            }
        }
    }
    for (int i = 1; i <= n; ++i) {
        priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> heap;
        for (int u = 1; u <= n; ++u) { dis[i][u] = INF; }
        dis[i][i] = 0;
        heap.push(make_pair(0, i));
        while (!heap.empty()) {
            pair<int, int> p = heap.top();
            heap.pop();
            int u = p.second;
            if (dis[i][u] < p.first) { continue; }
            for (auto e: E[u]) {
                int v = e.first;
                long long d = dis[i][u] + h[u] - h[v] + e.second;
                if (d < dis[i][v]) {
                    dis[i][v] = d;
                    heap.push(make_pair(d, e.first));
                }
            }
        }
        for (int u = 1; u <= n; ++u) {
            if (dis[i][u] != INF) { dis[i][u] += h[u] - h[i]; }
        }
    }
    return true;
}
