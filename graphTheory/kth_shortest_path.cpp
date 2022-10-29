namespace Kth_SSP {
    struct heap {
        int val, to, dis;
        heap *ch[2];
    };
    struct node {
        heap *h;
        int val;

        bool operator>(node const &_) const {
            return val > _.val;
        }
    };
    struct edge {
        int to, id, cost;
    };
    heap pool[17 * M], *allc, *H[N], *tmp[N];
    vector<edge> E[N], R[N];
    int A[N], dis[N], pre[N], deg[N], tree_edge[N];
    int edge_cnt;
    void add_edge(int u, int v, int c) {
        ++edge_cnt;
        E[u].push_back((edge){v, edge_cnt, c});
        R[v].push_back((edge){u, edge_cnt, c});
    }
    heap *merge(heap *p, heap *q) {
        if (p == nullptr) {
            return q;
        }
        if (q == nullptr) {
            return p;
        }
        if (q->val < p->val) {
            swap(p, q);
        }
        heap *r = allc++;
        *r = *p;
        r->ch[1] = merge(r->ch[1], q);
        if (r->ch[0] == nullptr || r->ch[0]->dis < r->ch[1]->dis) {
            swap(r->ch[0], r->ch[1]);
        }
        r->dis = r->ch[0]->dis + 1;
        return r;
    }
    int Dijkstra(int S, int T, int n) {
        for (int u = 1; u <= n; ++u) {
            dis[u] = INF;
            pre[u] = 0;
        }
        dis[T] = 0;
        priority_queue<pair<double, int>, vector<pair<double, int>>, greater<pair<double, int>>> Q;
        Q.push(make_pair(0, T));
        int sz = 0;
        while (!Q.empty()) {
            pair<double, int> p = Q.top();
            Q.pop();
            int u = p.second;
            if (dis[u] != p.first)
                continue;
            A[sz++] = u;
            for (auto e: R[u]) {
                int v = e.to, d = p.first + e.cost;
                if (d < dis[v]) {
                    dis[v] = d;
                    pre[v] = u;
                    tree_edge[v] = e.id;
                    Q.push(make_pair(d, v));
                }
            }
        }
        return sz;
    }
    int calc(int S, int T, int n, int k) {
        int m = Dijkstra(S, T, n);
        if (dis[S] == INF) {
            return -1;
        }
        if (S != T && --k == 0) {
            return dis[S];
        }
        allc = pool;
        for (int i = 0; i < m; ++i) {
            int u = A[i], sz = 0;
            for (auto e: E[u]) {
                if (dis[e.to] < INF && e.id != tree_edge[u]) {
                    tmp[++sz] = allc++;
                    *tmp[sz] = (heap){dis[e.to] + e.cost - dis[u], e.to, 1, {nullptr, nullptr}};
                }
            }
            make_heap(tmp + 1, tmp + sz + 1, [](heap *p, heap *q) {
                return p->val > q->val;
            });
            for (int j = sz; 0 < j; --j) {
                if ((j << 1) <= sz) {
                    tmp[j]->ch[0] = tmp[j << 1];
                }
                if ((j << 1 | 1) <= sz) {
                    tmp[j]->ch[1] = tmp[j << 1 | 1];
                }
                tmp[j]->dis = tmp[j]->ch[0] == nullptr ? 1 : tmp[j]->ch[0]->dis + 1;
            }
            H[u] = sz == 0 ? nullptr : tmp[1];
            if (pre[u] != 0) {
                H[u] = merge(H[u], H[pre[u]]);
            }
        }
        priority_queue<node, vector<node>, greater<node>> Q;
        if (H[S] != nullptr) {
            Q.push((node){H[S], dis[S] + H[S]->val});
        }
        while (!Q.empty()) {
            node p = Q.top();
            Q.pop();
            if (--k == 0) {
                return p.val;
            }
            if (p.h->ch[0] != nullptr) {
                Q.push((node){p.h->ch[0], p.val - p.h->val + p.h->ch[0]->val});
            }
            if (p.h->ch[1] != nullptr) {
                Q.push((node){p.h->ch[1], p.val - p.h->val + p.h->ch[1]->val});
            }
            int v = p.h->to;
            if (H[v] != nullptr) {
                Q.push((node){H[v], p.val + H[v]->val});
            }
        }
        return -1;
    }
}
