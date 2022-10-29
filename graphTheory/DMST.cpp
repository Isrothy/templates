namespace Edmonds {
    struct heap {
        int from, c, lazy;
        heap *ch[2];
        void push_down() {
            if (lazy != 0) {
                if (ch[0] != nullptr) {
                    ch[0]->lazy += lazy;
                    ch[0]->c += lazy;
                }
                if (ch[1] != nullptr) {
                    ch[1]->lazy += lazy;
                    ch[1]->c += lazy;
                }
            }
            lazy = 0;
        }
    };
    heap pool[M], *allc = pool, *H[M];
    int fa[M], mark[M];
    heap *merge(heap *p, heap *q) {
        if (p == nullptr) { return q; }
        if (q == nullptr) { return p; }
        if (q->c < p->c) { swap(p, q); }
        p->push_down();
        p->ch[1] = merge(p->ch[1], q);
        swap(p->ch[0], p->ch[1]);
        return p;
    }
    void pop(heap *&p) {
        p->push_down();
        p = merge(p->ch[0], p->ch[1]);
    }
    void add_edge(int u, int v, int c) {
        if (u == v) { return; }
        *allc = (heap) {u, c, 0, {nullptr, nullptr}};
        H[v] = merge(H[v], allc);
        ++allc;
    }
    int find(int u) {
        if (fa[u] == u) { return u; }
        fa[u] = find(fa[u]);
        return fa[u];
    }

    int DMST(int n, int root) {
        for (int i = 1; i <= n; ++i) { fa[i] = i; }
        mark[root] = root;
        int res = 0;
        for (int i = 1; i <= n; ++i) {
            if (mark[i] != 0 || i == root) continue;
            int u = i;
            mark[u] = i;
            for (;;) {
                while (H[u] != nullptr && find(H[u]->from) == find(u)) { pop(H[u]); }
                if (H[u] == nullptr) { return -1; }
                int v = find(H[u]->from);
                res += H[u]->c;
                if (mark[v] == i) {
                    H[u]->lazy -= H[u]->c;
                    pop(H[u]);
                    while (v != u) {
                        int w = find(H[v]->from);
                        H[v]->lazy -= H[v]->c;
                        pop(H[v]);
                        fa[find(v)] = u;
                        H[u] = merge(H[u], H[v]);
                        v = w;
                    }
                } else if (mark[v] == 0) {
                    u = v;
                    mark[u] = i;
                } else { break; }
            }
        }
        return res;
    }
}
