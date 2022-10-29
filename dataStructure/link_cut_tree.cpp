struct LCT {
    static const int M = 300005;
    int par[M], ch[M][2];
    bool rev[M];
    bool is_root(int p) { return ch[par[p]][0] != p && ch[par[p]][1] != p; }
    bool dir(int p) { return ch[par[p]][1] == p; }
    void update_rev(int p) {
        rev[p] ^= true;
        swap(ch[p][0], ch[p][1]);
    }
    void push_up(int p) {}
    void push_down(int p) {
        if (rev[p]) {
            if (ch[p][0] != 0) { update_rev(ch[p][0]); }
            if (ch[p][1] != 0) { update_rev(ch[p][1]); }
            rev[p] = false;
        }
    }
    void rotate(int p) {
        bool f = dir(p);
        int q = par[p], r = ch[p][!f];
        if (!is_root(q)) { ch[par[q]][dir(q)] = p; }
        par[p] = par[q];
        if (r != 0){ par[r] = q; }
        ch[q][f] = r;
        par[q] = p;
        ch[p][!f] = q;
        push_up(q);
    }
    void splay(int p) {
        static int stk[M];
        int top = 0;
        for (int q = p; !is_root(q); q = par[q]){ stk[top++] = par[q]; }
        while (top != 0){ push_down(stk[--top]); }
        push_down(p);
        while (!is_root(p)) {
            int q = par[p];
            if (!is_root(q)){ rotate(dir(p) == dir(q) ? q : p); }
            rotate(p);
        }
        push_up(p);
    }
    void access(int p) {
        int q = 0;
        while (p != 0) {
            splay(p);
            ch[p][1] = q;
            push_up(p);
            q = p;
            p = par[p];
        }
    }
    void make_root(int p) {
        access(p);
        splay(p);
        update_rev(p);
    }
    void link(int p, int q) {
        make_root(p);
        par[p] = q;
    }
    void cut(int p, int q) {
        make_root(p);
        access(q);
        splay(q);
        par[ch[q][0]] = 0;
        ch[q][0] = 0;
        push_up(q);
    }
    int find_root(int p) {
        access(p);
        splay(p);
        while (ch[p][0] != 0) {
            push_down(p);
            p = ch[p][0];
        }
        splay(p);
        return p;
    }
    bool cotree(int p, int q) {
        return find_root(p) == find_root(q);
    }
};
