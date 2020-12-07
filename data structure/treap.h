struct treap {
    int val, pri, sz;
    treap *ch[2];
    void push_up() {
        sz = 1;
        if (ch[0] != nullptr)
            sz += ch[0]->sz;
        if (ch[1] != nullptr)
            sz += ch[1]->sz;
    }
};
typedef pair<treap*,treap*> ptt;
treap pool[M], *allc = pool;
int n;

int Sz(treap *p) {
    return p == nullptr ? 0 : p->sz;
}

void rotate(treap *&p, bool f) {
    treap *q = p->ch[f];
    p->ch[f] = q->ch[!f];
    q->ch[!f] = p;
    p->push_up();
    p = q;
}

void insert(treap *&p, int x) {
    if (p == nullptr) {
        p = allc++;
        *p = (treap) {x, rand(), 1, {nullptr, nullptr}};
        return;
    }
    bool f = p->val < x;
    insert(p->ch[f], x);
    if (p->ch[f]->pri < p->pri)
        rotate(p, f);
    p->push_up();
}

void remove(treap *&p, int x) {
    if (x == p->val) {
        if (p->ch[0] == nullptr || p->ch[1] == nullptr) {
            p = p->ch[p->ch[0] == nullptr];
            return;
        }
        bool f = p->ch[1]->pri < p->ch[0]->pri;
        rotate(p, f);
        remove(p->ch[!f], x);
    } else
        remove(p->ch[p->val < x], x);
    p->push_up();
}

int index(treap *p, int x) {
    int res = 1;
    while (p != nullptr) {
        if (p->val < x) {
            res += Sz(p->ch[0]) + 1;
            p = p->ch[1];
        } else
            p = p->ch[0];
    }
    return res;
}

int kth(treap *p, int k) {
    for (;;) {
        int s = Sz(p->ch[0]);
        if (s + 1 == k)
            return p->val;
        if (k <= s)
            p = p->ch[0];
        else {
            k -= s + 1;
            p = p->ch[1];
        }
    }
}

int pre(treap *p, int x) {
    int res = -1;
    while (p != nullptr) {
        if (p->val < x) {
            res = p->val;
            p = p->ch[1];
        } else
            p = p->ch[0];
    }
    return res;
}

int nxt(treap *p, int x) {
    int res = -1;
    while (p != nullptr) {
        if (x < p->val) {
            res = p->val;
            p = p->ch[0];
        } else
            p = p->ch[1];
    }
    return res;
}
treap* merge(treap *p, treap *q) {
    if (p == nullptr)
        return q;
    if (q == nullptr)
        return p;
    p->push_down();
    q->push_down();
    if (p->pri < q->pri) {
        p->ch[1] = merge(p->ch[1], q);
        p->push_up();
        return p;
    } else {
        q->ch[0] = merge(p, q->ch[0]);
        q->push_up();
        return q;
    }
}

ptt split(treap *p, int k) {
    if (p == nullptr)
        return make_pair(nullptr, nullptr);
    p->push_down();
    if (k <= Sz(p->ch[0])) {
        ptt o = split(p->ch[0], k);
        p->ch[0] = o.second;
        p->push_up();
        return make_pair(o.first, p);
    } else {
        ptt o = split(p->ch[1], k - Sz(p->ch[0]) - 1);
        p->ch[1] = o.first;
        p->push_up();
        return make_pair(p, o.second);
    }
}
