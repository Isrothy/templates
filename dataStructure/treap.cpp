struct treap {
    int val, size;
    unsigned long long pri;
    treap *ch[2];
    void push_up() {
        size = 1;
        if (ch[0] != nullptr) {
            size += ch[0]->size;
        }
        if (ch[1] != nullptr) {
            size += ch[1]->size;
        }
    }
};
typedef pair<treap *, treap *> ptt;
treap pool[M], *allc = pool;
mt19937_64 mt_rand(time(NULL));
int Size(treap *p) {
    return p == nullptr ? 0 : p->size;
}
treap *rotate(treap *p, bool f) {
    treap *q = p->ch[f];
    p->ch[f] = q->ch[!f];
    q->ch[!f] = p;
    p->push_up();
    return q;
}
treap *insert(treap *p, int x) {
    if (p == nullptr) {
        *allc = (treap) {x, 1, mt_rand(), {nullptr, nullptr}};
        return allc++;
    }
    bool f = p->val < x;
    p->ch[f] = insert(p->ch[f], x);
    if (p->ch[f]->pri < p->pri) {
        p = rotate(p, f);
    }
    p->push_up();
    return p;
}
treap *erase(treap *p, int x) {
    if (x == p->val) {
        if (p->ch[0] == nullptr || p->ch[1] == nullptr) {
            return p->ch[p->ch[0] == nullptr];
        }
        bool f = p->ch[1]->pri < p->ch[0]->pri;
        p = rotate(p, f);
        p->ch[!f] = erase(p->ch[!f], x);
    } else {
        bool f = p->val < x;
        p->ch[f] = erase(p->ch[f], x);
    }
    p->push_up();
    return p;
}
int index(treap *p, int x) {
    int res = 1;
    while (p != nullptr) {
        if (p->val < x) {
            res += Size(p->ch[0]) + 1;
            p = p->ch[1];
        } else {
            p = p->ch[0];
        }
    }
    return res;
}
int kth(treap *p, int k) {
    for (;;) {
        int s = Size(p->ch[0]);
        if (s + 1 == k) {
            return p->val;
        }
        if (k <= s) {
            p = p->ch[0];
        }
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
        } else {
            p = p->ch[0];
        }
    }
    return res;
}
int nxt(treap *p, int x) {
    int res = -1;
    while (p != nullptr) {
        if (x < p->val) {
            res = p->val;
            p = p->ch[0];
        } else {
            p = p->ch[1];
        }
    }
    return res;
}
treap *merge(treap *p, treap *q) {
    if (p == nullptr) {
        return q;
    }
    if (q == nullptr) {
        return p;
    }
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
    if (p == nullptr) {
        return make_pair(nullptr, nullptr);
    }
    p->push_down();
    if (k <= Size(p->ch[0])) {
        ptt o = split(p->ch[0], k);
        p->ch[0] = o.second;
        p->push_up();
        return make_pair(o.first, p);
    } else {
        ptt o = split(p->ch[1], k - Size(p->ch[0]) - 1);
        p->ch[1] = o.first;
        p->push_up();
        return make_pair(p, o.second);
    }
}
ptt split_by_value(treap *p, int v) {
    if (p == nullptr) {
        return make_pair(nullptr, nullptr);
    }
    if (v < p->val) {
        ptt o = split_by_value(p->ch[0], v);
        p->ch[0] = o.second;
        p->push_up();
        return make_pair(o.first, p);
    } else {
        ptt o = split_by_value(p->ch[1], v);
        p->ch[1] = o.first;
        p->push_up();
        return make_pair(p, o.second);
    }
}
treap *heuristic_merge(treap *p, treap *q) {
    if (p == nullptr) {
        return q;
    }
    if (q == nullptr) {
        return p;
    }
    if (p->pri < q->pri) {
        swap(p, q);
    }
    ptt o = split_by_value(p, q->val);
    q->ch[0] = heuristic_merge(q->ch[0], o.first);
    q->ch[1] = heuristic_merge(q->ch[1], o.second);
    q->push_up();
    return q;
}
