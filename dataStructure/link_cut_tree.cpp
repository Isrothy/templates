struct lct {
    unsigned val;
    unsigned sum;
    bool rev;
    lct *ch[2], *fa;
    lct() = default;
    explicit lct(unsigned val)
        : val(val), sum(val), rev(false), ch{nullptr, nullptr}, fa(nullptr) {}
    lct *push_up() {
        sum = val;
        if (ch[0] != nullptr) {
            sum ^= ch[0]->sum;
        }
        if (ch[1] != nullptr) {
            sum ^= ch[1]->sum;
        }
        return this;
    }
    bool dir() {
        return fa->ch[1] == this;
    }
    bool is_root() {
        return fa == nullptr || (fa->ch[0] != this && fa->ch[1] != this);
    }
    lct *update_rev() {
        rev ^= 1;
        std::swap(ch[0], ch[1]);
        return this;
    }
    lct *push_down() {
        if (rev) {
            if (ch[0] != nullptr) {
                ch[0]->update_rev();
            }
            if (ch[1] != nullptr) {
                ch[1]->update_rev();
            }
            rev = false;
        }
        return this;
    }
    lct *rotate(bool f) {
        lct *q = ch[f];
        ch[f] = q->ch[!f];
        if (ch[f] != nullptr) {
            ch[f]->fa = this;
        }
        q->fa = fa;
        if (!is_root()) {
            fa->ch[dir()] = q;
        }
        q->ch[!f] = this;
        fa = q;
        push_up();
        return q;
    }
    lct *splay() {
        std::stack<lct *> stk;
        for (lct *q = this; !q->is_root(); q = q->fa) {
            stk.push(q->fa);
        }
        while (!stk.empty()) {
            stk.top()->push_down();
            stk.pop();
        }
        push_down();
        while (!is_root()) {
            bool f = dir();
            if (!fa->is_root()) {
                (fa->dir() == f ? fa->fa : fa)->rotate(f);
            }
            fa->rotate(dir());
        }
        return push_up();
    }
    lct *access() {
        lct *p = nullptr;
        for (lct *q = this; q != nullptr; q = q->fa) {
            q->splay();
            q->ch[1] = p;
            q->push_up();
            p = q;
        }
        return this;
    }
    lct *make_root() {
        return access()->splay()->update_rev();
    }
    lct *find_root() {
        access()->splay();
        lct *q = this;
        while (q->ch[0] != nullptr) {
            q = q->ch[0];
        }
        return q->splay();
    }
};

void link(lct *p, lct *q) {
    if (p->find_root() != q->find_root()) {
        p->make_root()->fa = q;
    }
}

void cut(lct *p, lct *q) {
    p->make_root();
    q->access()->splay();
    lct *r = q->ch[0];
    if (r == nullptr) {
        return;
    }
    while (r->ch[1] != nullptr) {
        r = r->ch[1];
    }
    if (r == p) {
        q->ch[0]->fa = nullptr;
        q->ch[0] = nullptr;
        q->push_up();
    }
    r->splay();
}

unsigned query(lct *p, lct *q) {
    p->make_root();
    q->access()->splay();
    return q->sum;
}

void update(lct *p, unsigned val) {
    p->make_root();
    p->val = val;
    p->push_up();
}
