struct avl {
    int val, size, height;
    avl *ch[2];
    void push_up() {
        size = height = 1;
        if (ch[0] != nullptr) {
            size += ch[0]->size;
            height = max(height, ch[0]->height + 1);
        }
        if (ch[1] != nullptr) {
            size += ch[1]->size;
            height = max(height, ch[1]->height + 1);
        }
    }
};
avl pool[M], *allc = pool;
int Size(avl *p) {
    return p == nullptr ? 0 : p->size;
}
int Height(avl *p) {
    return p == nullptr ? 0 : p->height;
}
avl *rotate(avl *p, bool f) {
    avl *q = p->ch[f];
    p->ch[f] = q->ch[!f];
    q->ch[!f] = p;
    p->push_up();
    return q;
}
avl *maintain(avl *p) {
    if (Height(p->ch[0]) - Height(p->ch[1]) == 2) {
        if (Height(p->ch[0]->ch[0]) >= Height(p->ch[0]->ch[1])) {
            p = rotate(p, 0);
        } else {
            p->ch[0] = rotate(p->ch[0], 1);
            p = rotate(p, 0);
        }
    } else if (Height(p->ch[1]) - Height(p->ch[0]) == 2) {
        if (Height(p->ch[1]->ch[1]) >= Height(p->ch[1]->ch[0])) {
            p = rotate(p, 1);
        } else {
            p->ch[1] = rotate(p->ch[1], 0);
            p = rotate(p, 1);
        }
    }
    p->push_up();
    return p;
}
avl *insert(avl *p, int x) {
    if (p == nullptr) {
        *allc = (avl) {x, 1, 1, {nullptr, nullptr}};
        return allc++;
    }
    bool f = p->val < x;
    p->ch[f] = insert(p->ch[f], x);
    p->push_up();
    return maintain(p);
}
avl *erase(avl *p, int x) {
    if (p->val == x) {
        if (p->ch[0] == nullptr || p->ch[1] == nullptr) {
            return p->ch[p->ch[0] == nullptr];
        }
        avl *q = p->ch[1];
        while (q->ch[0] != nullptr) {
            q = q->ch[0];
        }
        p->val = q->val;
        p->ch[1] = erase(p->ch[1], q->val);
    } else {
        bool f = p->val < x;
        p->ch[f] = erase(p->ch[f], x);
    }
    p->push_up();
    return maintain(p);
}
