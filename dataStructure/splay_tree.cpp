struct splay_node {
    int val, size;
    bool rev;
    splay_node *ch[2], *fa;
    bool dir() {
        return fa->ch[1] == this;
    }
    void update_rev() {
        rev ^= 1;
        std::swap(ch[0], ch[1]);
    }
    void push_up() {
        size = 1;
        if (ch[0] != nullptr) {
            size += ch[0]->size;
        }
        if (ch[1] != nullptr) {
            size += ch[1]->size;
        }
    }
    void push_down() {
        if (rev) {
            if (ch[0] != nullptr) {
                ch[0]->update_rev();
            }
            if (ch[1] != nullptr) {
                ch[1]->update_rev();
            }
            rev = false;
        }
    }
};

void rotate(splay_node *p, bool f) {
    splay_node *q = p->ch[f];
    p->ch[f] = q->ch[!f];
    if (q->ch[!f] != nullptr) {
        q->ch[!f]->fa = p;
    }
    q->fa = p->fa;
    if (p->fa != nullptr) {
        p->fa->ch[p->dir()] = q;
    }
    q->ch[!f] = p;
    p->fa = q;
    p->push_up();
}

splay_node *splay_to(splay_node *p, splay_node *t = nullptr) {
    std::stack<splay_node *> st;
    for (splay_node *x = p; x != t; x = x->fa) {
        st.push(x);
    }
    while (!st.empty()) {
        st.top()->push_down();
        st.pop();
    }
    while (p->fa != t) {
        if (p->fa->fa != t && p->dir() == p->fa->dir()) {
            rotate(p->fa->fa, p->dir());
        }
        rotate(p->fa, p->dir());
    }
    p->push_up();
    return p;
}
