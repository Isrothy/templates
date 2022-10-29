mt19937_64 mt_rand(0);
struct treap {
    double tag;
    int pos, size;
    unsigned long long pri;
    treap *ch[2];
    void push_up() {
        size = 1;
        if (ch[0] != nullptr) { size += ch[0]->size; }
        if (ch[1] != nullptr) { size += ch[1]->size; }
    }
};
double tag[M];
char S[M];
void rotate(treap *&p, bool f) {
    treap *q = p->ch[f];
    p->ch[f] = q->ch[!f];
    q->ch[!f] = p;
    p->push_up();
    p = q;
}
int Size(treap *p) { return p == nullptr ? 0 : p->size; }
bool comp(int i, int j) { return S[i] == S[j] ? tag[i - 1] < tag[j - 1] : S[i] < S[j]; }
bool comp(char *S, char *Q, int l1, int l2) {
    for (int i = 0; i < l1 && i < l2; ++i) {
        if (S[l1 - i] != Q[l2 - i]) { return S[l1 - i] < Q[l2 - i]; }
    }
    return l1 < l2;
}
void re_tag(treap *p, double l, double r) {
    if (p == nullptr) { return; }
    double mid = (l + r) * 0.5;
    tag[p->pos] = p->tag = mid;
    re_tag(p->ch[0], l, mid);
    re_tag(p->ch[1], mid, r);
}
void insert(treap *&p, int i, double l, double r) {
    double mid = (l + r) * 0.5;
    if (p == nullptr) {
        tag[i] = mid;
        *p = new (treap){mid, i, 1, mt_rand(), {nullptr, nullptr}};
        return;
    }
    bool f = comp(p->pos, i);
    if (f) { insert(p->ch[1], i, mid, r);
    } else { insert(p->ch[0], i, l, mid); }
    if (p->ch[f]->pri < p->pri) {
        rotate(p, f);
        re_tag(p, l, r);
    }
    p->push_up();
}
void remove(treap *&p, int i, double l, double r) {
    double mid = (l + r) * 0.5;
    if (p->pos == i) {
        if (p->ch[0] == nullptr || p->ch[1] == nullptr) {
            p = p->ch[p->ch[0] == nullptr];
            if (p != nullptr) { re_tag(p, l, r); }
            return;
        }
        bool f = p->ch[1]->pri < p->ch[0]->pri;
        rotate(p, f);
        if (f) {
            remove(p->ch[!f], i, l, mid);
            re_tag(p->ch[f], mid, r);
        } else {
            remove(p->ch[!f], i, mid, r);
            re_tag(p->ch[f], l, mid);
        }
        tag[p->pos] = p->tag = mid;
    } else {
        bool f = comp(p->pos, i);
        if (f) { remove(p->ch[f], i, mid, r);
        } else { remove(p->ch[f], i, l, mid); }
    }
    p->push_up();
}