#include <cstring>
#include <iostream>
#include <random>
std::mt19937_64 mt_rand(std::random_device{}());
struct treap;
size_t safe_size(treap *p);
struct treap {
    double tag;
    int pos;
    size_t size;
    unsigned long long pri;
    treap *ch[2]{};
    treap(double tag, int pos) : tag(tag), pos(pos), size(1), pri(mt_rand()) { ch[0] = ch[1] = nullptr; }
    treap *push_up() {
        size = safe_size(ch[0]) + safe_size(ch[1]) + 1;
        return this;
    }
    treap *rotate(bool f) {
        treap *q = ch[f];
        ch[f] = q->ch[!f];
        q->ch[!f] = this;
        push_up();
        return q;
    }
};
size_t safe_size(treap *p) { return p == nullptr ? 0 : p->size; }
bool suffix_comp(const char *S, const double *tag, int i, int j) {
    if (S[i] == S[j]) {
        return (i == 0 ? 0 : tag[i - 1]) < (j == 0 ? 0 : tag[j - 1]);
    } else {
        return S[i] < S[j];
    }
}
bool reverse_comp(const char *S, const char *Q, int l1, int l2) {
    for (int i = 0; i < l1 && i < l2; ++i) {
        if (S[l1 - i - 1] != Q[l2 - i - 1]) { return S[l1 - i - 1] < Q[l2 - i - 1]; }
    }
    return l1 < l2;
}
void re_tag(double *tags, treap *p, double l, double r) {
    if (p == nullptr) { return; }
    double mid = (l + r) * 0.5;
    tags[p->pos] = p->tag = mid;
    re_tag(tags, p->ch[0], l, mid);
    re_tag(tags, p->ch[1], mid, r);
}
treap *insert(char *S, double *tags, treap *p, int i, double l, double r) {
    double mid = (l + r) * 0.5;
    if (p == nullptr) {
        tags[i] = mid;
        return new treap(mid, i);
    }
    bool f = suffix_comp(S, tags, p->pos, i);
    if (f) {
        p->ch[1] = insert(S, tags, p->ch[1], i, mid, r);
    } else {
        p->ch[0] = insert(S, tags, p->ch[0], i, l, mid);
    }
    if (p->ch[f]->pri < p->pri) {
        p = p->rotate(f);
        re_tag(tags, p, l, r);
    }
    return p->push_up();
}
treap *remove(char *S, double *tags, treap *p, int i, double l, double r) {
    double mid = (l + r) * 0.5;
    if (p->pos == i) {
        if (p->ch[0] == nullptr || p->ch[1] == nullptr) {
            p = p->ch[p->ch[0] == nullptr];
            if (p != nullptr) { re_tag(tags, p, l, r); }
            return p;
        }
        bool f = p->ch[1]->pri < p->ch[0]->pri;
        p = p->rotate(f);
        if (f) {
            p->ch[!f] = remove(S, tags, p->ch[!f], i, l, mid);
            re_tag(tags, p->ch[f], mid, r);
        } else {
            p->ch[!f] = remove(S, tags, p->ch[!f], i, mid, r);
            re_tag(tags, p->ch[f], l, mid);
        }
        tags[p->pos] = p->tag = mid;
    } else {
        bool f = suffix_comp(S, tags, p->pos, i);
        if (f) {
            p->ch[f] = remove(S, tags, p->ch[f], i, mid, r);
        } else {
            p->ch[f] = remove(S, tags, p->ch[f], i, l, mid);
        }
    }
    p->push_up();
    return p;
}
size_t index(treap *p, char *S, char *Q, int len) {
    size_t res = 1;
    while (p != nullptr) {
        bool f = reverse_comp(S, Q, p->pos + 1, len);
        if (f) { res += safe_size(p->ch[0]) + 1; }
        p = p->ch[f];
    }
    return res;
}
