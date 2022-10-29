const double alpha = 0.7;

struct kd_tree {
    int x, y, val, sum, sz;
    int x_max, x_min, y_max, y_min;
    kd_tree *ch[2];
    void push_up() {
        sz = 1;
        sum = val;
        x_min = x_max = x;
        y_min = y_max = y;
        for (int k = 0; k < 2; ++k) {
            if (ch[k] != nullptr) {
                sz += ch[k]->sz;
                sum += ch[k]->sum;
                x_min = min(x_min, ch[k]->x_min);
                x_max = max(x_max, ch[k]->x_max);
                y_min = min(y_min, ch[k]->y_min);
                y_max = max(y_max, ch[k]->y_max);
            }
        }
    }
};
kd_tree pool[M], *allc = pool;
kd_tree *build(kd_tree **p, int l, int r, int D = 1) {
    if (r < l)
        return nullptr;
    int mid = (l + r) >> 1;
    nth_element(p + l, p + mid, p + r + 1, [=](kd_tree *p, kd_tree *q) {
        if (D == 1) {
            return p->x == q->x ? p->y < q->y : p->x < q->x;
        } else {
            return p->y == q->y ? p->x < q->x : p->y < q->y;
        }
    });
    kd_tree *root = p[mid];
    root->ch[0] = build(p, l, mid - 1, D ^ 1);
    root->ch[1] = build(p, mid + 1, r, D ^ 1);
    root->push_up();
    return root;
}
int query(kd_tree *p, int x1, int y1, int x2, int y2) {
    if (p == nullptr) {
        return 0;
    }
    if (x2 < p->x_min || x1 > p->x_max || y2 < p->y_min || y1 > p->y_max) {
        return 0;
    }
    if (x1 <= p->x_min && p->x_max <= x2 && y1 <= p->y_min && p->y_max <= y2) {
        return p->sum;
    }
    int s = 0;
    if (x1 <= p->x && p->x <= x2 && y1 <= p->y && p->y <= y2) {
        s += p->val;
    }
    return s + query(p->ch[0], x1, y1, x2, y2) + query(p->ch[1], x1, y1, x2, y2);
}
kd_tree *rebuild(kd_tree *p, int D) {
    static kd_tree *Q[M];
    int head = 0, tail = 0;
    Q[tail++] = p;
    while (head < tail) {
        kd_tree *q = Q[head++];
        if (q->ch[0] != nullptr) {
            Q[tail++] = q->ch[0];
        }
        if (q->ch[1] != nullptr) {
            Q[tail++] = q->ch[1];
        }
    }
    return build(Q, 0, tail - 1, D);
}
void insert(kd_tree *&p, int x, int y, int a, int D = 1) {
    if (p == nullptr) {
        p = allc++;
        *p = (kd_tree){x, y, a, a, 1, x, x, y, y, {nullptr, nullptr}};
        return;
    }
    if (p->x == x && p->y == y) {
        p->val += a;
        p->sum += a;
        return;
    }
    bool f;
    if (D == 1) {
        f = x == p->x ? p->y < y : p->x < x;
    } else {
        f = y == p->y ? p->x < x : p->y < y;
    }
    insert(p->ch[f], x, y, a, D ^ 1);
    p->push_up();
    if (p->sz * alpha < p->ch[f]->sz) {
        p = rebuild(p, D);
    }
}