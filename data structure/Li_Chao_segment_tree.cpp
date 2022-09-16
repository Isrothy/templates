struct segment {
    long long x1, y1, x2, y2;
    int id;
};

segment S[4 * M];

bool cmp(segment const &A, segment const &B, long long x0) {
    long long h1
        = (B.x2 - B.x1) * (x0 - A.x1) * (A.y2 - A.y1) + A.y1 * (A.x2 - A.x1) * (B.x2 - B.x1);
    long long h2
        = (A.x2 - A.x1) * (x0 - B.x1) * (B.y2 - B.y1) + B.y1 * (B.x2 - B.x1) * (A.x2 - A.x1);
    return h1 == h2 ? A.id > B.id : h1 < h2;
}

void update(int p, int l, int r, int a, int b, segment L) {
    int mid = (l + r) >> 1;
    if (l == a && r == b) {
        if (S[p].id == 0) {
            S[p] = L;
        } else {
            if (cmp(S[p], L, mid)) {
                swap(S[p], L);
            }
            if (l != r) {
                if (cmp(S[p], L, l)) {
                    update(p << 1, l, mid, a, mid, L);
                }
                if (cmp(S[p], L, r)) {
                    update(p << 1 | 1, mid + 1, r, mid + 1, b, L);
                }
            }
        }
        return;
    }
    if (b <= mid) {
        update(p << 1, l, mid, a, b, L);
    } else if (mid < a) {
        update(p << 1 | 1, mid + 1, r, a, b, L);
    } else {
        update(p << 1, l, mid, a, mid, L);
        update(p << 1 | 1, mid + 1, r, mid + 1, b, L);
    }
}

segment query(int p, int l, int r, long long x) {
    if (l == r) {
        return S[p];
    }
    int mid = (l + r) >> 1;
    segment res;
    if (x <= mid) {
        res = query(p << 1, l, mid, x);
    } else {
        res = query(p << 1 | 1, mid + 1, r, x);
    }
    if (S[p].id != 0 && (res.id == 0 || cmp(res, S[p], x))) {
        res = S[p];
    }
    return res;
}
