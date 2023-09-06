#include <algorithm>
constexpr int M = 100000 + 10;
constexpr int INF = 0x3f3f3f3f;
int mx[4 * M], mi[4 * M], second_max[4 * M], second_min[4 * M], cnt_max[4 * M], cnt_min[4 * M], len[4 * M];
int lazy_min[4 * M], lazy_max[4 * M], lazy_add[4 * M];
long long sum[4 * M];
void push_up(int p) {
    sum[p] = sum[p << 1] + sum[p << 1 | 1];
    mx[p] = std::max(mx[p << 1], mx[p << 1 | 1]);
    mi[p] = std::min(mi[p << 1], mi[p << 1 | 1]);
    second_max[p] = std::max(second_max[p << 1], second_max[p << 1 | 1]);
    second_min[p] = std::min(second_min[p << 1], second_min[p << 1 | 1]);
    cnt_max[p] = cnt_min[p] = 0;
    if (mx[p] == mx[p << 1]) {
        cnt_max[p] += cnt_max[p << 1];
    } else {
        second_max[p] = std::max(second_max[p], mx[p << 1]);
    }
    if (mx[p] == mx[p << 1 | 1]) {
        cnt_max[p] += cnt_max[p << 1 | 1];
    } else {
        second_max[p] = std::max(second_max[p], mx[p << 1 | 1]);
    }
    if (mi[p] == mi[p << 1]) {
        cnt_min[p] += cnt_min[p << 1];
    } else {
        second_min[p] = std::min(second_min[p], mi[p << 1]);
    }
    if (mi[p] == mi[p << 1 | 1]) {
        cnt_min[p] += cnt_min[p << 1 | 1];
    } else {
        second_min[p] = std::min(second_min[p], mi[p << 1 | 1]);
    }
}
void add(int p, int x) {
    mx[p] += x;
    mi[p] += x;
    sum[p] += (long long) x * len[p];
    if (second_max[p] != -INF) { second_max[p] += x; }
    if (second_min[p] != INF) { second_min[p] += x; }
    lazy_add[p] += x;
    if (lazy_min[p] != INF) { lazy_min[p] += x; }
    if (lazy_max[p] != -INF) { lazy_max[p] += x; }
}
void check_max(int p, int x) {
    if (x <= mi[p]) { return; }
    sum[p] += (long long) cnt_min[p] * (x - mi[p]);
    if (mx[p] == mi[p]) {
        mx[p] = lazy_min[p] = x;
    } else if (second_max[p] == mi[p]) {
        second_max[p] = x;
    }
    mi[p] = lazy_max[p] = x;
}
void check_min(int p, int x) {
    if (mx[p] <= x) { return; }
    sum[p] -= (long long) cnt_max[p] * (mx[p] - x);
    if (mx[p] == mi[p]) {
        mi[p] = lazy_max[p] = x;
    } else if (second_min[p] == mx[p]) {
        second_min[p] = x;
    }
    mx[p] = lazy_min[p] = x;
}
void push_down(int p) {
    if (lazy_add[p]) {
        add(p << 1, lazy_add[p]);
        add(p << 1 | 1, lazy_add[p]);
        lazy_add[p] = 0;
    }
    if (lazy_max[p] != -INF) {
        check_max(p << 1, lazy_max[p]);
        check_max(p << 1 | 1, lazy_max[p]);
        lazy_max[p] = -INF;
    }
    if (lazy_min[p] != INF) {
        check_min(p << 1, lazy_min[p]);
        check_min(p << 1 | 1, lazy_min[p]);
        lazy_min[p] = INF;
    }
}
void build(int p, int l, int r, int *A) {
    len[p] = r - l + 1;
    lazy_add[p] = 0;
    lazy_max[p] = -INF;
    lazy_min[p] = INF;
    if (l == r) {
        sum[p] = mx[p] = mi[p] = A[l];
        cnt_max[p] = cnt_min[p] = 1;
        lazy_add[p] = 0;
        second_max[p] = -INF;
        second_min[p] = INF;
        return;
    }
    int mid = (l + r) >> 1;
    build(p << 1, l, mid, A);
    build(p << 1 | 1, mid + 1, r, A);
    push_up(p);
}
void update_add(int p, int l, int r, int a, int b, int x) {
    if (l == a && r == b) {
        add(p, x);
        return;
    }
    int mid = (l + r) >> 1;
    push_down(p);
    if (b <= mid) {
        update_add(p << 1, l, mid, a, b, x);
    } else if (mid < a) {
        update_add(p << 1 | 1, mid + 1, r, a, b, x);
    } else {
        update_add(p << 1, l, mid, a, mid, x);
        update_add(p << 1 | 1, mid + 1, r, mid + 1, b, x);
    }
    push_up(p);
}
void check_max(int p, int l, int r, int a, int b, int x) {
    if (x <= mi[p]) { return; }
    if (l == a && r == b && x < second_min[p]) {
        check_max(p, x);
        return;
    }
    int mid = (l + r) >> 1;
    push_down(p);
    if (b <= mid) {
        check_max(p << 1, l, mid, a, b, x);
    } else if (mid < a) {
        check_max(p << 1 | 1, mid + 1, r, a, b, x);
    } else {
        check_max(p << 1, l, mid, a, mid, x);
        check_max(p << 1 | 1, mid + 1, r, mid + 1, b, x);
    }
    push_up(p);
}
void check_min(int p, int l, int r, int a, int b, int x) {
    if (mx[p] <= x) { return; }
    if (l == a && r == b && second_max[p] < x) {
        check_min(p, x);
        return;
    }
    int mid = (l + r) >> 1;
    push_down(p);
    if (b <= mid) {
        check_min(p << 1, l, mid, a, b, x);
    } else if (mid < a) {
        check_min(p << 1 | 1, mid + 1, r, a, b, x);
    } else {
        check_min(p << 1, l, mid, a, mid, x);
        check_min(p << 1 | 1, mid + 1, r, mid + 1, b, x);
    }
    push_up(p);
}
