bool cmp_x(Point const &A, Point const &B) {
    return A.x == B.x ? A.y < B.y : A.x < B.x;
}
bool cmp_y(Point const &A, Point const &B) {
    return A.y == B.y ? A.x < B.x : A.y < B.y;
}
void divide_and_conquer(Point *A, int l, int r, double &res) {
    if (r - l < 3) {
        for (int i = l; i < r; ++i) {
            for (int j = i + 1; j < r; ++j) {
                res = min(res, (A[i] - A[j]).len2());
            }
        }
        sort(A + l, A + r, cmp_y);
        return;
    }
    int mid = (l + r) >> 1;
    double xmid = A[mid].x;
    divide_and_conquer(A, l, mid, res);
    divide_and_conquer(A, mid, r, res);
    inplace_merge(A + l, A + mid, A + r, cmp_y);
    static Point B[M];
    int k = 0;
    for (int i = l; i < r; ++i) {
        if ((xmid - A[i].x) * (xmid - A[i].x) >= res) {
            continue;
        }
        for (int j = k - 1; 0 <= j && (A[i].y - B[j].y) * (A[i].y - B[j].y) < res; --j) {
            res = min(res, (A[i] - B[j]).len2());
        }
        B[k++] = A[i];
    }
}
double closest_pair_of_points(Point *A, int n) {
    sort(A, A + n, cmp_x);
    double res = 1e30;
    divide_and_conquer(A, 0, n, res);
    return sqrt(res);
}
