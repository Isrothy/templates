int half_planes_intersection(pair<Point, Point> *A, Point *B, int n) {
    static Point t[M];
    static pair<Point, Point> Q[M];
    sort(A, A + n, [](pair<Point, Point> l1, pair<Point, Point> l2) {
        int d = dcmp((l1.second - l1.first).angle(), (l2.second - l2.first).angle());
        return d == 0 ? dcmp(det(l2.first - l1.first, l2.second - l1.first)) > 0 : d < 0;
    });
    int l = 0, r = 0;
    Q[r++] = A[0];
    for (int i = 1; i < n; ++i) {
        if (dcmp(det(Q[r - 1].first - Q[r - 1].second, A[i].first - A[i].second)) == 0){
            continue;
        }
        while (r - l > 1 && dcmp(det(A[i].first - t[r - 1], A[i].second - t[r - 1])) <= 0) {
            --r;
        }
        while (r - l > 1 && dcmp(det(A[i].first - t[l + 1], A[i].second - t[l + 1])) <= 0) {
            ++l;
        }
        intersection(Q[r - 1].first, Q[r - 1].second, A[i].first, A[i].second, t[r]);
        Q[r++] = A[i];
    }
    while (r - l > 1 && dcmp(det(Q[l].first - t[r - 1], Q[l].second - t[r - 1])) <= 0) {
        --r;
    }
    intersection(Q[l].first, Q[l].second, Q[r - 1].first, Q[r - 1].second, t[l]);
    for (int i = l; i < r; ++i) {
        A[i - l] = Q[i];
        B[i - l] = t[i];
    }
    return r - l;
}
