int Andrew(Point *a, int n) {
    static Point b[M];
    static int stk[M];
    static bool used[M];
    for (int i = 0; i < n; ++i)
        b[i] = a[i];
    sort(b, b + n, [](Point const &A, Point const &B) {
        return A.x == B.x ? A.y < B.y : A.x < B.x;
    });
    int top = 0;
    stk[top++] = 0;
    for (int i = 1; i < n; ++i) {
        while (2 <= top && det(b[stk[top - 1]] - b[stk[top - 2]], b[i] - b[stk[top - 1]]) <= 0)
            used[stk[--top]] = false;
        stk[top++] = i;
        used[i] = true;
    }
    int tmp = top;
    for (int i = n - 2; 0 <= i; --i) {
        if (used[i])
            continue;
        while (tmp < top && det(b[stk[top - 1]] - b[stk[top - 2]], b[i] - b[stk[top - 1]]) <= 0)
            used[stk[--top]] = false;
        stk[top++] = i;
        used[i] = true;
    }
    for (int i = 0; i < top; ++i)
        a[i] = b[stk[i]];
    return top - 1;
}
