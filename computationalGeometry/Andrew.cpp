#include <algorithm>
int Andrew(Point *a, int n) {
    auto *b = new Point[n];
    std::sort(a, a + n, [](Point const &A, Point const &B) {
        return A.x == B.x ? A.y < B.y : A.x < B.x;
    });
    int top = 0;
    b[top++] = a[0];
    for (int i = 1; i < n; ++i) {
        while (2 <= top && det(b[top - 1] - b[top - 2], a[i] - b[top - 1]) <= 0) {
            --top;
        }
        b[top++] = a[i];
    }
    int tmp = top;
    for (int i = n - 2; i >= 0; --i) {
        while (tmp < top && det(b[top - 1] - b[top - 2], a[i] - b[top - 1]) <= 0) {
            --top;
        }
        b[top++] = a[i];
    }
    for (int i = 0; i < top; ++i) {
        a[i] = b[i];
    }
    return top - 1;
}
