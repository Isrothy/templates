#include <cmath>
const double EPS = 1e-10;
double Simpson(double l, double r, double (*F)(double)) {
    double mid = (l + r) / 2;
    return (r - l) * (F(l) + 4 * F(mid) + F(r)) / 6;
}
double ASR(double l, double r, double tmp, double (*F)(double)) {
    double mid = (l + r) * 0.5;
    double sl = Simpson(l, mid, F), sr = Simpson(mid, r, F);
    if (fabs(sl + sr - tmp) < EPS) {
        return sl + sr + (sl + sr - tmp);
    }
    return ASR(l, mid, sl, F) + ASR(mid, r, sr, F);
}
