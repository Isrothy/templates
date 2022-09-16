double Simpson(double l, double r) {
    double mid = (l + r) / 2;
    return (r - l) * (F(l) + 4 * F(mid) + F(r)) / 6;
}

double ASR(double l, double r, double tmp) {
    double mid = (l + r) * 0.5;
    double sl = Simpson(l, mid), sr = Simpson(mid, r);
    if (fabs(sl + sr - tmp) < EPS) {
        return sl + sr + (sl + sr - tmp);
    }
    return ASR(l, mid, sl) + ASR(mid, r, sr);
}
