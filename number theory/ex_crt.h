void ex_crt(long long a1, long long b1, long long a2, long long b2, long long &a, long long &b) {
    long long x, y;
    long long d = ex_gcd(a2, a1, x, y);
    long long _a = a1 / d;
    long long _b = multiply(x, (b1 - b2) / d, _a);
    a = a2 * _a;
    b = (multiply(a2, _b, a) + b2) % a;
}
