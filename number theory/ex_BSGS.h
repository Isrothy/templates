int ex_BSGS(long long A, long long B, long long mod) {
    A %= mod;
    B %= mod;
    if (B == 1 || mod == 1)
        return 0;
    int G = gcd(A, mod);
    if (B % G != 0)
        return -1;
    if (G == 1)
        return BSGS(A, B, mod);
    int g = inv(A / G, mod / G);
    int x = ex_BSGS(A, B / G * g, mod / G);
    return x == -1 ? -1 : x + 1;
}
