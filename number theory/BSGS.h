int BSGS(long long a, long long b,int mod) {
    if (a == 0)
        return 1;
    unordered_map<long long, int> Mp;
    long long w = 1, x = 1;
    int S = sqrt(mod) + 1;
    for (int k = 1; k <= S; ++k) {
        b = b * a % mod;
        w = w * a % mod;
        Mp[b] = k;
    }
    for (int k = 1;; ++k) {
        x = x * w % mod;
        if (Mp.count(x) != 0)
            return k * S - Mp[x];
    }
}
