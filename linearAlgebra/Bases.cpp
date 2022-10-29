struct Bases {
    unsigned long long A[K];
    void insert(unsigned long long x) {
        for (int k = K - 1; k >= 0; --k) {
            if (((x >> k) & 1) == 1) {
                if (A[k] == 0) { A[k] = x; break;
                } else { x ^= A[k]; }
            }
        }
    }
    unsigned long long maximum_xor_sum(unsigned long long res = 0) {
        for (int k = K - 1; k >= 0; --k) {
            res = max(res, res ^ A[k]);
        }
        return res;
    }
};
