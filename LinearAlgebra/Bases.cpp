#include <algorithm>
#include <cctype>
#include <cstdint>
struct Bases {
    static constexpr size_t K = 64;
    uint64_t A[K];
    void insert(uint64_t x) {
        for (int k = K - 1; k >= 0; --k) {
            if (((x >> k) & 1) == 1) {
                if (A[k] == 0) {
                    A[k] = x;
                    break;
                } else {
                    x ^= A[k];
                }
            }
        }
    }
    auto maximum_xor_sum(uint64_t res = 0) {
        for (int k = K - 1; k >= 0; --k) { res = std::max(res, res ^ A[k]); }
        return res;
    }
};
