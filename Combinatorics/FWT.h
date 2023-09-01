#include <bit>
#include <span>
void FWT_or(std::span<int> a, int mod) {
    auto n = std::bit_floor(a.size());
    for (int i = 1; i < 1 << n; i <<= 1)
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a.data() + j, *q = a.data() + i + j;
            for (int k = 0; k < i; ++k) { q[k] = (q[k] + p[k]) % mod; }
        }
}
void IFWT_or(std::span<int> a, int mod) {
    auto n = std::bit_floor(a.size());
    for (int i = 1; i < 1 << n; i <<= 1)
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a.data() + j, *q = a.data() + i + j;
            for (int k = 0; k < i; ++k) { q[k] = (q[k] - p[k]) % mod; }
        }
}
void FWT_and(std::span<int> a, int mod) {
    auto n = std::bit_floor(a.size());
    for (int i = 1; i < 1 << n; i <<= 1)
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a.data() + j, *q = a.data() + i + j;
            for (int k = 0; k < i; ++k) { p[k] = (p[k] + q[k]) % mod; }
        }
}
void IFWT_and(std::span<int> a, int mod) {
    auto n = std::bit_floor(a.size());
    for (int i = 1; i < 1 << n; i <<= 1)
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a.data() + j, *q = a.data() + i + j;
            for (int k = 0; k < i; ++k) { p[k] = (p[k] - q[k]) % mod; }
        }
}
void FWT_xor(std::span<int> a, int mod) {
    auto n = std::bit_floor(a.size());
    for (int i = 1; i < 1 << n; i <<= 1)
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a.data() + j, *q = a.data() + i + j;
            for (int k = 0; k < i; ++k) {
                int x = p[k], y = q[k];
                p[k] = (x + y) % mod;
                q[k] = (x - y) % mod;
            }
        }
}
void IFWT_xor(std::span<int> a, int mod) {
    auto n = std::bit_floor(a.size());
    int64_t w = (mod + 1) / 2;
    for (int i = 1; i < 1 << n; i <<= 1)
        for (int j = 0; j < 1 << n; j += i << 1) {
            int *p = a.data() + j, *q = a.data() + i + j;
            for (int k = 0; k < i; ++k) {
                auto x = p[k], y = q[k];
                p[k] = static_cast<int>((x + y) * w % mod);
                q[k] = static_cast<int>((x - y) * w % mod);
            }
        }
}
