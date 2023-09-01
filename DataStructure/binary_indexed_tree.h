#include <cctype>
#include <span>
template<size_t M>
struct binary_indexed_tree {
    long long b0[M], b1[M];
    int n;
    void update(size_t l, size_t r, int x) {
        for (auto i = l; i <= n; i += i & -i) {
            b0[i] -= (long long) (l - 1) * x;
            b1[i] += x;
        }
        for (auto i = r + 1; i <= n; i += i & -i) {
            b0[i] += (long long) r * x;
            b1[i] -= x;
        }
    }
    auto query(size_t i) {
        long long x = 0, y = 0;
        for (auto j = i; j != 0; j -= j & -j) {
            x += b0[j];
            y += b1[j];
        }
        return x + y * i;
    }
    auto query(size_t l, size_t r) { return query(r) - query(l - 1); }
    void build(std::span<int> a) {
        auto n = a.size();
        this->n = n;
        for (int i = 1; i <= n; ++i) { b0[i] = b1[i] = 0; }
        for (int i = 1; i <= n; ++i) {
            b0[i] += a[i];
            if (i + (i & -i) < M) { b0[i + (i & -i)] += b0[i]; }
        }
    }
};
