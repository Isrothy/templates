#include <cstring>
#include <numeric>
#include <string_view>
#include <vector>
auto suffix_sort(std::string_view s) {
    size_t n = s.size() + 1, m = 128;
    std::vector<size_t> tmp(n + 5), cnt(std::max(n, m) + 5), rank(n), sa(n), height(n);
    for (int i = 0; i < n; ++i) {
        rank[i] = i + 1 == n ? 0 : static_cast<unsigned char>(s[i]);
        ++cnt[rank[i]];
    }
    std::partial_sum(cnt.data(), cnt.data() + m, cnt.data());
    for (auto i = n; i-- > 0;) { sa[--cnt[rank[i]]] = i; }
    for (size_t k = 1; k <= n; k <<= 1) {
        std::iota(tmp.data(), tmp.data() + k, n - k);
        for (size_t i = 0, p = k; i < n; ++i) {
            if (sa[i] >= k) { tmp[p++] = sa[i] - k; }
        }
        std::fill(cnt.data(), cnt.data() + m, 0);
        for (size_t i = 0; i < n; ++i) { ++cnt[rank[i]]; }
        std::partial_sum(cnt.data(), cnt.data() + m, cnt.data());
        for (auto i = n; i-- > 0;) { sa[--cnt[rank[tmp[i]]]] = tmp[i]; }
        m = 1;
        tmp[sa[0]] = 0;
        for (size_t i = 1; i < n; ++i) {
            if (rank[sa[i]] != rank[sa[i - 1]] || rank[sa[i] + k] != rank[sa[i - 1] + k]) { ++m; }
            tmp[sa[i]] = m - 1;
        }
        std::swap(rank, tmp);
        if (n == m) { break; }
    }
    --n;
    for (size_t i = 0, h = 0; i < n; ++i) {
        if (h != 0) { --h; }
        auto j = sa[rank[i] - 1];
        while (i + h < n && j + h < n && s[i + h] == s[j + h]) { ++h; }
        height[rank[i] - 1] = h;
    }
    height[n] = 0;
    return std::make_pair(sa, height);
}
