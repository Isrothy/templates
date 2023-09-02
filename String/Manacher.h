#include <span>
#include <string_view>
#include <vector>
void Manacher(std::string_view s, std::span<size_t> p) {
    auto n = s.size();
    std::vector<char> t(2 * n + 3, '#');
    for (size_t i = 0; i < n; ++i) { t[2 * i + 1] = s[i]; }
    for (size_t i = 0, l = 0, r = 0; i <= 2 * n; ++i) {
        auto k = r < i ? 0 : std::min(p[l + r - i], r - i);
        while (0 <= i - k - 1 && i + k + 1 <= 2 * n && t[i - k - 1] == t[i + k + 1]) { ++k; }
        p[i] = k;
        if (r < i + k) {
            l = i - k;
            r = i + k;
        }
    }
}
