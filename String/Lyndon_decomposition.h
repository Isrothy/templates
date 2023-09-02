#include <string_view>
#include <vector>
auto lyndon_decomposition(std::string_view s) {
    auto n = s.size();
    size_t l = 0, r = 1, d = 1;
    std::vector<size_t> res;
    while (r <= n) {
        if (s[r] < s[r - d]) {
            while (l + d <= r) {
                l += d;
                res.push_back(l);
            }
            r = l;
            d = 1;
        } else if (s[r - d] < s[r]) {
            d = r - l + 1;
        }
        ++r;
    }
    return res;
}
