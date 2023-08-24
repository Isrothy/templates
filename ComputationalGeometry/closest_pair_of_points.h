#include "2D_computational_geometry.h"
namespace {
    auto cmp_y(const Point &A, const Point &B) {
        return A.y != B.y ? A.y < B.y : A.x < B.x;
    };
    auto closer(const Segment &s1, const Segment &s2) {
        return len2(s1) < len2(s2);
    }
    auto clostest_pair_of_points_helper(std::span<Point> a) -> Segment {
        auto n = a.size();
        if (n <= 3) {
            auto res = std::make_pair(a[0], a[1]);
            for (size_t i = 0; i < n; ++i) {
                for (size_t j = i + 1; j < n; ++j) {
                    res = std::min(res, std::make_pair(a[i], a[j]), closer);
                }
            }
            std::sort(a.begin(), a.end(), cmp_y);
            return res;
        }
        auto mid = n >> 1;
        auto xmid = a[mid].x;
        auto ls = clostest_pair_of_points_helper(a.subspan(0, mid));
        auto rs = clostest_pair_of_points_helper(a.subspan(mid, n - mid));
        auto res = std::min(ls, rs, closer);
        std::inplace_merge(a.begin(), a.begin() + static_cast<ptrdiff_t>(mid), a.end(), cmp_y);
        std::deque<Point> b;
        for (const auto &P: a) {
            if (xmid - P.x >= len(res)) {
                continue;
            }
            for (const auto &Q: b) {
                res = std::min(res, std::make_pair(P, Q), closer);
            }
            while (!b.empty() && P.y - b.front().y >= len(res)) {
                b.pop_front();
            }
            b.push_back(P);
        }
        return res;
    }
}// namespace
auto closest_pair_of_points(std::vector<Point> a) {
    assert(a.size() >= 2);
    std::sort(a.begin(), a.end(), [](const auto &P, const auto &Q) {
        return P.x != Q.x ? P.x < Q.x : P.y < Q.y;
    });
    return clostest_pair_of_points_helper(a);
}
