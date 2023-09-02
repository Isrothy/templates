#include <random>
namespace {
    std::mt19937_64 mt_rand(std::random_device{}());
}// namespace
class SuffixBalancedTree {
    struct treap {
        uint64_t pri{};
        double tag{};
        size_t size{}, ch[2]{};
        treap(double tag, size_t size) : pri(mt_rand()), tag(tag), size(size) {}
    };
    std::vector<treap> nodes;
    std::string s;
    size_t root{};
    using Interval = std::pair<double, double>;
    static auto middle(const Interval &inter) { return (inter.first + inter.second) * 0.5; }
    static auto half_inter(const Interval &interval, int f) {
        auto mid = middle(interval);
        return f ? Interval(mid, interval.second) : Interval(interval.first, mid);
    }
    bool suffix_comp(size_t i, size_t j) {
        if (auto diff = s[i] - s[j]) { return diff < 0; }
        return nodes[i - 1].tag < nodes[j - 1].tag;
    }
    auto &child(size_t p, int f) { return nodes[p].ch[f]; }
    void push_up(size_t p) { nodes[p].size = nodes[child(p, 0)].size + nodes[child(p, 1)].size + 1; }
    auto rotate(size_t p, bool f) {
        auto q = child(p, f);
        child(p, f) = child(q, !f);
        child(q, !f) = p;
        push_up(p);
        return q;
    }
    static bool reverse_comp(const std::string_view &s, const std::string_view &t) {
        auto l1 = s.size();
        auto l2 = t.size();
        for (size_t i = 0; i < l1 && i < l2; ++i) {
            if (auto diff = s[l1 - i - 1] - t[l2 - i - 1]) { return diff < 0; }
        }
        return l1 < l2;
    }
    void retag(size_t p, const Interval &inter) {
        if (!p) { return; }
        nodes[p].tag = middle(inter);
        retag(child(p, 0), half_inter(inter, 0));
        retag(child(p, 1), half_inter(inter, 1));
    }
    auto insert(size_t p, size_t i, const Interval &inter) -> size_t {
        if (!p) {
            nodes.emplace_back(middle(inter), 1);
            return nodes.size() - 1;
        }
        bool f = suffix_comp(p, i);
        child(p, f) = insert(child(p, f), i, half_inter(inter, f));
        if (nodes[child(p, f)].pri < nodes[p].pri) {
            p = rotate(p, f);
            retag(p, inter);
        }
        push_up(p);
        return p;
    }
    auto remove(size_t p, size_t i, const Interval &inter) {
        if (p == i) {
            if (!child(p, 1) || !child(p, 0)) {
                auto q = child(p, 0) | child(p, 1);
                retag(q, inter);
                nodes.pop_back();
                return q;
            }
            auto f = nodes[child(p, 1)].pri < nodes[child(p, 0)].pri;
            p = rotate(p, f);
            child(p, !f) = remove(child(p, !f), i, half_inter(inter, !f));
            retag(child(p, f), half_inter(inter, f));
            nodes[p].tag = middle(inter);
        } else {
            auto f = suffix_comp(p, i);
            child(p, f) = remove(child(p, f), i, half_inter(inter, f));
        }
        push_up(p);
        return p;
    }

  public:
    explicit SuffixBalancedTree(std::string s) : s('\0' + std::move(s)) {
        nodes.emplace_back(0, 0);
        for (size_t i = 1; i < this->s.size(); ++i) { root = insert(root, i, {0, 1}); }
    }
    SuffixBalancedTree() : SuffixBalancedTree("") {}
    auto index(const std::string_view &q) {
        size_t res = 1;
        for (auto p = root; p;) {
            auto f = reverse_comp({s.data() + 1, p}, q);
            if (f) { res += nodes[child(p, 0)].size + 1; }
            p = child(p, f);
        }
        return res;
    }
    void push_back(char ch) {
        s.push_back(ch);
        root = insert(root, s.size() - 1, {0, 1});
    }
    void pop_back() {
        root = remove(root, s.size() - 1, {0, 1});
        s.pop_back();
    }
};
