#include <random>
namespace {
    std::mt19937_64 mt_rand(std::random_device{}());
}// namespace
class SuffixBalancedTree {
    struct Treap {
        uint64_t pri;
        double tag;
        size_t size, ch[2];
        Treap(uint64_t pri, double tag, size_t size, size_t l, size_t r) : pri(pri), tag(tag), size(size), ch{l, r} {}
    };
    std::vector<Treap> nodes;
    std::string s;
    size_t root{};
    bool suffix_comp(size_t i, size_t j) {
        if (auto diff = s[i] - s[j]) { return diff < 0; }
        return nodes[i - 1].tag < nodes[j - 1].tag;
    }
    auto &child(size_t p, int f) { return nodes[p].ch[f]; }
    void push_up(size_t p) { nodes[p].size = nodes[child(p, 0)].size + nodes[child(p, 1)].size + 1; }
    static bool reverse_comp(const std::string_view &s, const std::string_view &t) {
        auto l1 = s.size();
        auto l2 = t.size();
        for (size_t i = 0; i < l1 && i < l2; ++i) {
            if (auto diff = s[l1 - i - 1] - t[l2 - i - 1]) { return diff < 0; }
        }
        return l1 < l2;
    }
    void retag(size_t p, double l, double r) {
        if (!p) { return; }
        auto mid = (l + r) * 0.5;
        nodes[p].tag = mid;
        retag(child(p, 0), l, mid);
        retag(child(p, 1), mid, r);
    }
    auto merge(size_t p, size_t q) {
        if (!p || !q) { return p | q; }
        if (nodes[p].pri < nodes[q].pri) {
            child(p, 1) = merge(child(p, 1), q);
            push_up(p);
            return p;
        } else {
            child(q, 0) = merge(p, child(q, 0));
            push_up(q);
            return q;
        }
    }
    auto split(size_t p, size_t x) -> std::pair<size_t, size_t> {
        if (!p) { return {0, 0}; }
        if (suffix_comp(p, x)) {
            auto [l, r] = split(child(p, 1), x);
            child(p, 1) = l;
            push_up(p);
            return {p, r};
        } else {
            auto [l, r] = split(child(p, 0), x);
            child(p, 0) = r;
            push_up(p);
            return {l, p};
        }
    }
    void insert() {
        auto *ptr = &root;
        auto pri = mt_rand();
        auto x = nodes.size();
        double l = 0, r = 1;
        while (*ptr && nodes[*ptr].pri < pri) {
            ++nodes[*ptr].size;
            auto f = suffix_comp(*ptr, x);
            ptr = &child(*ptr, f);
            (f ? l : r) = (l + r) * 0.5;
        }
        auto [ls, rs] = split(*ptr, x);
        *ptr = x;
        nodes.emplace_back(pri, (l + r) * 0.5, 1, ls, rs);
        retag(x, l, r);
        push_up(x);
    }
    auto remove() {
        auto *ptr = &root;
        auto x = nodes.size() - 1;
        double l = 0, r = 1;
        while (*ptr != x) {
            --nodes[*ptr].size;
            auto f = suffix_comp(*ptr, x);
            ptr = &child(*ptr, f);
            (f ? l : r) = (l + r) * 0.5;
        }
        *ptr = merge(child(x, 0), child(x, 1));
        retag(*ptr, l, r);
        nodes.pop_back();
    }
  public:
    explicit SuffixBalancedTree(std::string s) : s('\0' + std::move(s)) {
        nodes.emplace_back(0, 0, 0, 0, 0);
        for (size_t i = 1; i < this->s.size(); ++i) { insert(); }
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
        insert();
    }
    void pop_back() {
        remove();
        s.pop_back();
    }
};
