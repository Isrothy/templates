#include <deque>
#include <string_view>
template<size_t M, size_t Sigma> struct PalindromicAutomaton {
    size_t next[M][Sigma]{}, fail[M]{};
    int length[M]{};
    std::deque<char> str;
    int n, left, right;
    PalindromicAutomaton() : n(2), left(0), right(0) {
        length[1] = -1;
        fail[0] = 1;
    }
    void push_back(char c) {
        auto get_fail = [&](size_t p) {
            while (str.size() <= length[p] + 1 || str[str.size() - length[p] - 2] != str.back()) { p = fail[p]; }
            return p;
        };
        str.push_back(c);
        auto p = get_fail(right);
        auto &q = next[p][c - 'a'];
        if (!q) {
            auto r = ++n;
            length[r] = length[p] + 2;
            fail[r] = next[get_fail(fail[p])][c - 'a'];
            q = r;
        }
        p = q;
        right = p;
        if (length[p] == str.size()) { left = p; }
    }
    void push_front(char c) {
        auto get_fail = [&](size_t p) {
            while (str.size() <= length[p] + 1 || str[length[p] + 1] != str.front()) { p = fail[p]; }
            return p;
        };
        str.push_front(c);
        auto p = get_fail(left);
        auto &q = next[p][c - 'a'];
        if (!q) {
            auto r = ++n;
            length[r] = length[p] + 2;
            fail[r] = next[get_fail(fail[p])][c - 'a'];
            q = r;
        }
        p = q;
        left = p;
        if (length[p] == str.size()) { right = p; }
    }
};
