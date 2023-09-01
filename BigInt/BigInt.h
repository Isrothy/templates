#include <cassert>
#include <cstring>
#include <string>
#include <vector>
struct BigUnsigned : protected std::vector<int> {
    static constexpr int bit = 9;
    static constexpr int base = 1e9;
#define PLACE_HOLDER "%09d"
    BigUnsigned(unsigned long long x = 0) {
        while (x != 0) {
            push_back((int) (x % base));
            x /= base;
        }
    }
    explicit BigUnsigned(const char *s) {
        for (int i = (int) strlen(s); i > 0; i -= bit) {
            int x = 0;
            for (int j = std::max(i - bit, 0); j < i; ++j) { x = x * 10 + s[j] - '0'; }
            push_back(x);
        }
        while (!empty() && back() == 0) { pop_back(); }
    }
    explicit BigUnsigned(const std::string &s) : BigUnsigned(s.c_str()) {}
    bool is_zero() const { return empty(); }
    void print(FILE *fout = stdout) const {
        int n = int(size()) - 1;
        if (n == -1) {
            fprintf(fout, "0");
        } else {
            fprintf(fout, "%d", (*this)[n]);
            n--;
            while (n != -1) {
                fprintf(fout, PLACE_HOLDER, (*this)[n]);
                n--;
            }
        }
    }
#undef PLACE_HOLDER
    BigUnsigned operator+(const BigUnsigned &that) const {
        BigUnsigned res;
        res.resize(std::max(that.size(), size()));
        for (int i = 0; i < size(); ++i) { res[i] += (*this)[i]; }
        for (int i = 0; i < that.size(); ++i) { res[i] += that[i]; }
        int carry = 0;
        for (int i = 0; i < res.size(); ++i) {
            res[i] += carry;
            carry = 0;
            if (res[i] >= base) {
                res[i] -= base;
                ++carry;
            }
        }
        if (carry != 0) {
            res.push_back(carry);
            carry = 0;
        }
        return res;
    }
    BigUnsigned operator-(const BigUnsigned &that) const {
        BigUnsigned res;
        res.resize(std::max(that.size(), size()));
        for (int i = 0; i < size(); ++i) { res[i] += (*this)[i]; }
        for (int i = 0; i < that.size(); ++i) { res[i] -= that[i]; }
        int carry = 0;
        for (int i = 0; i < res.size(); ++i) {
            res[i] += carry;
            carry = 0;
            if (res[i] < 0) {
                res[i] += base;
                carry--;
            }
        }
        while (!res.empty() && res.back() == 0) { res.pop_back(); }
        return res;
    }
    int cmp(const BigUnsigned &that) const {
        if (size() != that.size()) { return (int) size() - (int) that.size(); }
        for (int i = (int) size() - 1; i >= 0; --i) {
            if ((*this)[i] != that[i]) { return (*this)[i] - that[i]; }
        }
        return 0;
    }
    bool operator<(const BigUnsigned &that) const { return cmp(that) < 0; }
    bool operator>(const BigUnsigned &that) const { return cmp(that) > 0; }
    bool operator<=(const BigUnsigned &that) const { return cmp(that) <= 0; }
    bool operator>=(const BigUnsigned &that) const { return cmp(that) >= 0; }
    bool operator==(const BigUnsigned &that) const { return cmp(that) == 0; }
    bool operator!=(const BigUnsigned &that) const { return cmp(that) != 0; }
    BigUnsigned operator*(const BigUnsigned &that) const {
        std::vector<long long> res;
        res.resize((size()) + that.size());
        for (int i = 0; i < (size()); ++i)
            for (int j = 0; j < that.size(); ++j) {
                res[i + j] += (long long) (*this)[i] * that[j];
                res[i + j + 1] += res[i + j] / base;
                res[i + j] %= base;
            }
        long long carry = 0;
        for (auto &digit: res) {
            digit += carry;
            carry = digit / base;
            digit %= base;
        }
        while (!res.empty() && res.back() == 0) { res.pop_back(); }
        BigUnsigned ans;
        ans.resize(res.size());
        for (int i = 0; i < res.size(); ++i) { ans[i] = (int) res[i]; }
        return ans;
    }
    static std::pair<BigUnsigned, BigUnsigned> divide(BigUnsigned A, BigUnsigned b) {
        BigUnsigned c;
        std::vector<long long> a;
        a.resize(A.size());
        std::copy(A.begin(), A.end(), a.begin());
        if (a.size() >= b.size()) { c.resize(a.size() - b.size() + 1); }
        for (int i = (int) a.size(); i >= b.size(); --i) {
            auto try_devide = [&a, b, &c, i](int t) {
                int carry = 0;
                c[i - b.size()] += t;
                for (int j = i - b.size(), k = 0; k < b.size(); ++j, ++k) {
                    long long val = a[j] + carry - (long long) t * b[k];
                    carry = val / base;
                    val %= base;
                    if (val < 0) {
                        --carry;
                        val += base;
                    }
                    a[j] = val;
                }
                a[i - 1] += (long long) carry * base;
            };
            while (a[i - 1] > b.back() + 1) { try_devide(a[i - 1] / (b.back() + 1)); }
            bool left = true;
            for (int j = i - 1, k = (int) b.size() - 1; ~k; --k, --j) {
                if (a[j] != b[k]) {
                    left = a[j] >= b[k];
                    break;
                }
            }
            if (left) { try_devide(1); }
            if (i > b.size()) {
                a[i - 2] += a[i - 1] * base;
                a[i - 1] = 0;
            }
        }
        auto fix = [](auto &c) {
            int carry = 0;
            for (int i = 0; i < c.size(); ++i) {
                c[i] += carry;
                carry = c[i] / base;
                c[i] %= base;
                if (c[i] < 0) {
                    c[i] += base;
                    --carry;
                }
            }
            while (!c.empty() && c.back() == 0) { c.pop_back(); }
        };
        fix(a);
        fix(c);
        A.resize(a.size());
        std::copy(a.begin(), a.end(), A.begin());
        return make_pair(c, A);
    }
    BigUnsigned operator/(const BigUnsigned &that) const { return divide(*this, that).first; }
    BigUnsigned operator%(const BigUnsigned &that) const { return divide(*this, that).second; }
    BigUnsigned operator+=(const BigUnsigned &b) {
        *this = *this + b;
        return *this;
    }
    BigUnsigned operator-=(const BigUnsigned &b) {
        *this = *this - b;
        return *this;
    }
    BigUnsigned operator*=(const BigUnsigned &b) {
        *this = *this * b;
        return *this;
    }
    BigUnsigned operator/=(const BigUnsigned &b) {
        *this = *this / b;
        return *this;
    }
    BigUnsigned operator%=(const BigUnsigned &b) {
        *this = *this % b;
        return *this;
    }
    BigUnsigned operator++() {
        *this += 1;
        return *this;
    }
    BigUnsigned operator--() {
        *this -= 1;
        return *this;
    }
    BigUnsigned operator++(int) {
        BigUnsigned ans = *this;
        *this += 1;
        return ans;
    }
    BigUnsigned operator--(int) {
        BigUnsigned ans = *this;
        *this -= 1;
        return ans;
    }
};
struct BigInt {
    bool _is_negative;
    BigUnsigned _digits;
    bool is_negative() const { return _is_negative; }
    BigUnsigned digits() const { return _digits; }
    [[nodiscard]] bool is_zero() const { return _digits.is_zero(); }
    int sign() const {
        if (is_zero()) {
            return 0;
        } else if (is_negative()) {
            return -1;
        } else {
            return 1;
        }
    }
    void print(FILE *fout = stdout) const {
        if (is_negative()) { fprintf(fout, "-"); }
        _digits.print(fout);
    }
    BigInt(long long x = 0) {
        _is_negative = false;
        if (x < 0) {
            x = -x;
            _is_negative = true;
        }
        _digits = x;
    }
    explicit BigInt(const std::string &s) : BigInt(s.c_str()) {}
    explicit BigInt(const char *s) {
        _is_negative = false;
        if (s[0] == '-') {
            s++;
            _is_negative = true;
        }
        _digits = BigUnsigned(s);
    }
    BigInt(bool _isNegative, BigUnsigned _digits) : _is_negative(_isNegative), _digits(std::move(_digits)) {}
    BigInt operator-() const {
        if (!is_zero()) { return {!is_negative(), digits()}; }
        return {};
    }
    BigInt operator+(const BigInt &that) const {
        if (is_negative() == that.is_negative()) {
            return {is_negative(), digits() + that.digits()};
        } else {
            if (is_negative()) {
                if (digits() <= that.digits()) return {false, that.digits() - digits()};
                else return {true, digits() - that.digits()};
            } else {
                if (that.digits() <= digits()) {
                    return {false, digits() - that.digits()};
                } else {
                    return {true, that.digits() - digits()};
                }
            }
        }
    }
    BigInt operator-(const BigInt &that) const { return (*this) + (-that); }
    BigInt operator*(const BigInt &that) const {
        if (is_zero() || that.is_zero()) { return {}; }
        return {static_cast<bool>(is_negative() ^ that.is_negative()), digits() * that.digits()};
    }
    BigInt operator/(const BigInt &that) const {
        assert(!that.is_zero());
        if (is_zero()) { return 0; }
        return {static_cast<bool>(is_negative() ^ that.is_negative()), digits() / that.digits()};
    }
    BigInt operator%(const BigInt &that) const {
        assert(!that.is_zero());
        if (is_zero()) { return 0; }
        return {is_negative(), digits() % that.digits()};
    }
    int cmp(const BigInt &that) const {
        if (sign() != that.sign()) { return sign() - that.sign(); }
        if (sign() > 0) {
            return digits().cmp(that.digits());
        } else if (sign() == 0) {
            return 0;
        } else {
            return -digits().cmp(that.digits());
        }
    }
    bool operator<(const BigInt &that) const { return cmp(that) < 0; }
    bool operator>(const BigInt &that) const { return cmp(that) > 0; }
    bool operator<=(const BigInt &that) const { return cmp(that) <= 0; }
    bool operator>=(const BigInt &that) const { return cmp(that) >= 0; }
    bool operator==(const BigInt &that) const { return cmp(that) == 0; }
    bool operator!=(const BigInt &that) const { return cmp(that) != 0; }
    BigInt &operator+=(const BigInt &that) {
        *this = *this + that;
        return *this;
    }
    BigInt &operator-=(const BigInt &that) {
        *this = *this - that;
        return *this;
    }
    BigInt &operator*=(const BigInt &that) {
        *this = *this * that;
        return *this;
    }
    BigInt &operator/=(const BigInt &that) {
        *this = *this / that;
        return *this;
    }
    BigInt &operator%=(const BigInt &that) {
        *this = *this % that;
        return *this;
    }
    BigInt &operator++() {
        *this += 1;
        return *this;
    }
    BigInt &operator--() {
        *this -= 1;
        return *this;
    }
    BigInt operator++(int) {
        BigInt ans = *this;
        *this += 1;
        return ans;
    }
    BigInt operator--(int) {
        BigInt ans = *this;
        *this -= 1;
        return ans;
    }
};
