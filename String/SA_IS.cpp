#include <numeric>
#include <vector>

#define L_TYPE false
#define S_TYPE true

bool is_lms_char(const std::vector<bool> &type, int x) {
    return x > 0 && type[x] == S_TYPE && type[x - 1] == L_TYPE;
}

bool substring_equal(const std::vector<int> &s, const std::vector<bool> &type, int x, int y) {
    do {
        if (s[x] != s[y]) {
            return false;
        }
        ++x;
        ++y;
    } while (!is_lms_char(type, x) && !is_lms_char(type, y));
    return s[x] == s[y];
}

void induce_sort(
    const std::vector<int> &s,
    std::vector<int> &sa,
    const std::vector<bool> &type,
    std::vector<int> &bucket,
    std::vector<int> &lbucket,
    std::vector<int> &sbucket,
    size_t n,
    size_t sigma
) {
    lbucket[0] = 0;
    sbucket[0] = 1;
    for (int i = 1; i < sigma; ++i) {
        lbucket[i] = bucket[i - 1];
        sbucket[i] = bucket[i];
    }
    for (int i = 0; i < n; ++i) {
        if (sa[i] > 0 && type[sa[i] - 1] == L_TYPE) {
            sa[lbucket[s[sa[i] - 1]]++] = sa[i] - 1;
        }
    }
    for (int i = (int) n - 1; i >= 0; --i) {
        if (sa[i] > 0 && type[sa[i] - 1] == S_TYPE) {
            sa[--sbucket[s[sa[i] - 1]]] = sa[i] - 1;
        }
    }
}

std::vector<int> SA_IS(const std::vector<int> &s, size_t sigma) {
    size_t n = s.size();
    std::vector<bool> type(n);
    std::vector<int> bucket(sigma), lbucket(sigma), sbucket(sigma);
    for (int i = 0; i < n; i++) {
        bucket[s[i]]++;
    }
    std::partial_sum(bucket.begin(), bucket.end(), bucket.begin());
    type[n - 1] = S_TYPE;
    for (int i = (int) n - 2; i >= 0; i--) {
        type[i] = s[i] < s[i + 1] || (s[i] == s[i + 1] && type[i + 1] == S_TYPE) ? S_TYPE : L_TYPE;
    }

    std::vector<int> position;
    for (int i = 1; i < n; ++i) {
        if (type[i] == S_TYPE && type[i - 1] == L_TYPE) {
            position.push_back(i);
        }
    }
    size_t m = position.size();

    sbucket[0] = 1;
    for (int i = 1; i < sigma; ++i) {
        sbucket[i] = bucket[i];
    }

    std::vector<int> sa(n, -1);
    for (int i = 0; i < m; i++) {
        sa[--sbucket[s[position[i]]]] = position[i];
    }
    induce_sort(s, sa, type, bucket, lbucket, sbucket, n, sigma);

    std::vector<int> name(n, -1);
    int last = -1, name_cnt = 1;
    bool duplicated = false;
    for (int i = 0; i < n; i++) {
        int x = sa[i];
        if (is_lms_char(type, x)) {
            if (last != -1) {
                if (!substring_equal(s, type, x, last)) {
                    ++name_cnt;
                } else {
                    duplicated = true;
                }
            }
            name[x] = name_cnt - 1;
            last = x;
        }
    }

    std::vector<int> s1;
    s1.reserve(m);
    for (int i = 0; i < n; ++i) {
        if (name[i] >= 0) {
            s1.push_back(name[i]);
        }
    }
    std::vector<int> sa1;
    if (!duplicated) {
        sa1.resize(m);
        for (int i = 0; i < m; ++i) {
            sa1[s1[i]] = i;
        }
    } else {
        sa1 = SA_IS(s1, name_cnt);
    }

    sbucket[0] = 1;
    for (int i = 1; i < sigma; ++i) {
        sbucket[i] = bucket[i];
    }
    fill(sa.begin(), sa.end(), -1);
    for (int i = (int) m - 1; i >= 0; --i) {
        sa[--sbucket[s[position[sa1[i]]]]] = position[sa1[i]];
    }
    induce_sort(s, sa, type, bucket, lbucket, sbucket, n, sigma);
    return sa;
}

std::pair<std::vector<int>, std::vector<int>> suffix_sort(char *S) {
    size_t n = strlen(S);
    std::vector<int> s(n + 1);
    for (int i = 0; i < n; ++i) {
        s[i] = (int) (unsigned) S[i];
    }
    auto sa = SA_IS(s, 256);
    std::vector<int> rank(n + 1), height(n);
    for (int i = 0; i <= n; ++i) {
        rank[sa[i]] = i;
    }
    for (int i = 0, h = 0; i < n; ++i) {
        if (h != 0) {
            --h;
        }
        int j = sa[rank[i] - 1];
        while (i + h < n && j + h < n && S[i + h] == S[j + h]) {
            ++h;
        }
        height[rank[i] - 1] = h;
    }
    return {sa, height};
}

#undef L_TYPE
#undef S_TYPE
