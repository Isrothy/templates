std::vector<int> Lyndon_decomposition(char *S) {
    int n = strlen(S);
    int l = 0, r = 1, d = 1;
    std::vector<int> res;
    while (r <= n) {
        if (S[r] < S[r - d]) {
            while (l + d <= r) { l += d; res.push_back(l); }
            r = l;
            d = 1;
        } else if (S[r - d] < S[r]) { d = r - l + 1; }
        ++r;
    }
    return res;
}
