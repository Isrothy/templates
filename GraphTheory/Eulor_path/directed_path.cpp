pair<int, int> stk[N];
int cur[M];
bool vis[N];
bool Euler(int S, int *ans, int m) {
    memset(cur, 0, sizeof cur);
    memset(vis, 0, sizeof vis);
    int sz = 0, top = 0;
    stk[top++] = make_pair(S, 0);
    while (top) {
        pair<int, int> p = stk[top - 1];
        int u = p.first, &i = cur[u];
        while (i < (int) E[u].size() && vis[abs(E[u][i].second)]) { ++i; }
        if (i < (int) E[u].size()) {
            stk[top++] = E[u][i];
            vis[abs(E[u][i].second)] = true;
            ++i;
        } else {
            ans[++sz] = p.second;
            --top;
        }
    }
    return sz - 1 == m;
}
