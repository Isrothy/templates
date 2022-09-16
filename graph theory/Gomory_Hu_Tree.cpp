void Gomory_Hu_Tree(int *p, int l, int r, int n) {
    if (l == r) {
        return;
    }
    Dinic::clear();
    int u = p[l], v = p[r];
    add_edge(u, v, Dinic::max_flow(u, v, n));
    Dinic::BFS(u, v, n);
    int i = l + 1, j = r - 1;
    while (i <= j) {
        while (i <= j && !Dinic::vis[p[i]]) {
            ++i;
        }
        while (i <= j && Dinic::vis[p[j]]) {
            --j;
        }
        if (i <= j) {
            swap(p[i], p[j]);
        }
    }
    Gomory_Hu_Tree(p, l, i - 1, n);
    Gomory_Hu_Tree(p, j + 1, r, n);
}
