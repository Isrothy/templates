#include <vector>
struct ACAutomaton {
    static const int OMEGA = 26;
    struct node {
        std::vector<int> next;
        int fail;
        int cnt;
        node() : next(OMEGA, -1), fail(-1), cnt(0) {}
    };
    std::vector<node> nodes;
    std::vector<int> que, position;
    ACAutomaton() : nodes(1) {}
    int node_count = 1;
    void insert(const char *S) {
        int p = 0;
        for (int i = 0; S[i] != '\0'; ++i) {
            int c = S[i] - 'a';
            if (nodes[p].next[c] == -1) {
                nodes[p].next[c] = node_count++;
                nodes.emplace_back();
            }
            p = nodes[p].next[c];
        }
        position.push_back(p);
    }
    void build() {
        que.clear();
        nodes[0].fail = 0;
        for (int i = 0; i < OMEGA; ++i) {
            if (nodes[0].next[i] == -1) {
                nodes[0].next[i] = 0;
            } else {
                nodes[nodes[0].next[i]].fail = 0;
                que.push_back(nodes[0].next[i]);
            }
        }
        for (int i = 0; i < (int) que.size(); ++i) {
            int p = que[i];
            auto &u = nodes[p];
            for (int j = 0; j < OMEGA; ++j) {
                if (u.next[j] == -1) {
                    u.next[j] = nodes[u.fail].next[j];
                } else {
                    nodes[u.next[j]].fail = nodes[u.fail].next[j];
                    que.push_back(u.next[j]);
                }
            }
        }
    }
};
