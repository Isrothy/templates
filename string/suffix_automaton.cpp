#include <vector>
struct SuffixAutomaton {
    static const int OMEGA = 26;
    struct Node {
        int len, link, cnt;
        std::vector<int> next;
        Node(int len, int link) : len(len), link(link), cnt(0), next(OMEGA, -1) {}
    };
    std::vector<Node> nodes;
    int last;
    SuffixAutomaton() : nodes(1, Node(0, -1)), last(0) {}
    void append(int c) {
        int p = (int) nodes.size();
        nodes.emplace_back(nodes[last].len + 1, -1);
        int q = last;
        while (q != -1 && nodes[q].next[c] == -1) {
            nodes[q].next[c] = p;
            q = nodes[q].link;
        }
        if (q == -1) {
            nodes[p].link = 0;
        } else {
            int r = nodes[q].next[c];
            if (nodes[q].len + 1 == nodes[r].len) {
                nodes[p].link = r;
            } else {
                int o = (int) nodes.size();
                nodes.emplace_back(nodes[q].len + 1, nodes[r].link);
                nodes[o].next = nodes[r].next;
                while (q != -1 && nodes[q].next[c] == r) {
                    nodes[q].next[c] = o;
                    q = nodes[q].link;
                }
                nodes[r].link = nodes[p].link = o;
            }
        }
        last = p;
    }
};
