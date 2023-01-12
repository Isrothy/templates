#include <vector>
struct SuffixAutomaton {
    static const int SIGMA = 26;
    struct Node {
        int len, link, cnt;
        std::vector<int> next;
        Node(int len, int link) : len(len), link(link), cnt(0), next(SIGMA, -1) {}
    };
    std::vector<Node> nodes;
    int last;
    SuffixAutomaton() : nodes(1, Node(0, -1)), last(0) {}
    void append(int c) {
        int p = last;
        int q = (int) nodes.size();
        nodes.emplace_back(nodes[p].len + 1, -1);
        while (p != -1 && nodes[p].next[c] == -1) {
            nodes[p].next[c] = q;
            p = nodes[p].link;
        }
        if (p == -1) {
            nodes[q].link = 0;
        } else {
            int r = nodes[p].next[c];
            if (nodes[p].len + 1 == nodes[r].len) {
                nodes[q].link = r;
            } else {
                int o = (int) nodes.size();
                nodes.emplace_back(nodes[p].len + 1, nodes[r].link);
                nodes[o].next = nodes[r].next;
                while (p != -1 && nodes[p].next[c] == r) {
                    nodes[p].next[c] = o;
                    p = nodes[p].link;
                }
                nodes[r].link = nodes[q].link = o;
            }
        }
        last = q;
    }
};
