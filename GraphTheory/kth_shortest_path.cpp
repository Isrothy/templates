#include <numeric>
#include <queue>
#include <span>
template<typename T> struct KthSsp {
    static constexpr T INF = std::numeric_limits<T>::max();
    template<typename C> using MinHeap = std::priority_queue<C, std::vector<C>, std::greater<>>;
    struct Heap {
        T value;
        size_t to, height, ch[2]{};
        Heap(T val, size_t to, size_t heigtht) : value(val), to(to), height(heigtht) {}
    };
    struct Edge {
        T cost;
        size_t to;
        Edge(size_t to, T cost) : cost(cost), to(to) {}
    };
    std::vector<Heap> heaps;
    std::vector<std::vector<Edge>> edges, reverseEdges;
    std::vector<T> dis;
    std::vector<size_t> a, options;
    MinHeap<std::pair<T, size_t>> q;
    size_t n;
    bool first{};
    std::optional<T> ssp;
    explicit KthSsp(size_t n) : edges(n + 1), reverseEdges(n + 1), dis(n + 1), options(n + 1), n(n) {}
    size_t new_heap(T val, size_t to) {
        heaps.emplace_back(val, to, 1);
        return heaps.size() - 1;
    }
    auto &child(size_t p, size_t i) { return heaps[p].ch[i]; }
    size_t merge(size_t p, size_t q) {
        if (!p || !q) { return p | q; }
        if (heaps[p].value > heaps[q].value) { std::swap(p, q); }
        auto r = new_heap(heaps[p].value, heaps[p].to);
        child(r, 0) = child(p, 0);
        child(r, 1) = merge(child(p, 1), q);
        if (heaps[child(r, 0)].height < heaps[child(r, 1)].height) { std::swap(child(r, 0), child(r, 1)); }
        heaps[r].height = heaps[child(r, 0)].height + 1;
        return r;
    }
    void add_edge(size_t u, size_t v, double cost) {
        edges[u].emplace_back(v, cost);
        reverseEdges[v].emplace_back(u, cost);
    }
    void dijkstra(int s, std::span<std::vector<Edge>> edges) {
        std::fill(dis.begin(), dis.end(), INF);
        MinHeap<std::pair<T, size_t>> q;
        dis[s] = 0;
        q.emplace(0, s);
        while (!q.empty()) {
            auto [d, u] = q.top();
            q.pop();
            if (d != dis[u]) { continue; }
            for (auto &e: edges[u]) {
                if (auto newDis = d + e.cost; newDis < dis[e.to]) {
                    dis[e.to] = newDis;
                    q.emplace(newDis, e.to);
                }
            }
        }
    }
    void build(size_t s, size_t t) {
        dijkstra(t, reverseEdges);
        heaps = std::vector<Heap>({{0, 0, 0}});
        q = MinHeap<std::pair<T, size_t>>();
        std::fill(options.begin(), options.end(), 0);
        first = true;
        ssp = dis[s] == INF ? std::nullopt : std::optional<T>(dis[s]);
        for (auto u: a) {
            Edge *tree_edge = nullptr;
            for (auto &e: edges[u]) {
                if (dis[u] == dis[e.to] + e.cost) { tree_edge = &e; }
            }
            std::vector<size_t> tmp;
            tmp.reserve(edges[u].size());
            for (auto &e: edges[u]) {
                if (&e != tree_edge && dis[e.to] != INF) { tmp.push_back(new_heap(dis[e.to] - dis[u] + e.cost, e.to)); }
            }
            std::make_heap(tmp.begin(), tmp.end(), [this](size_t p, size_t q) { return heaps[p].value > heaps[q].value; });
            for (size_t i = tmp.size(); i-- > 0;) {
                if (auto j = i * 2 + 1; j < tmp.size()) { child(tmp[i], 0) = tmp[j]; }
                if (auto j = i * 2 + 2; j < tmp.size()) { child(tmp[i], 1) = tmp[j]; }
                heaps[tmp[i]].height = heaps[child(tmp[i], 0)].height + 1;
            }
            options[u] = tmp.empty() ? 0 : tmp[0];
            if (tree_edge) { options[u] = merge(options[u], options[tree_edge->to]); }
        }
        if (options[s]) { q.emplace(heaps[options[s]].value + dis[s], options[s]); }
    }
    auto get_next() -> std::optional<T> {
        if (first) {
            first = false;
            return ssp;
        }
        if (q.empty()) { return std::nullopt; }
        auto [val, h] = q.top();
        q.pop();
        if (auto u = child(h, 0)) { q.emplace(val - heaps[h].value + heaps[u].value, u); }
        if (auto u = child(h, 1)) { q.emplace(val - heaps[h].value + heaps[u].value, u); }
        if (auto i = options[heaps[h].to]) { q.emplace(val + heaps[i].value, i); }
        return val;
    }
};
