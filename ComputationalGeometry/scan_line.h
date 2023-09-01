#include "2D_computational_geometry.h"
#include <algorithm>
#include <cmath>
#include <set>
#include <vector>
void build_tree(Point *o, long long *r, int n, int *par) {
    struct arc {
        int id;
        Point o;
        double r;
        bool convax;
        arc(int id, Point o, double r, bool convax) : id(id), o(o), r(r), convax(convax) {}
        double intersection(double x) const { return o.y + (convax ? -1 : 1) * sqrt(r * r - (x - o.x) * (x - o.x)); }
    };
    struct events {
        double x;
        int id;
        bool type;
        events(double x, int id, bool type) : x(x), id(id), type(type) {}
        bool operator<(const events &other) const { return x == other.x ? type < other.type : x < other.x; }
    };
    std::vector<events> events;
    for (int i = 0; i < n; i++) {
        events.emplace_back(o[i].x - r[i], i, 0);
        events.emplace_back(o[i].x + r[i], i, 1);
    }
    std::sort(events.begin(), events.end());
    double X = 0;
    auto comp = [&](const arc &a, const arc &b) {
        double y1 = a.intersection(X);
        double y2 = b.intersection(X);
        return y1 == y2 ? a.convax > b.convax : y1 < y2;
    };
    std::set<arc, decltype(comp)> s(comp);
    for (auto &e: events) {
        X = e.x;
        if (e.type == 0) {
            auto it = s.upper_bound(arc(e.id, o[e.id], r[e.id], false));
            par[e.id] = it != s.end() && !it->convax ? it->id : -1;
            s.emplace(e.id, o[e.id], r[e.id], true);
            s.emplace(e.id, o[e.id], r[e.id], false);
        } else {
            s.erase(arc(e.id, o[e.id], r[e.id], false));
            s.erase(arc(e.id, o[e.id], r[e.id], true));
        }
    }
}
