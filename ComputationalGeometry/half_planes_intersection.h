#include "2D-computational_geometry.h"
#include <queue>

Line *half_planes_intersection(Line *start, Line *end, Point *points) {
    auto n = end - start;
    std::deque<Line> q;
    std::deque<Point> t;
    std::sort(start, end, [](const Line &l1, const Line &l2) {
        int d = sign((l1.second - l1.first).angle() - (l2.second - l2.first).angle());
        return d == 0 ? sign(det(l2.first - l1.first, l2.second - l1.first)) > 0 : d < 0;
    });
    q.emplace_back(*start);
    for (auto line = start + 1; line != end; ++line) {
        if (sign(det(q.back().first - q.back().second, line->first - line->second)) == 0) {
            continue;
        }
        while (!t.empty() && sign(det(line->first - t.back(), line->second - t.back())) <= 0) {
            t.pop_back();
            q.pop_back();
        }
        while (!t.empty() && sign(det(line->first - t.front(), line->second - t.front())) <= 0) {
            t.pop_front();
            q.pop_front();
        }
        Point I;
        line_line_intersection(q.back(), *line, I);
        q.emplace_back(*line);
        t.emplace_back(I);
    }
    while (!t.empty() && sign(det(q.front().first - t.back(), q.front().second - t.back())) <= 0) {
        t.pop_back();
        q.pop_back();
    }
    Point I;
    line_line_intersection(q.front(), q.back(), I);
    t.emplace_front(I);
    std::copy(q.begin(), q.end(), start);
    std::copy(t.begin(), t.end(), points);
    return start + t.size();
}
