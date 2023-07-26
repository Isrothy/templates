#include <algorithm>
#include <array>
#include <cmath>
#include <optional>
#include <variant>
#include <vector>

constexpr double EPS = 1e-10;
constexpr int sign(double x) {
    if (x < -EPS) {
        return -1;
    }
    if (EPS < x) {
        return 1;
    }
    return 0;
}
struct Point {
    double x = 0, y = 0;
    Point() = default;
    Point(double x, double y) : x(x), y(y){};
    [[nodiscard]] double len2() const {
        return x * x + y * y;
    }
    [[nodiscard]] constexpr double len() const {
        return std::hypot(x, y);
    }
    Point operator*(double k) const {
        return {x * k, y * k};
    }
    Point operator/(double k) const {
        return {x / k, y / k};
    }
    [[nodiscard]] constexpr Point unit() const {
        return *this / len();
    }
    [[nodiscard]] Point normal() const {
        return {-y, x};
    }
    [[nodiscard]] double constexpr angle() const {
        return std::atan2(y, x);
    }
};
using Vector = Point;
using Line = std::pair<Point, Point>;
using Ray = Line;
using Segment = Line;
using Circle = std::pair<Point, double>;
using Polygon = std::vector<Point>;
using Triangle = std::tuple<Point, Point, Point>;

Vector operator+(const Vector &a, const Vector &b) {
    return {a.x + b.x, a.y + b.y};
}
Vector operator-(const Vector &a, const Vector &b) {
    return {a.x - b.x, a.y - b.y};
}
Vector operator*(double k, const Vector &a) {
    return {a.x * k, a.y * k};
}
constexpr bool operator==(const Point &A, const Point &B) {
    return sign((A - B).len()) == 0;
}
constexpr double dot(const Vector &A, const Vector &B) {
    return A.x * B.x + A.y * B.y;
}
constexpr double det(const Vector &A, const Vector &B) {
    return A.x * B.y - A.y * B.x;
}
constexpr Point middle(const Point &A, const Point &B) {
    return 0.5 * (A + B);
}
constexpr double angle(const Vector &a, const Vector &b) {
    auto tmp = a.len() * b.len();
    return sign(sqrt(tmp)) == 0 ? 0 : acos(dot(a, b) / tmp);
}
constexpr Point projection(const Point &P, const Line &l) {
    const auto &A = l.first;
    const auto &B = l.second;
    Vector v = B - A;
    return A + dot(v, P - A) * v / v.len2();
}
constexpr Point symmetry(const Point &P, const Line &l) {
    return 2 * projection(P, l) - P;
}
constexpr auto point_line_distance(const Point &P, const Line &l) {
    const auto &[A, B] = l;
    Vector v1 = B - A, v2 = P - A;
    return std::fabs(det(v1, v2) / v1.len());
}
bool point_on_segment(const Point &P, const Segment &s) {
    const auto &[A, B] = s;
    return sign(det(A - P, B - P)) == 0 && sign(dot(A - P, B - P)) <= 0;
}
double point_segment_distance(const Point &P, const Segment &s) {
    const auto &[A, B] = s;
    if (A == B) {
        return (P - A).len();
    }
    Vector v1 = B - A, v2 = P - A, v3 = P - B;
    if (sign(dot(v1, v2)) < 0) {
        return v2.len();
    }
    if (sign(dot(v1, v3)) > 0) {
        return v3.len();
    }
    return det(v1, v2) / v1.len();
}

enum class LineLineRelation {
    parallel,
    identical,
    intersecting,
};

std::pair<LineLineRelation, std::optional<Point>>
line_line_intersection(const Line &l1, const Line &l2) {
    const auto &[A, B] = l1;
    const auto &[C, D] = l2;
    if (sign(det(B - A, D - C)) == 0) {
        if (sign(det(B - A, C - A)) == 0) {
            return {LineLineRelation::identical, std::nullopt};
            ;
        }
        return {LineLineRelation::parallel, std::nullopt};
    }
    double s1 = det(D - A, C - A);
    double s2 = det(C - B, D - B);
    return {LineLineRelation::intersecting, A + (B - A) * (s1 / (s1 + s2))};
}
enum class SegmentSegmentRelation {
    disjoint,
    intersecting,
    touching,
};
std::pair<SegmentSegmentRelation, std::optional<Point>>
segment_segment_intersection(const Segment &s1, const Segment &s2) {
    const auto &[A, B] = s1;
    const auto &[C, D] = s2;
    const auto &[relation, p] = line_line_intersection(s1, s2);
    switch (relation) {
        case LineLineRelation::parallel:
            return {SegmentSegmentRelation::disjoint, std::nullopt};
        case LineLineRelation::identical: {
            if (sign(dot(C - A, C - B)) <= 0 || sign(dot(D - A, D - B)) <= 0
                || sign(dot(A - C, A - D)) <= 0 || sign(dot(B - C, B - D)) <= 0) {
                return {SegmentSegmentRelation::touching, std::nullopt};
            }
            return {SegmentSegmentRelation::disjoint, std::nullopt};
        }
        case LineLineRelation::intersecting: {
            auto O = p.value();
            if (sign(dot(O - A, O - A)) <= 0 && sign(dot(O - C, O - D)) <= 0) {
                return {SegmentSegmentRelation::intersecting, O};
            }
            return {SegmentSegmentRelation::disjoint, std::nullopt};
        }
    }
}
enum class PointShapeRelation {
    inside,
    on,
    outside,
};
auto point_in_polygon(const Point &P, const Polygon &poly) {
    int result = 0;
    auto n = poly.size();
    for (int i = 0; i < n; ++i) {
        Point A = poly[i];
        Point B = poly[(i + 1) % n];
        if (point_on_segment(P, {A, B})) {
            return PointShapeRelation::on;
        }
        if (sign(A.y - B.y) > 0) {
            std::swap(A, B);
        }
        if (sign(A.y - P.y) <= 0 && sign(P.y - B.y) < 0 && sign(det(A - P, B - P)) > 0) {
            result ^= 1;
        }
    }
    return result & 1 ? PointShapeRelation::inside : PointShapeRelation::outside;
}
double triangle_area(const Point &A, const Point &B, const Point &C) {
    return fabs(det(B - A, C - A)) * 0.5;
}

constexpr auto point_circle_relation(const Point &P, const Circle &c) {
    const auto &O = c.first;
    const auto &r = c.second;
    double d = (P - O).len();
    return std::array<PointShapeRelation, 3>{
        PointShapeRelation::inside,
        PointShapeRelation::on,
        PointShapeRelation::outside}[sign((d + r) * (d - r)) + 1];
}

enum class CircleCircleRelation {
    identital,
    disjoint,
    externally_tangent,
    internally_tangent_1_to_2,
    internally_tangent_2_to_1,
    circle1_contains_circle2,
    circle2_contains_circle1,
    intersecting,
};
constexpr auto circie_circle_relation(const Circle &c1, const Circle &c2) {
    const auto &[O1, r1] = c1;
    const auto &[O2, r2] = c2;
    auto d = (O2 - O1).len();
    if (sign(d) == 0 && sign(r1 - r2) == 0) {
        return CircleCircleRelation::identital;
    }
    switch (sign(d - r1 - r2)) {
        case 0:
            return CircleCircleRelation::externally_tangent;
        case 1:
            return CircleCircleRelation::disjoint;
        default:
            switch (sign(d - std::fabs(r1 - r2))) {
                case 0:
                    return r1 > r2 ? CircleCircleRelation::internally_tangent_2_to_1
                                   : CircleCircleRelation::internally_tangent_1_to_2;
                case -1:
                    return r1 > r2 ? CircleCircleRelation::circle1_contains_circle2
                                   : CircleCircleRelation::circle2_contains_circle1;
                default:
                    return CircleCircleRelation::intersecting;
            }
    }
}

std::pair<CircleCircleRelation, std::variant<std::monostate, Point, std::pair<Point, Point>>>
circle_circle_intersection(const Circle &c1, const Circle &c2) {
    const auto &[O1, r1] = c1;
    const auto &[O2, r2] = c2;
    auto relation = circie_circle_relation(c1, c2);
    auto d = (O2 - O1).len();
    auto d1 = 0.5 * (d + (r1 + r2) * (r1 - r2) / d);
    auto H = O1 + d1 * (O2 - O1) / d;
    switch (relation) {
        case CircleCircleRelation::disjoint:
        case CircleCircleRelation::identital:
        case CircleCircleRelation::circle1_contains_circle2:
        case CircleCircleRelation::circle2_contains_circle1:
            return {relation, {}};
        case CircleCircleRelation::externally_tangent:
        case CircleCircleRelation::internally_tangent_1_to_2:
        case CircleCircleRelation::internally_tangent_2_to_1:
            return {relation, {H}};
        case CircleCircleRelation::intersecting:
            auto v = (O2 - O1).unit().normal() * sqrt((r1 + d1) * (r1 - d1));
            return {CircleCircleRelation::intersecting, std::make_pair(H - v, H + v)};
    }
}

enum class CircleLineRelation {
    disjoint,
    tangent,
    intersecting,
};

constexpr auto circle_line_relation(const Circle &c, const Line &l) {
    const auto &[O, r] = c;
    auto d = point_line_distance(O, l);
    return std::array<CircleLineRelation, 3>{
        CircleLineRelation::disjoint,
        CircleLineRelation::tangent,
        CircleLineRelation::intersecting}[sign(d - r) + 1];
}

std::pair<CircleLineRelation, std::variant<std::monostate, Point, std::pair<Point, Point>>>
circle_line_intersection(const Circle &c, const Line &l) {
    const auto &[O, r] = c;
    const auto &[A, B] = l;
    Point H = projection(O, l);
    auto relation = circle_line_relation(c, l);
    switch (relation) {
        case CircleLineRelation::disjoint:
            return {relation, {}};
        case CircleLineRelation::tangent:
            return {relation, {H}};
        case CircleLineRelation::intersecting:
            double d = (H - O).len();
            auto v = (A - B).unit() * sqrt((r + d) * (r - d));
            return {CircleLineRelation::intersecting, std::make_pair(H - v, H + v)};
    }
}

std::pair<PointShapeRelation, std::variant<std::monostate, Point, std::pair<Point, Point>>>
circle_point_tangent(const Circle &c, const Point &P) {
    const auto &[O, r] = c;
    auto relation = point_circle_relation(P, c);
    switch (relation) {
        case PointShapeRelation::inside:
            return {relation, {}};
        case PointShapeRelation::on:
            return {relation, {P}};
        case PointShapeRelation::outside:
            auto d = (P - O).len();
            auto H = O + (P - O) * (r * r / d);
            auto v = (P - O).unit().normal() * sqrt((r + d) * (r - d));
            return {relation, std::make_pair(H - v, H + v)};
    }
}
std::optional<Circle> circumscribed_circle(const Triangle &t) {
    const auto &[A, B, C] = t;
    auto [relation, p] = line_line_intersection(
        {middle(A, B), middle(A, B) + (A - B).normal()},
        {middle(B, C), middle(B, C) + (B - C).normal()}
    );
    if (relation == LineLineRelation::intersecting) {
        return Circle{p.value(), (p.value() - A).len()};
    }
    return std::nullopt;
}
std::optional<Circle> inscribed_circle(const Triangle &t) {
    const auto &[A, B, C] = t;
    double a = (B - C).len(), b = (C - A).len(), c = (A - B).len();
    auto I = (A * a + B * b + C * c) / (a + b + c);
    double d = point_line_distance(I, {A, B});
    return sign(d) <= 0 ? std::nullopt : std::optional<Circle>{Circle{I, d}};
}
std::variant<std::monostate, Line, std::pair<Line, Line>>
external_co_tangent(const Circle &c1, const Circle &c2) {
    const auto &[O1, r1] = c1;
    const auto &[O2, r2] = c2;
    if (r1 < r2) {
        return external_co_tangent(c2, c1);
    }
    const auto &[relation, p] = circle_point_tangent({O1, r1 - r2}, O2);
    switch (relation) {
        case PointShapeRelation::inside:
            return {};
        case PointShapeRelation::on: {
            const auto P = std::get<1>(p);
            return Line{P, P + (O1 - P).normal()};
        }
        case PointShapeRelation::outside:
            const auto &[P1, P2] = std::get<2>(p);
            auto v1 = (P1 - O1).unit() * r2;
            auto v2 = (P2 - O1).unit() * r2;
            return std::make_pair(Line{P1 + v1, P2 + v2}, Line{O1 + v1, O2 + v2});
    }
}
std::variant<std::monostate, Line, std::pair<Line, Line>>
internal_co_tangent(const Circle &c1, const Circle &c2) {
    const auto &[O1, r1] = c1;
    const auto &[O2, r2] = c2;
    if (r1 < r2) {
        return internal_co_tangent(c2, c1);
    }
    const auto &[relation, p] = circle_point_tangent({O1, r1 + r2}, O2);
    switch (relation) {
        case PointShapeRelation::inside:
            return {};
        case PointShapeRelation::on: {
            auto P = std::get<1>(p);
            return {Line{P, P + (O1 - P).normal()}};
        }
        case PointShapeRelation::outside: {
            const auto &[P1, P2] = std::get<2>(p);
            auto v1 = (P1 - O1).unit() * r2;
            auto v2 = (P2 - O1).unit() * r2;
            return {std::make_pair(Line{P1 - v1, O2 + v1}, Line{P2 - v2, O2 + v1})};
        }
    }
}
