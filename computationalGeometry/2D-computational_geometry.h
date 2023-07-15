#include <algorithm>
#include <cmath>
#include <vector>

constexpr double EPS = 1e-10;
int sign(double x) {
    if (x < -EPS) {
        return -1;
    }
    if (EPS < x) {
        return 1;
    }
    return 0;
}
int sign(double x, double y) {
    return sign(x - y);
}
struct Point {
    double x, y;
    explicit Point(double _x = 0, double _y = 0) : x(_x), y(_y) {}
    [[nodiscard]] double len2() const {
        return x * x + y * y;
    }
    [[nodiscard]] double len() const {
        return sqrt(len2());
    }
    Point operator+(Point const &_) const {
        return Point(x + _.x, y + _.y);
    }
    Point operator-(Point const &_) const {
        return Point(x - _.x, y - _.y);
    }
    Point operator*(double p) const {
        return Point(x * p, y * p);
    }
    Point operator/(double p) const {
        return Point(x / p, y / p);
    }
    bool operator==(Point const &_) const {
        return sign(x, _.x) == 0 && sign(y, _.y) == 0;
    }
    [[nodiscard]] Point unit() const {
        return *this / len();
    }
    [[nodiscard]] constexpr double angle() const {
        return std::atan2(y, x);
    }
    [[nodiscard]] Point normal() const {
        return Point(-y, x);
    }
};
using Vector = Point;
using Line = std::pair<Point, Point>;
using Ray = Line;
using Segment = Line;
using Polygon = std::vector<Point>;

Point operator*(double p, Point const &_) {
    return Point(p * _.x, p * _.y);
}
constexpr double dot(const Vector &A, const Vector &B) {
    return A.x * B.x + A.y * B.y;
}
constexpr double det(const Vector &A, const Vector &B) {
    return A.x * B.y - A.y * B.x;
}
Point middle(const Point &A, const Point &B) {
    return 0.5 * (A + B);
}
double angle(const Vector &a, const Vector &b) {
    double tmp = a.len() * b.len();
    return sign(tmp) == 0 ? 0 : acos(dot(a, b) / tmp);
}
double point_line_distance(const Point &P, const Point &A, const Point &B) {
    Vector v1 = B - A, v2 = P - A;
    return std::fabs(det(v1, v2) / v1.len());
}
double point_line_distance(const Point &P, const Line &L) {
    return point_line_distance(P, L.first, L.second);
}
bool point_on_segment(const Point &P, const Point &A, const Point &B) {
    return sign(det(A - P, B - P)) == 0 && sign(dot(A - P, B - P)) <= 0;
}
bool point_on_segment(const Point &P, const Segment &s) {
    return point_on_segment(P, s.first, s.second);
}
double point_segment_distance(const Point &P, const Point &A, const Point &B) {
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
double point_segment_distance(const Point &P, const Segment &s) {
    return point_segment_distance(P, s.first, s.second);
}
Point projection(const Point &P, const Point &A, const Point &B) {
    Vector v = B - A;
    return A + v / v.len2() * dot(v, P - A);
}
Point projection(const Point &P, const Line &l) {
    return projection(P, l.first, l.second);
}
Point symmetry(const Point &P, const Point &A, const Point &B) {
    return 2 * projection(P, A, B) - P;
}
Point symmetry(const Point &P, const Line &l) {
    return symmetry(P, l.first, l.second);
}
int line_line_intersection(
    const Point &A, const Point &B, const Point &C, const Point &D, Point &O
) {
    if (sign(det(B - A, D - C)) == 0) {
        if (sign(det(B - A, C - A)) == 0) {
            return -1;
        }
        return 0;
    }
    double s1 = det(D - A, C - A);
    double s2 = det(C - B, D - B);
    O = A + (B - A) * (s1 / (s1 + s2));
    return 1;
}
int line_line_intersection(const Line &l1, const Line &l2, Point &O) {
    return line_line_intersection(l1.first, l1.second, l2.first, l2.second, O);
}
int segment_segment_intersection(
    const Point &A, const Point &B, const Point &C, const Point &D, Point &O
) {
    int ret = line_line_intersection(A, B, C, D, O);
    if (ret != 1) {
        return ret;
    }
    return sign(dot(A - O, B - O)) <= 0 && sign(dot(C - O, D - O)) <= 0 ? 1 : 0;
}
int segment_segment_intersection(const Segment &s1, const Segment &s2, Point &O) {
    return segment_segment_intersection(s1.first, s1.second, s2.first, s2.second, O);
}
int point_in_triangle(const Point &P, const Point &A, const Point &B, const Point &C) {
    if (point_on_segment(P, A, B) || point_on_segment(P, B, C) || point_on_segment(P, C, A)) {
        return 0;
    }
    double a = angle(A - P, B - P) + angle(B - P, C - P) + angle(C - P, A - P);
    return sign(a) == 0 ? 1 : -1;
}
int point_in_convex_polygon(const Point &P, const Polygon &poly) {
    auto n = poly.size();
    for (int i = 0; i < n; ++i) {
        if (point_on_segment(P, poly[i], poly[(i + 1) % n])) {
            return 0;
        }
        if (sign(det(poly[i] - P, poly[(i + 1) % n] - P)) < 0) {
            return 1;
        }
    }
    return -1;
}
int point_in_polygon(const Point &P, const Polygon &poly) {
    int result = 0;
    auto n = poly.size();
    for (int i = 0; i < n; ++i) {
        Point A = poly[i];
        Point B = poly[(i + 1) % n];
        if (point_on_segment(P, A, B)) {
            return 0;
        }
        if (sign(A.y - B.y) > 0) {
            std::swap(A, B);
        }
        if (sign(A.y - P.y) <= 0 && sign(P.y - B.y) < 0 && sign(det(A - P, B - P)) > 0) {
            result ^= 1;
        }
    }
    return result & 1 ? -1 : 1;
}
double triangle_area(const Point &A, const Point &B, const Point &C) {
    return fabs(det(B - A, C - A)) * 0.5;
}
int point_in_circle(Point const &P, Point const &O, double r) {
    return sign((P - O).len2(), r * r);
}
int circle_line_intersection(
    Point const &O, double r, Point const &A, Point const &B, Point &P1, Point &P2
) {
    Point H = projection(O, A, B);
    double tmp = r * r - (H - O).len2();
    int tmp1 = sign(tmp);
    if (tmp1 == -1) {
        return 0;
    } else if (tmp1 == 0) {
        P1 = H;
        return 1;
    } else {
        Vector v = (A - B) / (A - B).len() * sqrt(tmp);
        P1 = H + v;
        P2 = H - v;
        return 2;
    }
}

int circle_segment_intersection(
    const Point &O, double r, const Point &A, const Point &B, Point &P1, Point &P2
) {
    int tmp = circle_line_intersection(O, r, A, B, P1, P2);
    bool o1 = tmp >= 1 && sign(dot(P1 - A, P1 - B)) <= 0;
    bool o2 = tmp == 2 && sign(dot(P2 - A, P2 - B)) <= 0;
    if (o1 && o2) {
        return 2;
    }
    if (o1) {
        return 1;
    }
    if (o2) {
        P1 = P2;
        return 1;
    }
    return 0;
}
int circle_circle_intersection(
    const Point &O1, double r1, const Point &O2, double r2, Point &P1, Point &P2
) {
    if (O2 == O1) {
        return 0;
    }
    double tmp1 = ((O2 - O1).len2() + r1 * r1 - r2 * r2) / (2 * (O2 - O1).len());
    double tmp2 = r1 * r1 - tmp1 * tmp1;
    int tmp3 = sign(tmp2);
    if (tmp3 == -1) {
        return 0;
    } else if (tmp3 == 0) {
        P1 = O1 + (O2 - O1).unit() * tmp1;
        return 1;
    } else {
        Point H = O1 + (O2 - O1).unit() * tmp1;
        Vector v = (O2 - O1).unit().normal() * sqrt(tmp2);
        P1 = H + v;
        P2 = H - v;
        return 2;
    }
}
double circle_point_tangent(Point const &O, double r, Point const &A, Point &P1, Point &P2) {
    double tmp = (O - A).len2();
    if (sign(tmp, r * r) == -1) {
        return -1;
    }
    Point H = O + (A - O) * (r * r / tmp);
    tmp = r * r - (H - O).len2();
    int tmp1 = sign(tmp);
    if (tmp1 == -1) {
        return -1;
    } else if (tmp1 == 0) {
        P1 = H;
        return 0;
    } else {
        Vector v = (A - O).unit().normal() * sqrt(tmp);
        P1 = H + v;
        P2 = H - v;
        return (A - P1).len();
    }
}
double circumscribed_circle(const Point &A, const Point &B, const Point &C, Point &O) {
    if (line_line_intersection(
            middle(A, B),
            middle(A, B) + (A - B).normal(),
            middle(B, C),
            middle(B, C) + (B - C).normal(),
            O
        )
        <= 0) {
        return -1;
    }
    return (A - O).len();
}
double inscribed_circle(const Point &A, const Point &B, const Point &C, Point &I) {
    double a = (B - C).len(), b = (C - A).len(), c = (A - B).len();
    I = (A * a + B * b + C * c) / (a + b + c);
    return point_line_distance(I, A, B);
}
double external_co_tangent(
    const Point &O1,
    double r1,
    const Point &O2,
    double r2,
    Point &P1,
    Point &P2,
    Point &P3,
    Point &P4
) {
    if (r1 < r2) {
        return external_co_tangent(O2, r2, O1, r1, P3, P4, P1, P2);
    }
    double res = circle_point_tangent(O1, r1 - r2, O2, P1, P2);
    if (res <= 0) {
        return res;
    }
    Vector v1 = (P1 - O1).unit() * r2;
    Vector v2 = (P2 - O1).unit() * r2;
    P1 = P1 + v1;
    P2 = P2 + v2;
    P3 = O2 + v2;
    P4 = O2 + v1;
    return res;
}
double internal_co_tangent(
    const Point &O1,
    double r1,
    const Point &O2,
    double r2,
    Point &P1,
    Point &P2,
    Point &P3,
    Point &P4
) {
    if (r1 < r2) {
        return internal_co_tangent(O2, r2, O1, r1, P3, P4, P1, P2);
    }
    double res = circle_point_tangent(O1, r1 + r2, O2, P1, P2);
    if (res <= 0) {
        return res;
    }
    Vector v1 = (P1 - O1).unit() * r2;
    Vector v2 = (P2 - O1).unit() * r2;
    P1 = P1 - v1;
    P2 = P2 - v2;
    P3 = O2 - v2;
    P4 = O2 - v1;
    return res;
}
double sector_area(const Point &O, const double &r, const Point &A, const Point &B) {
    Vector v1 = A - O, v2 = B - O;
    double theta = acos(dot(v1, v2) / v1.len() / v2.len());
    return 0.5 * r * r * theta;
}

double circle_segment_area(const Point &O, const double &r, const Point &A, const Point &B) {
    bool i1 = point_in_circle(A, O, r) <= 0;
    bool i2 = point_in_circle(B, O, r) <= 0;
    if (i1 && i2) {
        return triangle_area(O, A, B);
    }
    Point P{}, Q{};
    int tmp = circle_line_intersection(O, r, A, B, P, Q);
    if (tmp < 2) {
        return sector_area(O, r, A, B);
    }
    if (i1) {
        return triangle_area(O, A, Q) + sector_area(O, r, Q, B);
    }
    if (i2) {
        return sector_area(O, r, A, P) + triangle_area(O, P, B);
    }
    if (dot(A - P, B - P) <= 0 && dot(A - Q, B - Q) <= 0) {
        return sector_area(O, r, A, P) + triangle_area(O, P, Q) + sector_area(O, r, Q, B);
    }
    return sector_area(O, r, A, B);
}
double circle_polygon_intersection_area(const Point &O, const double &r, Point const *p, int n) {
    double res = 0;
    for (int i = 0; i < n; ++i) {
        double area = circle_segment_area(O, r, p[i], p[(i + 1) % n]);
        res += (det(p[i] - O, p[(i + 1) % n] - O) > 0 ? 1 : -1) * area;
    }
    return res;
}
