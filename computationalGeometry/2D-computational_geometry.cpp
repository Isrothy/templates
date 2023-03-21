#include <cmath>
const double EPS = 1e-10;
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
    double len2() const {
        return x * x + y * y;
    }
    double len() const {
        return sqrt(len2());
    }
    Point operator+(Point const &_) const {
        return (Point){x + _.x, y + _.y};
    }
    Point operator-(Point const &_) const {
        return (Point){x - _.x, y - _.y};
    }
    Point operator*(double p) const {
        return (Point){x * p, y * p};
    }
    Point operator/(double p) const {
        return (Point){x / p, y / p};
    }
    bool operator==(Point const &_) const {
        return sign(x, _.x) == 0 && sign(y, _.y) == 0;
    }
    Point unit() const {
        return *this / len();
    }
    double angle() const {
        return atan2(y, x);
    }
    Point normal() const {
        return (Point){-y, x};
    }
};
typedef Point Vector;
Point operator*(double p, Point const &_) {
    return (Point){p * _.x, p * _.y};
}
double dot(Vector const &A, Vector const &B) {
    return A.x * B.x + A.y * B.y;
}
double det(Vector const &A, Vector const &B) {
    return A.x * B.y - A.y * B.x;
}
Point middle(Point const &A, Point const &B) {
    return 0.5 * (A + B);
}
double angle(Vector const &a, Vector const &b) {
    double tmp = a.len() * b.len();
    return sign(tmp) == 0 ? 0 : acos(dot(a, b) / tmp);
}
double point_line_distance(Point const &P, Point const &A, Point const &B) {
    Vector v1 = B - A, v2 = P - A;
    return fabs(det(v1, v2) / v1.len());
}
int point_on_segment(Point const &P, Point const &A, Point const &B) {
    return (int) (sign(det(A - P, B - P)) == 0 && sign(dot(A - P, B - P)) <= 0);
}
double point_segment_distance(Point const &P, Point const &A, Point const &B) {
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
Point projection(Point const &P, Point const &A, Point const &B) {
    Vector v = B - A;
    return A + v * (dot(v, P - A) / v.len2());
}
Point symmetry(Point const &P, Point const &A, Point const &B) {
    return 2 * projection(P, A, B) - P;
}
int intersection(Point const &A, Point const &B, Point const &C, Point const &D, Point &O) {
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
int segment_segment_intersection(
    Point const &A, Point const &B, Point const &C, Point const &D, Point &O
) {
    int ret = intersection(A, B, C, D, O);
    if (ret != 1) {
        return ret;
    }
    return (int) (point_on_segment(O, A, B) == 1 && point_on_segment(O, C, D) == 1);
}
int point_in_triangle(Point const &P, Point const &A, Point const &B, Point const &C) {
    double s1 = det(A - P, B - P) + det(B - P, C - P) + det(C - P, A - P);
    double s2 = det(A - B, C - B);
    return (int) (sign(s1, s2) == 0);
}
int point_in_convex_polygon(Point const &P, Point const *polygon, int n) {
    for (int i = 0; i < n; ++i) {
        if (sign(det(polygon[i] - P, polygon[(i + 1) % n] - P)) < 0) {
            return 0;
        }
    }
    return 1;
}
int point_in_polygon(Point const &P, Point const *polygon, int n) {
    int result = 0;
    for (int i = 0; i < n; ++i) {
        Point A = polygon[i];
        Point B = polygon[(i + 1) % n];
        if (point_on_segment(P, A, B)) {
            return 1;
        }
        if (sign(A.y - B.y) > 0) {
            std::swap(A, B);
        }
        if (sign(A.y - P.y) <= 0 && sign(P.y - B.y) < 0 && sign(det(A - P, B - P)) > 0) {
            result ^= 1;
        }
    }
    return result;
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
    Point const &O, double r, Point const &A, Point const &B, Point &P1, Point &P2
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
    Point const &O1, double r1, Point const &O2, double r2, Point &P1, Point &P2
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
double circumscribed_circle(Point const &A, Point const &B, Point const &C, Point &O) {
    if (intersection(
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
double inscribed_circle(Point const &A, Point const &B, Point const &C, Point &I) {
    double a = (B - C).len(), b = (C - A).len(), c = (A - B).len();
    I = (A * a + B * b + C * c) / (a + b + c);
    return point_line_distance(I, A, B);
}
double external_co_tangent(
    Point const &O1,
    double r1,
    Point const &O2,
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
    Point const &O1,
    double r1,
    Point const &O2,
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
