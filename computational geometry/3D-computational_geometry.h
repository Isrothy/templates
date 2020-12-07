struct Point {
    double x, y, z;

    Point operator+(const Point &_) const {
        return (Point) {x + _.x, y + _.y, z + _.z};
    }

    Point operator-(const Point &_) const {
        return (Point) {x - _.x, y - _.y, z - _.z};
    }

    Point operator*(const double &p) const {
        return (Point) {x * p, y * p, z * p};
    }

    Point operator/(const double &p) const {
        return (Point) {x / p, y / p, z / p};
    }

    bool operator==(const Point &_) const {
        return dcmp(x, _.x) == 0 && dcmp(y, _.y) == 0 && dcmp(z, _.z) == 0;
    }

    double operator/(const Point &_) const {
        if (unit() == _.unit())
            return len() / _.len();
        return -len() / _.len();
    }

    double len2() const {
        return x * x + y * y + z * z;
    }

    double len() const {
        return sqrtl(len2());
    }

    Point unit() const {
        return *this / len();
    }

    void read() {
        scanf("%lf%lf%lf", &x, &y, &z);
    }

    void write() const {
        printf("%lf %lf %lf\n", x, y, z);
    }

};

typedef Point Vector;

Point operator*(const double &p, const Point &A) {
    return (Point) {A.x * p, A.y * p, A.z * p};
}

double dot(const Vector &A, const Vector &B) {
    return A.x * B.x + A.y * B.y + A.z * B.z;
}

Vector det(const Vector &A, const Vector &B) {
    return (Vector) {A.y * B.z - A.z * B.y, A.z * B.x - A.x * B.z, A.x * B.y - A.y * B.x};
}

double point_line_distance(const Point &P, const Point &A, const Point &B) {
    Vector v1 = B - A, v2 = P - A;
    return det(v1, v2).len() / v1.len();
}

int point_on_segment(const Point &P, const Point &A, const Point &B) {
    return (int) (dcmp(det(A - P, B - P).len2()) == 0 && dcmp(dot(A - P, B - P)) <= 0);
}

double point_segment_distance(const Point &P, const Point &A, const Point &B) {
    if (A == B)
        return (P - A).len();
    Vector v1 = B - A, v2 = P - A, v3 = P - B;
    if (dcmp(dot(v1, v2)) < 0)
        return v2.len();
    if (dcmp(dot(v1, v3)) > 0)
        return v3.len();
    return det(v1, v2).len() / v1.len();
}

double point_plane_diatance(const Point &P, const Point &P0, const Vector &n) {
    return fabs(dot(P - P0, n));
}

Point point_plane_projection(const Point &P, const Point &P0, const Vector &n) {
    return P - n * dot(P - P0, n);
}

int coplaner(const Point &A, const Point &B, const Point &C, const Point &D) {
    return (int) (dcmp(det(det(C - A, D - A), det(C - B, D - B)).len2()) == 0);
}

int line_line_intersection(const Point &A, const Point &B, const Point &C, const Point &D, Point &O) {
    if (!coplaner(A, B, C, D) || dcmp(det(B - A, D - C).len2()) == 0)
        return 0;
    Vector s1 = det(D - A, C - A);
    Vector s2 = det(C - B, D - B);
    O = A + (B - A) * (s1 / (s1 + s2));
    return 1;
}

int line_segment_line_segment_intersection(const Point &A, const Point &B, const Point &C, const Point &D, Point &O) {
    if (!line_line_intersection(A, B, C, D, O))
        return 0;
    return (int) (point_on_segment(O, A, B) == 1 && point_on_segment(O, C, D) == 1);
}

double line_segment_line_segment_distance(const Point &P1, const Point &P2, const Point &Q1, const Point &Q2) {
    double A1 = (P1 - P2).len2(), B1 = -dot(P2 - P1, Q2 - Q1), C1 = dot(P1 - P2, P1 - Q1);
    double A2 = -dot(P2 - P1, Q2 - Q1), B2 = (Q1 - Q2).len2(), C2 = dot(P1 - Q1, Q2 - Q1);
    double x1 = B2 * C1 - B1 * C2, y1 = A1 * B2 - A2 * B1;
    double x2 = A2 * C1 - A1 * C2, y2 = A2 * B1 - A1 * B2;

    double s = dcmp(x1) == 0 && dcmp(y1) == 0 ? 0 : x1 / y1;
    double t = dcmp(x2) == 0 && dcmp(y2) == 0 ? 0 : x2 / y2;

    if (0 <= dcmp(s) && dcmp(s, 1) <= 0 && 0 <= dcmp(t) && dcmp(t, 1) <= 0)
        return ((P1 + s * (P2 - P1)) - (Q1 + t * (Q2 - Q1))).len();

    double a = point_segment_distance(P1, Q1, Q2);
    double b = point_segment_distance(P2, Q1, Q2);
    double c = point_segment_distance(Q1, P1, P2);
    double d = point_segment_distance(Q2, P1, P2);
    return min(min(a, b), min(c, d));
}

int line_plane_intersection(const Point &P1, const Point &P2, const Point P0, const Vector &n, Point &P) {
    double s = dot(n, P2 - P1);
    if (dcmp(s) == 0)
        return 0;
    Vector v = P2 - P1;
    double t = dot(n, P0 - P1) / s;
    P = P1 + v * t;
    return 1;
}

int point_in_triangle(const Point &P, const Point &A, const Point &B, const Point &C) {
    double S1 = det(A - P, B - P).len() + det(B - P, C - P).len() + det(C - P, A - P).len();
    double S2 = det(A - B, C - B).len();
    return (int) (dcmp(S1, S2) == 0);
}

double point_triangle_distance(const Point &P, const Point &A, const Point &B, const Point &C) {
    Point P0 = point_plane_projection(P, A, det(B - A, C - A).unit());
    if (point_in_triangle(P0, A, B, C))
        return (P - P0).len();
    return min(min(point_segment_distance(P, A, B), point_segment_distance(P, A, C)), point_segment_distance(P, B, C));
}

bool
line_segment_triangle_intersection(const Point &A, const Point &B, const Point &P0, const Point &P1, const Point &P2,
                                   Point &P) {
    if (line_segment_line_segment_intersection(A, B, P0, P1, P))
        return true;
    if (line_segment_line_segment_intersection(A, B, P1, P2, P))
        return true;
    if (line_segment_line_segment_intersection(A, B, P2, P0, P))
        return true;
    Vector n = det(P1 - P0, P2 - P0);
    if (dcmp(dot(n, B - A)) == 0)
        return false;
    double t = dot(n, P0 - A) / dot(n, B - A);
    if (dcmp(t) < 0 || 0 < dcmp(t, 1))
        return false;
    P = A + (B - A) * t;
    return point_in_triangle(P, P0, P1, P2);
}

bool triangle_triangle_intersection(Point *T1, Point *T2) {
    Point P;
    for (int i = 0; i < 3; ++i) {
        if (line_segment_triangle_intersection(T1[i], T1[(i + 1) % 3], T2[0], T2[1], T2[2], P))
            return true;
        if (line_segment_triangle_intersection(T2[i], T2[(i + 1) % 3], T1[0], T1[1], T1[2], P))
            return true;
    }
    return false;
}

double triangle_triangle_distrance(Point *T1, Point *T2) {
    if (triangle_triangle_intersection(T1, T2))
        return 0;
    double res = inf;
    for (int i = 0; i < 3; ++i) {
        res = min(res, point_triangle_distance(T1[i], T2[0], T2[1], T2[2]));
        res = min(res, point_triangle_distance(T2[i], T1[0], T1[1], T1[2]));
        for (int j = 0; j < 3; ++j)
            res = min(res, line_segment_line_segment_distance(T1[i], T1[(i + 1) % 3], T2[j], T2[(j + 1) % 3]));
    }
    return res;
}
