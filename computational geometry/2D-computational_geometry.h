const double pi = acos(-1);

int dcmp(double x) {
    if (x < -EPS)
        return -1;
    if (EPS < x)
        return 1;
    return 0;
}

int dcmp(double x, double y) {
    if (x - y < -EPS)
        return -1;
    if (EPS < x - y)
        return 1;
    return 0;
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
        return (Point) {x + _.x, y + _.y};
    }

    Point operator-(Point const &_) const {
        return (Point) {x - _.x, y - _.y};
    }

    Point operator*(double p) const {
        return (Point) {x * p, y * p};
    }

    Point operator/(double p) const {
        return (Point) {x / p, y / p};
    }

    bool operator==(Point const &_) const {
        return dcmp(x, _.x) == 0 && dcmp(y, _.y) == 0;
    }

    Point unit() const {
        return *this / len();
    }

    double angle() const {
        return atan2(y, x);
    }

    Point normal() const {
        return (Point) {-y, x};
    }

    void read() {
        scanf("%lf%lf", &x, &y);
    }

    void write() const {
        printf("%lf %lf\n", x, y);
    }
};

typedef Point Vector;

Point operator*(double p, Point const &_) {
    return (Point) {p * _.x, p * _.y};
}

double dot(Point const &A, Point const &B) {
    return A.x * B.x + A.y * B.y;
}

double det(Point const &A, Point const &B) {
    return A.x * B.y - A.y * B.x;
}

Point middle(Point const &A, Point const &B) {
    return 0.5 * (A + B);
}

double point_line_distance(Point const &P,  Point const &A, Point const &B) {
    Vector v1 = B - A, v2 = P - A;
    return fabs(det(v1, v2) / v1.len());
}

int point_on_line_segment(Point const &P, Point const &A, Point const &B) {
    return (int) (dcmp(det(A - P, B - P)) == 0 && dcmp(dot(A - P, B - P)) <= 0);
}

double point_line_segment_distance(Point const &P, Point const &A, Point const &B) {
    if (A == B)
        return (P - A).len();
    Vector v1 = B - A, v2 = P - A, v3 = P - B;
    if (dcmp(dot(v1, v2)) < 0)
        return v2.len();
    if (dcmp(dot(v1, v3)) > 0)
        return v3.len();
    return det(v1, v2) / v1.len();
}

Point projection(Point const &P, Point const &A, Point const &B) {
    Vector v = B - A;
    return A + v * (dot(v, P - A) / v.len2());
}

Point symmetry(Point const &P, Point const &A, Point const &B) {
    return 2 * projection(P, A, B) - P;
}

int intersection(Point const&A,  Point const &B, Point const &C, Point const &D, Point &O) {
    if (dcmp(det(B - A, D - C)) == 0)
        return 0;
    double s1 = det(D - A, C - A);
    double s2 = det(C - B, D - B);
    O = A + (B - A) * (s1 / (s1 + s2));
    return 1;
}

int line_segment_line_segment_intersection(Point const &A, Point const &B, Point const &C, Point const &D, Point &O) {
    if (!intersection(A, B, C, D, O))
        return 0;
    return (int) (point_on_line_segment(O, A, B) == 1 && point_on_line_segment(O, C, D) == 1);
}

int point_in_triangle(Point const &P, Point const &A, Point const &B, Point const &C) {
    double s1 = det(A - P, B - P) + det(B - P, C - P) + det(C - P, A - P);
    double s2 = det(A - B, C - B);
    return (int) (dcmp(s1, s2) == 0);
}

int circle_line_intersection(Point const &O, double r, Point const &A, Point const &B, Point &P1, Point &P2) {
    Point H = projection(O, A, B);
    double tmp = r * r - (H - O).len2();
    int tmp1 = dcmp(tmp);
    if (tmp1 == -1)
        return 0;
    else if (tmp1 == 0) {
        P1 = H;
        return 1;
    } else {
        Vector v = (A - B) / (A - B).len() * sqrt(tmp);
        P1 = H + v;
        P2 = H - v;
        return 2;
    }
}

int circle_line_segment_intersection(Point const &O, double r, Point const &A, Point const &B, Point &P1, Point &P2) {
    int tmp = circle_line_intersection(O, r, A, B, P1, P2);
    bool o1 = tmp >= 1 && dcmp(dot(P1 - A, P1 - B)) <= 0;
    bool o2 = tmp == 2 && dcmp(dot(P2 - A, P2 - B)) <= 0;
    if (o1 && o2)
        return 2;
    if (o1)
        return 1;
    if (o2) {
        P1 = P2;
        return 1;
    }
    return 0;
}

int circle_circle_intersection(Point const &O1, double r1, Point const &O2, double r2, Point &P1, Point &P2) {
    if (O2 == O1)
        return 0;
    double tmp1 = ((O2 - O1).len2() + r1 * r1 - r2 * r2) / (2 * (O2 - O1).len());
    double tmp2 = r1 * r1 - tmp1 * tmp1;
    int tmp3 = dcmp(tmp2);
    if (tmp3 == -1)
        return 0;
    else if (tmp3 == 0) {
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
    if (dcmp(tmp, r * r) == -1)
        return -1;
    Point H = O + (A - O) * (r * r / tmp);
    tmp = r * r - (H - O).len2();
    int tmp1 = dcmp(tmp);
    if (tmp1 == -1)
        return -1;
    else if (tmp1 == 0) {
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
    if (intersection(middle(A, B), middle(A, B) + (A - B).normal(), middle(B, C), middle(B, C) + (B - C).normal(), O) ==
        0)
        return -1;
    return (A - O).len();
}

double inscribed_circle(Point const &A, Point const &B, Point const &C, Point &I) {
    double a = (B - C).len(), b = (C - A).len(), c = (A - B).len();
    I = (A * a + B * b + C * c) / (a + b + c);
    return point_line_distance(I, A, B);
}

