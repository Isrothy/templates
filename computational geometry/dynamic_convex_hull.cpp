struct dynamic_convex_hull {

    struct cmp1 {
        bool operator()(Point A, Point B) {
            return A.x == B.x ? A.y < B.y : A.x < B.x;
        }
    };

    struct cmp2 {
        bool operator()(Point const &A, Point const &B) {
            return A.x == B.x ? A.y > B.y : A.x > B.x;
        }
    };

    set<Point, cmp1> L;
    set<Point, cmp2> U;

    template<typename T>
    bool contain(T &S, Point const &P) {
        if (S.size() < 2) {
            return false;
        }
        typename T::iterator i = S.lower_bound(P);
        return *i == P || (i != S.end() && i != S.begin() && dcmp(det(*prev(i) - P, *i - P)) >= 0);
    }

    template<typename T>
    void insert(T &S, Point const &P) {
        if (contain(S, P)) {
            return;
        }
        S.insert(P);
        auto p = S.lower_bound(P), l_bound = S.begin(), r_bound = prev(S.end());
        if (p != l_bound) {
            for (typename T::iterator i = prev(p), j; i != l_bound; i = j) {
                j = prev(i);
                if (dcmp(det(P - *j, *i - *j)) < 0) {
                    break;
                }
                S.erase(i);
                i = j;
            }
        }
        if (p != r_bound) {
            for (typename T::iterator i = next(p), j; i != r_bound; i = j) {
                j = next(i);
                if (dcmp(det(P - *j, *i - *j)) > 0) {
                    break;
                }
                S.erase(i);
                i = j;
            }
        }
    }

    bool contain(Point const &P) {
        return contain(L, P) && contain(U, P);
    }

    void insert(Point const &P) {
        insert(L, P);
        insert(U, P);
    }
};
