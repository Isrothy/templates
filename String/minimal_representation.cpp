#include <algorithm>
#include <cstring>
int minimal_representation(char *S) {
    int n = (int) strlen(S);
    int i = 0, j = 1, k = 0;
    while (i < n && j < n && k < n) {
        if (S[(i + k) % n] == S[(j + k) % n]) {
            ++k;
        } else {
            if (S[(i + k) % n] < S[(j + k) % n]) {
                j += k + 1;
            } else {
                i += k + 1;
            }
            if (i == j) { ++i; }
            k = 0;
        }
    }
    return std::min(i, j);
}
