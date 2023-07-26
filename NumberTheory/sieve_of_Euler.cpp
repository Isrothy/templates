void sieve_of_Euler(int n, int *prime, int &prime_cnt, bool *is_composite) {
    for (int i = 2; i <= n; ++i) {
        if (!is_composite[i]) {
            prime[++prime_cnt] = i;
        }
        for (int j = 1; j <= prime_cnt && i * prime[j] <= n; ++j) {
            is_composite[i * prime[j]] = true;
            if (i % prime[j] == 0) {
                break;
            }
        }
    }
}
