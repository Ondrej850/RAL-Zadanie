#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <inttypes.h>

typedef uint64_t u64;
typedef int64_t  i64;

/* p must be < 2^32 so that a*b < 2^64, no overflow */
#define MAXK 3

typedef struct {
    u64 a[MAXK][MAXK];
    int k;
} Mat;

u64 mulmod(u64 a, u64 b, u64 mod) {
    /* safe only when a,b < 2^32 and mod < 2^32 */
    return (a * b) % mod;
}

u64 addmod(u64 a, u64 b, u64 mod) {
    return (a + b) % mod;
}

Mat mat_mul(Mat A, Mat B, u64 mod) {
    Mat C;
    C.k = A.k;
    memset(C.a, 0, sizeof(C.a));
    for (int i = 0; i < A.k; i++)
        for (int j = 0; j < A.k; j++)
            for (int l = 0; l < A.k; l++)
                C.a[i][j] = addmod(C.a[i][j], mulmod(A.a[i][l], B.a[l][j], mod), mod);
    return C;
}

Mat mat_pow(Mat M, u64 exp, u64 mod) {
    Mat result;
    result.k = M.k;
    memset(result.a, 0, sizeof(result.a));
    for (int i = 0; i < M.k; i++)
        result.a[i][i] = 1;

    while (exp > 0) {
        if (exp & 1) result = mat_mul(result, M, mod);
        M = mat_mul(M, M, mod);
        exp >>= 1;
    }
    return result;
}

u64 linear_rec(u64 *c, u64 *init, int k, u64 n, u64 mod) {
    if (n < (u64)k) return init[n] % mod;

    Mat M;
    M.k = k;
    memset(M.a, 0, sizeof(M.a));
    for (int j = 0; j < k; j++)
        M.a[0][j] = c[j];
    for (int i = 1; i < k; i++)
        M.a[i][i-1] = 1;

    Mat Mn = mat_pow(M, n - (k - 1), mod);

    u64 result = 0;
    for (int j = 0; j < k; j++)
        result = addmod(result, mulmod(Mn.a[0][j], init[k-1-j] % mod, mod), mod);

    return result;
}

int main(void) {
    int k;
    printf("=== Linearна rekurentna postupnost mod p ===\n\n");

    printf("Rad postupnosti (2 alebo 3): ");
    if (scanf("%d", &k) != 1 || (k != 2 && k != 3)) {
        fprintf(stderr, "Chyba: rad musi byt 2 alebo 3.\n");
        return 1;
    }

    u64 p;
    printf("Prvocislo p, mod (max 4294967291): ");
    if (scanf("%" SCNu64, &p) != 1 || p < 2 || p > 4294967291ULL) {
        fprintf(stderr, "Chyba: p musi byt prvocislo, aspon 2 a najviac 4294967291.\n");
        return 1;
    }

    u64 n;
    printf("Index n (ktory clen hladate): ");
    if (scanf("%" SCNu64, &n) != 1) {
        fprintf(stderr, "Chyba: neplatne n.\n");
        return 1;
    }

    printf("\nZadajte %d koeficienty charakteristickeho polynomu\n", k + 1);
    if (k == 2)
        printf("(napr. pre x^2 - x - 1 zadajte: 1, potom -1, potom -1)\n");
    else
        printf("(napr. pre x^3 - x^2 - x - 1 zadajte: 1, potom -1, potom -1, potom -1)\n");

    i64 leading;
    printf("  Veduci koeficient (musi byt 1): ");
    if (scanf("%" SCNd64, &leading) != 1 || leading != 1) {
        fprintf(stderr, "Chyba: veduci koeficient musi byt 1.\n");
        return 1;
    }

    u64 c[MAXK];
    for (int i = 0; i < k; i++) {
        i64 coeff;
        printf("  Koeficient pri x^%d: ", k - 1 - i);
        if (scanf("%" SCNd64, &coeff) != 1) {
            fprintf(stderr, "Chyba: neplatny koeficient.\n");
            return 1;
        }
        /* rekurencny koeficient = -coeff mod p */
        i64 rec = (-coeff) % (i64)p;
        if (rec < 0) rec += (i64)p;
        c[i] = (u64)rec;
    }

    printf("\nZadajte %d pociatocnych clenov postupnosti:\n", k);
    u64 init[MAXK];
    for (int i = 0; i < k; i++) {
        printf("  a(%d): ", i);
        if (scanf("%" SCNu64, &init[i]) != 1) {
            fprintf(stderr, "Chyba: neplatny pociatocny clen.\n");
            return 1;
        }
    }

    u64 result = linear_rec(c, init, k, n, p);
    printf("\na(%" PRIu64 ") mod %" PRIu64 " = %" PRIu64 "\n", n, p, result);

    return 0;
}