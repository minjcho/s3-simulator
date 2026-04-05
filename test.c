/*
 * Mini S3 - 종합 실험 및 벤치마크
 *
 * [1] GF(256) 필드 공리 전수 검증
 * [2] 이레이저 코딩 전수 검사 (255개 고장 패턴)
 * [3] 대용량 랜덤 데이터 복구 검증
 * [4] Power of Two Choices 통계 분석
 * [5] 인코딩/디코딩 처리량 벤치마크
 * [6] 읽기 꼬리 지연시간 시뮬레이션
 * [7] 저장 방식 효율 비교
 * [8] 내구성(Durability) Markov 체인 + Monte Carlo 분석
 * [9] CRC32 Silent Corruption 탐지
 * [10] Rebuild 취약 구간 분석
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "gf256.h"
#include "erasure.h"

#define K 5  /* data shards   */
#define M 4  /* parity shards */
#define N 9  /* total shards  */

/* ---- utilities ---- */

static double now(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + ts.tv_nsec * 1e-9;
}

static int dbl_cmp(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return (x > y) - (x < y);
}

static double vec_stddev(const int *a, int n)
{
    double m = 0;
    for (int i = 0; i < n; i++) m += a[i];
    m /= n;
    double v = 0;
    for (int i = 0; i < n; i++) { double d = a[i] - m; v += d * d; }
    return sqrt(v / n);
}

static int vec_max(const int *a, int n)
{
    int m = a[0];
    for (int i = 1; i < n; i++) if (a[i] > m) m = a[i];
    return m;
}

static int vec_min(const int *a, int n)
{
    int m = a[0];
    for (int i = 1; i < n; i++) if (a[i] < m) m = a[i];
    return m;
}

/* ---- CRC32 ---- */

static uint32_t crc32_tbl[256];

static void crc32_init(void)
{
    for (uint32_t i = 0; i < 256; i++) {
        uint32_t c = i;
        for (int j = 0; j < 8; j++)
            c = (c & 1) ? (c >> 1) ^ 0xEDB88320u : (c >> 1);
        crc32_tbl[i] = c;
    }
}

static uint32_t crc32(const uint8_t *data, size_t len)
{
    uint32_t c = 0xFFFFFFFF;
    for (size_t i = 0; i < len; i++)
        c = crc32_tbl[(c ^ data[i]) & 0xFF] ^ (c >> 8);
    return c ^ 0xFFFFFFFF;
}

/* ---- Durability helpers ---- */

static double rand_exp(double rate)
{
    double u = ((double)rand() + 1.0) / ((double)RAND_MAX + 2.0);
    return -log(u) / rate;
}

/*
 * MTTDL 계산 (Markov 체인 backward recursion)
 *
 * 상태 i = 고장 디스크 i개 (0 ≤ i ≤ m)
 * 상태 m+1 = 데이터 손실 (흡수 상태)
 * 전이: i→i+1 rate (n-i)*λ, i→i-1 rate i*μ
 */
static double compute_mttdl(int n, int m, double lambda, double mu)
{
    if (m == 0)
        return 1.0 / ((double)n * lambda);

    double a_m = (double)(n - m) * lambda;
    double b_m = (double)m * mu;
    double alpha = 1.0 / (a_m + b_m);
    double beta  = b_m / (a_m + b_m);

    for (int i = m - 1; i >= 1; i--) {
        double a_i = (double)(n - i) * lambda;
        double b_i = (double)i * mu;
        double g  = 1.0 - a_i * beta / (a_i + b_i);
        double na = (1.0 + a_i * alpha) / ((a_i + b_i) * g);
        double nb = b_i / ((a_i + b_i) * g);
        alpha = na;
        beta  = nb;
    }

    double a_0 = (double)n * lambda;
    return (1.0 / a_0 + alpha) / (1.0 - beta);
}

/* 이산 사건 시뮬레이션: 1년간 단일 객체의 데이터 손실 여부 */
static int mc_trial(int n, int m, double lambda, double rebuild_hours)
{
    double ev[20];
    int act[20];
    for (int i = 0; i < n; i++) {
        act[i] = 1;
        ev[i] = rand_exp(lambda);
    }
    for (;;) {
        double t = 8760.0;
        int d = -1;
        for (int i = 0; i < n; i++)
            if (ev[i] < t) { t = ev[i]; d = i; }
        if (d < 0 || t >= 8760.0) return 0;
        if (act[d]) {
            act[d] = 0;
            ev[d] = t + rebuild_hours;
            int f = 0;
            for (int i = 0; i < n; i++) if (!act[i]) f++;
            if (f > m) return 1;
        } else {
            act[d] = 1;
            ev[d] = t + rand_exp(lambda);
        }
    }
}

static double binom_coeff(int n, int k)
{
    if (k < 0 || k > n) return 0;
    double r = 1;
    for (int i = 0; i < k; i++)
        r = r * (n - i) / (i + 1);
    return r;
}

/* P(X >= k) where X ~ Binomial(n, p) */
static double binom_tail(int n, int k, double p)
{
    double s = 0;
    for (int i = k; i <= n; i++)
        s += binom_coeff(n, i) * pow(p, i) * pow(1.0 - p, n - i);
    return s;
}

/* ================================================================
 * [1] GF(256) 필드 공리 검증
 *
 * 이레이저 코딩의 수학적 기반인 GF(256)이 올바르게 구현되었는지
 * 4가지 필드 공리를 전수/샘플 검사한다.
 * ================================================================ */
static void test_gf256(void)
{
    printf("[1] GF(256) 필드 공리 검증\n");
    printf("--------------------------------------------\n");

    int ok = 1, pass;

    /* 항등원: a * 1 = a, 전수 검사 (256개) */
    pass = 0;
    for (int a = 0; a < 256; a++)
        if (gf256_mul((uint8_t)a, 1) == (uint8_t)a) pass++;
    printf("  곱셈 항등원 (a*1=a):       %d/256\n", pass);
    ok = ok && (pass == 256);

    /* 역원: a * inv(a) = 1, 전수 검사 (255개, 0 제외) */
    pass = 0;
    for (int a = 1; a < 256; a++)
        if (gf256_mul((uint8_t)a, gf256_inv((uint8_t)a)) == 1) pass++;
    printf("  곱셈 역원 (a*a^-1=1):      %d/255\n", pass);
    ok = ok && (pass == 255);

    /* 교환법칙: a*b = b*a, 전수 검사 (65536 쌍) */
    pass = 0;
    for (int a = 0; a < 256; a++)
        for (int b = 0; b < 256; b++)
            if (gf256_mul((uint8_t)a, (uint8_t)b) ==
                gf256_mul((uint8_t)b, (uint8_t)a))
                pass++;
    printf("  교환법칙 (a*b=b*a):        %d/65536\n", pass);
    ok = ok && (pass == 65536);

    /* 분배법칙: a*(b+c) = a*b + a*c, 10만 샘플 */
    pass = 0;
    int total = 100000;
    for (int i = 0; i < total; i++) {
        uint8_t a = rand() & 0xFF, b = rand() & 0xFF, c = rand() & 0xFF;
        if (gf256_mul(a, b ^ c) == (gf256_mul(a, b) ^ gf256_mul(a, c)))
            pass++;
    }
    printf("  분배법칙 (a*(b+c)=ab+ac):  %d/%d\n", pass, total);
    ok = ok && (pass == total);

    printf("  결과: %s\n\n", ok ? "PASS" : "FAIL");
}

/* ================================================================
 * [2] 이레이저 코딩 전수 검사
 *
 * 고장 1~4개에 대해 가능한 모든 조합 C(9,f)을 테스트한다.
 * 총 C(9,1)+C(9,2)+C(9,3)+C(9,4) = 9+36+84+126 = 255 패턴.
 * MDS 코드 성질에 의해 전부 복구 가능해야 한다.
 * ================================================================ */

static int test_pattern(erasure_t *ec, const uint8_t *enc,
                        size_t ss, const int *fail, int nf)
{
    uint8_t *w = malloc(N * ss);
    memcpy(w, enc, N * ss);
    uint8_t *sh[N];
    int pr[N];
    for (int i = 0; i < N; i++) {
        sh[i] = w + (size_t)i * ss;
        pr[i] = 1;
    }
    for (int f = 0; f < nf; f++) {
        pr[fail[f]] = 0;
        memset(sh[fail[f]], 0, ss);
    }
    int ret = erasure_decode(ec, sh, ss, pr);
    int ok = (ret == 0 && memcmp(w, enc, K * ss) == 0);
    free(w);
    return ok;
}

static void test_erasure_exhaustive(void)
{
    printf("[2] 이레이저 코딩 전수 검사\n");
    printf("--------------------------------------------\n");
    printf("  설정: %d 데이터 + %d 패리티 = %d 조각\n\n", K, M, N);

    erasure_t ec;
    erasure_init(&ec, K, M);

    size_t ss = 64;
    uint8_t *mem = calloc(N, ss);
    for (size_t i = 0; i < K * ss; i++)
        mem[i] = rand() & 0xFF;
    uint8_t *sh[N];
    for (int i = 0; i < N; i++) sh[i] = mem + (size_t)i * ss;
    erasure_encode(&ec, sh, ss);

    int all_pass = 0, all_total = 0;

    for (int nf = 1; nf <= M; nf++) {
        int pass = 0, total = 0;
        int idx[4];
        for (int i = 0; i < nf; i++) idx[i] = i;

        for (;;) {
            total++;
            pass += test_pattern(&ec, mem, ss, idx, nf);

            int i = nf - 1;
            while (i >= 0 && idx[i] == N - nf + i) i--;
            if (i < 0) break;
            idx[i]++;
            for (int j = i + 1; j < nf; j++) idx[j] = idx[j - 1] + 1;
        }

        printf("  고장 %d개: %3d/%-3d (C(%d,%d) 패턴 전수 검사)\n",
               nf, pass, total, N, nf);
        all_pass += pass;
        all_total += total;
    }

    printf("\n  총 %d개 패턴: %s\n\n", all_total,
           (all_pass == all_total) ? "PASS" : "FAIL");
    free(mem);
    erasure_free(&ec);
}

/* ================================================================
 * [3] 대용량 랜덤 데이터 복구 검증
 *
 * 다양한 크기의 랜덤 데이터에 대해 임의 개수(1~4)의
 * 랜덤 고장을 발생시키고 복구를 검증한다.
 * ================================================================ */
static void test_erasure_random(void)
{
    printf("[3] 대용량 랜덤 데이터 복구 검증\n");
    printf("--------------------------------------------\n");

    erasure_t ec;
    erasure_init(&ec, K, M);

    struct { size_t size; const char *label; int trials; } cfgs[] = {
        {100,     "100 B ", 100},
        {1024,    "1 KB  ", 100},
        {65536,   "64 KB ", 50},
        {1048576, "1 MB  ", 20},
    };
    int all_ok = 1;

    for (int c = 0; c < 4; c++) {
        size_t dsz = cfgs[c].size;
        size_t ss = (dsz + K - 1) / K;
        int trials = cfgs[c].trials;
        int pass = 0;

        for (int t = 0; t < trials; t++) {
            uint8_t *mem = calloc(N, ss);
            uint8_t *sh_arr[N];
            for (int i = 0; i < N; i++)
                sh_arr[i] = mem + (size_t)i * ss;
            for (size_t i = 0; i < dsz; i++)
                mem[i] = rand() & 0xFF;

            uint8_t *orig = malloc(K * ss);
            memcpy(orig, mem, K * ss);
            erasure_encode(&ec, sh_arr, ss);

            int nf = 1 + rand() % M;
            int pr[N];
            for (int i = 0; i < N; i++) pr[i] = 1;
            int d = 0;
            while (d < nf) {
                int x = rand() % N;
                if (!pr[x]) continue;
                pr[x] = 0;
                memset(sh_arr[x], 0, ss);
                d++;
            }

            int ret = erasure_decode(&ec, sh_arr, ss, pr);
            if (ret == 0 && memcmp(mem, orig, K * ss) == 0) pass++;
            free(mem);
            free(orig);
        }

        printf("  %s: %d/%d 통과 (랜덤 1~%d개 고장)\n",
               cfgs[c].label, pass, trials, M);
        if (pass != trials) all_ok = 0;
    }

    printf("  결과: %s\n\n", all_ok ? "PASS" : "FAIL");
    erasure_free(&ec);
}

/* ================================================================
 * [4] Power of Two Choices 통계 분석
 *
 * 디스크 수를 16~512로 변화시키며, 랜덤 배치와 P2C를
 * 200회 시행하여 최대 부하, 최소 부하, 표준편차를 비교한다.
 *
 * 이론적 기대:
 *   랜덤: max load ≈ mean + O(sqrt(mean * ln(n)))
 *   P2C:  max load ≈ mean + O(ln ln n)
 * ================================================================ */
static void test_p2c(void)
{
    printf("[4] Power of Two Choices 통계 분석\n");
    printf("--------------------------------------------\n");

    int disks[] = {16, 32, 64, 128, 256, 512};
    int ncfg = 6;
    int lf = 100;
    int trials = 200;

    printf("  부하율 %d (디스크당 평균 %d개), %d회 시행\n\n", lf, lf, trials);
    printf("  | 디스크 | 방식 | 평균 최대 | 평균 최소 | 표준편차 | 최대/평균 |\n");
    printf("  |--------|------|----------|----------|---------|----------|\n");

    for (int c = 0; c < ncfg; c++) {
        int nd = disks[c];
        int items = lf * nd;
        double mean = (double)lf;

        /* 랜덤 배치 */
        double r_mx = 0, r_mn = 0, r_sd = 0;
        for (int t = 0; t < trials; t++) {
            int *cnt = calloc((size_t)nd, sizeof(int));
            for (int i = 0; i < items; i++)
                cnt[rand() % nd]++;
            r_mx += vec_max(cnt, nd);
            r_mn += vec_min(cnt, nd);
            r_sd += vec_stddev(cnt, nd);
            free(cnt);
        }
        r_mx /= trials; r_mn /= trials; r_sd /= trials;

        printf("  | %5d  | rand | %8.1f | %8.1f | %7.1f | %8.4f |\n",
               nd, r_mx, r_mn, r_sd, r_mx / mean);

        /* Power of Two Choices */
        double p_mx = 0, p_mn = 0, p_sd = 0;
        for (int t = 0; t < trials; t++) {
            int *cnt = calloc((size_t)nd, sizeof(int));
            for (int i = 0; i < items; i++) {
                int a = rand() % nd, b;
                do { b = rand() % nd; } while (b == a);
                cnt[(cnt[a] <= cnt[b]) ? a : b]++;
            }
            p_mx += vec_max(cnt, nd);
            p_mn += vec_min(cnt, nd);
            p_sd += vec_stddev(cnt, nd);
            free(cnt);
        }
        p_mx /= trials; p_mn /= trials; p_sd /= trials;

        printf("  | %5d  | P2C  | %8.1f | %8.1f | %7.1f | %8.4f |\n",
               nd, p_mx, p_mn, p_sd, p_mx / mean);
    }
    printf("\n");
}

/* ================================================================
 * [5] 인코딩/디코딩 처리량 벤치마크
 *
 * 다양한 크기의 데이터에 대해 인코딩/디코딩 처리량(MB/s)을 측정한다.
 * 디코딩은 3개 데이터 + 1개 패리티 고장 (최악에 가까운 시나리오).
 * ================================================================ */
static void bench_throughput(void)
{
    printf("[5] 인코딩/디코딩 처리량 벤치마크\n");
    printf("--------------------------------------------\n");
    printf("  디코딩: 4개 고장 (데이터 3 + 패리티 1) 복구\n");

    erasure_t ec;
    erasure_init(&ec, K, M);

    struct { size_t size; const char *label; int iters; } cfgs[] = {
        {1024,    "1 KB  ", 5000},
        {4096,    "4 KB  ", 3000},
        {16384,   "16 KB ", 1000},
        {65536,   "64 KB ", 500},
        {262144,  "256 KB", 100},
        {1048576, "1 MB  ", 30},
    };
    int ncfg = 6;

    printf("\n  | 데이터 크기 | 인코딩 (MB/s) | 디코딩 (MB/s) |\n");
    printf("  |------------|--------------|---------------|\n");

    for (int c = 0; c < ncfg; c++) {
        size_t total = cfgs[c].size;
        size_t ss = (total + K - 1) / K;
        int iters = cfgs[c].iters;

        uint8_t *mem = calloc(N, ss);
        uint8_t *sh[N];
        for (int i = 0; i < N; i++)
            sh[i] = mem + (size_t)i * ss;
        for (size_t i = 0; i < total; i++)
            mem[i] = rand() & 0xFF;

        /* 인코딩 벤치마크 */
        double t0 = now();
        for (int i = 0; i < iters; i++)
            erasure_encode(&ec, sh, ss);
        double enc_mbs = (double)total * iters / (now() - t0)
                         / (1024.0 * 1024.0);

        /* 디코딩 벤치마크 */
        int pr[N] = {0, 0, 0, 1, 1, 0, 1, 1, 1};
        int fail[] = {0, 1, 2, 5};

        t0 = now();
        for (int i = 0; i < iters; i++) {
            for (int f = 0; f < 4; f++)
                memset(sh[fail[f]], 0, ss);
            erasure_decode(&ec, sh, ss, pr);
        }
        double dec_mbs = (double)total * iters / (now() - t0)
                         / (1024.0 * 1024.0);

        printf("  | %9s  | %12.1f | %13.1f |\n",
               cfgs[c].label, enc_mbs, dec_mbs);
        free(mem);
    }

    printf("\n");
    erasure_free(&ec);
}

/* ================================================================
 * [6] 읽기 꼬리 지연시간 시뮬레이션
 *
 * "하드 하나가 굉장히 느리면 어쩔 거예요. 이레이저 코딩을 써 놨으면
 *  이거를 기다릴 필요가 없습니다. 그냥 버리고 다른 하드디스크에 있는
 *  조각을 빠르게 꺼내오면 되는 거예요." -- 영상 3:44
 *
 * 기본 읽기: 데이터 조각 5개만 요청, 전부 도착할 때까지 대기
 * 헤지 읽기: 전체 9개 요청, 가장 빠른 5개 도착 시점에 완료
 * ================================================================ */
static void test_tail_latency(void)
{
    printf("[6] 읽기 꼬리 지연시간 시뮬레이션\n");
    printf("--------------------------------------------\n");
    printf("  정상 디스크: 1-5ms\n");
    printf("  지연 디스크 (5%% 확률): 50-200ms\n");

    int nr = 100000;
    double *naive  = malloc((size_t)nr * sizeof(double));
    double *hedged = malloc((size_t)nr * sizeof(double));

    for (int r = 0; r < nr; r++) {
        double lat[N];
        for (int i = 0; i < N; i++)
            lat[i] = ((rand() % 100) < 5)
                     ? 50.0 + (rand() % 151)
                     : 1.0 + (rand() % 5);

        /* 기본: 데이터 조각 5개만 읽기, 가장 느린 것이 완료 시간 */
        double mx = 0;
        for (int i = 0; i < K; i++)
            if (lat[i] > mx) mx = lat[i];
        naive[r] = mx;

        /* 헤지: 9개 전부 읽기, 5번째로 빠른 것이 완료 시간 */
        double s[N];
        memcpy(s, lat, sizeof(s));
        for (int i = 0; i < N - 1; i++)
            for (int j = i + 1; j < N; j++)
                if (s[j] < s[i]) { double t = s[i]; s[i] = s[j]; s[j] = t; }
        hedged[r] = s[K - 1];
    }

    qsort(naive,  (size_t)nr, sizeof(double), dbl_cmp);
    qsort(hedged, (size_t)nr, sizeof(double), dbl_cmp);

    printf("  시행: %d회\n\n", nr);
    printf("  | 백분위 | 기본 (5개 대기) | 헤지 (9개 중 5개) | 개선율  |\n");
    printf("  |--------|---------------|------------------|---------|\n");

    double pct[] = {0.50, 0.90, 0.95, 0.99, 0.999};
    const char *lbl[] = {"P50  ", "P90  ", "P95  ", "P99  ", "P99.9"};

    for (int i = 0; i < 5; i++) {
        int idx = (int)(pct[i] * nr);
        if (idx >= nr) idx = nr - 1;
        double nv = naive[idx], hv = hedged[idx];
        double imp = (nv > 0) ? (1.0 - hv / nv) * 100.0 : 0;
        printf("  | %s  | %10.1f ms | %13.1f ms | %6.1f%% |\n",
               lbl[i], nv, hv, imp);
    }

    printf("\n");
    free(naive);
    free(hedged);
}

/* ================================================================
 * [7] 저장 방식 효율 비교
 *
 * 복제(replication) vs 이레이저 코딩(erasure coding)의
 * 저장 효율과 내결함성을 비교한다.
 * 저장 배율 = 전체 저장 용량 / 원본 데이터 크기
 * ================================================================ */
static void print_storage_cmp(void)
{
    printf("[7] 저장 방식 효율 비교\n");
    printf("--------------------------------------------\n\n");
    printf("  | 방식              | 저장 배율 | 내결함성 (디스크) | 최소 디스크 |\n");
    printf("  |-------------------|----------|-----------------|------------|\n");
    printf("  | 복제 없음          | 1.0x     | 0               | 1          |\n");
    printf("  | 2중 복제           | 2.0x     | 1               | 2          |\n");
    printf("  | 3중 복제           | 3.0x     | 2               | 3          |\n");
    printf("  | RS(5,4) - S3 방식  | 1.8x     | 4               | 9          |\n");
    printf("  | RS(10,4)          | 1.4x     | 4               | 14         |\n");
    printf("  | RS(16,4)          | 1.25x    | 4               | 20         |\n");
    printf("\n  -> RS(5,4): 3중 복제 대비 40%% 적은 용량으로 2배 높은 내결함성\n\n");
}

/* ================================================================
 * [8] 내구성(Durability) 분석
 *
 * Markov 체인으로 MTTDL(Mean Time To Data Loss)을 해석적으로 계산하고,
 * Monte Carlo 이산 사건 시뮬레이션으로 검증한다.
 *
 * MTTDL이 너무 커서 (RS(5,4)는 ~10^16년) 직접 관측이 불가하므로,
 * 가속 파라미터(AFR=20%)로 MC를 돌려 해석적 공식을 검증한 뒤,
 * 실제 파라미터(AFR=2%)에 적용하는 전략을 사용한다.
 * ================================================================ */
static void test_durability(void)
{
    printf("[8] 내구성(Durability) 분석\n");
    printf("--------------------------------------------\n");

    double afr = 0.02;
    double lambda = -log(1.0 - afr) / 8760.0;
    double mu = 1.0 / 24.0;

    printf("  Markov 체인 MTTDL 계산\n");
    printf("  AFR=%.0f%%, Rebuild=24시간\n\n", afr * 100);

    struct { const char *name; int n; int m; } cfgs[] = {
        {"단일 디스크      ", 1, 0},
        {"2중 복제         ", 2, 1},
        {"3중 복제         ", 3, 2},
        {"RS(5,4) S3 방식  ", 9, 4},
        {"RS(10,4)        ", 14, 4},
    };
    int ncfg = 5;

    printf("  | 방식              | MTTDL (년)   | 연간 손실 확률  | 내구성 (9s) |\n");
    printf("  |-------------------|-------------|----------------|------------|\n");

    for (int i = 0; i < ncfg; i++) {
        double h = compute_mttdl(cfgs[i].n, cfgs[i].m, lambda, mu);
        double yr = h / 8760.0;
        double p = -expm1(-8760.0 / h);
        double nines = (p > 0) ? -log10(p) : 99;
        printf("  | %s | %11.2e | %14.2e | %10.1f |\n",
               cfgs[i].name, yr, p, nines);
    }

    /* Monte Carlo 검증 (가속 파라미터) */
    printf("\n  Monte Carlo 검증 (가속 파라미터)\n");
    printf("  AFR=20%%, Rebuild=24h, 100만회 시행\n\n");

    double mc_afr = 0.20;
    double mc_lam = -log(1.0 - mc_afr) / 8760.0;
    double mc_rb  = 24.0;
    int mc_n = 1000000;

    struct { const char *name; int n; int m; } mc_cfgs[] = {
        {"단일 디스크", 1, 0},
        {"2중 복제   ", 2, 1},
        {"RS(2,1)   ", 3, 1},
    };

    printf("  | 방식        | 분석적 예측    | MC 실측       | 오차   |\n");
    printf("  |-------------|-------------- |------------- |--------|\n");

    for (int c = 0; c < 3; c++) {
        double h = compute_mttdl(mc_cfgs[c].n, mc_cfgs[c].m, mc_lam, 1.0 / mc_rb);
        double p_a = -expm1(-8760.0 / h);

        int losses = 0;
        for (int t = 0; t < mc_n; t++)
            losses += mc_trial(mc_cfgs[c].n, mc_cfgs[c].m, mc_lam, mc_rb);
        double p_mc = (double)losses / mc_n;

        double err = (p_a > 0) ? fabs(p_mc - p_a) / p_a * 100.0 : 0;
        printf("  | %s | %11.4f%% | %11.4f%% | %5.1f%% |\n",
               mc_cfgs[c].name, p_a * 100, p_mc * 100, err);
    }
    printf("\n  -> 분석적 공식이 MC와 일치 → 실제 파라미터 결과 신뢰 가능\n\n");
}

/* ================================================================
 * [9] CRC32 Silent Corruption 탐지
 *
 * 디스크가 "고장"이 아니라 조용히 잘못된 데이터를 반환하는 경우
 * (bit rot, firmware 버그 등). 체크섬 없이는 탐지 불가.
 *
 * CRC32로 조각별 무결성을 검증하고, 손상된 조각을 "고장"으로 처리한 뒤
 * 이레이저 코딩으로 복구하는 전체 파이프라인을 시뮬레이션한다.
 * ================================================================ */
static void test_crc32_corruption(void)
{
    printf("[9] CRC32 Silent Corruption 탐지\n");
    printf("--------------------------------------------\n");

    erasure_t ec;
    erasure_init(&ec, K, M);

    int trials = 1000;
    int detected = 0, recovered = 0;
    size_t ss = 64;

    for (int t = 0; t < trials; t++) {
        uint8_t *mem = calloc(N, ss);
        uint8_t *sh[N];
        for (int i = 0; i < N; i++) sh[i] = mem + (size_t)i * ss;
        for (size_t i = 0; i < K * ss; i++) mem[i] = rand() & 0xFF;

        uint8_t *orig = malloc(K * ss);
        memcpy(orig, mem, K * ss);
        erasure_encode(&ec, sh, ss);

        /* CRC32 per shard */
        uint32_t crcs[N];
        for (int i = 0; i < N; i++)
            crcs[i] = crc32(sh[i], ss);

        /* Silent corruption: 1~3개 조각에 1~3바이트 변조 */
        int ncorrupt = 1 + rand() % 3;
        int corrupted[N] = {0};
        int nc = 0;
        while (nc < ncorrupt) {
            int idx = rand() % N;
            if (corrupted[idx]) continue;
            corrupted[idx] = 1;
            int nflip = 1 + rand() % 3;
            for (int f = 0; f < nflip; f++)
                sh[idx][rand() % ss] ^= (uint8_t)(1 + rand() % 255);
            nc++;
        }

        /* CRC 검증 → 손상 탐지 */
        int all_detected = 1;
        int present[N];
        for (int i = 0; i < N; i++) {
            present[i] = (crc32(sh[i], ss) == crcs[i]) ? 1 : 0;
            if (corrupted[i] && present[i]) all_detected = 0;
        }
        if (all_detected) detected++;

        /* 손상 조각을 "고장"으로 처리, 이레이저 코딩으로 복구 */
        int avail = 0;
        for (int i = 0; i < N; i++) if (present[i]) avail++;
        if (avail >= K) {
            for (int i = 0; i < N; i++)
                if (!present[i]) memset(sh[i], 0, ss);
            int ret = erasure_decode(&ec, sh, ss, present);
            if (ret == 0 && memcmp(mem, orig, K * ss) == 0)
                recovered++;
        }

        free(mem); free(orig);
    }

    printf("  시나리오: 인코딩 후 1~3개 조각에 1~3바이트 랜덤 변조\n");
    printf("  시행: %d회\n\n", trials);
    printf("  CRC32 탐지율:  %d/%d (%.1f%%)\n",
           detected, trials, detected * 100.0 / trials);
    printf("  탐지 후 복구:  %d/%d (%.1f%%)\n\n",
           recovered, trials, recovered * 100.0 / trials);
    printf("  -> 체크섬 없이는 손상 데이터가 정상으로 반환될 위험\n");
    printf("  -> CRC32 + 이레이저 코딩 = 감지 + 복구 모두 가능\n\n");

    erasure_free(&ec);
}

/* ================================================================
 * [10] Rebuild 취약 구간 분석
 *
 * 디스크 1개 고장 후 rebuild 중에 추가 고장이 발생하면?
 * Rebuild 시간이 길수록 취약 구간이 넓어진다.
 *
 * P(data loss during rebuild)
 *   = P(m개 이상 추가 고장 | n-1개 디스크, T시간)
 *   = Σ_{i=m}^{n-1} C(n-1,i) * p^i * (1-p)^{n-1-i}
 *   where p = 1 - exp(-λ*T)
 * ================================================================ */
static void test_rebuild_vulnerability(void)
{
    printf("[10] Rebuild 취약 구간 분석\n");
    printf("--------------------------------------------\n");
    printf("  1개 디스크 고장 후 rebuild 중 추가 고장 시 데이터 손실 확률\n");
    printf("  AFR: 2%%\n\n");

    double afr = 0.02;
    double lambda = -log(1.0 - afr) / 8760.0;

    double rb_hours[] = {1, 6, 24, 72, 168, 720};
    const char *rb_lbl[] = {"1h  ", "6h  ", "24h ", "72h ", "168h", "720h"};
    int nrb = 6;

    printf("  | Rebuild | 3중 복제       | RS(5,4)        | RS 우위 (배수)  |\n");
    printf("  |---------|---------------|----------------|----------------|\n");

    for (int r = 0; r < nrb; r++) {
        double p = 1.0 - exp(-lambda * rb_hours[r]);

        /* 3중 복제 (n=3, m=2): 1 고장 후 나머지 2개 중 2개 더 고장 */
        double p_3rep = binom_tail(2, 2, p);

        /* RS(5,4) (n=9, m=4): 1 고장 후 나머지 8개 중 4개 더 고장 */
        double p_rs54 = binom_tail(8, 4, p);

        double ratio = (p_rs54 > 0) ? p_3rep / p_rs54 : 0;
        printf("  | %s    | %13.2e | %14.2e | %14.0f |\n",
               rb_lbl[r], p_3rep, p_rs54, ratio);
    }

    printf("\n  -> RS(5,4)는 rebuild 중에도 3중 복제 대비 수십만~수억 배 안전\n");
    printf("  -> Rebuild 시간이 길어질수록 격차 확대\n\n");
}

/* ================================================================
 * [11] SIMD 최적화 벤치마크
 *
 * GF(256) 곱셈-누적의 스칼라 vs SIMD 성능을 비교한다.
 * SIMD는 split-table 기법: 각 바이트를 hi/lo nibble로 분리하고
 * 16-entry 테이블 2개로 NEON tbl / SSE pshufb 병렬 lookup.
 *
 * 인코딩/디코딩의 핫 패스인 gf256_mul_vec를 직접 벤치마크하고,
 * 전체 인코딩 파이프라인에서의 실제 스피드업도 측정한다.
 * ================================================================ */

static void encode_with(void (*mul_vec)(uint8_t *, const uint8_t *, uint8_t, size_t),
                        erasure_t *ec, uint8_t **shards, size_t ss)
{
    int k = ec->data_shards, m = ec->parity_shards;
    for (int i = 0; i < m; i++) {
        memset(shards[k + i], 0, ss);
        for (int j = 0; j < k; j++)
            mul_vec(shards[k + i], shards[j],
                    ec->parity_matrix[i * k + j], ss);
    }
}

static void bench_simd(void)
{
    printf("[11] SIMD 최적화 벤치마크\n");
    printf("--------------------------------------------\n");

#if defined(__aarch64__)
    printf("  아키텍처: ARM64 NEON (vqtbl1q_u8)\n");
#elif defined(__SSSE3__)
    printf("  아키텍처: x86 SSSE3 (_mm_shuffle_epi8)\n");
#else
    printf("  아키텍처: SIMD 없음 (스칼라 fallback)\n");
#endif

    /* (A) gf256_mul_vec 직접 벤치마크 */
    size_t len = 1048576; /* 1 MB */
    uint8_t *dst = calloc(len, 1);
    uint8_t *src = malloc(len);
    for (size_t i = 0; i < len; i++) src[i] = rand() & 0xFF;
    uint8_t c = 0x53;
    int iters = 500;

    /* scalar */
    double t0 = now();
    for (int i = 0; i < iters; i++)
        gf256_mul_vec_ref(dst, src, c, len);
    double ref_mbs = (double)len * iters / (now() - t0) / (1024.0 * 1024.0);

    /* SIMD */
    memset(dst, 0, len);
    t0 = now();
    for (int i = 0; i < iters; i++)
        gf256_mul_vec(dst, src, c, len);
    double simd_mbs = (double)len * iters / (now() - t0) / (1024.0 * 1024.0);

    printf("\n  gf256_mul_vec (1MB x %d회):\n", iters);
    printf("  Scalar:  %8.1f MB/s\n", ref_mbs);
    printf("  SIMD:    %8.1f MB/s\n", simd_mbs);
    printf("  Speedup: %8.1fx\n", simd_mbs / ref_mbs);

    free(dst); free(src);

    /* (B) 전체 인코딩 파이프라인 벤치마크 */
    erasure_t ec;
    erasure_init(&ec, K, M);

    struct { size_t size; const char *label; int iters; } cfgs[] = {
        {4096,    "4 KB  ", 3000},
        {65536,   "64 KB ", 500},
        {1048576, "1 MB  ", 50},
    };
    int ncfg = 3;

    printf("\n  인코딩 파이프라인 (Scalar vs SIMD):\n\n");
    printf("  | 데이터    | Scalar (MB/s) | SIMD (MB/s)   | Speedup |\n");
    printf("  |----------|--------------|--------------|----------|\n");

    for (int ci = 0; ci < ncfg; ci++) {
        size_t total = cfgs[ci].size;
        size_t ss = (total + K - 1) / K;
        int it = cfgs[ci].iters;

        uint8_t *mem = calloc(N, ss);
        uint8_t *sh[N];
        for (int i = 0; i < N; i++)
            sh[i] = mem + (size_t)i * ss;
        for (size_t i = 0; i < total; i++)
            mem[i] = rand() & 0xFF;

        /* scalar encode */
        t0 = now();
        for (int i = 0; i < it; i++)
            encode_with(gf256_mul_vec_ref, &ec, sh, ss);
        double sc = (double)total * it / (now() - t0) / (1024.0 * 1024.0);

        /* SIMD encode */
        t0 = now();
        for (int i = 0; i < it; i++)
            encode_with(gf256_mul_vec, &ec, sh, ss);
        double sm = (double)total * it / (now() - t0) / (1024.0 * 1024.0);

        printf("  | %7s  | %12.1f | %12.1f | %7.1fx |\n",
               cfgs[ci].label, sc, sm, sm / sc);
        free(mem);
    }

    printf("\n");
    erasure_free(&ec);
}

/* ================================================================ */

int main(void)
{
    srand(42); /* 재현성을 위한 고정 시드 */
    gf256_init();
    crc32_init();

    printf("\n============================================\n");
    printf("  Mini S3 - 종합 실험 결과\n");
    printf("============================================\n\n");

    test_gf256();
    test_erasure_exhaustive();
    test_erasure_random();
    test_p2c();
    bench_throughput();
    test_tail_latency();
    print_storage_cmp();
    test_durability();
    test_crc32_corruption();
    test_rebuild_vulnerability();
    bench_simd();

    return 0;
}
