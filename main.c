/*
 * Mini S3 시뮬레이터
 * AWS S3의 핵심 기술 3가지를 C로 구현
 *
 * [1] 데이터 분산 저장 (Data Striping)
 * [2] Power of Two Choices 로드 밸런싱
 * [3] 이레이저 코딩 (Erasure Coding)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "gf256.h"
#include "erasure.h"

#define DATA_SHARDS   5
#define PARITY_SHARDS 4
#define TOTAL_SHARDS  (DATA_SHARDS + PARITY_SHARDS)
#define NUM_DISKS     16

/* ---- 유틸리티 ---- */

static void print_hex(const uint8_t *data, size_t len, size_t max)
{
    size_t n = len < max ? len : max;
    for (size_t i = 0; i < n; i++)
        printf("%02x ", data[i]);
    if (len > max)
        printf("...");
}

static void print_ascii(const uint8_t *data, size_t len, size_t max)
{
    size_t n = len < max ? len : max;
    putchar('"');
    for (size_t i = 0; i < n; i++)
        putchar((data[i] >= 32 && data[i] < 127) ? data[i] : '.');
    if (len > max)
        printf("...");
    putchar('"');
}

static double stddev(const int *arr, int n)
{
    double mean = 0;
    for (int i = 0; i < n; i++) mean += arr[i];
    mean /= n;
    double var = 0;
    for (int i = 0; i < n; i++) {
        double d = arr[i] - mean;
        var += d * d;
    }
    return sqrt(var / n);
}

/* ================================================================
 * [1] 데이터 분산 저장
 * "데이터를 받으면 아주 잘게 쪼개 버리고, 여러 하드에 분산시켜서
 *  저장해 두면 끝이에요. 5개로 쪼개서 각각 다른 하드에 저장하면
 *  저장 속도가 5배 가까이 증가할 수 있는 거고요."
 * ================================================================ */
static void demo_striping(void)
{
    printf("[1] 데이터 분산 저장 (Data Striping)\n");
    printf("--------------------------------------------\n");

    const char *msg =
        "Hello AWS S3! Data is split across multiple disks for speed.";
    size_t len = strlen(msg);
    size_t shard_sz = (len + DATA_SHARDS - 1) / DATA_SHARDS;

    printf("원본: \"%s\" (%zu bytes)\n", msg, len);
    printf("%d개 조각으로 분할 (조각 크기: %zu bytes)\n\n", DATA_SHARDS, shard_sz);

    int usage[NUM_DISKS] = {0};

    for (int i = 0; i < DATA_SHARDS; i++) {
        size_t off = (size_t)i * shard_sz;
        size_t sz  = (off + shard_sz <= len) ? shard_sz : len - off;

        /* Power of Two Choices로 디스크 선택 */
        int a = rand() % NUM_DISKS;
        int b;
        do { b = rand() % NUM_DISKS; } while (b == a);
        int disk = (usage[a] <= usage[b]) ? a : b;
        usage[disk] += (int)sz;

        printf("  조각 %d -> 디스크 %2d: ", i, disk);
        print_hex((const uint8_t *)msg + off, sz, 8);
        printf("  ");
        print_ascii((const uint8_t *)msg + off, sz, 14);
        printf("\n");
    }

    printf("\n-> %d개 디스크에서 병렬 읽기 = 속도 %dx!\n\n", DATA_SHARDS, DATA_SHARDS);
}

/* ================================================================
 * [2] Power of Two Choices
 * "데이터를 저장할 때 랜덤하게 하나를 고르는게 아니라 두 개를 고르는 거예요.
 *  그다음에 두 개를 비교해서 사용률이 적은 하드에 저장해 버리는 거예요.
 *  이러면 놀랍게도 모든 하드들의 사용률이 이런 식으로 찍히게 됩니다."
 * ================================================================ */
static void demo_power_of_two(void)
{
    printf("[2] Power of Two Choices 로드 밸런싱\n");
    printf("--------------------------------------------\n");

    const int ndisks = NUM_DISKS;
    const int nitems = 10000;

    /* (A) 랜덤 배치 */
    int rand_cnt[NUM_DISKS] = {0};
    for (int i = 0; i < nitems; i++)
        rand_cnt[rand() % ndisks]++;

    /* (B) Power of Two Choices */
    int p2c_cnt[NUM_DISKS] = {0};
    for (int i = 0; i < nitems; i++) {
        int a = rand() % ndisks;
        int b;
        do { b = rand() % ndisks; } while (b == a);
        p2c_cnt[(p2c_cnt[a] <= p2c_cnt[b]) ? a : b]++;
    }

    /* 바 차트 스케일링 */
    int mx = 0;
    for (int i = 0; i < ndisks; i++) {
        if (rand_cnt[i] > mx) mx = rand_cnt[i];
        if (p2c_cnt[i] > mx)  mx = p2c_cnt[i];
    }

    printf("%d개 데이터를 %d개 디스크에 분배:\n\n", nitems, ndisks);

    /* 랜덤 */
    printf("랜덤 배치:\n");
    for (int i = 0; i < ndisks; i++) {
        printf("  디스크 %2d: ", i);
        int w = (int)((double)rand_cnt[i] / mx * 35);
        for (int j = 0; j < w; j++) putchar('#');
        for (int j = w; j < 35; j++) putchar(' ');
        printf(" %d\n", rand_cnt[i]);
    }
    printf("  표준편차: %.1f\n\n", stddev(rand_cnt, ndisks));

    /* P2C */
    printf("Power of Two Choices:\n");
    for (int i = 0; i < ndisks; i++) {
        printf("  디스크 %2d: ", i);
        int w = (int)((double)p2c_cnt[i] / mx * 35);
        for (int j = 0; j < w; j++) putchar('#');
        for (int j = w; j < 35; j++) putchar(' ');
        printf(" %d\n", p2c_cnt[i]);
    }
    double p2c_sd  = stddev(p2c_cnt, ndisks);
    double rand_sd = stddev(rand_cnt, ndisks);
    printf("  표준편차: %.1f\n\n", p2c_sd);

    double pct = (1.0 - p2c_sd / rand_sd) * 100.0;
    printf("-> 표준편차 %.0f%% 감소! 훨씬 균등한 분배!\n\n", pct);
}

/* ================================================================
 * [3] 이레이저 코딩 (Erasure Coding)
 * "이레이저 코딩을 사용합니다. 데이터를 5개로 쪼갤 때 거기 있는
 *  정보들을 이용해서 4개의 덩어리를 추가로 만들 수가 있습니다.
 *  요걸 패리티라고 부르고요. 그러면 이 덩어리들 중에서 아무거나
 *  5개를 뽑아도 데이터를 원상 복구가 가능하게 만들 수가 있습니다."
 * ================================================================ */
static void demo_erasure(void)
{
    printf("[3] 이레이저 코딩 (Erasure Coding)\n");
    printf("--------------------------------------------\n");

    erasure_t ec;
    erasure_init(&ec, DATA_SHARDS, PARITY_SHARDS);

    printf("설정: %d 데이터 + %d 패리티 = %d 조각\n",
           DATA_SHARDS, PARITY_SHARDS, TOTAL_SHARDS);
    printf("-> 최대 %d개 디스크 고장까지 복구 가능\n\n", PARITY_SHARDS);

    /* 패리티 행렬 출력 */
    printf("패리티 행렬 (Cauchy matrix over GF(256)):\n");
    for (int i = 0; i < PARITY_SHARDS; i++) {
        printf("  P[%d]: ", i);
        for (int j = 0; j < DATA_SHARDS; j++)
            printf("%02x ", ec.parity_matrix[i * DATA_SHARDS + j]);
        printf("\n");
    }
    printf("\n");

    /* 데이터 준비 */
    const char *msg =
        "Hello AWS S3! Erasure coding restores data even when disks fail!";
    size_t msg_len    = strlen(msg);
    size_t shard_size = (msg_len + DATA_SHARDS - 1) / DATA_SHARDS;

    uint8_t *shard_mem = calloc(TOTAL_SHARDS, shard_size);
    uint8_t *shards[TOTAL_SHARDS];
    for (int i = 0; i < TOTAL_SHARDS; i++)
        shards[i] = shard_mem + (size_t)i * shard_size;

    /* 데이터를 조각으로 분할 */
    memcpy(shard_mem, msg, msg_len);

    printf("원본: \"%s\" (%zu bytes)\n", msg, msg_len);
    printf("조각 크기: %zu bytes\n\n", shard_size);

    /* 인코딩 */
    erasure_encode(&ec, shards, shard_size);

    printf("인코딩 결과:\n");
    for (int i = 0; i < TOTAL_SHARDS; i++) {
        const char *type = (i < DATA_SHARDS) ? "DATA  " : "PARITY";
        int idx = (i < DATA_SHARDS) ? i : i - DATA_SHARDS;
        printf("  [%s %d] ", type, idx);
        print_hex(shards[i], shard_size, 10);
        if (i < DATA_SHARDS) {
            printf("  ");
            print_ascii(shards[i], shard_size, 14);
        }
        printf("\n");
    }
    printf("\n");

    /* 원본 백업 (검증용) */
    uint8_t *original = malloc(DATA_SHARDS * shard_size);
    memcpy(original, shard_mem, DATA_SHARDS * shard_size);

    /* 디스크 고장 시뮬레이션 */
    int present[TOTAL_SHARDS];
    for (int i = 0; i < TOTAL_SHARDS; i++)
        present[i] = 1;

    printf("디스크 고장 시뮬레이션 (%d개 파괴):\n", PARITY_SHARDS);
    int destroyed = 0;
    while (destroyed < PARITY_SHARDS) {
        int idx = rand() % TOTAL_SHARDS;
        if (!present[idx]) continue;
        present[idx] = 0;
        memset(shards[idx], 0, shard_size);

        const char *type = (idx < DATA_SHARDS) ? "데이터" : "패리티";
        int num = (idx < DATA_SHARDS) ? idx : idx - DATA_SHARDS;
        printf("  X 조각 %d (%s %d) 파괴!\n", idx, type, num);
        destroyed++;
    }

    printf("\n사용 가능한 조각: ");
    int avail = 0;
    for (int i = 0; i < TOTAL_SHARDS; i++) {
        if (present[i]) {
            printf("%d ", i);
            avail++;
        }
    }
    printf("(%d개)\n", avail);

    /* 복구 */
    printf("\n복구 중... ");
    int ret = erasure_decode(&ec, shards, shard_size, present);

    if (ret == 0) {
        int match = (memcmp(shard_mem, original, DATA_SHARDS * shard_size) == 0);
        printf("%s\n\n", match ? "성공!" : "실패 (데이터 불일치)");

        if (match) {
            char recovered[256] = {0};
            memcpy(recovered, shard_mem, msg_len);
            printf("복구된 데이터: \"%s\"\n", recovered);
        }
    } else {
        printf("실패! 조각이 부족합니다.\n");
    }

    free(shard_mem);
    free(original);
    erasure_free(&ec);
}

/* ================================================================ */

int main(void)
{
    srand((unsigned)time(NULL));
    gf256_init();

    printf("\n============================================\n");
    printf("  Mini S3 - AWS S3 핵심 기술 구현 (C)\n");
    printf("============================================\n\n");

    demo_striping();
    demo_power_of_two();
    demo_erasure();

    printf("\n");
    return 0;
}
