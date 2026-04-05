/*
 * Reed-Solomon 이레이저 코딩 (Cauchy 행렬 기반)
 *
 * 영상 설명: "데이터를 5개로 쪼갤 때 거기 있는 정보들을 이용해서
 *            4개의 덩어리를 추가로 만들 수가 있습니다. 요걸 패리티라고 부르고요.
 *            그러면 이 덩어리들 중에서 아무거나 5개를 뽑아도
 *            데이터를 원상 복구가 가능하게 만들 수가 있습니다."
 *
 * 인코딩 행렬 E (n x k):
 *   [ I_k ]   <- 데이터 조각은 그대로 저장
 *   [  P  ]   <- 패리티 = Cauchy 행렬로 계산
 *
 * Cauchy 행렬: P[i][j] = 1/(x_i XOR y_j) in GF(256)
 * 이 행렬의 모든 정방 부분행렬이 역행렬을 가짐 → MDS 코드 보장
 */
#include "erasure.h"
#include "gf256.h"
#include <stdlib.h>
#include <string.h>

int erasure_init(erasure_t *ec, int data_shards, int parity_shards)
{
    ec->data_shards   = data_shards;
    ec->parity_shards = parity_shards;
    ec->total_shards  = data_shards + parity_shards;

    int k = data_shards, m = parity_shards;
    ec->parity_matrix = malloc((size_t)(m * k));
    if (!ec->parity_matrix)
        return -1;

    /*
     * Cauchy 행렬 구성
     * x = {0, 1, ..., m-1},  y = {m, m+1, ..., m+k-1}
     * P[i][j] = gf256_inv(i ^ (m + j))
     *
     * x와 y 집합이 겹치지 않으므로 i^(m+j) != 0 보장
     */
    for (int i = 0; i < m; i++)
        for (int j = 0; j < k; j++)
            ec->parity_matrix[i * k + j] =
                gf256_inv((uint8_t)(i ^ (m + j)));

    return 0;
}

void erasure_free(erasure_t *ec)
{
    free(ec->parity_matrix);
    ec->parity_matrix = NULL;
}

void erasure_encode(erasure_t *ec, uint8_t **shards, size_t shard_size)
{
    int k = ec->data_shards;
    int m = ec->parity_shards;

    /*
     * parity[i][b] = SUM_j( P[i][j] * data[j][b] )  in GF(256)
     *
     * 영상 비유: "다항식에 있는 미지수 5개를 찾아내는 문제"
     * 점(조각)을 여러 개 만들어 놓으면 아무 k개로 복구 가능
     */
    for (int i = 0; i < m; i++) {
        memset(shards[k + i], 0, shard_size);
        for (int j = 0; j < k; j++)
            gf256_mul_vec(shards[k + i], shards[j],
                          ec->parity_matrix[i * k + j], shard_size);
    }
}

/* GF(256) 위의 k x k 행렬 역원 (가우스 소거법) */
static int gf256_mat_inv(uint8_t *mat, uint8_t *inv, int k)
{
    for (int col = 0; col < k; col++) {
        /* 피벗 탐색 */
        int pivot = -1;
        for (int r = col; r < k; r++) {
            if (mat[r * k + col]) { pivot = r; break; }
        }
        if (pivot < 0)
            return -1;

        /* 행 교환 */
        if (pivot != col) {
            for (int j = 0; j < k; j++) {
                uint8_t t;
                t = mat[col*k+j]; mat[col*k+j] = mat[pivot*k+j]; mat[pivot*k+j] = t;
                t = inv[col*k+j]; inv[col*k+j] = inv[pivot*k+j]; inv[pivot*k+j] = t;
            }
        }

        /* 피벗 행 스케일링 */
        uint8_t s = gf256_inv(mat[col * k + col]);
        for (int j = 0; j < k; j++) {
            mat[col*k+j] = gf256_mul(mat[col*k+j], s);
            inv[col*k+j] = gf256_mul(inv[col*k+j], s);
        }

        /* 소거 */
        for (int r = 0; r < k; r++) {
            if (r == col) continue;
            uint8_t f = mat[r * k + col];
            if (!f) continue;
            for (int j = 0; j < k; j++) {
                mat[r*k+j] ^= gf256_mul(f, mat[col*k+j]);
                inv[r*k+j] ^= gf256_mul(f, inv[col*k+j]);
            }
        }
    }
    return 0;
}

int erasure_decode(erasure_t *ec, uint8_t **shards, size_t shard_size,
                   const int *present)
{
    int k = ec->data_shards;
    int n = ec->total_shards;

    /* 사용 가능한 조각 수 확인 */
    int avail = 0;
    for (int i = 0; i < n; i++)
        if (present[i]) avail++;
    if (avail < k)
        return -1;

    /* 모든 데이터 조각이 살아있으면 복구 불필요 */
    int need = 0;
    for (int i = 0; i < k; i++)
        if (!present[i]) { need = 1; break; }
    if (!need)
        return 0;

    /* 사용할 k개 조각 선택 */
    int *sel = malloc((size_t)k * sizeof(int));
    int s = 0;
    for (int i = 0; i < n && s < k; i++)
        if (present[i]) sel[s++] = i;

    /*
     * 인코딩 행렬 E = [I_k ; P] 에서 선택된 행으로 부분행렬 구성
     * sel[i] < k  -> 단위행렬의 sel[i]번째 행
     * sel[i] >= k -> 패리티 행렬의 (sel[i]-k)번째 행
     */
    uint8_t *submat = calloc((size_t)(k * k), 1);
    for (int i = 0; i < k; i++) {
        int row = sel[i];
        if (row < k)
            submat[i * k + row] = 1;
        else
            memcpy(&submat[i * k],
                   &ec->parity_matrix[(row - k) * k], (size_t)k);
    }

    /* 부분행렬의 역행렬 계산 */
    uint8_t *inv_mat = calloc((size_t)(k * k), 1);
    for (int i = 0; i < k; i++)
        inv_mat[i * k + i] = 1;

    int ret = gf256_mat_inv(submat, inv_mat, k);
    if (ret != 0) {
        free(sel); free(submat); free(inv_mat);
        return -1;
    }

    /*
     * 손실된 데이터 조각 복원:
     * data[j] = SUM_i( inv[j][i] * available[i] )
     */
    for (int j = 0; j < k; j++) {
        if (present[j]) continue;
        memset(shards[j], 0, shard_size);
        for (int i = 0; i < k; i++)
            gf256_mul_vec(shards[j], shards[sel[i]],
                          inv_mat[j * k + i], shard_size);
    }

    free(sel);
    free(submat);
    free(inv_mat);
    return 0;
}
