#ifndef ERASURE_H
#define ERASURE_H

#include <stdint.h>
#include <stddef.h>

typedef struct {
    int      data_shards;   /* k: 원본 데이터 조각 수 */
    int      parity_shards; /* m: 패리티(복구용) 조각 수 */
    int      total_shards;  /* n = k + m */
    uint8_t *parity_matrix; /* m x k Cauchy 행렬 */
} erasure_t;

int  erasure_init(erasure_t *ec, int data_shards, int parity_shards);
void erasure_free(erasure_t *ec);

/* 패리티 조각 생성: shards[k..n-1]에 패리티 기록 */
void erasure_encode(erasure_t *ec, uint8_t **shards, size_t shard_size);

/* 손실된 조각 복구: present[i]=0인 조각을 복원 */
int  erasure_decode(erasure_t *ec, uint8_t **shards, size_t shard_size,
                    const int *present);

#endif
