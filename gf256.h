#ifndef GF256_H
#define GF256_H

#include <stdint.h>
#include <stddef.h>

void    gf256_init(void);
uint8_t gf256_mul(uint8_t a, uint8_t b);
uint8_t gf256_inv(uint8_t a);

/*
 * 벡터 곱셈-누적: dst[i] ^= c * src[i]  (GF(256))
 *
 * 이레이저 코딩의 핫 패스. 인코딩/디코딩 모두 이 연산이 지배.
 * SIMD 버전은 split-table 기법을 사용:
 *   c * b = lo_tbl[b & 0x0F] ^ hi_tbl[b >> 4]
 *   16-byte 테이블 × pshufb/tbl 명령 → 16바이트 동시 처리
 */
void gf256_mul_vec(uint8_t *dst, const uint8_t *src, uint8_t c, size_t len);

/* 스칼라 레퍼런스 (벤치마크 비교용) */
void gf256_mul_vec_ref(uint8_t *dst, const uint8_t *src, uint8_t c, size_t len);

#endif
