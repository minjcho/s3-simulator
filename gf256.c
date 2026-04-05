/*
 * GF(256) 갈루아 필드 연산
 * 원시 다항식: x^8 + x^4 + x^3 + x^2 + 1 (0x11D)
 * 이레이저 코딩의 수학적 기반
 */
#include "gf256.h"
#include <string.h>

static uint8_t exp_tbl[512];
static uint8_t log_tbl[256];

void gf256_init(void)
{
    uint16_t x = 1;
    for (int i = 0; i < 255; i++) {
        exp_tbl[i]       = (uint8_t)x;
        exp_tbl[i + 255] = (uint8_t)x;
        log_tbl[(uint8_t)x] = (uint8_t)i;
        x <<= 1;
        if (x & 0x100)
            x ^= 0x11D;
    }
    log_tbl[0] = 0;
}

uint8_t gf256_mul(uint8_t a, uint8_t b)
{
    if (a == 0 || b == 0)
        return 0;
    return exp_tbl[log_tbl[a] + log_tbl[b]];
}

uint8_t gf256_inv(uint8_t a)
{
    return exp_tbl[255 - log_tbl[a]];
}

/* ---- 벡터 곱셈-누적: dst[i] ^= c * src[i] ---- */

void gf256_mul_vec_ref(uint8_t *dst, const uint8_t *src, uint8_t c, size_t len)
{
    if (c == 0) return;
    if (c == 1) {
        for (size_t i = 0; i < len; i++) dst[i] ^= src[i];
        return;
    }
    for (size_t i = 0; i < len; i++)
        dst[i] ^= gf256_mul(c, src[i]);
}

/*
 * SIMD split-table 기법:
 *
 * GF(256)에서 c * b 를 계산할 때, b를 두 nibble로 분해:
 *   b = (b >> 4) << 4 | (b & 0x0F)
 *
 * 분배법칙에 의해:
 *   c * b = c * (hi_part + lo_part) = c * hi_part ^ c * lo_part
 *
 * 16-entry 테이블 2개를 미리 계산:
 *   lo_tbl[i] = c * i       (i = 0..15)
 *   hi_tbl[i] = c * (i<<4)  (i = 0..15)
 *
 * NEON vqtbl1q / SSE pshufb 로 16바이트 동시 테이블 lookup
 */

#if defined(__aarch64__)
/* ---- ARM NEON (Apple Silicon, AWS Graviton, etc.) ---- */
#include <arm_neon.h>

void gf256_mul_vec(uint8_t *dst, const uint8_t *src, uint8_t c, size_t len)
{
    if (c == 0) return;
    if (c == 1) {
        size_t i = 0;
        for (; i + 16 <= len; i += 16) {
            uint8x16_t d = vld1q_u8(dst + i);
            uint8x16_t s = vld1q_u8(src + i);
            vst1q_u8(dst + i, veorq_u8(d, s));
        }
        for (; i < len; i++) dst[i] ^= src[i];
        return;
    }

    /* split-table 구성 */
    uint8_t lo_arr[16], hi_arr[16];
    for (int i = 0; i < 16; i++) {
        lo_arr[i] = gf256_mul(c, (uint8_t)i);
        hi_arr[i] = gf256_mul(c, (uint8_t)(i << 4));
    }

    uint8x16_t lo_tbl = vld1q_u8(lo_arr);
    uint8x16_t hi_tbl = vld1q_u8(hi_arr);
    uint8x16_t mask   = vdupq_n_u8(0x0F);

    size_t i = 0;
    for (; i + 16 <= len; i += 16) {
        uint8x16_t s = vld1q_u8(src + i);
        uint8x16_t d = vld1q_u8(dst + i);

        uint8x16_t lo_val = vqtbl1q_u8(lo_tbl, vandq_u8(s, mask));
        uint8x16_t hi_val = vqtbl1q_u8(hi_tbl, vshrq_n_u8(s, 4));

        vst1q_u8(dst + i, veorq_u8(d, veorq_u8(lo_val, hi_val)));
    }
    for (; i < len; i++)
        dst[i] ^= gf256_mul(c, src[i]);
}

#elif defined(__SSSE3__)
/* ---- x86 SSSE3 (Intel/AMD) ---- */
#include <tmmintrin.h>

void gf256_mul_vec(uint8_t *dst, const uint8_t *src, uint8_t c, size_t len)
{
    if (c == 0) return;
    if (c == 1) {
        size_t i = 0;
        for (; i + 16 <= len; i += 16) {
            __m128i d = _mm_loadu_si128((__m128i *)(dst + i));
            __m128i s = _mm_loadu_si128((__m128i *)(src + i));
            _mm_storeu_si128((__m128i *)(dst + i), _mm_xor_si128(d, s));
        }
        for (; i < len; i++) dst[i] ^= src[i];
        return;
    }

    uint8_t lo_arr[16], hi_arr[16];
    for (int i = 0; i < 16; i++) {
        lo_arr[i] = gf256_mul(c, (uint8_t)i);
        hi_arr[i] = gf256_mul(c, (uint8_t)(i << 4));
    }

    __m128i lo_tbl = _mm_loadu_si128((__m128i *)lo_arr);
    __m128i hi_tbl = _mm_loadu_si128((__m128i *)hi_arr);
    __m128i mask   = _mm_set1_epi8(0x0F);

    size_t i = 0;
    for (; i + 16 <= len; i += 16) {
        __m128i s = _mm_loadu_si128((__m128i *)(src + i));
        __m128i d = _mm_loadu_si128((__m128i *)(dst + i));

        __m128i lo_val = _mm_shuffle_epi8(lo_tbl, _mm_and_si128(s, mask));
        __m128i hi_val = _mm_shuffle_epi8(hi_tbl,
                             _mm_and_si128(_mm_srli_epi16(s, 4), mask));

        _mm_storeu_si128((__m128i *)(dst + i),
                         _mm_xor_si128(d, _mm_xor_si128(lo_val, hi_val)));
    }
    for (; i < len; i++)
        dst[i] ^= gf256_mul(c, src[i]);
}

#else
/* ---- fallback: no SIMD ---- */
void gf256_mul_vec(uint8_t *dst, const uint8_t *src, uint8_t c, size_t len)
{
    gf256_mul_vec_ref(dst, src, c, len);
}
#endif
