/*
 * GF(256) 갈루아 필드 연산
 * 원시 다항식: x^8 + x^4 + x^3 + x^2 + 1 (0x11D)
 * 이레이저 코딩의 수학적 기반
 */
#include "gf256.h"

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
