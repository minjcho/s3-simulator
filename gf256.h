#ifndef GF256_H
#define GF256_H

#include <stdint.h>

void    gf256_init(void);
uint8_t gf256_mul(uint8_t a, uint8_t b);
uint8_t gf256_inv(uint8_t a);

#endif
