#ifndef TIIMAT3_UTIL_H
#define TIIMAT3_UTIL_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

void *tiimat3_util_alloc(size_t len, size_t size);
void  tiimat3_util_dealloc(void *ptr);

/* https://en.wikipedia.org/wiki/Hacker%27s_Delight */
uint64_t tiimat3_util_bitrev(uint64_t a, size_t log2);
uint64_t tiimat3_util_invmod(uint64_t a, uint64_t mod);
size_t   tiimat3_util_msb64(uint64_t a);
uint32_t tiimat3_util_popcount32(uint32_t a);
uint64_t tiimat3_util_rol64(uint64_t a, size_t r);

#endif /* TIIMAT3_UTIL_H */
