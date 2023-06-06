#ifndef TIIMAT3_RANDOM_H
#define TIIMAT3_RANDOM_H

#include <openssl/rand.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "util.h"

typedef struct tiimat3_seed tiimat3_Seed;
struct tiimat3_seed {
	uint64_t state[25];
	unsigned char value[16];
	size_t pos;
	int init;
};

void tiimat3_random_bytes(unsigned char *out, size_t len, tiimat3_Seed *seed);
void tiimat3_random_cbd1_i8(int8_t *out, size_t len, tiimat3_Seed *seed);
void tiimat3_random_cbd1_i64(int64_t *out, size_t len, tiimat3_Seed *seed);
void tiimat3_random_cbd21_i8(int8_t *out, size_t len, tiimat3_Seed *seed);
void tiimat3_random_cbd21_i64(int64_t *out, size_t len, tiimat3_Seed *seed);
void tiimat3_random_seed(tiimat3_Seed *seed);
void tiimat3_random_uniform_i64(int64_t *out, size_t len, uint64_t mod, tiimat3_Seed *seed);
void tiimat3_random_uniform_u64(uint64_t *out, size_t len, uint64_t mod, tiimat3_Seed *seed);

#endif /* TIIMAT3_RANDOM_H */
