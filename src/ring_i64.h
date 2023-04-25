#ifndef TIIMAT3_RING_H
#define TIIMAT3_RING_H

#include <assert.h>
#include <gmp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>

#include "util.h"
#include "random.h"

#define TIIMAT3_D     32768
#define TIIMAT3_QLEN  10
#define TIIMAT3_PLEN  4
#define TIIMAT3_QPLEN (TIIMAT3_QLEN + TIIMAT3_PLEN)

typedef int64_t             tiimat3_Digit;
typedef struct tiimat3_mod  tiimat3_Mod;
typedef struct tiimat3_poly tiimat3_Poly;

struct tiimat3_mod {
	int64_t value;
	int64_t barrett;
	int64_t inv;
	int64_t pow;
	int64_t root;
	int64_t t;
	int64_t P;
	size_t idx;
	unsigned drop;
};

struct tiimat3_poly {
	int64_t value[TIIMAT3_D];
};

extern tiimat3_Mod tiimat3_q[TIIMAT3_QPLEN];
extern const tiimat3_Digit tiimat3_rns[TIIMAT3_QPLEN];

void          tiimat3_mod_deinit(size_t idx);
void          tiimat3_mod_digits(tiimat3_Digit *q);
void          tiimat3_mod_drop(size_t idx, unsigned count);
void          tiimat3_mod_init_mpz(size_t idx, mpz_t t);
void          tiimat3_mod_init_u64(size_t idx, uint64_t t);
tiimat3_Digit tiimat3_mod_inv(size_t idx, const size_t *invs, size_t len);
void          tiimat3_mod_mpz(mpz_t q);
tiimat3_Digit tiimat3_mod_negtinv(size_t idx);

void tiimat3_poly_add(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2);
void tiimat3_poly_addmul(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2, const tiimat3_Poly *op3, const tiimat3_Poly *op4);
void tiimat3_poly_cmod(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p);
void tiimat3_poly_deinit(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p);
void tiimat3_poly_ext(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly **p, const size_t *base, size_t len);
void tiimat3_poly_init(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p);
void tiimat3_poly_init_error(size_t idx, tiimat3_Poly *p, const tiimat3_Poly *e);
void tiimat3_poly_intt(size_t idx, tiimat3_Poly *p);
void tiimat3_poly_mod_mpz(size_t idx, tiimat3_Poly *rop, const mpz_t *p);
void tiimat3_poly_mod_u64(size_t idx, tiimat3_Poly *rop, const uint64_t *p);
void tiimat3_poly_mul(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2);
void tiimat3_poly_muladd(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2, const tiimat3_Poly *op3);
void tiimat3_poly_mulc(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Digit op2);
void tiimat3_poly_mulcadd(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Digit op2, const tiimat3_Poly *op3);
void tiimat3_poly_neg(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p);
void tiimat3_poly_ntt(size_t idx, tiimat3_Poly *p);
void tiimat3_poly_rot(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p, size_t steps);
void tiimat3_poly_sample_error(tiimat3_Poly *e, tiimat3_Seed *seed);
void tiimat3_poly_sample_secret(tiimat3_Poly *s, tiimat3_Seed *seed);
void tiimat3_poly_sample_uniform(size_t idx, tiimat3_Poly *p, tiimat3_Seed *seed);
void tiimat3_poly_sub(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2);

#endif /* TIIMAT3_RING_H */
