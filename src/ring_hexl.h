#ifndef BGV_RING_H
#define BGV_RING_H

#include <assert.h>
#include <gmp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "random.h"
#include "util.h"

#define BGV_D     32768
#define BGV_QLEN  10
#define BGV_PLEN  4
#define BGV_QPLEN (BGV_QLEN + BGV_PLEN)

typedef uint64_t           bgv_Digit;
typedef struct bgv_mod  bgv_Mod;
typedef struct bgv_poly bgv_Poly;

struct bgv_mod {
	uint64_t value;
	uint64_t t;
	uint64_t P;
	void *ntt;
	size_t idx;
	unsigned drop;
};

struct bgv_poly {
	uint64_t value[BGV_D];
};

#ifdef BGV_RING_IMPL
const bgv_Digit bgv_rns[BGV_QPLEN] = {
	/* bgv_q */
	576460752308273153,
	576460752312401921,
	576460752313712641,
	576460752314368001,
	576460752315482113,
	576460752315678721,
	576460752318824449,
	576460752319021057,
	576460752319414273,
	576460752320790529,
	/* bgv_P */
	576460752321642497,
	576460752325705729,
	576460752328130561,
	576460752328327169
};
#else
extern bgv_Digit bgv_rns[BGV_QPLEN];
#endif

extern bgv_Mod bgv_q[BGV_QPLEN];

void         bgv_mod_deinit(size_t idx);
void         bgv_mod_digits(bgv_Digit *q);
void         bgv_mod_drop(size_t idx, unsigned count);
void         bgv_mod_init_mpz(size_t idx, mpz_t t);
void         bgv_mod_init_u64(size_t idx, uint64_t t);
bgv_Digit bgv_mod_inv(size_t idx, const size_t *invs, size_t len);
void         bgv_mod_mpz(mpz_t q);
bgv_Digit bgv_mod_negtinv(size_t idx);

void bgv_poly_add(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2);
void bgv_poly_addmul(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2, const bgv_Poly *op3, const bgv_Poly *op4);
void bgv_poly_cmod(size_t idx, bgv_Poly *rop, const bgv_Poly *p);
void bgv_poly_deinit(size_t idx, bgv_Poly *rop, const bgv_Poly *p);
void bgv_poly_ext(size_t idx, bgv_Poly *rop, const bgv_Poly **p, const size_t *base, size_t len);
void bgv_poly_init(size_t idx, bgv_Poly *rop, const bgv_Poly *p);
void bgv_poly_init_error(size_t idx, bgv_Poly *p, const bgv_Poly *e);
void bgv_poly_intt(size_t idx, bgv_Poly *p);
void bgv_poly_mod_mpz(size_t idx, bgv_Poly *rop, const mpz_t *p);
void bgv_poly_mod_u64(size_t idx, bgv_Poly *rop, const uint64_t *p);
void bgv_poly_mul(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2);
void bgv_poly_muladd(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2, const bgv_Poly *op3);
void bgv_poly_mulc(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Digit op2);
void bgv_poly_mulcadd(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Digit op2, const bgv_Poly *op3);
void bgv_poly_neg(size_t idx, bgv_Poly *rop, const bgv_Poly *p);
void bgv_poly_ntt(size_t idx, bgv_Poly *p);
void bgv_poly_rot(size_t idx, bgv_Poly *rop, const bgv_Poly *p, size_t steps);
void bgv_poly_sample_error(bgv_Poly *e, bgv_Seed *seed);
void bgv_poly_sample_secret(bgv_Poly *s, bgv_Seed *seed);
void bgv_poly_sample_uniform(size_t idx, bgv_Poly *p, bgv_Seed *seed);
void bgv_poly_sub(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2);

#endif /* BGV_RING_H */
