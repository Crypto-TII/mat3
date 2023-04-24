#ifndef BGV_RING_H
#define BGV_RING_H

#include <assert.h>
#include <gmp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <x86intrin.h>

#include "random.h"
#include "util.h"

#define BGV_D     32768
#define BGV_QLEN  10
#define BGV_PLEN  4
#define BGV_QPLEN (BGV_QLEN + BGV_PLEN)

typedef int64_t            bgv_Digit;
typedef struct bgv_mod  bgv_Mod;
typedef struct bgv_poly bgv_Poly;

struct bgv_mod {
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

struct bgv_poly {
	int64_t value[BGV_D];
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


#ifdef BGV_RING_IMPL
bgv_Mod   bgv_q[BGV_QPLEN];
bgv_Digit bgv_rns_mods[BGV_QPLEN][BGV_QPLEN];
bgv_Digit bgv_rns_invs[BGV_QPLEN][BGV_QPLEN];

__extension__ typedef          __int128  int128_t;
__extension__ typedef unsigned __int128 uint128_t;

#define BSR64(x) __bsrq(x)

/* https://eprint.iacr.org/2018/039 */
static int64_t barrett_precomp(int64_t mod);
static int64_t barrett_reduce(int64_t a, int64_t mod, int64_t precomp);

/* idx-based functions */
static bgv_Digit digit_add(size_t idx, bgv_Digit a, bgv_Digit b);
static bgv_Digit digit_deinit(size_t idx, bgv_Digit a);
static bgv_Digit digit_init(size_t idx, bgv_Digit a);
static bgv_Digit digit_inv(size_t idx, bgv_Digit a);
static bgv_Digit digit_mul(size_t idx, bgv_Digit a, bgv_Digit b);
static bgv_Digit digit_neg(size_t idx, bgv_Digit a);
static bgv_Digit digit_sub(size_t idx, bgv_Digit a, bgv_Digit b);

/* https://eprint.iacr.org/2018/039 */
static int64_t montgomery_deinit(int64_t a, int64_t mod, int64_t inv);
static int64_t montgomery_init(int64_t a, int64_t mod, int64_t inv, int64_t pow);
static int64_t montgomery_mulred(int64_t a, int64_t b, int64_t mod, int64_t inv);
static int64_t montgomery_powred(int64_t a, size_t exp, int64_t mod, int64_t inv, int64_t pow);
static int64_t montgomery_precomp_inv(int64_t mod);
static int64_t montgomery_precomp_pow(int64_t mod);
static int64_t montgomery_reduce(int128_t a, int64_t mod, int64_t inv);

/* https://eprint.iacr.org/2018/039 */
static void    seiler_intt(int64_t *poly, size_t degree, int64_t mod, int64_t root, int64_t barrett, int64_t inv, int64_t pow);
static void    seiler_ntt(int64_t *poly, size_t degree, int64_t mod, int64_t root, int64_t barrett, int64_t inv, int64_t pow);
static int64_t seiler_precomp(size_t degree, int64_t mod, int64_t barrett, int64_t inv, int64_t pow);

static int64_t
barrett_precomp(int64_t mod)
{
	uint128_t precomp;

	precomp = (uint128_t)1 << (63 + BSR64(mod));
	precomp = (precomp + (mod >> 1)) / mod;

	return -(int64_t)precomp;
}

static int64_t
barrett_reduce(int64_t a, int64_t mod, int64_t precomp)
{
	int128_t tmp128;
	int64_t tmp64;

	/* check for arithmetic shift */
	assert(-1 >> 1 == -1);

	tmp128 = a;
	tmp64 = (tmp128 * precomp) >> 64;
	tmp128 = tmp64 >> (BSR64(mod) - 1);
	tmp64 = (int64_t)tmp128 * mod + a;

	return tmp64 + ((tmp64 >> 63) & mod);
}

bgv_Digit
digit_add(size_t idx, bgv_Digit a, bgv_Digit b)
{

	return barrett_reduce(a + b, bgv_q[idx].value, bgv_q[idx].barrett);
}

bgv_Digit
digit_deinit(size_t idx, bgv_Digit a)
{

	return montgomery_deinit(a, bgv_q[idx].value, bgv_q[idx].inv);
}

bgv_Digit
digit_init(size_t idx, bgv_Digit a)
{

	return montgomery_init(a, bgv_q[idx].value, bgv_q[idx].inv, bgv_q[idx].pow);
}

bgv_Digit
digit_inv(size_t idx, bgv_Digit a)
{
	int64_t inv;

	inv = bgv_util_invmod(digit_deinit(idx, a), bgv_q[idx].value);
	return digit_init(idx, inv);
}

bgv_Digit
digit_mul(size_t idx, bgv_Digit a, bgv_Digit b)
{

	return montgomery_mulred(a, b, bgv_q[idx].value, bgv_q[idx].inv);
}

bgv_Digit
digit_neg(size_t idx, bgv_Digit a)
{

	(void)idx;

	return -a;
}

bgv_Digit
digit_sub(size_t idx, bgv_Digit a, bgv_Digit b)
{

	return barrett_reduce(a - b, bgv_q[idx].value, bgv_q[idx].barrett);
}

static int64_t
montgomery_deinit(int64_t a, int64_t mod, int64_t inv)
{
	int64_t tmp;

	/* check for arithmetic shift */
	assert(-1 >> 1 == -1);

	tmp = montgomery_mulred(a, 1, mod, inv);

	return tmp + ((tmp >> 63) & mod);
}

static int64_t
montgomery_init(int64_t a, int64_t mod, int64_t inv, int64_t pow)
{

	return montgomery_mulred(a, pow, mod, inv);
}

static int64_t
montgomery_mulred(int64_t a, int64_t b, int64_t mod, int64_t inv)
{

	return montgomery_reduce((int128_t)a * b, mod, inv);
}

static int64_t
montgomery_powred(int64_t a, size_t exp, int64_t mod, int64_t inv, int64_t pow)
{
	int64_t ret;

	ret = montgomery_init(1, mod, inv, pow);
	while (exp > 0) {
		if (exp & 1)
			ret = montgomery_mulred(a, ret, mod, inv);
		a = montgomery_mulred(a, a, mod, inv);

		exp >>= 1;
	}

	return ret;
}

static int64_t
montgomery_precomp_inv(int64_t mod)
{
	uint64_t inv;
	int i;

	/* https://crypto.stackexchange.com/questions/47493 */
	inv = mod & 1;
	for (i = 0; i < 6; ++i)
		inv *= 2 - mod * inv;

	return inv;
}

static int64_t
montgomery_precomp_pow(int64_t mod)
{
	int64_t pow;
	int i;

	pow = 1;
	for (i = 0; i < 128; ++i)
		pow <<= 1, pow %= mod;

	return pow;
}

static int64_t
montgomery_reduce(int128_t a, int64_t mod, int64_t inv)
{
	int64_t tmp;

	/* check for arithmetic shift */
	assert(-1 >> 1 == -1);

	tmp = (int128_t)(int64_t)a * inv;
	tmp = (a - (int128_t)tmp * mod) >> 64;

	return tmp;
}

static void
seiler_intt(int64_t *poly, size_t degree, int64_t mod, int64_t root, int64_t barrett, int64_t inv, int64_t pow)
{
	size_t dlog, dinv, j, k, l, s;

	dlog = BSR64(degree);
	dinv = montgomery_init(bgv_util_invmod(degree, mod), mod, inv, pow);

	k = 0;
	for (l = 1; l < degree; l <<= 1) {
		for (s = 0; s < degree; s = j + l) {
			int64_t r;
			size_t exp;

			exp = 1 + bgv_util_bitrev(k++, dlog);
			r = montgomery_powred(root, exp, mod, inv, pow);
			r = montgomery_deinit(r, mod, inv);
			r = montgomery_init(bgv_util_invmod(r, mod), mod, inv, pow);

			for (j = s; j < s + l; ++j) {
				int64_t tmp;

				tmp = poly[j];
				poly[j] = barrett_reduce(tmp + poly[j + l], mod, barrett);
				poly[j + l] = barrett_reduce(tmp - poly[j + l], mod, barrett);
				poly[j + l] = montgomery_mulred(r, poly[j + l], mod, inv);
			}
		}
	}

	for (j = 0; j < degree; ++j)
		poly[j] = montgomery_mulred(poly[j], dinv, mod, inv);
}

static void
seiler_ntt(int64_t *poly, size_t degree, int64_t mod, int64_t root, int64_t barrett, int64_t inv, int64_t pow)
{
	size_t dlog, j, k, l, s;

	dlog = BSR64(degree);

	k = 1;
	for (l = degree >> 1; l > 0; l >>= 1) {
		for (s = 0; s < degree; s = j + l) {
			int64_t r;
			size_t exp;

			exp = bgv_util_bitrev(k++, dlog);
			r = montgomery_powred(root, exp, mod, inv, pow);

			for (j = s; j < s + l; ++j) {
				int64_t tmp;

				tmp = montgomery_mulred(r, poly[j + l], mod, inv);
				poly[j + l] = barrett_reduce(poly[j] - tmp, mod, barrett);
				poly[j] = barrett_reduce(poly[j] + tmp, mod, barrett);
			}
		}
	}
}

static int64_t
seiler_precomp(size_t degree, int64_t mod, int64_t barrett, int64_t inv, int64_t pow)
{
	int64_t root, cnt, one;
	size_t ord, exp;

	ord = mod - 1;
	exp = ord >> (BSR64(degree) + 1);
	one = montgomery_init(1, mod, inv, pow);

	cnt = one;
	for (;;) {
		cnt = barrett_reduce(cnt + one, mod, barrett);
		if (montgomery_powred(cnt, ord, mod, inv, pow) != one)
			continue;

		root = montgomery_powred(cnt, exp, mod, inv, pow);
		if (montgomery_powred(root, degree, mod, inv, pow) != one)
			break;
	}

	return root;
}

void
bgv_mod_deinit(size_t idx)
{

	(void)idx;
}

void
bgv_mod_digits(bgv_Digit *q)
{
	size_t idx, i;

	idx = 0;
	for (i = 0; i < BGV_QLEN; ++i) {
		if (bgv_q[i].drop != 0)
			continue;

		q[idx++] = bgv_q[i].value;
	}
}

void
bgv_mod_drop(size_t idx, unsigned count)
{

	bgv_q[idx].drop = count;
}

void
bgv_mod_init_mpz(size_t idx, mpz_t t)
{
	mpz_t tmp;

	mpz_init(tmp);

	mpz_mod_ui(tmp, t, bgv_rns[idx]);
	bgv_mod_init_u64(idx, mpz_get_ui(tmp));

	mpz_clear(tmp);
}

void
bgv_mod_init_u64(size_t idx, uint64_t t)
{
	int64_t barrett, mod, inv, P, pow;
	size_t i;

	mod = bgv_rns[idx];
	barrett = barrett_precomp(mod);
	inv = montgomery_precomp_inv(mod);
	pow = montgomery_precomp_pow(mod);

	bgv_q[idx].value = mod;
	bgv_q[idx].barrett = barrett;
	bgv_q[idx].inv = inv;
	bgv_q[idx].pow = pow;
	bgv_q[idx].root = seiler_precomp(BGV_D, mod, barrett, inv, pow);
	bgv_q[idx].t = montgomery_init(t, mod, inv, pow);
	bgv_q[idx].idx = idx;
	bgv_q[idx].drop = 0;

	for (i = 0; i < BGV_QPLEN; ++i) {
		bgv_Digit tmp;

		tmp = digit_init(idx, bgv_rns[i]);

		bgv_rns_mods[idx][i] = tmp;
		if (idx == i)
			bgv_rns_invs[idx][i] = digit_init(idx, 1);
		else
			bgv_rns_invs[idx][i] = digit_inv(idx, tmp);
	}

	P = digit_init(idx, 1);
	for (i = BGV_QLEN; i < BGV_QPLEN; ++i)
		P = digit_mul(idx, P, bgv_rns_mods[idx][i]);
	bgv_q[idx].P = P;
}

bgv_Digit
bgv_mod_inv(size_t idx, const size_t *invs, size_t len)
{
	bgv_Digit inv;
	size_t i;

	inv = digit_init(idx, 1);
	for (i = 0; i < len; ++i)
		inv = digit_mul(idx, inv, bgv_rns_invs[idx][invs[i]]);

	return inv;
}

void
bgv_mod_mpz(mpz_t q)
{
	size_t i;

	mpz_set_ui(q, 1);
	for (i = 0; i < BGV_QLEN; ++i) {
		if (bgv_q[i].drop != 0)
			continue;

		mpz_mul_ui(q, q, bgv_q[i].value);
	}
}

bgv_Digit
bgv_mod_negtinv(size_t idx)
{
	bgv_Digit inv;

	inv = digit_inv(idx, bgv_q[idx].t);

	return digit_neg(idx, inv);
}

void
bgv_poly_add(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j)
		rop->value[j] = digit_add(idx, op1->value[j], op2->value[j]);
}

void
bgv_poly_addmul(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2, const bgv_Poly *op3, const bgv_Poly *op4)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j) {
		int64_t add12, add34;

		add12 = digit_add(idx, op1->value[j], op2->value[j]);
		add34 = digit_add(idx, op3->value[j], op4->value[j]);
		rop->value[j] = digit_mul(idx, add12, add34);
	}
}

void
bgv_poly_cmod(size_t idx, bgv_Poly *rop, const bgv_Poly *p)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j) {
		rop->value[j] = p->value[j];
		if (p->value[j] > bgv_q[idx].value / 2)
			rop->value[j] -= bgv_q[idx].value;
		if (p->value[j] < -bgv_q[idx].value / 2)
			rop->value[j] += bgv_q[idx].value;
	}
}

void
bgv_poly_deinit(size_t idx, bgv_Poly *rop, const bgv_Poly *p)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j)
		rop->value[j] = digit_deinit(idx, p->value[j]);
}

void
bgv_poly_ext(size_t idx, bgv_Poly *rop, const bgv_Poly **p, const size_t *base, size_t len)
{
	bgv_Digit inv[BGV_QPLEN], div[BGV_QPLEN];
	size_t i, j;

	for (i = 0; i < len; ++i) {
		size_t ii;

		inv[i] = digit_init(base[i], 1);
		for (ii = 0; ii < len; ++ii) {
			if (bgv_q[base[ii]].drop != 0)
				continue;
			inv[i] = digit_mul(base[i], inv[i], bgv_rns_invs[base[i]][base[ii]]);
		}

		div[i] = bgv_rns_invs[idx][base[i]];
		for (ii = 0; ii < len; ++ii) {
			if (idx == base[i] && idx == base[ii])
				continue;
			if (bgv_q[base[ii]].drop != 0)
				continue;
			div[i] = digit_mul(idx, div[i], bgv_rns_mods[idx][base[ii]]);
		}
	}

	for (j = 0; j < BGV_D; ++j) {
		bgv_Digit sum;

		sum = 0;
		for (i = 0; i < len; ++i) {
			bgv_Digit tmp;

			if (bgv_q[base[i]].drop != 0)
				continue;

			tmp = digit_mul(base[i], p[i]->value[j], inv[i]);
			tmp = digit_deinit(base[i], tmp);
			tmp = digit_init(idx, tmp);
			tmp = digit_mul(idx, tmp, div[i]);
			sum = digit_add(idx, sum, tmp);
		}
		rop->value[j] = sum;
	}
}

void
bgv_poly_init(size_t idx, bgv_Poly *rop, const bgv_Poly *p)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j)
		rop->value[j] = digit_init(idx, p->value[j]);
}

void
bgv_poly_init_error(size_t idx, bgv_Poly *p, const bgv_Poly *e)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j) {
		bgv_Digit init = digit_init(idx, e->value[j]);
		p->value[j] = digit_mul(idx, init, bgv_q[idx].t);
	}
	bgv_poly_ntt(idx, p);
}

void
bgv_poly_intt(size_t idx, bgv_Poly *p)
{

	seiler_intt(p->value, BGV_D, bgv_q[idx].value, bgv_q[idx].root, bgv_q[idx].barrett, bgv_q[idx].inv, bgv_q[idx].pow);
}

void
bgv_poly_mod_mpz(size_t idx, bgv_Poly *rop, const mpz_t *p)
{
	mpz_t tmp;
	size_t j;

	mpz_init(tmp);

	for (j = 0; j < BGV_D; ++j) {
		mpz_mod_ui(tmp, p[j], bgv_q[idx].value);
		rop->value[j] = digit_init(idx, mpz_get_ui(tmp));
	}

	mpz_clear(tmp);
}

void
bgv_poly_mod_u64(size_t idx, bgv_Poly *rop, const uint64_t *p)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j)
		rop->value[j] = digit_init(idx, p[j] % bgv_q[idx].value);
}

void
bgv_poly_mul(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j)
		rop->value[j] = digit_mul(idx, op1->value[j], op2->value[j]);
}

void
bgv_poly_muladd(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2, const bgv_Poly *op3)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j)
		rop->value[j] = digit_add(idx, digit_mul(idx, op1->value[j], op2->value[j]), op3->value[j]);
}

void
bgv_poly_mulc(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Digit op2)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j)
		rop->value[j] = digit_mul(idx, op1->value[j], op2);
}

void
bgv_poly_mulcadd(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Digit op2, const bgv_Poly *op3)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j)
		rop->value[j] = digit_add(idx, digit_mul(idx, op1->value[j], op2), op3->value[j]);
}

void
bgv_poly_neg(size_t idx, bgv_Poly *rop, const bgv_Poly *p)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j)
		rop->value[j] = digit_neg(idx, p->value[j]);
}

void
bgv_poly_ntt(size_t idx, bgv_Poly *p)
{

	seiler_ntt(p->value, BGV_D, bgv_q[idx].value, bgv_q[idx].root, bgv_q[idx].barrett, bgv_q[idx].inv, bgv_q[idx].pow);
}

void
bgv_poly_rot(size_t idx, bgv_Poly *rop, const bgv_Poly *p, size_t steps)
{
	size_t step, i, j;

	assert(rop != p);
	(void)idx;

	steps %= BGV_D / 2;
	if (steps == 0) {
		memcpy(rop, p, sizeof *rop);
		return;
	}

	step = 3;
	for (i = 1; i < steps; ++i) {
		step *= 3;
		step %= 2 * BGV_D;
	}

	i = 0;
	for (j = 0; j < BGV_D; ++j) {
		bgv_Digit tmp = p->value[j];

		if (i / BGV_D & 1)
			tmp = -tmp;
		rop->value[i % BGV_D] = tmp;

		i += step;
		i %= 2 * BGV_D;
	}
}

void
bgv_poly_sample_error(bgv_Poly *e, bgv_Seed *seed)
{

	bgv_sample_cbd21_i64(e->value, BGV_D, seed);
}

void
bgv_poly_sample_secret(bgv_Poly *s, bgv_Seed *seed)
{

	bgv_sample_cbd1_i64(s->value, BGV_D, seed);
}

void
bgv_poly_sample_uniform(size_t idx, bgv_Poly *p, bgv_Seed *seed)
{

	bgv_sample_uniform_i64(p->value, BGV_D, bgv_q[idx].value, seed);
}

void
bgv_poly_sub(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j)
		rop->value[j] = digit_sub(idx, op1->value[j], op2->value[j]);
}

#endif /* BGV_RING_IMPL */
