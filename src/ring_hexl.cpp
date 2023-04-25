#include <gmp.h>
extern "C" {
	#include "ring_hexl.h"
}

#include "hexl/hexl.hpp"

tiimat3_Mod tiimat3_q[TIIMAT3_QPLEN];
const tiimat3_Digit tiimat3_rns[TIIMAT3_QPLEN] = {
	/* tiimat3_q */
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
	/* tiimat3_P */
	576460752321642497,
	576460752325705729,
	576460752328130561,
	576460752328327169
};

void
tiimat3_mod_deinit(size_t idx)
{

	delete (intel::hexl::NTT *)tiimat3_q[idx].ntt;
}

void
tiimat3_mod_digits(tiimat3_Digit *q)
{
	size_t idx, i;

	idx = 0;
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		if (tiimat3_q[i].drop != 0)
			continue;

		q[idx++] = tiimat3_rns[i];
	}
}

void
tiimat3_mod_drop(size_t idx, unsigned count)
{

	tiimat3_q[idx].drop = count;
}

void
tiimat3_mod_init_mpz(size_t idx, mpz_t t)
{
	mpz_t tmp;

	mpz_init(tmp);

	mpz_mod_ui(tmp, t, tiimat3_rns[idx]);
	tiimat3_mod_init_u64(idx, mpz_get_ui(tmp));

	mpz_clear(tmp);
}

void
tiimat3_mod_init_u64(size_t idx, uint64_t t)
{
	int64_t mod, P;
	size_t i;

	mod = tiimat3_rns[idx];
	tiimat3_q[idx].value = mod;
	tiimat3_q[idx].t = t;

	P = 1;
	for (i = TIIMAT3_QLEN; i < TIIMAT3_QPLEN; ++i)
		P = intel::hexl::MultiplyMod(P, tiimat3_rns[i], mod);
	tiimat3_q[idx].P = P;

	tiimat3_q[idx].ntt = (void *)(new intel::hexl::NTT(TIIMAT3_D, mod));
	tiimat3_q[idx].idx = idx;
	tiimat3_q[idx].drop = 0;
}

tiimat3_Digit
tiimat3_mod_inv(size_t idx, const size_t *invs, size_t len)
{
	tiimat3_Digit inv;
	size_t i;

	inv = 1;
	for (i = 0; i < len; ++i) {
		tiimat3_Digit tmp;

		if (tiimat3_q[invs[i]].drop != 0)
			continue;

		tmp = intel::hexl::InverseMod(tiimat3_rns[invs[i]], tiimat3_rns[idx]);
		inv = intel::hexl::MultiplyMod(inv, tmp, tiimat3_rns[idx]);
	}

	return inv;
}

void
tiimat3_mod_mpz(mpz_t q)
{
	size_t i;

	mpz_set_ui(q, 1);
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		if (tiimat3_q[i].drop != 0)
			continue;

		mpz_mul_ui(q, q, tiimat3_rns[i]);
	}
}

tiimat3_Digit
tiimat3_mod_negtinv(size_t idx)
{
	tiimat3_Digit inv;

	inv = intel::hexl::InverseMod(tiimat3_q[idx].t, tiimat3_rns[idx]);

	return tiimat3_rns[idx] - inv;
}

void
tiimat3_poly_add(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2)
{

	intel::hexl::EltwiseAddMod(rop->value, op1->value, op2->value, TIIMAT3_D, tiimat3_rns[idx]);
}

void
tiimat3_poly_addmul(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2, const tiimat3_Poly *op3, const tiimat3_Poly *op4)
{
	tiimat3_Poly tmp;

	intel::hexl::EltwiseAddMod(tmp.value, op1->value, op2->value, TIIMAT3_D, tiimat3_rns[idx]);
	intel::hexl::EltwiseAddMod(rop->value, op3->value, op4->value, TIIMAT3_D, tiimat3_rns[idx]);
	intel::hexl::EltwiseMultMod(rop->value, rop->value, tmp.value, TIIMAT3_D, tiimat3_rns[idx], 1);
}

void
tiimat3_poly_cmod(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p)
{

	(void)idx;
	memcpy(rop, p, sizeof *rop);
}

void
tiimat3_poly_deinit(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p)
{

	(void)idx;
	memcpy(rop, p, sizeof *rop);
}

void
tiimat3_poly_ext(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly **p, const size_t *base, size_t len)
{
	tiimat3_Digit inv[TIIMAT3_QPLEN], div[TIIMAT3_QPLEN];
	size_t i;

	for (i = 0; i < len; ++i) {
		size_t ii;

		inv[i] = tiimat3_mod_inv(base[i], base, len);
		div[i] = intel::hexl::InverseMod(tiimat3_rns[base[i]], tiimat3_rns[idx]);
		for (ii = 0; ii < len; ++ii) {
			if (idx == base[i] && idx == base[ii])
				continue;

			if (tiimat3_q[base[ii]].drop != 0)
				continue;

			div[i] = intel::hexl::MultiplyMod(div[i], tiimat3_rns[base[ii]], tiimat3_rns[idx]);
		}
	}

	memset(rop, 0, sizeof *rop);
	for (i = 0; i < len; ++i) {
		tiimat3_Poly tmp;

		if (tiimat3_q[base[i]].drop != 0)
			continue;

		tiimat3_poly_mulc(base[i], &tmp, p[i], inv[i]);
		tiimat3_poly_mulcadd(idx, rop, &tmp, div[i], rop);
	}
}

void
tiimat3_poly_init(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p)
{
	size_t j;

	/* check for arithmetic shift */
	assert(-1 >> 1 == -1);

	for (j = 0; j < TIIMAT3_D; ++j) {
		int64_t val = p->value[j];
		rop->value[j] = val + ((val >> 63) & tiimat3_rns[idx]);
	}
}

void
tiimat3_poly_init_error(size_t idx, tiimat3_Poly *p, const tiimat3_Poly *e)
{

	tiimat3_poly_init(idx, p, e);
	tiimat3_poly_mulc(idx, p, p, tiimat3_q[idx].t);
	tiimat3_poly_ntt(idx, p);
}

void
tiimat3_poly_intt(size_t idx, tiimat3_Poly *p)
{
	intel::hexl::NTT *ntt;

	ntt = (intel::hexl::NTT *)tiimat3_q[idx].ntt;
	ntt->ComputeInverse(p->value, p->value, 1, 1);
}

void
tiimat3_poly_mod_mpz(size_t idx, tiimat3_Poly *rop, const mpz_t *p)
{
	mpz_t tmp;
	size_t j;

	mpz_init(tmp);

	for (j = 0; j < TIIMAT3_D; ++j) {
		mpz_mod_ui(tmp, p[j], tiimat3_rns[idx]);
		rop->value[j] = mpz_get_ui(tmp);
	}

	mpz_clear(tmp);
}

void
tiimat3_poly_mod_u64(size_t idx, tiimat3_Poly *rop, const uint64_t *p)
{

	intel::hexl::EltwiseReduceMod(rop->value, p, TIIMAT3_D, tiimat3_rns[idx], 1, 1);
}

void
tiimat3_poly_mul(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2)
{

	intel::hexl::EltwiseMultMod(rop->value, op1->value, op2->value, TIIMAT3_D, tiimat3_rns[idx], 1);
}

void
tiimat3_poly_muladd(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2, const tiimat3_Poly *op3)
{
	tiimat3_Poly tmp;

	intel::hexl::EltwiseMultMod(tmp.value, op1->value, op2->value, TIIMAT3_D, tiimat3_rns[idx], 1);
	intel::hexl::EltwiseAddMod(rop->value, tmp.value, op3->value, TIIMAT3_D, tiimat3_rns[idx]);
}

void
tiimat3_poly_mulc(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Digit op2)
{

	intel::hexl::EltwiseFMAMod(rop->value, op1->value, op2, 0, TIIMAT3_D, tiimat3_rns[idx], 1);
}

void
tiimat3_poly_mulcadd(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Digit op2, const tiimat3_Poly *op3)
{

	intel::hexl::EltwiseFMAMod(rop->value, op1->value, op2, op3->value, TIIMAT3_D, tiimat3_rns[idx], 1);
}

void
tiimat3_poly_neg(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		rop->value[j] = tiimat3_rns[idx] - p->value[j];
}

void
tiimat3_poly_ntt(size_t idx, tiimat3_Poly *p)
{
	intel::hexl::NTT *ntt;

	ntt = (intel::hexl::NTT *)tiimat3_q[idx].ntt;
	ntt->ComputeForward(p->value, p->value, 1, 1);
}

void
tiimat3_poly_rot(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p, size_t steps)
{
	size_t step, i, j;

	assert(rop != p);

	steps %= TIIMAT3_D / 2;
	if (steps == 0) {
		memcpy(rop, p, sizeof *rop);
		return;
	}

	step = 3;
	for (i = 1; i < steps; ++i) {
		step *= 3;
		step %= 2 * TIIMAT3_D;
	}

	i = 0;
	for (j = 0; j < TIIMAT3_D; ++j) {
		tiimat3_Digit tmp = p->value[j];

		if (i / TIIMAT3_D & 1)
			tmp = tiimat3_rns[idx] - tmp;
		rop->value[i % TIIMAT3_D] = tmp;

		i += step;
		i %= 2 * TIIMAT3_D;
	}
}

void
tiimat3_poly_sample_error(tiimat3_Poly *e, tiimat3_Seed *seed)
{

	tiimat3_random_cbd21_i64((int64_t *)e->value, TIIMAT3_D, seed);
}

void
tiimat3_poly_sample_secret(tiimat3_Poly *s, tiimat3_Seed *seed)
{

	tiimat3_random_cbd1_i64((int64_t *)s->value, TIIMAT3_D, seed);
}

void
tiimat3_poly_sample_uniform(size_t idx, tiimat3_Poly *p, tiimat3_Seed *seed)
{

	tiimat3_random_uniform_u64(p->value, TIIMAT3_D, tiimat3_rns[idx], seed);
}

void
tiimat3_poly_sub(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2)
{

	intel::hexl::EltwiseSubMod(rop->value, op1->value, op2->value, TIIMAT3_D, tiimat3_rns[idx]);
}
