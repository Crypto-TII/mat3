#include <gmp.h>

extern "C" {
	#include "ring_hexl.h"
}

#include "hexl/hexl.hpp"

bgv_Mod bgv_q[BGV_QPLEN];

void
bgv_mod_deinit(size_t idx)
{

	delete (intel::hexl::NTT *)bgv_q[idx].ntt;
}

void
bgv_mod_digits(bgv_Digit *q)
{
	size_t idx, i;

	idx = 0;
	for (i = 0; i < BGV_QLEN; ++i) {
		if (bgv_q[i].drop != 0)
			continue;

		q[idx++] = bgv_rns[i];
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
	int64_t mod, P;
	size_t i;

	mod = bgv_rns[idx];
	bgv_q[idx].value = mod;
	bgv_q[idx].t = t;

	P = 1;
	for (i = BGV_QLEN; i < BGV_QPLEN; ++i)
		P = intel::hexl::MultiplyMod(P, bgv_rns[i], mod);
	bgv_q[idx].P = P;

	bgv_q[idx].ntt = (void *)(new intel::hexl::NTT(BGV_D, mod));
	bgv_q[idx].idx = idx;
	bgv_q[idx].drop = 0;
}

bgv_Digit
bgv_mod_inv(size_t idx, const size_t *invs, size_t len)
{
	bgv_Digit inv;
	size_t i;

	inv = 1;
	for (i = 0; i < len; ++i) {
		bgv_Digit tmp;

		if (bgv_q[invs[i]].drop != 0)
			continue;

		tmp = intel::hexl::InverseMod(bgv_rns[invs[i]], bgv_rns[idx]);
		inv = intel::hexl::MultiplyMod(inv, tmp, bgv_rns[idx]);
	}

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

		mpz_mul_ui(q, q, bgv_rns[i]);
	}
}

bgv_Digit
bgv_mod_negtinv(size_t idx)
{
	bgv_Digit inv;

	inv = intel::hexl::InverseMod(bgv_q[idx].t, bgv_rns[idx]);

	return bgv_rns[idx] - inv;
}

void
bgv_poly_add(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2)
{

	intel::hexl::EltwiseAddMod(rop->value, op1->value, op2->value, BGV_D, bgv_rns[idx]);
}

void
bgv_poly_addmul(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2, const bgv_Poly *op3, const bgv_Poly *op4)
{
	bgv_Poly tmp;

	intel::hexl::EltwiseAddMod(tmp.value, op1->value, op2->value, BGV_D, bgv_rns[idx]);
	intel::hexl::EltwiseAddMod(rop->value, op3->value, op4->value, BGV_D, bgv_rns[idx]);
	intel::hexl::EltwiseMultMod(rop->value, rop->value, tmp.value, BGV_D, bgv_rns[idx], 1);
}

void
bgv_poly_cmod(size_t idx, bgv_Poly *rop, const bgv_Poly *p)
{

	(void)idx;
	memcpy(rop, p, sizeof *rop);
}

void
bgv_poly_deinit(size_t idx, bgv_Poly *rop, const bgv_Poly *p)
{

	(void)idx;
	memcpy(rop, p, sizeof *rop);
}

void
bgv_poly_ext(size_t idx, bgv_Poly *rop, const bgv_Poly **p, const size_t *base, size_t len)
{
	bgv_Digit inv[BGV_QPLEN], div[BGV_QPLEN];
	size_t i;

	for (i = 0; i < len; ++i) {
		size_t ii;

		inv[i] = bgv_mod_inv(base[i], base, len);
		div[i] = intel::hexl::InverseMod(bgv_rns[base[i]], bgv_rns[idx]);
		for (ii = 0; ii < len; ++ii) {
			if (idx == base[i] && idx == base[ii])
				continue;

			if (bgv_q[base[ii]].drop != 0)
				continue;

			div[i] = intel::hexl::MultiplyMod(div[i], bgv_rns[base[ii]], bgv_rns[idx]);
		}
	}

	memset(rop, 0, sizeof *rop);
	for (i = 0; i < len; ++i) {
		bgv_Poly tmp;

		if (bgv_q[base[i]].drop != 0)
			continue;

		bgv_poly_mulc(base[i], &tmp, p[i], inv[i]);
		bgv_poly_mulcadd(idx, rop, &tmp, div[i], rop);
	}
}

void
bgv_poly_init(size_t idx, bgv_Poly *rop, const bgv_Poly *p)
{
	size_t j;

	/* check for arithmetic shift */
	assert(-1 >> 1 == -1);

	for (j = 0; j < BGV_D; ++j) {
		int64_t val = p->value[j];
		rop->value[j] = val + ((val >> 63) & bgv_rns[idx]);
	}
}

void
bgv_poly_init_error(size_t idx, bgv_Poly *p, const bgv_Poly *e)
{

	bgv_poly_init(idx, p, e);
	bgv_poly_mulc(idx, p, p, bgv_q[idx].t);
	bgv_poly_ntt(idx, p);
}

void
bgv_poly_intt(size_t idx, bgv_Poly *p)
{
	intel::hexl::NTT *ntt;

	ntt = (intel::hexl::NTT *)bgv_q[idx].ntt;
	ntt->ComputeInverse(p->value, p->value, 1, 1);
}

void
bgv_poly_mod_mpz(size_t idx, bgv_Poly *rop, const mpz_t *p)
{
	mpz_t tmp;
	size_t j;

	mpz_init(tmp);

	for (j = 0; j < BGV_D; ++j) {
		mpz_mod_ui(tmp, p[j], bgv_rns[idx]);
		rop->value[j] = mpz_get_ui(tmp);
	}

	mpz_clear(tmp);
}

void
bgv_poly_mod_u64(size_t idx, bgv_Poly *rop, const uint64_t *p)
{

	intel::hexl::EltwiseReduceMod(rop->value, p, BGV_D, bgv_rns[idx], 1, 1);
}

void
bgv_poly_mul(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2)
{

	intel::hexl::EltwiseMultMod(rop->value, op1->value, op2->value, BGV_D, bgv_rns[idx], 1);
}

void
bgv_poly_muladd(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2, const bgv_Poly *op3)
{
	bgv_Poly tmp;

	intel::hexl::EltwiseMultMod(tmp.value, op1->value, op2->value, BGV_D, bgv_rns[idx], 1);
	intel::hexl::EltwiseAddMod(rop->value, tmp.value, op3->value, BGV_D, bgv_rns[idx]);
}

void
bgv_poly_mulc(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Digit op2)
{

	intel::hexl::EltwiseFMAMod(rop->value, op1->value, op2, 0, BGV_D, bgv_rns[idx], 1);
}

void
bgv_poly_mulcadd(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Digit op2, const bgv_Poly *op3)
{

	intel::hexl::EltwiseFMAMod(rop->value, op1->value, op2, op3->value, BGV_D, bgv_rns[idx], 1);
}

void
bgv_poly_neg(size_t idx, bgv_Poly *rop, const bgv_Poly *p)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j)
		rop->value[j] = bgv_rns[idx] - p->value[j];
}

void
bgv_poly_ntt(size_t idx, bgv_Poly *p)
{
	intel::hexl::NTT *ntt;

	ntt = (intel::hexl::NTT *)bgv_q[idx].ntt;
	ntt->ComputeForward(p->value, p->value, 1, 1);
}

void
bgv_poly_rot(size_t idx, bgv_Poly *rop, const bgv_Poly *p, size_t steps)
{
	size_t step, i, j;

	assert(rop != p);

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
			tmp = bgv_rns[idx] - tmp;
		rop->value[i % BGV_D] = tmp;

		i += step;
		i %= 2 * BGV_D;
	}
}

void
bgv_poly_sample_error(bgv_Poly *e, tiimat3_Seed *seed)
{

	tiimat3_random_cbd21_i64((int64_t *)e->value, BGV_D, seed);
}

void
bgv_poly_sample_secret(bgv_Poly *s, tiimat3_Seed *seed)
{

	tiimat3_random_cbd1_i64((int64_t *)s->value, BGV_D, seed);
}

void
bgv_poly_sample_uniform(size_t idx, bgv_Poly *p, tiimat3_Seed *seed)
{

	tiimat3_random_uniform_u64(p->value, BGV_D, bgv_rns[idx], seed);
}

void
bgv_poly_sub(size_t idx, bgv_Poly *rop, const bgv_Poly *op1, const bgv_Poly *op2)
{

	intel::hexl::EltwiseSubMod(rop->value, op1->value, op2->value, BGV_D, bgv_rns[idx]);
}
