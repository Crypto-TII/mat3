#ifndef BGV_BGV_H
#define BGV_BGV_H

#define BGV_OMEGA    3
#define BGV_KMAX     4

#if BGV_USE_HEXL
	#include "ring_hexl.h"
	#define bgv_mod_init bgv_mod_init_mpz
	#define bgv_msg_crt  bgv_msg_crt_u64
	#define bgv_poly_mod bgv_poly_mod_mpz
#else
	#include "ring_i64.h"
	#define bgv_mod_init bgv_mod_init_mpz
	#define bgv_msg_crt  bgv_msg_crt_i64
	#define bgv_poly_mod bgv_poly_mod_mpz
#endif /* BGV_USE_HEXL */

#include "msg_mpz.h"

typedef struct bgv_ciphertext        bgv_Ciphertext;
typedef struct bgv_ciphertext_switch bgv_CiphertextSwitch;
typedef struct bgv_delta             bgv_Delta;
typedef struct bgv_key_public        bgv_KeyPublic;
typedef struct bgv_key_secret        bgv_KeySecret;
typedef struct bgv_key_switch        bgv_KeySwitch;

struct bgv_ciphertext {
	bgv_Poly poly[3];
	size_t degree;
	int ntt;
};

struct bgv_ciphertext_switch {
	bgv_Ciphertext ct;
	bgv_Poly ext[BGV_OMEGA];
	size_t poly;
};

struct bgv_delta {
	bgv_Ciphertext ct;
	bgv_Digit inv;
	size_t idx;
};

struct bgv_key_public {
	bgv_Seed seed;
	bgv_Poly poly;
};

struct bgv_key_secret {
	bgv_Poly poly;
};

struct bgv_key_switch {
	bgv_Seed seed;
	bgv_Poly poly[BGV_OMEGA];
};

extern size_t bgv_chunk_idx[BGV_OMEGA][BGV_KMAX];
extern size_t bgv_chunk_len[BGV_OMEGA];

void  bgv_add(size_t idx, bgv_Ciphertext *rop, bgv_Ciphertext *op1, bgv_Ciphertext *op2);
void *bgv_alloc(size_t len, size_t size);
void  bgv_dealloc(void *ptr);
void  bgv_decode(bgv_Message *m, const bgv_Poly *p);
void  bgv_decrypt(size_t idx, bgv_Poly *p, const bgv_KeySecret *sk, bgv_Ciphertext *ct);
void  bgv_deinit(void);
void  bgv_encode(size_t idx, bgv_Poly *p, const bgv_Message *m);
void  bgv_encrypt(size_t idx, bgv_Ciphertext *ct, bgv_KeyPublic *pk, const bgv_Poly *p, bgv_Seed *seed);
void  bgv_init(void);
void  bgv_keygen_public(size_t idx, bgv_KeyPublic *pk, const bgv_KeySecret *sk, bgv_Seed *seed);
void  bgv_keygen_secret(bgv_KeySecret *sk);
void  bgv_keygen_switch(size_t idx, bgv_KeySwitch *ksw, const bgv_KeySecret *sk, const bgv_KeySecret *sk2, bgv_Seed *seed);
void  bgv_keygen_switch2(size_t idx, bgv_KeySwitch *ksw, const bgv_KeySecret *sk, bgv_Seed *seed);
void  bgv_keygen_switchr(size_t idx, bgv_KeySwitch *ksw, const bgv_KeySecret *sk, size_t steps, bgv_Seed *seed);
void  bgv_keyswitch(bgv_Ciphertext *ct, bgv_KeySwitch *ksw, size_t dsw);
void  bgv_keyswitch_delta(size_t idx, bgv_Delta *delta, bgv_Ciphertext *ct);
void  bgv_keyswitch_dot(size_t idx, bgv_Ciphertext *ct, bgv_CiphertextSwitch *csw, bgv_KeySwitch *ksw);
void  bgv_keyswitch_ext(bgv_CiphertextSwitch *csw, bgv_Ciphertext *ct, size_t poly);
void  bgv_keyswitch_switch(size_t idx, bgv_Ciphertext *ct, bgv_Delta *delta);
void  bgv_modswitch(size_t idx, bgv_Ciphertext *ct, bgv_Delta *delta);
void  bgv_modswitch_delta(size_t idx, bgv_Delta *delta, bgv_Ciphertext *ct);
void  bgv_modswitch_ext(size_t idx, bgv_Ciphertext *ct, bgv_Delta *delta);
void  bgv_mulc(size_t idx, bgv_Ciphertext *rop, bgv_Ciphertext *ct, const bgv_Poly *p);
void  bgv_mul(size_t idx, bgv_Ciphertext *rop, bgv_Ciphertext *op1, bgv_Ciphertext *op2);
void  bgv_rot(size_t idx, bgv_Ciphertext *rop, bgv_Ciphertext *ct, size_t steps);
void  bgv_rot_csw(size_t idx, bgv_CiphertextSwitch *rop, bgv_CiphertextSwitch *csw, size_t steps);
void  bgv_rot_csw_inplace(size_t idx, bgv_CiphertextSwitch *csw, size_t steps);
void  bgv_rot_inplace(size_t idx, bgv_Ciphertext *ct, size_t steps);

#endif /* BGV_BGV_H */


#ifdef BGV_BGV_IMPL
size_t bgv_chunk_idx[BGV_OMEGA][BGV_KMAX];
size_t bgv_chunk_len[BGV_OMEGA];

#define MAX(a, b)   ((a) > (b) ? (a) : (b))

static void ciphertext_intt(size_t idx, bgv_Ciphertext *ct);
static void ciphertext_ntt(size_t idx, bgv_Ciphertext *ct);
static void ciphertext_switch_intt(size_t idx, bgv_CiphertextSwitch *csw);
static void ciphertext_switch_ntt(size_t idx, bgv_CiphertextSwitch *csw);
static void delta_ext(size_t idx, bgv_Delta *rop, bgv_Delta *delta, const size_t *base, size_t len);

static void
ciphertext_intt(size_t idx, bgv_Ciphertext *ct)
{

	if (ct->ntt == 0)
		return;
	ct->ntt = 0;

	bgv_poly_intt(idx, &ct->poly[0]);
	bgv_poly_intt(idx, &ct->poly[1]);

	if (ct->degree == 2)
		return;

	bgv_poly_intt(idx, &ct->poly[2]);
}

static void
ciphertext_ntt(size_t idx, bgv_Ciphertext *ct)
{

	if (ct->ntt == 1)
		return;
	ct->ntt = 1;

	bgv_poly_ntt(idx, &ct->poly[0]);
	bgv_poly_ntt(idx, &ct->poly[1]);

	if (ct->degree == 2)
		return;

	bgv_poly_ntt(idx, &ct->poly[2]);
}

static void
ciphertext_switch_intt(size_t idx, bgv_CiphertextSwitch *csw)
{
	size_t k;

	if (csw->ct.ntt == 0)
		return;

	ciphertext_intt(idx, &csw->ct);
	for (k = 0; k < BGV_OMEGA; ++k)
		bgv_poly_intt(idx, &csw->ext[k]);
}

static void
ciphertext_switch_ntt(size_t idx, bgv_CiphertextSwitch *csw)
{
	size_t k;

	if (csw->ct.ntt == 1)
		return;

	ciphertext_ntt(idx, &csw->ct);
	for (k = 0; k < BGV_OMEGA; ++k)
		bgv_poly_ntt(idx, &csw->ext[k]);
}

void
delta_ext(size_t idx, bgv_Delta *rop, bgv_Delta *delta, const size_t *base, size_t len)
{
	const bgv_Poly *p[BGV_QPLEN];
	size_t i;

	ciphertext_intt(idx, &delta->ct);

	rop->ct.degree = delta->ct.degree;
	rop->ct.ntt = 0;
	rop->inv = bgv_mod_inv(idx, base, len);

	for (i = 0; i < len; ++i)
		p[i] = &delta[i].ct.poly[0];
	bgv_poly_ext(idx, &rop->ct.poly[0], p, base, len);

	for (i = 0; i < len; ++i)
		p[i] = &delta[i].ct.poly[1];
	bgv_poly_ext(idx, &rop->ct.poly[1], p, base, len);

	if (rop->ct.degree == 2)
		return;

	for (i = 0; i < len; ++i)
		p[i] = &delta[i].ct.poly[2];
	bgv_poly_ext(idx, &rop->ct.poly[2], p, base, len);
}

void
bgv_add(size_t idx, bgv_Ciphertext *rop, bgv_Ciphertext *op1, bgv_Ciphertext *op2)
{

	if (op1->ntt != op2->ntt) {
		ciphertext_ntt(idx, op1);
		ciphertext_ntt(idx, op2);
	}

	rop->degree = MAX(op1->degree, op2->degree);
	rop->ntt = op1->ntt;

	bgv_poly_add(idx, &rop->poly[0], &op1->poly[0], &op2->poly[0]);
	bgv_poly_add(idx, &rop->poly[1], &op1->poly[1], &op2->poly[1]);

	if (rop->degree == 2)
		return;

	bgv_poly_add(idx, &rop->poly[2], &op1->poly[2], &op2->poly[2]);
}

void *
bgv_alloc(size_t len, size_t size)
{

	return bgv_util_alloc(len, size);
}

void
bgv_dealloc(void *ptr)
{

	bgv_util_dealloc(ptr);
}

void
bgv_decode(bgv_Message *m, const bgv_Poly *p)
{
	bgv_Poly *cpy;
	const bgv_Digit *rns[BGV_QLEN];
	bgv_Digit mods[BGV_QLEN];
	mpz_t mod;
	size_t len, i, ii;

	cpy = bgv_alloc(BGV_QLEN, sizeof *cpy);
	mpz_init(mod);

	len = 0;
	for (i = 0; i < BGV_QLEN; ++i) {
		if (bgv_q[i].drop != 0)
			continue;

		memcpy(&cpy[i], &p[i], sizeof cpy[i]);
		bgv_poly_intt(i, &cpy[i]);
		bgv_poly_deinit(i, &cpy[i], &cpy[i]);
		rns[len++] = cpy[i].value;
	}

	bgv_mod_digits(mods);
	bgv_msg_crt(m, rns, mods, len);

	bgv_mod_mpz(mod);
	bgv_msg_cmod(m, m, mod);

	for (i = 0; i < BGV_QLEN; ++i)
		for (ii = 0; ii < bgv_q[i].drop; ++ii)
			bgv_msg_mulc(m, m, bgv_q[i].value);

	bgv_dealloc(cpy);
	mpz_clear(mod);
}

void
bgv_decrypt(size_t idx, bgv_Poly *p, const bgv_KeySecret *sk, bgv_Ciphertext *ct)
{
	bgv_Poly *s;

	s = bgv_alloc(1, sizeof *s);

	ciphertext_ntt(idx, ct);

	bgv_poly_init(idx, s, &sk->poly);
	bgv_poly_ntt(idx, s);
	bgv_poly_muladd(idx, p, &ct->poly[1], s, &ct->poly[0]);

	if (ct->degree == 3) {
		bgv_poly_mul(idx, s, s, s);
		bgv_poly_muladd(idx, p, &ct->poly[2], s, p);
	}

	bgv_dealloc(s);
}

void
bgv_deinit(void)
{
	size_t i;

	bgv_msgmod_deinit();
	for (i = 0; i < BGV_QPLEN; ++i)
		bgv_mod_deinit(i);
}

void
bgv_encode(size_t idx, bgv_Poly *p, const bgv_Message *m)
{
	
	bgv_poly_mod(idx, p, m->value);
	bgv_poly_ntt(idx, p);
}

void
bgv_encrypt(size_t idx, bgv_Ciphertext *ct, bgv_KeyPublic *pk, const bgv_Poly *p, bgv_Seed *seed)
{
	bgv_Poly *coins;

	coins = bgv_alloc(4, sizeof *coins);

	/* sample and init coins: a, u, e0, e1 */
	bgv_poly_sample_uniform(idx, &coins[0], &pk->seed);
	pk->seed.init = 1;

	bgv_poly_sample_secret(&coins[1], seed);
	bgv_poly_init(idx, &coins[1], &coins[1]);
	bgv_poly_ntt(idx, &coins[1]);

	bgv_poly_sample_error(&coins[2], seed);
	bgv_poly_init_error(idx, &coins[2], &coins[2]);

	bgv_poly_sample_error(&coins[3], seed);
	bgv_poly_init_error(idx, &coins[3], &coins[3]);
	seed->init = 1;

	/* compute ciphertext */
	bgv_poly_muladd(idx, &ct->poly[0], &pk->poly, &coins[1], p);
	bgv_poly_muladd(idx, &ct->poly[1], &coins[0], &coins[1], &coins[3]);
	memset(&ct->poly[2], 0, sizeof ct->poly[2]);
	ct->degree = 2;
	ct->ntt = 1;

	bgv_dealloc(coins);
}

void
bgv_init(void)
{
	mpz_t tmp;
	size_t div, rest, len, i, k;

	mpz_init(tmp);

	bgv_msgmod_init();

	for (i = 0; i < BGV_QPLEN; ++i)
		bgv_mod_init(i, bgv_t.value);

	/* initialize chunks */
	div = BGV_QLEN / BGV_OMEGA;
	rest = BGV_QLEN % BGV_OMEGA;

	for (k = 0; k < BGV_OMEGA; ++k) {
		len = div;
		if (rest > 0)
			++len, --rest;
		bgv_chunk_len[k] = len;
	}

	for (i = 0; i < BGV_QLEN; ++i)
		bgv_chunk_idx[i % BGV_OMEGA][i / BGV_OMEGA] = i;

	mpz_clear(tmp);
}

void
bgv_keygen_public(size_t idx, bgv_KeyPublic *pk, const bgv_KeySecret *sk, bgv_Seed *seed)
{
	bgv_Poly *p;

	p = bgv_alloc(3, sizeof *p);

	/* sample and init polynomials */
	bgv_sample_seed(&pk->seed);
	bgv_poly_sample_uniform(idx, &p[0], &pk->seed);
	pk->seed.init = 1;

	bgv_poly_init(idx, &p[1], &sk->poly);
	bgv_poly_ntt(idx, &p[1]);

	bgv_poly_sample_error(&p[2], seed);
	bgv_poly_init_error(idx, &p[2], &p[2]);
	seed->init = 1;

	/* compute public key */
	bgv_poly_neg(idx, &p[0], &p[0]);
	bgv_poly_muladd(idx, &pk->poly, &p[0], &p[1], &p[2]);

	bgv_dealloc(p);
}

void
bgv_keygen_secret(bgv_KeySecret *sk)
{
	bgv_Seed seed;

	bgv_sample_seed(&seed);
	bgv_poly_sample_secret(&sk->poly, &seed);
}

void
bgv_keygen_switch(size_t idx, bgv_KeySwitch *ksw, const bgv_KeySecret *sk, const bgv_KeySecret *sk2, bgv_Seed *seed)
{
	bgv_Poly *a, *e, *s;
	size_t k;

	a = bgv_alloc(BGV_OMEGA, sizeof *a);
	e = bgv_alloc(BGV_OMEGA, sizeof *e);
	s = bgv_alloc(2, sizeof *s);

	bgv_sample_seed(&ksw->seed);
	for (k = 0; k < BGV_OMEGA; ++k) {
		bgv_poly_sample_uniform(idx, &a[k], &ksw->seed);
		bgv_poly_neg(idx, &a[k], &a[k]);
	}
	ksw->seed.init = 1;

	for (k = 0; k < BGV_OMEGA; ++k) {
		bgv_poly_sample_error(&e[k], &seed[k]);
		bgv_poly_init_error(idx, &e[k], &e[k]);
		seed[k].init = 1;
	}

	bgv_poly_init(idx, &s[0], &sk->poly);
	bgv_poly_ntt(idx, &s[0]);

	bgv_poly_init(idx, &s[1], &sk2->poly);
	bgv_poly_ntt(idx, &s[1]);

	for (k = 0; k < BGV_OMEGA; ++k) {
		bgv_poly_muladd(idx, &ksw->poly[k], &a[k], &s[0], &e[k]);

		if (idx % BGV_OMEGA == k)
			bgv_poly_mulcadd(idx, &ksw->poly[k], &s[1], bgv_q[idx].P, &ksw->poly[k]);
	}

	bgv_dealloc(a);
	bgv_dealloc(e);
	bgv_dealloc(s);
}

void
bgv_keygen_switch2(size_t idx, bgv_KeySwitch *ksw, const bgv_KeySecret *sk, bgv_Seed *seed)
{
	bgv_KeySecret *s;

	s = bgv_alloc(1, sizeof *s);

	bgv_poly_init(idx, &s->poly, &sk->poly);
	bgv_poly_ntt(idx, &s->poly);
	bgv_poly_mul(idx, &s->poly, &s->poly, &s->poly);
	bgv_poly_intt(idx, &s->poly);
	bgv_poly_deinit(idx, &s->poly, &s->poly);
	bgv_poly_cmod(idx, &s->poly, &s->poly);

	bgv_keygen_switch(idx, ksw, sk, s, seed);

	bgv_dealloc(s);
}

void
bgv_keygen_switchr(size_t idx, bgv_KeySwitch *ksw, const bgv_KeySecret *sk, size_t steps, bgv_Seed *seed)
{
	bgv_KeySecret *s;

	s = bgv_alloc(1, sizeof *s);

	bgv_poly_rot(idx, &s->poly, &sk->poly, steps);
	bgv_keygen_switch(idx, ksw, sk, s, seed);

	bgv_dealloc(s);
}

void
bgv_keyswitch(bgv_Ciphertext *ct, bgv_KeySwitch *ksw, size_t dsw)
{
	bgv_CiphertextSwitch *csw;
	bgv_Delta *delta;
	size_t i;

	csw = bgv_alloc(BGV_QPLEN, sizeof *csw);
	delta = bgv_alloc(BGV_PLEN, sizeof *delta);

	bgv_keyswitch_ext(csw, ct, dsw);

	for (i = 0; i < BGV_QPLEN; ++i) {
		if (bgv_q[i].drop != 0)
			continue;

		bgv_keyswitch_dot(i, &ct[i], &csw[i], &ksw[i]);
	}

	for (i = BGV_QLEN; i < BGV_QPLEN; ++i)
		bgv_keyswitch_delta(i, &delta[i - BGV_QLEN], &ct[i]);

	for (i = 0; i < BGV_QLEN; ++i) {
		if (bgv_q[i].drop != 0)
			continue;

		bgv_keyswitch_switch(i, &ct[i], delta);
	}

	bgv_dealloc(delta);
	bgv_dealloc(csw);
}
void
bgv_keyswitch_delta(size_t idx, bgv_Delta *delta, bgv_Ciphertext *ct)
{

	bgv_modswitch_delta(idx, delta, ct);
}

void
bgv_keyswitch_dot(size_t idx, bgv_Ciphertext *ct, bgv_CiphertextSwitch *csw, bgv_KeySwitch *ksw)
{
	bgv_Poly *a;
	size_t k;

	a = bgv_alloc(1, sizeof *a);

	ciphertext_switch_ntt(idx, csw);

	memset(ct, 0, sizeof *ct);
	ct->degree = 2;
	ct->ntt = 1;

	for (k = 0; k < BGV_OMEGA; ++k)
		bgv_poly_muladd(idx, &ct->poly[0], &csw->ext[k], &ksw->poly[k], &ct->poly[0]);
	bgv_poly_mulcadd(idx, &ct->poly[0], &csw->ct.poly[0], bgv_q[idx].P, &ct->poly[0]);

	for (k = 0; k < BGV_OMEGA; ++k) {
		bgv_poly_sample_uniform(idx, a, &ksw->seed);
		bgv_poly_muladd(idx, &ct->poly[1], &csw->ext[k], a, &ct->poly[1]);
	}
	bgv_poly_mulcadd(idx, &ct->poly[1], &csw->ct.poly[1], csw->poly == 2 ? bgv_q[idx].P : 0, &ct->poly[1]);
	ksw->seed.init = 1;

	bgv_dealloc(a);
}

void
bgv_keyswitch_ext(bgv_CiphertextSwitch *csw, bgv_Ciphertext *ct, size_t poly)
{
	const bgv_Poly *p[BGV_KMAX];
	size_t idx, len, i, k;

	for (idx = 0; idx < BGV_QLEN; ++idx) {
		if (bgv_q[idx].drop != 0)
			continue;

		ciphertext_intt(idx, &ct[idx]);
		memcpy(&csw[idx].ct, &ct[idx], sizeof csw[idx].ct);
	}

	for (idx = 0; idx < BGV_QPLEN; ++idx) {
		if (bgv_q[idx].drop != 0)
			continue;

		for (k = 0; k < BGV_OMEGA; ++k) {
			len = bgv_chunk_len[k];
			for (i = 0; i < len; ++i)
				p[i] = &ct[bgv_chunk_idx[k][i]].poly[poly];
			bgv_poly_ext(idx, &csw[idx].ext[k], p, bgv_chunk_idx[k], len);
		}
		csw[idx].poly = poly;
	}
}

void
bgv_keyswitch_switch(size_t idx, bgv_Ciphertext *ct, bgv_Delta *delta)
{
	bgv_Delta *tmp;
	size_t P[BGV_PLEN], i;

	tmp = bgv_alloc(1, sizeof *tmp);

	for (i = 0; i < BGV_PLEN; ++i)
		P[i] = BGV_QLEN + i;

	delta_ext(idx, tmp, delta, P, BGV_PLEN);
	bgv_modswitch(idx, ct, tmp);

	bgv_dealloc(tmp);
}

void
bgv_modswitch(size_t idx, bgv_Ciphertext *ct, bgv_Delta *delta)
{

	assert(ct->degree == delta->ct.degree);

	ciphertext_ntt(idx, ct);
	ciphertext_ntt(idx, &delta->ct);

	bgv_poly_mulcadd(idx, &ct->poly[0], &delta->ct.poly[0], bgv_q[idx].t, &ct->poly[0]);
	bgv_poly_mulc(idx, &ct->poly[0], &ct->poly[0], delta->inv);

	bgv_poly_mulcadd(idx, &ct->poly[1], &delta->ct.poly[1], bgv_q[idx].t, &ct->poly[1]);
	bgv_poly_mulc(idx, &ct->poly[1], &ct->poly[1], delta->inv);

	if (ct->degree == 2)
		return;

	bgv_poly_mulcadd(idx, &ct->poly[2], &delta->ct.poly[2], bgv_q[idx].t, &ct->poly[2]);
	bgv_poly_mulc(idx, &ct->poly[2], &ct->poly[2], delta->inv);
}


void
bgv_modswitch_delta(size_t idx, bgv_Delta *delta, bgv_Ciphertext *ct)
{
	bgv_Digit inv;

	ciphertext_ntt(idx, ct);

	delta->ct.degree = ct->degree;
	delta->ct.ntt = 0;
	delta->inv = 0;
	delta->idx = idx;

	inv = bgv_mod_negtinv(idx);

	bgv_poly_mulc(idx, &delta->ct.poly[0], &ct->poly[0], inv);
	bgv_poly_intt(idx, &delta->ct.poly[0]);

	bgv_poly_mulc(idx, &delta->ct.poly[1], &ct->poly[1], inv);
	bgv_poly_intt(idx, &delta->ct.poly[1]);

	if (delta->ct.degree == 2)
		return;

	bgv_poly_mulc(idx, &delta->ct.poly[2], &ct->poly[2], inv);
	bgv_poly_intt(idx, &delta->ct.poly[2]);
}

void
bgv_modswitch_ext(size_t idx, bgv_Ciphertext *ct, bgv_Delta *delta)
{
	bgv_Delta *cpy;

	cpy = bgv_alloc(1, sizeof *cpy);

	delta_ext(idx, cpy, delta, &delta->idx, 1);
	bgv_modswitch(idx, ct, cpy);

	bgv_dealloc(cpy);
}

void
bgv_mul(size_t idx, bgv_Ciphertext *rop, bgv_Ciphertext *op1, bgv_Ciphertext *op2)
{

	assert(op1->degree == 2);
	assert(op2->degree == 2);

	ciphertext_ntt(idx, op1);
	ciphertext_ntt(idx, op2);

	rop->degree = 3;
	rop->ntt = 1;

	bgv_poly_mul(idx, &rop->poly[2], &op1->poly[1], &op2->poly[1]);
	bgv_poly_addmul(idx, &rop->poly[1], &op1->poly[0], &op1->poly[1], &op2->poly[0], &op2->poly[1]);
	bgv_poly_mul(idx, &rop->poly[0], &op1->poly[0], &op2->poly[0]);
	bgv_poly_sub(idx, &rop->poly[1], &rop->poly[1], &rop->poly[0]);
	bgv_poly_sub(idx, &rop->poly[1], &rop->poly[1], &rop->poly[2]);
}

void
bgv_mulc(size_t idx, bgv_Ciphertext *rop, bgv_Ciphertext *ct, const bgv_Poly *p)
{

	ciphertext_ntt(idx, ct);

	rop->degree = ct->degree;
	rop->ntt = ct->ntt;

	bgv_poly_mul(idx, &rop->poly[0], &ct->poly[0], p);
	bgv_poly_mul(idx, &rop->poly[1], &ct->poly[1], p);

	if (rop->degree == 2)
		return;

	bgv_poly_mul(idx, &rop->poly[2], &ct->poly[2], p);
}

void
bgv_rot(size_t idx, bgv_Ciphertext *rop, bgv_Ciphertext *ct, size_t steps)
{

	assert(rop != ct);

	ciphertext_intt(idx, ct);

	rop->degree = ct->degree;
	rop->ntt = ct->ntt;

	bgv_poly_rot(idx, &rop->poly[0], &ct->poly[0], steps);
	bgv_poly_rot(idx, &rop->poly[1], &ct->poly[1], steps);

	if (rop->degree == 2)
		return;

	bgv_poly_rot(idx, &rop->poly[2], &ct->poly[2], steps);
}

void
bgv_rot_csw(size_t idx, bgv_CiphertextSwitch *rop, bgv_CiphertextSwitch *csw, size_t steps)
{
	size_t k;

	assert(rop != csw);

	ciphertext_switch_intt(idx, csw);

	bgv_rot(idx, &rop->ct, &csw->ct, steps);
	for (k = 0; k < BGV_OMEGA; ++k)
		bgv_poly_rot(idx, &rop->ext[k], &csw->ext[k], steps);
}

void
bgv_rot_csw_inplace(size_t idx, bgv_CiphertextSwitch *csw, size_t steps)
{
	bgv_CiphertextSwitch *cpy;

	cpy = bgv_alloc(1, sizeof *cpy);
	memcpy(cpy, csw, sizeof *cpy);

	bgv_rot_csw(idx, csw, cpy, steps);

	bgv_dealloc(cpy);
}

void
bgv_rot_inplace(size_t idx, bgv_Ciphertext *ct, size_t steps)
{
	bgv_Ciphertext *cpy;

	cpy = bgv_alloc(1, sizeof *cpy);
	memcpy(cpy, ct, sizeof *cpy);

	bgv_rot(idx, ct, cpy, steps);

	bgv_dealloc(cpy);
}

#endif /* BGV_BGV_IMPL */
