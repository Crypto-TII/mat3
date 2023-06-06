#include "bgv.h"

size_t tiimat3_chunk_idx[TIIMAT3_OMEGA][TIIMAT3_KMAX];
size_t tiimat3_chunk_len[TIIMAT3_OMEGA];

#define MAX(a, b) ((a) > (b) ? (a) : (b))

static void ciphertext_intt(size_t idx, tiimat3_Ciphertext *ct);
static void ciphertext_ntt(size_t idx, tiimat3_Ciphertext *ct);
static void ciphertext_switch_intt(size_t idx, tiimat3_CiphertextSwitch *csw);
static void ciphertext_switch_ntt(size_t idx, tiimat3_CiphertextSwitch *csw);
static void delta_ext(size_t idx, tiimat3_Delta *rop, tiimat3_Delta *delta, const size_t *base, size_t len);

static void
ciphertext_intt(size_t idx, tiimat3_Ciphertext *ct)
{

	if (ct->ntt == 0)
		return;
	ct->ntt = 0;

	tiimat3_poly_intt(idx, &ct->poly[0]);
	tiimat3_poly_intt(idx, &ct->poly[1]);

	if (ct->degree == 2)
		return;

	tiimat3_poly_intt(idx, &ct->poly[2]);
}

static void
ciphertext_ntt(size_t idx, tiimat3_Ciphertext *ct)
{

	if (ct->ntt == 1)
		return;
	ct->ntt = 1;

	tiimat3_poly_ntt(idx, &ct->poly[0]);
	tiimat3_poly_ntt(idx, &ct->poly[1]);

	if (ct->degree == 2)
		return;

	tiimat3_poly_ntt(idx, &ct->poly[2]);
}

static void
ciphertext_switch_intt(size_t idx, tiimat3_CiphertextSwitch *csw)
{
	size_t k;

	if (csw->ct.ntt == 0)
		return;

	ciphertext_intt(idx, &csw->ct);
	for (k = 0; k < TIIMAT3_OMEGA; ++k)
		tiimat3_poly_intt(idx, &csw->ext[k]);
}

static void
ciphertext_switch_ntt(size_t idx, tiimat3_CiphertextSwitch *csw)
{
	size_t k;

	if (csw->ct.ntt == 1)
		return;

	ciphertext_ntt(idx, &csw->ct);
	for (k = 0; k < TIIMAT3_OMEGA; ++k)
		tiimat3_poly_ntt(idx, &csw->ext[k]);
}

void
delta_ext(size_t idx, tiimat3_Delta *rop, tiimat3_Delta *delta, const size_t *base, size_t len)
{
	const tiimat3_Poly *p[TIIMAT3_QPLEN];
	size_t i;

	ciphertext_intt(idx, &delta->ct);

	rop->ct.degree = delta->ct.degree;
	rop->ct.ntt = 0;
	rop->inv = tiimat3_mod_inv(idx, base, len);

	for (i = 0; i < len; ++i)
		p[i] = &delta[i].ct.poly[0];
	tiimat3_poly_ext(idx, &rop->ct.poly[0], p, base, len);

	for (i = 0; i < len; ++i)
		p[i] = &delta[i].ct.poly[1];
	tiimat3_poly_ext(idx, &rop->ct.poly[1], p, base, len);

	if (rop->ct.degree == 2)
		return;

	for (i = 0; i < len; ++i)
		p[i] = &delta[i].ct.poly[2];
	tiimat3_poly_ext(idx, &rop->ct.poly[2], p, base, len);
}

void
tiimat3_bgv_add(size_t idx, tiimat3_Ciphertext *rop, tiimat3_Ciphertext *op1, tiimat3_Ciphertext *op2)
{

	if (op1->ntt != op2->ntt) {
		ciphertext_ntt(idx, op1);
		ciphertext_ntt(idx, op2);
	}

	rop->degree = MAX(op1->degree, op2->degree);
	rop->ntt = op1->ntt;

	tiimat3_poly_add(idx, &rop->poly[0], &op1->poly[0], &op2->poly[0]);
	tiimat3_poly_add(idx, &rop->poly[1], &op1->poly[1], &op2->poly[1]);

	if (rop->degree == 2)
		return;

	tiimat3_poly_add(idx, &rop->poly[2], &op1->poly[2], &op2->poly[2]);
}

void
tiimat3_bgv_decode(tiimat3_Message *m, const tiimat3_Poly *p)
{
	tiimat3_Poly *cpy;
	const tiimat3_Digit *rns[TIIMAT3_QLEN];
	tiimat3_Digit mods[TIIMAT3_QLEN];
	mpz_t mod;
	size_t len, i, ii;

	cpy = tiimat3_util_alloc(TIIMAT3_QLEN, sizeof *cpy);
	mpz_init(mod);

	len = 0;
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		if (tiimat3_q[i].drop != 0)
			continue;

		memcpy(&cpy[i], &p[i], sizeof cpy[i]);
		tiimat3_poly_intt(i, &cpy[i]);
		tiimat3_poly_deinit(i, &cpy[i], &cpy[i]);
		rns[len++] = cpy[i].value;
	}

	tiimat3_mod_digits(mods);
	tiimat3_msg_crt(m, rns, mods, len);

	tiimat3_mod_mpz(mod);
	tiimat3_msg_cmod(m, m, mod);

	for (i = 0; i < TIIMAT3_QLEN; ++i)
		for (ii = 0; ii < tiimat3_q[i].drop; ++ii)
			tiimat3_msg_mulc(m, m, tiimat3_q[i].value);

	tiimat3_util_dealloc(cpy);
	mpz_clear(mod);
}

void
tiimat3_bgv_decrypt(size_t idx, tiimat3_Poly *p, const tiimat3_KeySecret *sk, tiimat3_Ciphertext *ct)
{
	tiimat3_Poly *s;

	s = tiimat3_util_alloc(1, sizeof *s);

	ciphertext_ntt(idx, ct);

	tiimat3_poly_init(idx, s, &sk->poly);
	tiimat3_poly_ntt(idx, s);
	tiimat3_poly_muladd(idx, p, &ct->poly[1], s, &ct->poly[0]);

	if (ct->degree == 3) {
		tiimat3_poly_mul(idx, s, s, s);
		tiimat3_poly_muladd(idx, p, &ct->poly[2], s, p);
	}

	tiimat3_util_dealloc(s);
}

void
tiimat3_bgv_deinit(void)
{
	size_t i;

	tiimat3_msgmod_deinit();
	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_mod_deinit(i);
}

void
tiimat3_bgv_encode(size_t idx, tiimat3_Poly *p, const tiimat3_Message *m)
{
	
	tiimat3_poly_mod(idx, p, m->value);
	tiimat3_poly_ntt(idx, p);
}

void
tiimat3_bgv_encrypt(size_t idx, tiimat3_Ciphertext *ct, tiimat3_KeyPublic *pk, const tiimat3_Poly *p, tiimat3_Seed *seed)
{
	tiimat3_Poly *coins;

	coins = tiimat3_util_alloc(4, sizeof *coins);

	/* sample and init coins: a, u, e0, e1 */
	tiimat3_poly_sample_uniform(idx, &coins[0], &pk->seed);
	pk->seed.init = 1;

	tiimat3_poly_sample_secret(&coins[1], seed);
	tiimat3_poly_init(idx, &coins[1], &coins[1]);
	tiimat3_poly_ntt(idx, &coins[1]);

	tiimat3_poly_sample_error(&coins[2], seed);
	tiimat3_poly_init_error(idx, &coins[2], &coins[2]);

	tiimat3_poly_sample_error(&coins[3], seed);
	tiimat3_poly_init_error(idx, &coins[3], &coins[3]);
	seed->init = 1;

	/* compute ciphertext */
	tiimat3_poly_muladd(idx, &ct->poly[0], &pk->poly, &coins[1], p);
	tiimat3_poly_muladd(idx, &ct->poly[1], &coins[0], &coins[1], &coins[3]);
	memset(&ct->poly[2], 0, sizeof ct->poly[2]);
	ct->degree = 2;
	ct->ntt = 1;

	tiimat3_util_dealloc(coins);
}

void
tiimat3_bgv_init(void)
{
	mpz_t tmp;
	size_t div, rest, len, i, k;

	mpz_init(tmp);

	tiimat3_msgmod_init();

	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_mod_init(i, tiimat3_t.value);

	/* initialize chunks */
	div = TIIMAT3_QLEN / TIIMAT3_OMEGA;
	rest = TIIMAT3_QLEN % TIIMAT3_OMEGA;

	for (k = 0; k < TIIMAT3_OMEGA; ++k) {
		len = div;
		if (rest > 0)
			++len, --rest;
		tiimat3_chunk_len[k] = len;
	}

	for (i = 0; i < TIIMAT3_QLEN; ++i)
		tiimat3_chunk_idx[i % TIIMAT3_OMEGA][i / TIIMAT3_OMEGA] = i;

	mpz_clear(tmp);
}

void
tiimat3_bgv_keygen_public(size_t idx, tiimat3_KeyPublic *pk, const tiimat3_KeySecret *sk, tiimat3_Seed *seed)
{
	tiimat3_Poly *p;

	p = tiimat3_util_alloc(3, sizeof *p);

	/* sample and init polynomials */
	tiimat3_random_seed(&pk->seed);
	tiimat3_poly_sample_uniform(idx, &p[0], &pk->seed);
	pk->seed.init = 1;

	tiimat3_poly_init(idx, &p[1], &sk->poly);
	tiimat3_poly_ntt(idx, &p[1]);

	tiimat3_poly_sample_error(&p[2], seed);
	tiimat3_poly_init_error(idx, &p[2], &p[2]);
	seed->init = 1;

	/* compute public key */
	tiimat3_poly_neg(idx, &p[0], &p[0]);
	tiimat3_poly_muladd(idx, &pk->poly, &p[0], &p[1], &p[2]);

	tiimat3_util_dealloc(p);
}

void
tiimat3_bgv_keygen_secret(tiimat3_KeySecret *sk)
{
	tiimat3_Seed seed;

	tiimat3_random_seed(&seed);
	tiimat3_poly_sample_secret(&sk->poly, &seed);
}

void
tiimat3_bgv_keygen_switch(size_t idx, tiimat3_KeySwitch *ksw, const tiimat3_KeySecret *sk, const tiimat3_KeySecret *sk2, tiimat3_Seed *seed)
{
	tiimat3_Poly *a, *e, *s;
	size_t k;

	a = tiimat3_util_alloc(TIIMAT3_OMEGA, sizeof *a);
	e = tiimat3_util_alloc(TIIMAT3_OMEGA, sizeof *e);
	s = tiimat3_util_alloc(2, sizeof *s);

	tiimat3_random_seed(&ksw->seed);
	for (k = 0; k < TIIMAT3_OMEGA; ++k) {
		tiimat3_poly_sample_uniform(idx, &a[k], &ksw->seed);
		tiimat3_poly_neg(idx, &a[k], &a[k]);
	}
	ksw->seed.init = 1;

	for (k = 0; k < TIIMAT3_OMEGA; ++k) {
		tiimat3_poly_sample_error(&e[k], &seed[k]);
		tiimat3_poly_init_error(idx, &e[k], &e[k]);
		seed[k].init = 1;
	}

	tiimat3_poly_init(idx, &s[0], &sk->poly);
	tiimat3_poly_ntt(idx, &s[0]);

	tiimat3_poly_init(idx, &s[1], &sk2->poly);
	tiimat3_poly_ntt(idx, &s[1]);

	for (k = 0; k < TIIMAT3_OMEGA; ++k) {
		tiimat3_poly_muladd(idx, &ksw->poly[k], &a[k], &s[0], &e[k]);

		if (idx % TIIMAT3_OMEGA == k)
			tiimat3_poly_mulcadd(idx, &ksw->poly[k], &s[1], tiimat3_q[idx].P, &ksw->poly[k]);
	}

	tiimat3_util_dealloc(a);
	tiimat3_util_dealloc(e);
	tiimat3_util_dealloc(s);
}

void
tiimat3_bgv_keygen_switch2(size_t idx, tiimat3_KeySwitch *ksw, const tiimat3_KeySecret *sk, tiimat3_Seed *seed)
{
	tiimat3_KeySecret *s;

	s = tiimat3_util_alloc(1, sizeof *s);

	tiimat3_poly_init(idx, &s->poly, &sk->poly);
	tiimat3_poly_ntt(idx, &s->poly);
	tiimat3_poly_mul(idx, &s->poly, &s->poly, &s->poly);
	tiimat3_poly_intt(idx, &s->poly);
	tiimat3_poly_deinit(idx, &s->poly, &s->poly);
	tiimat3_poly_cmod(idx, &s->poly, &s->poly);

	tiimat3_bgv_keygen_switch(idx, ksw, sk, s, seed);

	tiimat3_util_dealloc(s);
}

void
tiimat3_bgv_keygen_switchr(size_t idx, tiimat3_KeySwitch *ksw, const tiimat3_KeySecret *sk, size_t steps, tiimat3_Seed *seed)
{
	tiimat3_KeySecret *s;

	s = tiimat3_util_alloc(1, sizeof *s);

	tiimat3_poly_rot(idx, &s->poly, &sk->poly, steps);
	tiimat3_bgv_keygen_switch(idx, ksw, sk, s, seed);

	tiimat3_util_dealloc(s);
}

void
tiimat3_bgv_keyswitch(tiimat3_Ciphertext *ct, tiimat3_KeySwitch *ksw, size_t dsw)
{
	tiimat3_CiphertextSwitch *csw;
	tiimat3_Delta *delta;
	size_t i;

	csw = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *csw);
	delta = tiimat3_util_alloc(TIIMAT3_PLEN, sizeof *delta);

	tiimat3_bgv_keyswitch_ext(csw, ct, dsw);

	for (i = 0; i < TIIMAT3_QPLEN; ++i) {
		if (tiimat3_q[i].drop != 0)
			continue;

		tiimat3_bgv_keyswitch_dot(i, &ct[i], &csw[i], &ksw[i]);
	}

	for (i = TIIMAT3_QLEN; i < TIIMAT3_QPLEN; ++i)
		tiimat3_bgv_keyswitch_delta(i, &delta[i - TIIMAT3_QLEN], &ct[i]);

	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		if (tiimat3_q[i].drop != 0)
			continue;

		tiimat3_bgv_keyswitch_switch(i, &ct[i], delta);
	}

	tiimat3_util_dealloc(delta);
	tiimat3_util_dealloc(csw);
}
void
tiimat3_bgv_keyswitch_delta(size_t idx, tiimat3_Delta *delta, tiimat3_Ciphertext *ct)
{

	tiimat3_bgv_modswitch_delta(idx, delta, ct);
}

void
tiimat3_bgv_keyswitch_dot(size_t idx, tiimat3_Ciphertext *ct, tiimat3_CiphertextSwitch *csw, tiimat3_KeySwitch *ksw)
{
	tiimat3_Poly *a;
	size_t k;

	a = tiimat3_util_alloc(1, sizeof *a);

	ciphertext_switch_ntt(idx, csw);

	memset(ct, 0, sizeof *ct);
	ct->degree = 2;
	ct->ntt = 1;

	for (k = 0; k < TIIMAT3_OMEGA; ++k)
		tiimat3_poly_muladd(idx, &ct->poly[0], &csw->ext[k], &ksw->poly[k], &ct->poly[0]);
	tiimat3_poly_mulcadd(idx, &ct->poly[0], &csw->ct.poly[0], tiimat3_q[idx].P, &ct->poly[0]);

	for (k = 0; k < TIIMAT3_OMEGA; ++k) {
		tiimat3_poly_sample_uniform(idx, a, &ksw->seed);
		tiimat3_poly_muladd(idx, &ct->poly[1], &csw->ext[k], a, &ct->poly[1]);
	}
	tiimat3_poly_mulcadd(idx, &ct->poly[1], &csw->ct.poly[1], csw->poly == 2 ? tiimat3_q[idx].P : 0, &ct->poly[1]);
	ksw->seed.init = 1;

	tiimat3_util_dealloc(a);
}

void
tiimat3_bgv_keyswitch_ext(tiimat3_CiphertextSwitch *csw, tiimat3_Ciphertext *ct, size_t poly)
{
	const tiimat3_Poly *p[TIIMAT3_KMAX];
	size_t idx, len, i, k;

	for (idx = 0; idx < TIIMAT3_QLEN; ++idx) {
		if (tiimat3_q[idx].drop != 0)
			continue;

		ciphertext_intt(idx, &ct[idx]);
		memcpy(&csw[idx].ct, &ct[idx], sizeof csw[idx].ct);
	}

	for (idx = 0; idx < TIIMAT3_QPLEN; ++idx) {
		if (tiimat3_q[idx].drop != 0)
			continue;

		for (k = 0; k < TIIMAT3_OMEGA; ++k) {
			len = tiimat3_chunk_len[k];
			for (i = 0; i < len; ++i)
				p[i] = &ct[tiimat3_chunk_idx[k][i]].poly[poly];
			tiimat3_poly_ext(idx, &csw[idx].ext[k], p, tiimat3_chunk_idx[k], len);
		}
		csw[idx].poly = poly;
	}
}

void
tiimat3_bgv_keyswitch_switch(size_t idx, tiimat3_Ciphertext *ct, tiimat3_Delta *delta)
{
	tiimat3_Delta *tmp;
	size_t P[TIIMAT3_PLEN], i;

	tmp = tiimat3_util_alloc(1, sizeof *tmp);

	for (i = 0; i < TIIMAT3_PLEN; ++i)
		P[i] = TIIMAT3_QLEN + i;

	delta_ext(idx, tmp, delta, P, TIIMAT3_PLEN);
	tiimat3_bgv_modswitch(idx, ct, tmp);

	tiimat3_util_dealloc(tmp);
}

void
tiimat3_bgv_modswitch(size_t idx, tiimat3_Ciphertext *ct, tiimat3_Delta *delta)
{

	assert(ct->degree == delta->ct.degree);

	ciphertext_ntt(idx, ct);
	ciphertext_ntt(idx, &delta->ct);

	tiimat3_poly_mulcadd(idx, &ct->poly[0], &delta->ct.poly[0], tiimat3_q[idx].t, &ct->poly[0]);
	tiimat3_poly_mulc(idx, &ct->poly[0], &ct->poly[0], delta->inv);

	tiimat3_poly_mulcadd(idx, &ct->poly[1], &delta->ct.poly[1], tiimat3_q[idx].t, &ct->poly[1]);
	tiimat3_poly_mulc(idx, &ct->poly[1], &ct->poly[1], delta->inv);

	if (ct->degree == 2)
		return;

	tiimat3_poly_mulcadd(idx, &ct->poly[2], &delta->ct.poly[2], tiimat3_q[idx].t, &ct->poly[2]);
	tiimat3_poly_mulc(idx, &ct->poly[2], &ct->poly[2], delta->inv);
}


void
tiimat3_bgv_modswitch_delta(size_t idx, tiimat3_Delta *delta, tiimat3_Ciphertext *ct)
{
	tiimat3_Digit inv;

	ciphertext_ntt(idx, ct);

	delta->ct.degree = ct->degree;
	delta->ct.ntt = 0;
	delta->inv = 0;
	delta->idx = idx;

	inv = tiimat3_mod_negtinv(idx);

	tiimat3_poly_mulc(idx, &delta->ct.poly[0], &ct->poly[0], inv);
	tiimat3_poly_intt(idx, &delta->ct.poly[0]);

	tiimat3_poly_mulc(idx, &delta->ct.poly[1], &ct->poly[1], inv);
	tiimat3_poly_intt(idx, &delta->ct.poly[1]);

	if (delta->ct.degree == 2)
		return;

	tiimat3_poly_mulc(idx, &delta->ct.poly[2], &ct->poly[2], inv);
	tiimat3_poly_intt(idx, &delta->ct.poly[2]);
}

void
tiimat3_bgv_modswitch_ext(size_t idx, tiimat3_Ciphertext *ct, tiimat3_Delta *delta)
{
	tiimat3_Delta *cpy;

	cpy = tiimat3_util_alloc(1, sizeof *cpy);

	delta_ext(idx, cpy, delta, &delta->idx, 1);
	tiimat3_bgv_modswitch(idx, ct, cpy);

	tiimat3_util_dealloc(cpy);
}

void
tiimat3_bgv_mul(size_t idx, tiimat3_Ciphertext *rop, tiimat3_Ciphertext *op1, tiimat3_Ciphertext *op2)
{

	assert(op1->degree == 2);
	assert(op2->degree == 2);

	ciphertext_ntt(idx, op1);
	ciphertext_ntt(idx, op2);

	rop->degree = 3;
	rop->ntt = 1;

	tiimat3_poly_mul(idx, &rop->poly[2], &op1->poly[1], &op2->poly[1]);
	tiimat3_poly_addmul(idx, &rop->poly[1], &op1->poly[0], &op1->poly[1], &op2->poly[0], &op2->poly[1]);
	tiimat3_poly_mul(idx, &rop->poly[0], &op1->poly[0], &op2->poly[0]);
	tiimat3_poly_sub(idx, &rop->poly[1], &rop->poly[1], &rop->poly[0]);
	tiimat3_poly_sub(idx, &rop->poly[1], &rop->poly[1], &rop->poly[2]);
}

void
tiimat3_bgv_mulc(size_t idx, tiimat3_Ciphertext *rop, tiimat3_Ciphertext *ct, const tiimat3_Poly *p)
{

	ciphertext_ntt(idx, ct);

	rop->degree = ct->degree;
	rop->ntt = ct->ntt;

	tiimat3_poly_mul(idx, &rop->poly[0], &ct->poly[0], p);
	tiimat3_poly_mul(idx, &rop->poly[1], &ct->poly[1], p);

	if (rop->degree == 2)
		return;

	tiimat3_poly_mul(idx, &rop->poly[2], &ct->poly[2], p);
}

void
tiimat3_bgv_rot(size_t idx, tiimat3_Ciphertext *rop, tiimat3_Ciphertext *ct, size_t steps)
{

	assert(rop != ct);

	ciphertext_intt(idx, ct);

	rop->degree = ct->degree;
	rop->ntt = ct->ntt;

	tiimat3_poly_rot(idx, &rop->poly[0], &ct->poly[0], steps);
	tiimat3_poly_rot(idx, &rop->poly[1], &ct->poly[1], steps);

	if (rop->degree == 2)
		return;

	tiimat3_poly_rot(idx, &rop->poly[2], &ct->poly[2], steps);
}

void
tiimat3_bgv_rot_csw(size_t idx, tiimat3_CiphertextSwitch *rop, tiimat3_CiphertextSwitch *csw, size_t steps)
{
	size_t k;

	assert(rop != csw);

	ciphertext_switch_intt(idx, csw);

	tiimat3_bgv_rot(idx, &rop->ct, &csw->ct, steps);
	for (k = 0; k < TIIMAT3_OMEGA; ++k)
		tiimat3_poly_rot(idx, &rop->ext[k], &csw->ext[k], steps);
}

void
tiimat3_bgv_rot_csw_inplace(size_t idx, tiimat3_CiphertextSwitch *csw, size_t steps)
{
	tiimat3_CiphertextSwitch *cpy;

	cpy = tiimat3_util_alloc(1, sizeof *cpy);
	memcpy(cpy, csw, sizeof *cpy);

	tiimat3_bgv_rot_csw(idx, csw, cpy, steps);

	tiimat3_util_dealloc(cpy);
}

void
tiimat3_bgv_rot_inplace(size_t idx, tiimat3_Ciphertext *ct, size_t steps)
{
	tiimat3_Ciphertext *cpy;

	cpy = tiimat3_util_alloc(1, sizeof *cpy);
	memcpy(cpy, ct, sizeof *cpy);

	tiimat3_bgv_rot(idx, ct, cpy, steps);

	tiimat3_util_dealloc(cpy);
}
