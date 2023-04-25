#include "ring_i64.h"

tiimat3_Mod   tiimat3_q[TIIMAT3_QPLEN];
tiimat3_Digit tiimat3_rns_mods[TIIMAT3_QPLEN][TIIMAT3_QPLEN];
tiimat3_Digit tiimat3_rns_invs[TIIMAT3_QPLEN][TIIMAT3_QPLEN];
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

__extension__ typedef          __int128  int128_t;
__extension__ typedef unsigned __int128 uint128_t;

#define BSR64(x) __bsrq(x)

/* https://eprint.iacr.org/2018/039 */
static int64_t barrett_precomp(int64_t mod);
static int64_t barrett_reduce(int64_t a, int64_t mod, int64_t precomp);

/* idx-based functions */
static tiimat3_Digit digit_add(size_t idx, tiimat3_Digit a, tiimat3_Digit b);
static tiimat3_Digit digit_deinit(size_t idx, tiimat3_Digit a);
static tiimat3_Digit digit_init(size_t idx, tiimat3_Digit a);
static tiimat3_Digit digit_inv(size_t idx, tiimat3_Digit a);
static tiimat3_Digit digit_mul(size_t idx, tiimat3_Digit a, tiimat3_Digit b);
static tiimat3_Digit digit_neg(size_t idx, tiimat3_Digit a);
static tiimat3_Digit digit_sub(size_t idx, tiimat3_Digit a, tiimat3_Digit b);

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

tiimat3_Digit
digit_add(size_t idx, tiimat3_Digit a, tiimat3_Digit b)
{

	return barrett_reduce(a + b, tiimat3_q[idx].value, tiimat3_q[idx].barrett);
}

tiimat3_Digit
digit_deinit(size_t idx, tiimat3_Digit a)
{

	return montgomery_deinit(a, tiimat3_q[idx].value, tiimat3_q[idx].inv);
}

tiimat3_Digit
digit_init(size_t idx, tiimat3_Digit a)
{

	return montgomery_init(a, tiimat3_q[idx].value, tiimat3_q[idx].inv, tiimat3_q[idx].pow);
}

tiimat3_Digit
digit_inv(size_t idx, tiimat3_Digit a)
{
	int64_t inv;

	inv = tiimat3_util_invmod(digit_deinit(idx, a), tiimat3_q[idx].value);
	return digit_init(idx, inv);
}

tiimat3_Digit
digit_mul(size_t idx, tiimat3_Digit a, tiimat3_Digit b)
{

	return montgomery_mulred(a, b, tiimat3_q[idx].value, tiimat3_q[idx].inv);
}

tiimat3_Digit
digit_neg(size_t idx, tiimat3_Digit a)
{

	(void)idx;

	return -a;
}

tiimat3_Digit
digit_sub(size_t idx, tiimat3_Digit a, tiimat3_Digit b)
{

	return barrett_reduce(a - b, tiimat3_q[idx].value, tiimat3_q[idx].barrett);
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
	dinv = montgomery_init(tiimat3_util_invmod(degree, mod), mod, inv, pow);

	k = 0;
	for (l = 1; l < degree; l <<= 1) {
		for (s = 0; s < degree; s = j + l) {
			int64_t r;
			size_t exp;

			exp = 1 + tiimat3_util_bitrev(k++, dlog);
			r = montgomery_powred(root, exp, mod, inv, pow);
			r = montgomery_deinit(r, mod, inv);
			r = montgomery_init(tiimat3_util_invmod(r, mod), mod, inv, pow);

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

			exp = tiimat3_util_bitrev(k++, dlog);
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
tiimat3_mod_deinit(size_t idx)
{

	(void)idx;
}

void
tiimat3_mod_digits(tiimat3_Digit *q)
{
	size_t idx, i;

	idx = 0;
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		if (tiimat3_q[i].drop != 0)
			continue;

		q[idx++] = tiimat3_q[i].value;
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
	int64_t barrett, mod, inv, P, pow;
	size_t i;

	mod = tiimat3_rns[idx];
	barrett = barrett_precomp(mod);
	inv = montgomery_precomp_inv(mod);
	pow = montgomery_precomp_pow(mod);

	tiimat3_q[idx].value = mod;
	tiimat3_q[idx].barrett = barrett;
	tiimat3_q[idx].inv = inv;
	tiimat3_q[idx].pow = pow;
	tiimat3_q[idx].root = seiler_precomp(TIIMAT3_D, mod, barrett, inv, pow);
	tiimat3_q[idx].t = montgomery_init(t, mod, inv, pow);
	tiimat3_q[idx].idx = idx;
	tiimat3_q[idx].drop = 0;

	for (i = 0; i < TIIMAT3_QPLEN; ++i) {
		tiimat3_Digit tmp;

		tmp = digit_init(idx, tiimat3_rns[i]);

		tiimat3_rns_mods[idx][i] = tmp;
		if (idx == i)
			tiimat3_rns_invs[idx][i] = digit_init(idx, 1);
		else
			tiimat3_rns_invs[idx][i] = digit_inv(idx, tmp);
	}

	P = digit_init(idx, 1);
	for (i = TIIMAT3_QLEN; i < TIIMAT3_QPLEN; ++i)
		P = digit_mul(idx, P, tiimat3_rns_mods[idx][i]);
	tiimat3_q[idx].P = P;
}

tiimat3_Digit
tiimat3_mod_inv(size_t idx, const size_t *invs, size_t len)
{
	tiimat3_Digit inv;
	size_t i;

	inv = digit_init(idx, 1);
	for (i = 0; i < len; ++i)
		inv = digit_mul(idx, inv, tiimat3_rns_invs[idx][invs[i]]);

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

		mpz_mul_ui(q, q, tiimat3_q[i].value);
	}
}

tiimat3_Digit
tiimat3_mod_negtinv(size_t idx)
{
	tiimat3_Digit inv;

	inv = digit_inv(idx, tiimat3_q[idx].t);

	return digit_neg(idx, inv);
}

void
tiimat3_poly_add(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		rop->value[j] = digit_add(idx, op1->value[j], op2->value[j]);
}

void
tiimat3_poly_addmul(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2, const tiimat3_Poly *op3, const tiimat3_Poly *op4)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j) {
		int64_t add12, add34;

		add12 = digit_add(idx, op1->value[j], op2->value[j]);
		add34 = digit_add(idx, op3->value[j], op4->value[j]);
		rop->value[j] = digit_mul(idx, add12, add34);
	}
}

void
tiimat3_poly_cmod(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j) {
		rop->value[j] = p->value[j];
		if (p->value[j] > tiimat3_q[idx].value / 2)
			rop->value[j] -= tiimat3_q[idx].value;
		if (p->value[j] < -tiimat3_q[idx].value / 2)
			rop->value[j] += tiimat3_q[idx].value;
	}
}

void
tiimat3_poly_deinit(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		rop->value[j] = digit_deinit(idx, p->value[j]);
}

void
tiimat3_poly_ext(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly **p, const size_t *base, size_t len)
{
	tiimat3_Digit inv[TIIMAT3_QPLEN], div[TIIMAT3_QPLEN];
	size_t i, j;

	for (i = 0; i < len; ++i) {
		size_t ii;

		inv[i] = digit_init(base[i], 1);
		for (ii = 0; ii < len; ++ii) {
			if (tiimat3_q[base[ii]].drop != 0)
				continue;
			inv[i] = digit_mul(base[i], inv[i], tiimat3_rns_invs[base[i]][base[ii]]);
		}

		div[i] = tiimat3_rns_invs[idx][base[i]];
		for (ii = 0; ii < len; ++ii) {
			if (idx == base[i] && idx == base[ii])
				continue;
			if (tiimat3_q[base[ii]].drop != 0)
				continue;
			div[i] = digit_mul(idx, div[i], tiimat3_rns_mods[idx][base[ii]]);
		}
	}

	for (j = 0; j < TIIMAT3_D; ++j) {
		tiimat3_Digit sum;

		sum = 0;
		for (i = 0; i < len; ++i) {
			tiimat3_Digit tmp;

			if (tiimat3_q[base[i]].drop != 0)
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
tiimat3_poly_init(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		rop->value[j] = digit_init(idx, p->value[j]);
}

void
tiimat3_poly_init_error(size_t idx, tiimat3_Poly *p, const tiimat3_Poly *e)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j) {
		tiimat3_Digit init = digit_init(idx, e->value[j]);
		p->value[j] = digit_mul(idx, init, tiimat3_q[idx].t);
	}
	tiimat3_poly_ntt(idx, p);
}

void
tiimat3_poly_intt(size_t idx, tiimat3_Poly *p)
{

	seiler_intt(p->value, TIIMAT3_D, tiimat3_q[idx].value, tiimat3_q[idx].root, tiimat3_q[idx].barrett, tiimat3_q[idx].inv, tiimat3_q[idx].pow);
}

void
tiimat3_poly_mod_mpz(size_t idx, tiimat3_Poly *rop, const mpz_t *p)
{
	mpz_t tmp;
	size_t j;

	mpz_init(tmp);

	for (j = 0; j < TIIMAT3_D; ++j) {
		mpz_mod_ui(tmp, p[j], tiimat3_q[idx].value);
		rop->value[j] = digit_init(idx, mpz_get_ui(tmp));
	}

	mpz_clear(tmp);
}

void
tiimat3_poly_mod_u64(size_t idx, tiimat3_Poly *rop, const uint64_t *p)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		rop->value[j] = digit_init(idx, p[j] % tiimat3_q[idx].value);
}

void
tiimat3_poly_mul(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		rop->value[j] = digit_mul(idx, op1->value[j], op2->value[j]);
}

void
tiimat3_poly_muladd(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2, const tiimat3_Poly *op3)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		rop->value[j] = digit_add(idx, digit_mul(idx, op1->value[j], op2->value[j]), op3->value[j]);
}

void
tiimat3_poly_mulc(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Digit op2)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		rop->value[j] = digit_mul(idx, op1->value[j], op2);
}

void
tiimat3_poly_mulcadd(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Digit op2, const tiimat3_Poly *op3)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		rop->value[j] = digit_add(idx, digit_mul(idx, op1->value[j], op2), op3->value[j]);
}

void
tiimat3_poly_neg(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		rop->value[j] = digit_neg(idx, p->value[j]);
}

void
tiimat3_poly_ntt(size_t idx, tiimat3_Poly *p)
{

	seiler_ntt(p->value, TIIMAT3_D, tiimat3_q[idx].value, tiimat3_q[idx].root, tiimat3_q[idx].barrett, tiimat3_q[idx].inv, tiimat3_q[idx].pow);
}

void
tiimat3_poly_rot(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *p, size_t steps)
{
	size_t step, i, j;

	assert(rop != p);
	(void)idx;

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
			tmp = -tmp;
		rop->value[i % TIIMAT3_D] = tmp;

		i += step;
		i %= 2 * TIIMAT3_D;
	}
}

void
tiimat3_poly_sample_error(tiimat3_Poly *e, tiimat3_Seed *seed)
{

	tiimat3_random_cbd21_i64(e->value, TIIMAT3_D, seed);
}

void
tiimat3_poly_sample_secret(tiimat3_Poly *s, tiimat3_Seed *seed)
{

	tiimat3_random_cbd1_i64(s->value, TIIMAT3_D, seed);
}

void
tiimat3_poly_sample_uniform(size_t idx, tiimat3_Poly *p, tiimat3_Seed *seed)
{

	tiimat3_random_uniform_i64(p->value, TIIMAT3_D, tiimat3_q[idx].value, seed);
}

void
tiimat3_poly_sub(size_t idx, tiimat3_Poly *rop, const tiimat3_Poly *op1, const tiimat3_Poly *op2)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		rop->value[j] = digit_sub(idx, op1->value[j], op2->value[j]);
}
