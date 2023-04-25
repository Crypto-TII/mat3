#include "msg_mpz.h"

tiimat3_MessageMod tiimat3_t;

#define BSR64(x) __bsrq(x)

/* https://en.wikipedia.org/wiki/Chinese_remainder_theorem */
static void mpz_cmod(mpz_t rop, const mpz_t op, const mpz_t mod);
static void mpz_crt(mpz_t *poly, size_t degree, const uint64_t **rns, const uint64_t *mods, size_t len);
static void mpz_precomp(mpz_t root, size_t degree, const mpz_t mod);
static void mpz_ntt(mpz_t *poly, size_t degree, const mpz_t mod, const mpz_t root);
static void mpz_intt(mpz_t *poly, size_t degree, const mpz_t mod, const mpz_t root);

static void
mpz_cmod(mpz_t rop, const mpz_t op, const mpz_t mod)
{
	mpz_t bound;

	mpz_init_set(bound, mod);
	mpz_tdiv_q_2exp(bound, mod, 1);

	mpz_mod(rop, op, mod);
	if (mpz_cmp_ui(op, 0) < 0)
		mpz_add(rop, rop, mod);
	if (mpz_cmp(rop, bound) >= 0)
		mpz_sub(rop, rop, mod);

	mpz_clear(bound);
}

static void
mpz_crt(mpz_t *poly, size_t degree, const uint64_t **rns, const uint64_t *mods, size_t len)
{
	mpz_t mod, inv, tmp;
	size_t i, j;

	mpz_inits(mod, inv, tmp, 0);

	mpz_set_ui(mod, mods[0]);
	for (j = 0; j < degree; ++j)
		mpz_set_ui(poly[j], rns[0][j]);

	for (i = 1; i < len; ++i) {
		mpz_set_ui(inv, mods[i]);
		mpz_invert(inv, inv, mod);
		mpz_mul_ui(mod, mod, mods[i]);

		for (j = 0; j < degree; ++j) {
			mpz_sub_ui(poly[j], poly[j], rns[i][j]);
			mpz_mul_ui(poly[j], poly[j], mods[i]);
			mpz_mod(poly[j], poly[j], mod);

			mpz_mul(poly[j], poly[j], inv);
			mpz_add_ui(poly[j], poly[j], rns[i][j]);
			mpz_mod(poly[j], poly[j], mod);
		}
	}

	mpz_clears(mod, inv, tmp, 0);
}

static void
mpz_precomp(mpz_t root, size_t degree, const mpz_t mod)
{
	mpz_t ord, exp, pow, tmp;

	mpz_inits(ord, exp, pow, tmp, 0);

	mpz_sub_ui(ord, mod, 1);
	mpz_tdiv_q_2exp(exp, ord, BSR64(degree) + 1);
	mpz_set_ui(tmp, 1);

	for (;;) {
		mpz_add_ui(tmp, tmp, 1);
		mpz_mod(tmp, tmp, mod);

		mpz_powm(pow, tmp, ord, mod);
		if (mpz_cmp_ui(pow, 1) != 0)
			continue;

		mpz_powm(root, tmp, exp, mod);
		mpz_powm_ui(pow, root, degree, mod);
		if (mpz_cmp_ui(pow, 1) != 0)
			break;
	}

	mpz_clears(ord, exp, pow, tmp, 0);
}

static void
mpz_ntt(mpz_t *poly, size_t degree, const mpz_t mod, const mpz_t root)
{
	mpz_t r, tmp;
	size_t dlog, exp, j, k, l, s;

	mpz_inits(r, tmp, 0);
	dlog = BSR64(degree);

	k = 1;
	for (l = degree >> 1; l > 0; l >>= 1) {
		for (s = 0; s < degree; s = j + l) {
			exp = tiimat3_util_bitrev(k++, dlog);
			mpz_powm_ui(r, root, exp, mod);

			for (j = s; j < s + l; ++j) {
				mpz_mul(tmp, r, poly[j + l]);
				mpz_mod(tmp, tmp, mod);

				mpz_sub(poly[j + l], poly[j], tmp);
				mpz_mod(poly[j + l], poly[j + l], mod);

				mpz_add(poly[j], poly[j], tmp);
				mpz_mod(poly[j], poly[j], mod);
			}
		}
	}

	mpz_clears(r, tmp, 0);
}

static void
mpz_intt(mpz_t *poly, size_t degree, const mpz_t mod, const mpz_t root)
{
	mpz_t dinv, r, tmp;
	size_t dlog, exp, j, k, l, s;

	mpz_inits(dinv, r, tmp, 0);
	mpz_set_ui(dinv, degree);
	mpz_invert(dinv, dinv, mod);
	dlog = BSR64(degree);

	k = 0;
	for (l = 1; l < degree; l <<= 1) {
		for (s = 0; s < degree; s = j + l) {
			exp = 1 + tiimat3_util_bitrev(k++, dlog);
			mpz_powm_ui(r, root, exp, mod);
			mpz_invert(r, r, mod);

			for (j = s; j < s + l; ++j) {
				mpz_set(tmp, poly[j]);

				mpz_add(poly[j], poly[j + l], tmp);
				mpz_mod(poly[j], poly[j], mod);

				mpz_sub(poly[j + l], tmp, poly[j + l]);
				mpz_mod(poly[j + l], poly[j + l], mod);

				mpz_mul(poly[j + l], poly[j + l], r);
				mpz_mod(poly[j + l], poly[j + l], mod);
			}
		}
	}

	for (j = 0; j < degree; ++j) {
		mpz_mul(poly[j], poly[j], dinv);
		mpz_mod(poly[j], poly[j], mod);
	}

	mpz_clears(dinv, r, tmp, 0);
}


tiimat3_Message *
tiimat3_alloc_msg(size_t len)
{
	tiimat3_Message *m;
	size_t i, j;

	m = tiimat3_util_alloc(len, sizeof *m);
	for (i = 0; i < len; ++i)
		for (j = 0; j < TIIMAT3_D; ++j)
			mpz_init(m[i].value[j]);

	return m;
}

void
tiimat3_dealloc_msg(tiimat3_Message *m, size_t len)
{
	size_t i, j;

	for (i = 0; i < len; ++i)
		for (j = 0; j < TIIMAT3_D; ++j)
			mpz_clear(m[i].value[j]);
	free(m);
}

void
tiimat3_msg_cmod(tiimat3_Message *rop, const tiimat3_Message *m, const mpz_t mod)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		mpz_cmod(rop->value[j], m->value[j], mod);
}

void
tiimat3_msg_crt_i64(tiimat3_Message *m, const int64_t **rns, const int64_t *mods, size_t len)
{
	const uint64_t **rns2;
	uint64_t *mods2;
	size_t i;

	mods2 = tiimat3_util_alloc(len, sizeof *mods2);
	rns2 = tiimat3_util_alloc(len, sizeof *rns2);

	for (i = 0; i < len; ++i) {
		mods2[i] = mods[i];
		rns2[i] = (uint64_t *)rns[i];
	}
	tiimat3_msg_crt_u64(m, rns2, mods2, len);

	tiimat3_util_dealloc(mods2);
	tiimat3_util_dealloc(rns2);
}

void
tiimat3_msg_crt_u64(tiimat3_Message *m, const uint64_t **rns, const uint64_t *mods, size_t len)
{

	mpz_crt(m->value, TIIMAT3_D, rns, mods, len);
}

void
tiimat3_msg_mulc(tiimat3_Message *rop, const tiimat3_Message *m, uint64_t c)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j) {
		mpz_mul_ui(rop->value[j], m->value[j], c);
		mpz_mod(rop->value[j], rop->value[j], tiimat3_t.value);
	}
}

void
tiimat3_msg_pack(tiimat3_Message *rop, const tiimat3_Message *m)
{
	tiimat3_Message *cpy;
	size_t dlog, idx, j;

	cpy = tiimat3_alloc_msg(1);

	/* reorder slots */
	idx = 1;
	dlog = BSR64(TIIMAT3_D);
	for (j = 0; j < TIIMAT3_D / 2; ++j) {
		size_t idx1, idx2;

		idx1 = tiimat3_util_bitrev((idx - 1) / 2, dlog);
		idx2 = tiimat3_util_bitrev((2 * TIIMAT3_D - idx - 1) / 2, dlog);

		mpz_set(cpy->value[idx1], m->value[j]);
		mpz_set(cpy->value[idx2], m->value[j | TIIMAT3_D / 2]);

		idx *= 3;
		idx %= 2 * TIIMAT3_D;
	}

	/* actual packing */
	mpz_intt(cpy->value, TIIMAT3_D, tiimat3_t.value, tiimat3_t.root);
	for (j = 0; j < TIIMAT3_D; ++j)
		mpz_set(rop->value[j], cpy->value[j]);

	tiimat3_dealloc_msg(cpy, 1);
}

void
tiimat3_msg_unpack(tiimat3_Message *rop, const tiimat3_Message *m)
{
	tiimat3_Message *cpy;
	size_t dlog, idx, j;

	cpy = tiimat3_alloc_msg(1);

	#if TIIMAT3_ERROR_CSV
	do {
		static FILE *f = 0;
		mpz_t tmp;

		if (f == 0) {
			f = fopen("error.csv", "w+");
			if (f == 0)
				perror("fopen"), exit(1);
		}

		mpz_init(tmp);
		for (j = 0; j < TIIMAT3_D; ++j) {
			mpz_tdiv_q(tmp, tmp, tiimat3_t.value);
			gmp_fprintf(f, j < TIIMAT3_D - 1 ? "%Zd, " : "%Zd\n", tmp);
		}
		mpz_clear(tmp);
	} while (0);
	#endif /* TIIMAT3_ERROR_CSV */

	for (j = 0; j < TIIMAT3_D; ++j)
		mpz_mod(cpy->value[j], m->value[j], tiimat3_t.value);

	/* actual unpacking */
	mpz_ntt(cpy->value, TIIMAT3_D, tiimat3_t.value, tiimat3_t.root);

	/* reorder slots */
	idx = 1;
	dlog = BSR64(TIIMAT3_D);
	for (j = 0; j < TIIMAT3_D / 2; ++j) {
		size_t idx1, idx2;

		idx1 = tiimat3_util_bitrev((idx - 1) / 2, dlog);
		idx2 = tiimat3_util_bitrev((2 * TIIMAT3_D - idx - 1) / 2, dlog);

		mpz_set(rop->value[j], cpy->value[idx1]);
		mpz_set(rop->value[j | TIIMAT3_D / 2], cpy->value[idx2]);

		idx *= 3;
		idx %= 2 * TIIMAT3_D;
	}

	tiimat3_dealloc_msg(cpy, 1);
}

void
tiimat3_msgmod_deinit(void)
{

	mpz_clears(tiimat3_t.value, tiimat3_t.root, 0);
}

void
tiimat3_msgmod_init(void)
{

	mpz_inits(tiimat3_t.value, tiimat3_t.root, 0);
	mpz_set_str(tiimat3_t.value, TIIMAT3_T, 0);
	mpz_precomp(tiimat3_t.root, TIIMAT3_D, tiimat3_t.value);
}
