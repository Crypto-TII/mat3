#include "triples.h"

static gmp_randstate_t randstate;
static tiimat3_KeySecret sk;
static tiimat3_KeyPublic pk[TIIMAT3_QLEN];
static tiimat3_KeySwitch ksw2[TIIMAT3_QPLEN];
static tiimat3_KeySwitch *kswA;
static tiimat3_KeySwitch *kswAneg;
static tiimat3_KeySwitch *kswB;

static mpz_t MAC[TIIMAT3_SHARES + 1];
static tiimat3_Ciphertext ctMAC[TIIMAT3_SHARES + 1][TIIMAT3_QLEN];
static tiimat3_Message Uphi[TIIMAT3_DIM][2];

static void block_add(tiimat3_Block *rop, const tiimat3_Block *op1, const tiimat3_Block *op2, size_t len);
static tiimat3_Block *block_alloc(size_t len);
static int  block_cmp(const tiimat3_Block *op1, const tiimat3_Block *op2, size_t len);
static void block_cpy(tiimat3_Block *rop, const tiimat3_Block *op);
static void block_dealloc(tiimat3_Block *op, size_t len);
static void block_decrypt(tiimat3_Block *rop, tiimat3_BlockEnc *op);
static void block_encrypt(tiimat3_BlockEnc *rop, const tiimat3_Block *blocks, size_t rot, void (*transform)(tiimat3_Block *, const tiimat3_Block *));
static void block_init_rand(tiimat3_Block *rop, size_t len);
static void block_init_zero(tiimat3_Block *rop, size_t len);
static void block_mul(tiimat3_Block *rop, const tiimat3_Block *op1, const tiimat3_Block *op2, size_t len);
static void block_mulc(tiimat3_Block *rop, const tiimat3_Block *op1, const mpz_t op2, size_t len);
static void block_pack(tiimat3_Message *m, const tiimat3_Block *op, size_t rot, void (*transform)(tiimat3_Block *, const tiimat3_Block *));
static void block_sigma(tiimat3_Block *rop, const tiimat3_Block *op);
static void block_tau(tiimat3_Block *rop, const tiimat3_Block *op);
static void block_unpack(tiimat3_Block *vec, const tiimat3_Message *m);

static void ciphertext_cpy(tiimat3_Ciphertext *rop, tiimat3_Ciphertext *ct, size_t len);
static void ciphertext_rotA(tiimat3_Ciphertext *rop, tiimat3_Ciphertext *ct);
static void ciphertext_rotA_hoisted(tiimat3_Ciphertext *rop, tiimat3_CiphertextSwitch *csw, int k);
static void ciphertext_rotB(tiimat3_Ciphertext *rop, tiimat3_Ciphertext *ct);
static void ciphertext_rotB_hoisted(tiimat3_Ciphertext *rop, tiimat3_CiphertextSwitch *csw, int k);

static void msg_cpy(tiimat3_Message *rop, const tiimat3_Message *m);
static void msg_deinit(tiimat3_Message *m);
static void msg_init(tiimat3_Message *m);
static void msg_rot(tiimat3_Message *rop, const tiimat3_Message *m, size_t steps);

static void mpz_addmod(mpz_t rop, const mpz_t op1, const mpz_t op2);
static void mpz_mulmod(mpz_t rop, const mpz_t op1, const mpz_t op2);

static void poly_flood(size_t idx, tiimat3_Poly *p, tiimat3_Seed *seed);

static void
block_add(tiimat3_Block *rop, const tiimat3_Block *op1, const tiimat3_Block *op2, size_t len)
{
	size_t i, row, col;

	for (i = 0; i < len; ++i)
		for (row = 0; row < TIIMAT3_DIM; ++row)
			for (col = 0; col < TIIMAT3_DIM; ++col)
				mpz_addmod(rop[i].value[row][col], op1[i].value[row][col], op2[i].value[row][col]);
}

static tiimat3_Block *
block_alloc(size_t len)
{
	tiimat3_Block *block;
	size_t i, row, col;

	block = tiimat3_util_alloc(len, sizeof *block);
	for (i = 0; i < len; ++i)
		for (row = 0; row < TIIMAT3_DIM; ++row)
			for (col = 0; col < TIIMAT3_DIM; ++col)
				mpz_init(block[i].value[row][col]);

	return block;
}

static int
block_cmp(const tiimat3_Block *op1, const tiimat3_Block *op2, size_t len)
{
	size_t i, row, col;

	for (i = 0; i < len; ++i)
		for (row = 0; row < TIIMAT3_DIM; ++row)
			for (col = 0; col < TIIMAT3_DIM; ++col)
				if (mpz_cmp(op1[i].value[row][col], op2[i].value[row][col]) != 0)
					return 1;

	return 0;
}

static void
block_cpy(tiimat3_Block *rop, const tiimat3_Block *op)
{
	size_t row, col;

	for (row = 0; row < TIIMAT3_DIM; ++row)
		for (col = 0; col < TIIMAT3_DIM; ++col)
			mpz_set(rop->value[row][col], op->value[row][col]);
}

static void
block_dealloc(tiimat3_Block *op, size_t len)
{
	size_t i, row, col;

	for (i = 0; i < len; ++i)
		for (row = 0; row < TIIMAT3_DIM; ++row)
			for (col = 0; col < TIIMAT3_DIM; ++col)
				mpz_clear(op[i].value[row][col]);
	tiimat3_util_dealloc(op);
}

static void
block_decrypt(tiimat3_Block *rop, tiimat3_BlockEnc *op)
{
	tiimat3_Seed seed;
	tiimat3_Message *m;
	tiimat3_Poly *p;
	size_t i;

	m = tiimat3_msg_alloc(1);
	p = tiimat3_util_alloc(TIIMAT3_QLEN, sizeof *p);
	tiimat3_random_seed(&seed);

	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		if (tiimat3_q[i].drop)
			continue;

		tiimat3_bgv_decrypt(i, &p[i], &sk, &op->value[i]);
		poly_flood(i, &p[i], &seed);
		seed.init = 1;
	}
	tiimat3_bgv_decode(m, p);
	block_unpack(rop, m);

	tiimat3_msg_dealloc(m, 1);
	tiimat3_util_dealloc(p);
}

static void
block_encrypt(tiimat3_BlockEnc *rop, const tiimat3_Block *blocks, size_t rot, void (*transform)(tiimat3_Block *, const tiimat3_Block *))
{
	tiimat3_Seed seed;
	tiimat3_Message *m;
	tiimat3_Poly *p;
	size_t i;

	m = tiimat3_msg_alloc(1);
	p = tiimat3_util_alloc(1, sizeof *p);
	tiimat3_random_seed(&seed);

	block_pack(m, blocks, rot, transform);
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		tiimat3_bgv_encode(i, p, m);
		tiimat3_bgv_encrypt(i, &rop->value[i], &pk[i], p, &seed);
	}

	tiimat3_msg_dealloc(m, 1);
	tiimat3_util_dealloc(p);
}

static void
block_init_rand(tiimat3_Block *rop, size_t len)
{
	size_t i, row, col;

	for (i = 0; i < len; ++i)
		for (row = 0; row < TIIMAT3_DIM; ++row)
			for (col = 0; col < TIIMAT3_DIM; ++col)
				mpz_urandomm(rop[i].value[row][col], randstate, tiimat3_t.value);
}

static void
block_init_zero(tiimat3_Block *rop, size_t len)
{
	size_t i, row, col;

	for (i = 0; i < len; ++i)
		for (row = 0; row < TIIMAT3_DIM; ++row)
			for (col = 0; col < TIIMAT3_DIM; ++col)
				mpz_set_ui(rop[i].value[row][col], 0);
}

void
block_mul(tiimat3_Block *rop, const tiimat3_Block *op1, const tiimat3_Block *op2, size_t len)
{
	tiimat3_Block *cpy;
	size_t i, row, col, k;
	mpz_t tmp;

	cpy = block_alloc(1);
	mpz_init(tmp);

	for (i = 0; i < len; ++i) {
		for (row = 0; row < TIIMAT3_DIM; ++row) {
			for (col = 0; col < TIIMAT3_DIM; ++col) {
				mpz_set_ui(cpy->value[row][col], 0);
				for (k = 0; k < TIIMAT3_DIM; ++k) {
					mpz_mulmod(tmp, op1[i].value[row][k], op2[i].value[k][col]);
					mpz_addmod(cpy->value[row][col], cpy->value[row][col], tmp);
				}
			}
		}

		block_cpy(&rop[i], cpy);
	}

	block_dealloc(cpy, 1);
	mpz_clear(tmp);
}

static void
block_mulc(tiimat3_Block *rop, const tiimat3_Block *op1, const mpz_t op2, size_t len)
{
	size_t i, row, col;

	for (i = 0; i < len; ++i)
		for (row = 0; row < TIIMAT3_DIM; ++row)
			for (col = 0; col < TIIMAT3_DIM; ++col)
				mpz_mulmod(rop[i].value[row][col], op1[i].value[row][col], op2);
}

static void
block_pack(tiimat3_Message *m, const tiimat3_Block *op, size_t rot, void (*transform)(tiimat3_Block *, const tiimat3_Block *))
{
	tiimat3_Block *cpy;
	size_t row, col, i;

	cpy = block_alloc(1);

	for (i = 0; i < TIIMAT3_PACK / 2; ++i) {
		transform(cpy, &op[i]);

		for (row = 0; row < TIIMAT3_DIM; ++row) {
			for (col = 0; col < TIIMAT3_DIM; ++col) {
				size_t j = (row * TIIMAT3_DIM + col) * TIIMAT3_PACK / 2 + i;
				mpz_set(m->value[j], cpy->value[row][col]);
			}
		}
	}

	for (i = 0; i < TIIMAT3_PACK / 2; ++i) {
		transform(cpy, &op[TIIMAT3_PACK / 2 + i]);

		for (row = 0; row < TIIMAT3_DIM; ++row) {
			for (col = 0; col < TIIMAT3_DIM; ++col) {
				size_t j = TIIMAT3_D / 2 + (row * TIIMAT3_DIM + col) * TIIMAT3_PACK / 2 + i;
				mpz_set(m->value[j], cpy->value[row][col]);
			}
		}
	}

	msg_rot(m, m, rot);
	tiimat3_msg_pack(m, m);

	block_dealloc(cpy, 1);
}

static void
block_sigma(tiimat3_Block *rop, const tiimat3_Block *op)
{
	tiimat3_Block *cpy;
	size_t row, col;

	cpy = block_alloc(1);

	for (row = 0; row < TIIMAT3_DIM; ++row)
		for (col = 0; col < TIIMAT3_DIM; ++col)
			mpz_set(cpy->value[row][col], op->value[row][(row + col) % TIIMAT3_DIM]);
	block_cpy(rop, cpy);

	block_dealloc(cpy, 1);
}

static void
block_tau(tiimat3_Block *rop, const tiimat3_Block *op)
{
	tiimat3_Block *cpy;
	size_t row, col;

	cpy = block_alloc(1);

	for (row = 0; row < TIIMAT3_DIM; ++row)
		for (col = 0; col < TIIMAT3_DIM; ++col)
			mpz_set(cpy->value[row][col], op->value[(row + col) % TIIMAT3_DIM][col]);
	block_cpy(rop, cpy);

	block_dealloc(cpy, 1);
}

static void
block_unpack(tiimat3_Block *vec, const tiimat3_Message *m)
{
	tiimat3_Message *cpy;
	size_t i, row, col;

	cpy = tiimat3_msg_alloc(1);
	tiimat3_msg_unpack(cpy, m);

	for (i = 0; i < TIIMAT3_PACK / 2; ++i) {
		for (row = 0; row < TIIMAT3_DIM; ++row) {
			for (col = 0; col < TIIMAT3_DIM; ++col) {
				size_t j = (row * TIIMAT3_DIM + col) * TIIMAT3_PACK / 2 + i;
				mpz_set(vec[i].value[row][col], cpy->value[j]);
			}
		}
	}

	for (i = 0; i < TIIMAT3_PACK / 2; ++i) {
		for (row = 0; row < TIIMAT3_DIM; ++row) {
			for (col = 0; col < TIIMAT3_DIM; ++col) {
				size_t j = TIIMAT3_D / 2 + (row * TIIMAT3_DIM + col) * TIIMAT3_PACK / 2 + i;
				mpz_set(vec[TIIMAT3_PACK / 2 + i].value[row][col], cpy->value[j]);
			}
		}
	}

	tiimat3_msg_dealloc(cpy, 1);
}

static void
ciphertext_cpy(tiimat3_Ciphertext *rop, tiimat3_Ciphertext *ct, size_t len)
{

	memcpy(rop, ct, len * sizeof *ct);
}

static void
ciphertext_rotA(tiimat3_Ciphertext *rop, tiimat3_Ciphertext *ct)
{
	size_t i;

	if (rop == ct) {
		for (i = 0; i < TIIMAT3_QLEN; ++i)
			tiimat3_bgv_rot_inplace(i, &ct[i], TIIMAT3_ROTA);
	} else {
		for (i = 0; i < TIIMAT3_QLEN; ++i)
			tiimat3_bgv_rot(i, &rop[i], &ct[i], TIIMAT3_ROTA);
	}

	tiimat3_bgv_keyswitch(rop, &kswA[TIIMAT3_QPLEN], 1);
}

static void
ciphertext_rotA_hoisted(tiimat3_Ciphertext *rop, tiimat3_CiphertextSwitch *csw, int k)
{
	tiimat3_CiphertextSwitch *cswr;
	tiimat3_Delta *delta;
	size_t i;

	cswr = tiimat3_util_alloc(1, sizeof *cswr);
	delta = tiimat3_util_alloc(TIIMAT3_PLEN, sizeof *delta);

	if (k < 0) {
		k = -k;
		for (i = 0; i < TIIMAT3_QPLEN; ++i) {
			tiimat3_bgv_rot_csw(i, cswr, &csw[i], TIIMAT3_D / 2 - TIIMAT3_ROTB + k * TIIMAT3_ROTA);
			tiimat3_bgv_keyswitch_dot(i, &rop[i], cswr, &kswAneg[k * TIIMAT3_QPLEN + i]);
		}
	} else {
		for (i = 0; i < TIIMAT3_QPLEN; ++i) {
			tiimat3_bgv_rot_csw(i, cswr, &csw[i], k * TIIMAT3_ROTA);
			tiimat3_bgv_keyswitch_dot(i, &rop[i], cswr, &kswA[k * TIIMAT3_QPLEN + i]);
		}
	}

	for (i = TIIMAT3_QLEN; i < TIIMAT3_QPLEN; ++i)
		tiimat3_bgv_keyswitch_delta(i, &delta[i - TIIMAT3_QLEN], &rop[i]);
	for (i = 0; i < TIIMAT3_QLEN; ++i)
		tiimat3_bgv_keyswitch_switch(i, &rop[i], delta);

	tiimat3_util_dealloc(cswr);
	tiimat3_util_dealloc(delta);
}

static void
ciphertext_rotB(tiimat3_Ciphertext *rop, tiimat3_Ciphertext *ct)
{
	size_t i;

	if (rop == ct) {
		for (i = 0; i < TIIMAT3_QLEN; ++i)
			tiimat3_bgv_rot_inplace(i, &ct[i], TIIMAT3_ROTB);
	} else {
		for (i = 0; i < TIIMAT3_QLEN; ++i)
			tiimat3_bgv_rot(i, &rop[i], &ct[i], TIIMAT3_ROTB);
	}

	tiimat3_bgv_keyswitch(rop, &kswB[TIIMAT3_QPLEN], 1);
}

static void
ciphertext_rotB_hoisted(tiimat3_Ciphertext *rop, tiimat3_CiphertextSwitch *csw, int k)
{
	tiimat3_CiphertextSwitch *cswr;
	tiimat3_Delta *delta;
	size_t i;

	cswr = tiimat3_util_alloc(1, sizeof *cswr);
	delta = tiimat3_util_alloc(TIIMAT3_PLEN, sizeof *delta);


	for (i = 0; i < TIIMAT3_QPLEN; ++i) {
		tiimat3_bgv_rot_csw(i, cswr, &csw[i], k * TIIMAT3_ROTB);
		tiimat3_bgv_keyswitch_dot(i, &rop[i], cswr, &kswB[k * TIIMAT3_QPLEN + i]);
	}

	for (i = TIIMAT3_QLEN; i < TIIMAT3_QPLEN; ++i)
		tiimat3_bgv_keyswitch_delta(i, &delta[i - TIIMAT3_QLEN], &rop[i]);
	for (i = 0; i < TIIMAT3_QLEN; ++i)
		tiimat3_bgv_keyswitch_switch(i, &rop[i], delta);

	tiimat3_util_dealloc(cswr);
	tiimat3_util_dealloc(delta);
}

static void
msg_cpy(tiimat3_Message *rop, const tiimat3_Message *m)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		mpz_set(rop->value[j], m->value[j]);
}

static void
msg_deinit(tiimat3_Message *m)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		mpz_clear(m->value[j]);
}

static void
msg_init(tiimat3_Message *m)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j)
		mpz_init(m->value[j]);
}

static void
msg_rot(tiimat3_Message *rop, const tiimat3_Message *m, size_t steps)
{
	tiimat3_Message *cpy;
	size_t j;

	cpy = tiimat3_msg_alloc(1);
	msg_cpy(cpy, m);

	for (j = 0; j < TIIMAT3_D / 2; ++j) {
		size_t idx = (j + steps) % (TIIMAT3_D / 2);
		mpz_set(rop->value[j], cpy->value[idx]);
	}
	for (; j < TIIMAT3_D; ++j) {
		size_t idx = TIIMAT3_D / 2 + ((j + steps) % (TIIMAT3_D / 2));
		mpz_set(rop->value[j], cpy->value[idx]);
	}

	tiimat3_msg_dealloc(cpy, 1);
}

static void
mpz_addmod(mpz_t rop, const mpz_t op1, const mpz_t op2)
{

	mpz_add(rop, op1, op2);
	mpz_mod(rop, rop, tiimat3_t.value);
}

static void
mpz_mulmod(mpz_t rop, const mpz_t op1, const mpz_t op2)
{

	mpz_mul(rop, op1, op2);
	mpz_mod(rop, rop, tiimat3_t.value);
}

static void
poly_flood(size_t idx, tiimat3_Poly *p, tiimat3_Seed *seed)
{
	unsigned char rand[TIIMAT3_FLOOD_BYTES];
	tiimat3_Poly *e;
	mpz_t tmp;
	size_t j;

	e = tiimat3_util_alloc(1, sizeof *e);
	mpz_init(tmp);

	for (j = 0; j < TIIMAT3_D; ++j) {
		tiimat3_random_bytes(rand, sizeof rand, seed);
		rand[0] >>= (CHAR_BIT - (TIIMAT3_FLOOD_BITS % CHAR_BIT)) % CHAR_BIT;
		mpz_import(tmp, sizeof rand, 1, 1, 1, 0, rand);

		mpz_mod_ui(tmp, tmp, tiimat3_q[idx].value);
		e->value[j] = mpz_get_ui(tmp);
	}
	tiimat3_poly_init_error(idx, e, e);
	tiimat3_poly_add(idx, p, p, e);

	tiimat3_util_dealloc(e);
	mpz_clear(tmp);
}

void
tiimat3_block_add_enc(tiimat3_BlockEnc *rop, tiimat3_BlockEnc *op1, tiimat3_BlockEnc *op2)
{
	size_t i;

	for (i = 0; i < TIIMAT3_QLEN; ++i)
		tiimat3_bgv_add(i, &rop->value[i], &op1->value[i], &op2->value[i]);
}
void
tiimat3_block_mul_enc_hoisted(tiimat3_BlockEnc *rop, tiimat3_BlockEnc *A, tiimat3_BlockEnc *B)
{
	tiimat3_CiphertextSwitch *cswA, *cswB;
	tiimat3_Ciphertext *rotA, *rotAneg, *rotB, *ct;
	tiimat3_Delta *delta;
	int i, k;

	cswA = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *cswA);
	cswB = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *cswB);
	rotA = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *rotA);
	rotAneg = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *rotAneg);
	rotB = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *rotB);
	ct = tiimat3_util_alloc(2, sizeof *ct);
	delta = tiimat3_util_alloc(2, sizeof *delta);

	/* reset dropped moduli */
	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	/* ciphertext extension */
	tiimat3_bgv_keyswitch_ext(cswA, A->value, 1);
	tiimat3_bgv_keyswitch_ext(cswB, B->value, 1);

	/* first multiplication */
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		ciphertext_cpy(&ct[0], &A->value[i], 1);
		ciphertext_cpy(&ct[1], &B->value[i], 1);

		if (i == 0) {
			tiimat3_bgv_modswitch_delta(0, &delta[0], &ct[0]);
			tiimat3_bgv_modswitch_delta(0, &delta[1], &ct[1]);
			continue;
		} else {
			tiimat3_bgv_modswitch_ext(i, &ct[0], &delta[0]);
			tiimat3_bgv_modswitch_ext(i, &ct[1], &delta[1]);
		}

		tiimat3_bgv_mul(i, &rop->value[i], &ct[0], &ct[1]);
	}

	/* squared block multiplication */
	for (k = 1; k < TIIMAT3_DIM; ++k) {
		ciphertext_rotA_hoisted(rotA, cswA, k);
		ciphertext_rotA_hoisted(rotAneg, cswA, -k);
		ciphertext_rotB_hoisted(rotB, cswB, k);

		/* multiplication */
		for (i = 0; i < TIIMAT3_QLEN; ++i) {
			tiimat3_Poly p;

			tiimat3_bgv_encode(i, &p, &Uphi[k][0]);
			tiimat3_bgv_mulc(i, &rotA[i], &rotA[i], &p);

			tiimat3_bgv_encode(i, &p, &Uphi[k][1]);
			tiimat3_bgv_mulc(i, &rotAneg[i], &rotAneg[i], &p);

			tiimat3_bgv_add(i, &rotA[i], &rotA[i], &rotAneg[i]);

			if (i == 0) {
				tiimat3_bgv_modswitch_delta(0, &delta[0], &rotA[0]);
				tiimat3_bgv_modswitch_delta(0, &delta[1], &rotB[0]);
				continue;
			} else {
				tiimat3_bgv_modswitch_ext(i, &rotA[i], &delta[0]);
				tiimat3_bgv_modswitch_ext(i, &rotB[i], &delta[1]);
			}

			tiimat3_bgv_mul(i, ct, &rotA[i], &rotB[i]);
			tiimat3_bgv_add(i, &rop->value[i], &rop->value[i], ct);
		}
	}
	tiimat3_mod_drop(0, 2);

	tiimat3_bgv_keyswitch(rop->value, ksw2, 2);

	for (i = 1; i < TIIMAT3_QLEN; ++i) {
		tiimat3_Delta delta;

		if (i == 1) {
			tiimat3_bgv_modswitch_delta(1, &delta, &rop->value[1]);
			continue;
		} else
			tiimat3_bgv_modswitch_ext(i, &rop->value[i], &delta);
	}
	tiimat3_mod_drop(1, 1);

	tiimat3_util_dealloc(cswA);
	tiimat3_util_dealloc(cswB);
	tiimat3_util_dealloc(rotA);
	tiimat3_util_dealloc(rotAneg);
	tiimat3_util_dealloc(rotB);
	tiimat3_util_dealloc(ct);
	tiimat3_util_dealloc(delta);
}

void
tiimat3_block_mul_enc_prerot(tiimat3_BlockEnc *rop, tiimat3_BlockEnc *A1, tiimat3_BlockEnc *A2, tiimat3_BlockEnc *B)
{
	tiimat3_BlockEnc *mul;
	size_t i, k;

	mul = tiimat3_util_alloc(TIIMAT3_DIM, sizeof *mul);

	/* reset dropped moduli */
	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	/* squared block multiplication */
	#pragma omp parallel for
	for (k = 0; k < TIIMAT3_DIM; ++k) {
		tiimat3_Ciphertext *ct;
		tiimat3_Delta *delta;
		size_t i;

		ct = tiimat3_util_alloc(2, sizeof *ct);
		delta = tiimat3_util_alloc(2, sizeof *delta);

		if (k == 0) {
			for (i = 0; i < TIIMAT3_QLEN; ++i) {
				ciphertext_cpy(&ct[0], &A2[0].value[i], 1);
				ciphertext_cpy(&ct[1], &B[0].value[i], 1);

				if (i == 0) {
					tiimat3_bgv_modswitch_delta(0, &delta[0], &ct[0]);
					tiimat3_bgv_modswitch_delta(0, &delta[1], &ct[1]);
				} else {
					tiimat3_bgv_modswitch_ext(i, &ct[0], &delta[0]);
					tiimat3_bgv_modswitch_ext(i, &ct[1], &delta[1]);
					tiimat3_bgv_mul(i, &mul[k].value[i], &ct[0], &ct[1]);
				}
			}
		} else {
			for (i = 0; i < TIIMAT3_QLEN; ++i) {
				tiimat3_Poly p;

				tiimat3_bgv_encode(i, &p, &Uphi[k][0]);
				tiimat3_bgv_mulc(i, &ct[0], &A2[k].value[i], &p);

				tiimat3_bgv_encode(i, &p, &Uphi[k][1]);
				tiimat3_bgv_mulc(i, &ct[1], &A1[k].value[i], &p);

				tiimat3_bgv_add(i, &ct[0], &ct[0], &ct[1]);
				ciphertext_cpy(&ct[1], &B[k].value[i], 1);

				if (i == 0) {
					tiimat3_bgv_modswitch_delta(0, &delta[0], &ct[0]);
					tiimat3_bgv_modswitch_delta(0, &delta[1], &ct[1]);
				} else {
					tiimat3_bgv_modswitch_ext(i, &ct[0], &delta[0]);
					tiimat3_bgv_modswitch_ext(i, &ct[1], &delta[1]);
					tiimat3_bgv_mul(i, &mul[k].value[i], &ct[0], &ct[1]);
				}
			}
		}

		tiimat3_util_dealloc(ct);
		tiimat3_util_dealloc(delta);
	}

	for (k = 1; k < TIIMAT3_DIM; ++k) {
		for (i = 0; i < TIIMAT3_QLEN; ++i) {
			if (k == 1)
				tiimat3_bgv_add(i, &rop->value[i], &mul[0].value[i], &mul[1].value[i]);
			else
				tiimat3_bgv_add(i, &rop->value[i], &rop->value[i], &mul[k].value[i]);
		}
	}
	tiimat3_mod_drop(0, 2);
	tiimat3_bgv_keyswitch(rop->value, ksw2, 2);

	for (i = 1; i < TIIMAT3_QLEN; ++i) {
		tiimat3_Delta delta;

		if (i == 1) {
			tiimat3_bgv_modswitch_delta(1, &delta, &rop->value[1]);
			continue;
		} else
			tiimat3_bgv_modswitch_ext(i, &rop->value[i], &delta);
	}
	tiimat3_mod_drop(1, 1);

	tiimat3_util_dealloc(mul);
}

void
tiimat3_block_mul_enc_reuse(tiimat3_BlockEnc *rop, tiimat3_BlockEnc *A, tiimat3_BlockEnc *B)
{
	tiimat3_Ciphertext *rotA1, *rotA2, *rotB, *ct;
	tiimat3_Delta *delta;
	size_t i, k;

	delta = tiimat3_util_alloc(2, sizeof *delta);
	rotA1 = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *rotA1);
	rotA2 = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *rotA2);
	rotB = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *rotB);
	ct = tiimat3_util_alloc(2, sizeof *ct);

	/* reset dropped moduli */
	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	/* initial rotation of A2 */
	ciphertext_rotB(rotA2, A->value);
	ciphertext_cpy(rotA1, A->value, TIIMAT3_QPLEN);
	ciphertext_cpy(rotB, B->value, TIIMAT3_QPLEN);

	/* first multiplication */
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		ciphertext_cpy(&ct[0], &rotA2[i], 1);
		ciphertext_cpy(&ct[1], &rotB[i], 1);

		if (i == 0) {
			tiimat3_bgv_modswitch_delta(0, &delta[0], &ct[0]);
			tiimat3_bgv_modswitch_delta(0, &delta[1], &ct[1]);
			continue;
		} else {
			tiimat3_bgv_modswitch_ext(i, &ct[0], &delta[0]);
			tiimat3_bgv_modswitch_ext(i, &ct[1], &delta[1]);
		}

		tiimat3_bgv_mul(i, &rop->value[i], &ct[0], &ct[1]);
	}

	/* squared block multiplication */
	for (k = 1; k < TIIMAT3_DIM; ++k) {
		ciphertext_rotA(rotA1, rotA1);
		ciphertext_rotA(rotA2, rotA2);
		ciphertext_rotB(rotB, rotB);

		/* multiplication */
		for (i = 0; i < TIIMAT3_QLEN; ++i) {
			tiimat3_Poly p;

			tiimat3_bgv_encode(i, &p, &Uphi[k][0]);
			tiimat3_bgv_mulc(i, &ct[0], &rotA2[i], &p);

			tiimat3_bgv_encode(i, &p, &Uphi[k][1]);
			tiimat3_bgv_mulc(i, &ct[1], &rotA1[i], &p);

			tiimat3_bgv_add(i, &ct[0], &ct[0], &ct[1]);
			ciphertext_cpy(&ct[1], &rotB[i], 1);

			if (i == 0) {
				tiimat3_bgv_modswitch_delta(0, &delta[0], &ct[0]);
				tiimat3_bgv_modswitch_delta(0, &delta[1], &ct[1]);
				continue;
			} else {
				tiimat3_bgv_modswitch_ext(i, &ct[0], &delta[0]);
				tiimat3_bgv_modswitch_ext(i, &ct[1], &delta[1]);
			}

			tiimat3_bgv_mul(i, &ct[0], &ct[0], &ct[1]);
			tiimat3_bgv_add(i, &rop->value[i], &rop->value[i], &ct[0]);
		}
	}
	tiimat3_mod_drop(0, 2);

	tiimat3_bgv_keyswitch(rop->value, ksw2, 2);

	for (i = 1; i < TIIMAT3_QLEN; ++i) {
		tiimat3_Delta delta;

		if (i == 1) {
			tiimat3_bgv_modswitch_delta(1, &delta, &rop->value[1]);
			continue;
		} else
			tiimat3_bgv_modswitch_ext(i, &rop->value[i], &delta);
	}
	tiimat3_mod_drop(1, 1);

	tiimat3_util_dealloc(delta);
	tiimat3_util_dealloc(rotA1);
	tiimat3_util_dealloc(rotA2);
	tiimat3_util_dealloc(rotB);
	tiimat3_util_dealloc(ct);
}


void
tiimat3_block_prerotA_enc_hoisted(tiimat3_BlockEnc *rotA1, tiimat3_BlockEnc *rotA2, tiimat3_BlockEnc *A)
{
	tiimat3_CiphertextSwitch *csw;
	size_t k;

	csw = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *csw);
	tiimat3_bgv_keyswitch_ext(csw, A->value, 1);

	#pragma omp parallel for
	for (k = 0; k < TIIMAT3_DIM; ++k) {
		if (k == 0) {
			ciphertext_cpy(rotA2[0].value, A->value, TIIMAT3_QPLEN);
		} else {
			ciphertext_rotA_hoisted(rotA1[k].value, csw, -k);
			ciphertext_rotA_hoisted(rotA2[k].value, csw, k);
		}
	}

	tiimat3_util_dealloc(csw);
}

void
tiimat3_block_prerotA_enc_mixed(tiimat3_BlockEnc *rotA1, tiimat3_BlockEnc *rotA2, tiimat3_BlockEnc *A)
{
	tiimat3_CiphertextSwitch *csw;
	size_t k1, k2;

	csw = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *csw);
	ciphertext_cpy(rotA2[0].value, A->value, TIIMAT3_QPLEN);

	for (k1 = 0; k1 < TIIMAT3_DIM / TIIMAT3_THREADS; ++k1) {
		tiimat3_bgv_keyswitch_ext(csw, rotA2[k1 * TIIMAT3_THREADS].value, 1);

		#pragma omp parallel for
		for (k2 = 1; k2 <= TIIMAT3_THREADS; ++k2) {
			size_t k = k1 * TIIMAT3_THREADS + k2;

			if (k < TIIMAT3_DIM) {
				ciphertext_rotA_hoisted(rotA1[k].value, csw, -k2);
				ciphertext_rotA_hoisted(rotA2[k].value, csw, k2);
			}
		}
	}

	tiimat3_util_dealloc(csw);
}


void
tiimat3_block_prerotA_enc_reuse(tiimat3_BlockEnc *rotA1, tiimat3_BlockEnc *rotA2, tiimat3_BlockEnc *A)
{
	size_t k;

	ciphertext_rotB(rotA2[0].value, A->value);
	ciphertext_cpy(rotA1[0].value, A->value, TIIMAT3_QPLEN);

	for (k = 1; k < TIIMAT3_DIM; ++k) {
		ciphertext_rotA(rotA1[k].value, rotA1[k - 1].value);
		ciphertext_rotA(rotA2[k].value, rotA2[k - 1].value);
	}
}

void
tiimat3_block_prerotB_enc_hoisted(tiimat3_BlockEnc *rotB, tiimat3_BlockEnc *B)
{
	tiimat3_CiphertextSwitch *csw;
	size_t k;

	csw = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *csw);
	tiimat3_bgv_keyswitch_ext(csw, B->value, 1);

	#pragma omp parallel for
	for (k = 0; k < TIIMAT3_DIM; ++k)
		if (k == 0)
			ciphertext_cpy(rotB[0].value, B->value, TIIMAT3_QPLEN);
		else
			ciphertext_rotB_hoisted(rotB[k].value, csw, k);

	tiimat3_util_dealloc(csw);
}

void
tiimat3_block_prerotB_enc_mixed(tiimat3_BlockEnc *rotB, tiimat3_BlockEnc *B)
{
	tiimat3_CiphertextSwitch *csw;
	size_t k1, k2;

	csw = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *csw);
	ciphertext_cpy(rotB[0].value, B->value, TIIMAT3_QPLEN);

	for (k1 = 0; k1 < TIIMAT3_DIM / TIIMAT3_THREADS; ++k1) {
		tiimat3_bgv_keyswitch_ext(csw, rotB[k1 * TIIMAT3_THREADS].value, 1);

		#pragma omp parallel for
		for (k2 = 1; k2 <= TIIMAT3_THREADS; ++k2) {
			size_t k = k1 * TIIMAT3_THREADS + k2;

			if (k < TIIMAT3_DIM) {
				ciphertext_rotB_hoisted(rotB[k].value, csw, k2);
			}
		}
	}

	tiimat3_util_dealloc(csw);
}


void
tiimat3_block_prerotB_enc_reuse(tiimat3_BlockEnc *rotB, tiimat3_BlockEnc *B)
{
	size_t k;

	ciphertext_cpy(rotB[0].value, B->value, TIIMAT3_QPLEN);

	for (k = 1; k < TIIMAT3_DIM; ++k)
		ciphertext_rotB(rotB[k].value, rotB[k - 1].value);
}


void
tiimat3_deinit(void)
{
	size_t i, k;

	tiimat3_bgv_deinit();
	gmp_randclear(randstate);

	tiimat3_util_dealloc(kswA);
	tiimat3_util_dealloc(kswAneg);
	tiimat3_util_dealloc(kswB);

	for (i = 0; i <= TIIMAT3_SHARES; ++i)
		mpz_clear(MAC[i]);
	for (k = 1; k < TIIMAT3_DIM; ++k) {
		msg_deinit(&Uphi[k][0]);
		msg_deinit(&Uphi[k][1]);
	}
}

void
tiimat3_deinit_hoisted(void)
{

	;
}

void
tiimat3_deinit_mixed(void)
{

	;
}

void
tiimat3_init(void)
{
	tiimat3_Seed seed[TIIMAT3_QPLEN];
	tiimat3_Message *m;
	tiimat3_Poly p;
	size_t i, j, k;

	tiimat3_bgv_init();
	gmp_randinit_default(randstate);

	#if TIIMAT3_HOIST
	kswA = tiimat3_util_alloc(TIIMAT3_DIM * TIIMAT3_QPLEN, sizeof *kswA);
	kswAneg = tiimat3_util_alloc(TIIMAT3_DIM * TIIMAT3_QPLEN, sizeof *kswAneg);
	kswB = tiimat3_util_alloc(TIIMAT3_DIM * TIIMAT3_QPLEN, sizeof *kswB);
	#elif TIIMAT3_MIXED
	kswA = tiimat3_util_alloc((TIIMAT3_THREADS + 1) * TIIMAT3_QPLEN, sizeof *kswA);
	kswAneg = tiimat3_util_alloc((TIIMAT3_THREADS + 1) * TIIMAT3_QPLEN, sizeof *kswAneg);
	kswB = tiimat3_util_alloc((TIIMAT3_THREADS + 1) * TIIMAT3_QPLEN, sizeof *kswB);
	#else
	kswA = tiimat3_util_alloc(2 * TIIMAT3_QPLEN, sizeof *kswA);
	kswAneg = tiimat3_util_alloc(2 * TIIMAT3_QPLEN, sizeof *kswAneg);
	kswB = tiimat3_util_alloc(2 * TIIMAT3_QPLEN, sizeof *kswB);
	#endif /* TIIMAT3_HOIST */

	for (i = 0; i <= TIIMAT3_SHARES; ++i)
		mpz_init(MAC[i]);
	for (k = 1; k < TIIMAT3_DIM; ++k) {
		msg_init(&Uphi[k][0]);
		msg_init(&Uphi[k][1]);
	}

	/* BGV keys */
	tiimat3_bgv_keygen_secret(&sk);
	tiimat3_random_seed(seed);
	for (i = 0; i < TIIMAT3_QLEN; ++i)
		tiimat3_bgv_keygen_public(i, &pk[i], &sk, seed);

	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_random_seed(&seed[i]);
	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_bgv_keygen_switchr(i, &kswA[TIIMAT3_QPLEN + i], &sk, TIIMAT3_ROTA, seed);

	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_random_seed(&seed[i]);
	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_bgv_keygen_switchr(i, &kswB[TIIMAT3_QPLEN + i], &sk, TIIMAT3_ROTB, seed);

	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_random_seed(&seed[i]);
	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_bgv_keygen_switch2(i, &ksw2[i], &sk, seed);

	/* MAC key */
	m = tiimat3_msg_alloc(1);
	mpz_set_ui(MAC[0], 0);
	memset(ctMAC[0], 0, sizeof ctMAC[0]);
	for (i = 1; i <= TIIMAT3_SHARES; ++i) {
		mpz_urandomm(MAC[i], randstate, tiimat3_t.value);
		mpz_addmod(MAC[0], MAC[0], MAC[i]);

		for (j = 0; j < TIIMAT3_D; ++j)
			mpz_set(m->value[j], MAC[i]);
		tiimat3_msg_pack(m, m);

		tiimat3_random_seed(seed);
		for (j = 0; j < TIIMAT3_QLEN; ++j) {
			tiimat3_bgv_encode(j, &p, m);
			tiimat3_bgv_encrypt(j, &ctMAC[i][j], &pk[j], &p, seed);
			tiimat3_bgv_add(j, &ctMAC[0][j], &ctMAC[0][j], &ctMAC[i][j]);
		}
	}
	tiimat3_msg_dealloc(m, 1);

	/* Uphi */
	for (k = 1; k < TIIMAT3_DIM; ++k) {
		for (j = 0; j < TIIMAT3_D; ++j)
			if (j % (TIIMAT3_PACK / 2 * TIIMAT3_DIM) < TIIMAT3_PACK / 2 * (TIIMAT3_DIM - k))
				mpz_set_ui(Uphi[k][0].value[j], 1);
			else
				mpz_set_ui(Uphi[k][0].value[j], 0);
		tiimat3_msg_pack(&Uphi[k][0], &Uphi[k][0]);

		for (j = 0; j < TIIMAT3_D; ++j)
			if (j % (TIIMAT3_PACK / 2 * TIIMAT3_DIM) < TIIMAT3_PACK / 2 * (TIIMAT3_DIM - k))
				mpz_set_ui(Uphi[k][1].value[j], 0);
			else
				mpz_set_ui(Uphi[k][1].value[j], 1);
		tiimat3_msg_pack(&Uphi[k][1], &Uphi[k][1]);
	}
}

void
tiimat3_init_hoisted(void)
{
	tiimat3_Seed seed[TIIMAT3_QPLEN];
	size_t i, k;

	for (k = 2; k < TIIMAT3_DIM; ++k) {
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_random_seed(&seed[i]);
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_bgv_keygen_switchr(i, &kswA[k * TIIMAT3_QPLEN + i], &sk, k * TIIMAT3_ROTA, seed);
	}

	for (k = 1; k < TIIMAT3_DIM; ++k) {
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_random_seed(&seed[i]);
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_bgv_keygen_switchr(i, &kswAneg[k * TIIMAT3_QPLEN + i], &sk, TIIMAT3_D / 2 - TIIMAT3_ROTB + k * TIIMAT3_ROTA, seed);
	}

	for (k = 2; k < TIIMAT3_DIM; ++k) {
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_random_seed(&seed[i]);
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_bgv_keygen_switchr(i, &kswB[k * TIIMAT3_QPLEN + i], &sk, k * TIIMAT3_ROTB, seed);
	}
}

void
tiimat3_init_mixed(void)
{
	tiimat3_Seed seed[TIIMAT3_QPLEN];
	size_t i, k;

	for (k = 2; k <= TIIMAT3_THREADS; ++k) {
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_random_seed(&seed[i]);
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_bgv_keygen_switchr(i, &kswA[k * TIIMAT3_QPLEN + i], &sk, k * TIIMAT3_ROTA, seed);
	}

	for (k = 1; k <= TIIMAT3_THREADS; ++k) {
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_random_seed(&seed[i]);
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_bgv_keygen_switchr(i, &kswAneg[k * TIIMAT3_QPLEN + i], &sk, TIIMAT3_D / 2 - TIIMAT3_ROTB + k * TIIMAT3_ROTA, seed);
	}

	for (k = 2; k <= TIIMAT3_THREADS; ++k) {
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_random_seed(&seed[i]);
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_bgv_keygen_switchr(i, &kswB[k * TIIMAT3_QPLEN + i], &sk, k * TIIMAT3_ROTB, seed);
	}
}

void
tiimat3_matrix_add(tiimat3_Matrix *rop, const tiimat3_Matrix *op1, const tiimat3_Matrix *op2)
{
	size_t row, col;

	assert(op1->drow == op2->drow);
	assert(op1->dcol == op2->dcol);

	for (row = 0; row < op1->drow; ++row) {
		for (col = 0; col < op1->dcol; ++col) {
			size_t idx = (row * op1->dcol + col) * TIIMAT3_PACK;
			block_add(&rop->value[idx], &op1->value[idx], &op2->value[idx], TIIMAT3_PACK);
		}
	}

}

void
tiimat3_matrix_add_enc(tiimat3_MatrixEnc *rop, tiimat3_MatrixEnc *op1, tiimat3_MatrixEnc *op2)
{
	size_t row, col;

	assert(op1->drow == op2->drow);
	assert(op1->dcol == op2->dcol);

	for (row = 0; row < op1->drow; ++row) {
		for (col = 0; col < op1->dcol; ++col) {
			size_t idx = row * op1->dcol + col;
			tiimat3_block_add_enc(&rop->value[idx], &op1->value[idx], &op2->value[idx]);
		}
	}

}

tiimat3_Matrix *
tiimat3_matrix_alloc(size_t drow, size_t dcol, size_t len)
{
	tiimat3_Matrix *mat;
	size_t i;

	mat = tiimat3_util_alloc(len, sizeof *mat);
	for (i = 0; i < len; ++i) {
		mat[i].value = block_alloc(drow * dcol * TIIMAT3_PACK);
		mat[i].drow = drow;
		mat[i].dcol = dcol;
	}

	return mat;
}

tiimat3_MatrixEnc *
tiimat3_matrix_alloc_enc(size_t drow, size_t dcol, size_t len)
{
	tiimat3_MatrixEnc *mat;
	size_t i;

	mat = tiimat3_util_alloc(len, sizeof *mat);
	for (i = 0; i < len; ++i) {
		mat[i].value = tiimat3_util_alloc(drow * dcol, sizeof *mat[i].value);
		mat[i].drow = drow;
		mat[i].dcol = dcol;
	}

	return mat;
}

tiimat3_MatrixEnc *
tiimat3_matrix_alloc_encdim(size_t drow, size_t dcol, size_t len)
{
	tiimat3_MatrixEnc *mat;
	size_t i;

	mat = tiimat3_util_alloc(len, sizeof *mat);
	for (i = 0; i < len; ++i) {
		mat[i].value = tiimat3_util_alloc(drow * dcol * TIIMAT3_DIM, sizeof *mat[i].value);
		mat[i].drow = drow;
		mat[i].dcol = dcol;
	}

	return mat;
}

int
tiimat3_matrix_cmp(const tiimat3_Matrix *op1, const tiimat3_Matrix *op2)
{
	size_t row, col;
	int cmp;

	assert(op1->drow == op2->drow);
	assert(op1->dcol == op2->dcol);

	cmp = 0;
	for (row = 0; row < op1->drow; ++row) {
		for (col = 0; col < op1->dcol; ++col) {
			size_t idx = (row * op1->dcol + col) * TIIMAT3_PACK;
			cmp += block_cmp(&op1->value[idx], &op2->value[idx], TIIMAT3_PACK);
		}
	}

	return cmp;
}

void
tiimat3_matrix_dealloc(tiimat3_Matrix *op, size_t len)
{
	size_t i;

	for (i = 0; i < len; ++i)
		block_dealloc(op[i].value, op->drow * op->dcol * TIIMAT3_PACK);
	tiimat3_util_dealloc(op);
}

void
tiimat3_matrix_dealloc_enc(tiimat3_MatrixEnc *op, size_t len)
{
	size_t i;

	for (i = 0; i < len; ++i)
		tiimat3_util_dealloc(op[i].value);
	tiimat3_util_dealloc(op);
}

void
tiimat3_matrix_decrypt(tiimat3_Matrix *rop, tiimat3_MatrixEnc *op)
{
	size_t row, col;

	for (row = 0; row < op->drow; ++row) {
		for (col = 0; col < op->dcol; ++col) {
			size_t idx = row * op->dcol + col;
			block_decrypt(&rop->value[idx * TIIMAT3_PACK], &op->value[idx]);
		}
	}
}

void
tiimat3_matrix_encrypt(tiimat3_MatrixEnc *rop, const tiimat3_Matrix *op)
{
	size_t row, col;

	for (row = 0; row < op->drow; ++row) {
		for (col = 0; col < op->dcol; ++col) {
			size_t idx = row * op->dcol + col;
			block_encrypt(&rop->value[idx], &op->value[idx * TIIMAT3_PACK], 0, block_cpy);
		}
	}
}

void
tiimat3_matrix_encryptA(tiimat3_MatrixEnc *rop, const tiimat3_Matrix *op)
{
	size_t row, col;

	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	for (row = 0; row < op->drow; ++row) {
		for (col = 0; col < op->dcol; ++col) {
			size_t idx = row * op->dcol + col;
			#if TIIMAT3_HOIST || TIIMAT3_MIXED
			block_encrypt(&rop->value[idx], &op->value[idx * TIIMAT3_PACK], 0, block_sigma);
			#else
			block_encrypt(&rop->value[idx], &op->value[idx * TIIMAT3_PACK], TIIMAT3_D / 2 - TIIMAT3_ROTB, block_sigma);
			#endif /* TIIMAT3_HOIST */
		}
	}
}

void
tiimat3_matrix_encryptB(tiimat3_MatrixEnc *rop, const tiimat3_Matrix *op)
{
	size_t row, col;

	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	for (row = 0; row < op->drow; ++row) {
		for (col = 0; col < op->dcol; ++col) {
			size_t idx = row * op->dcol + col;
			block_encrypt(&rop->value[idx], &op->value[idx * TIIMAT3_PACK], 0, block_tau);
		}
	}
}

void
tiimat3_matrix_init_asc(tiimat3_Matrix *op, size_t len)
{
	size_t mrow, mcol, brow, bcol, i, ii;

	for (i = 0; i < len; ++i) {
		for (mrow = 0; mrow < op[i].drow; ++mrow) {
			for (mcol = 0; mcol < op[i].dcol; ++mcol) {
				for (ii = 0; ii < TIIMAT3_PACK; ++ii) {
					for (brow = 0; brow < TIIMAT3_DIM; ++brow) {
						for (bcol = 0; bcol < TIIMAT3_DIM; ++bcol) {
							size_t idx, val;

							idx = (mrow * op[i].dcol + mcol) * TIIMAT3_PACK + ii;
							val = (mrow * op[i].dcol + brow) * op[i].dcol * TIIMAT3_DIM + mcol * TIIMAT3_DIM + bcol;
							mpz_set_ui(op[i].value[idx].value[brow][bcol], val);
						}
					}
				}
			}
		}
	}
}

void
tiimat3_matrix_init_rand(tiimat3_Matrix *op, size_t len)
{
	size_t i;

	for (i = 0; i < len; ++i)
		block_init_rand(op[i].value, op[i].drow * op[i].dcol * TIIMAT3_PACK);
}

void
tiimat3_matrix_init_zero(tiimat3_Matrix *op, size_t len)
{
	size_t i;

	for (i = 0; i < len; ++i)
		block_init_zero(op[i].value, op[i].drow * op[i].dcol * TIIMAT3_PACK);
}

void
tiimat3_matrix_init_zero_enc(tiimat3_MatrixEnc *op, size_t len)
{
	size_t row, col, i;

	for (i = 0; i < len; ++i) {
		for (row = 0; row < op[i].drow; ++row) {
			for (col = 0; col < op[i].dcol; ++col) {
				size_t idx;

				idx = row * op[i].dcol + col;
				memset(&op[i].value[idx], 0, sizeof op[i].value[idx]);
			}
		}
	}
}

void
tiimat3_matrix_mac(tiimat3_MatrixEnc *rop, tiimat3_MatrixEnc *op)
{
	size_t row, col, i;

	assert(rop->drow == op->drow);
	assert(rop->dcol == op->dcol);

	for (row = 0; row < rop->drow; ++row) {
		for (col = 0; col < rop->dcol; ++col) {
			for (i = 0; i < TIIMAT3_QLEN; ++i) {
				size_t idx;

				if (tiimat3_q[i].drop)
					continue;

				idx = row * rop->dcol + col;
				tiimat3_bgv_mul(i, &rop->value[idx].value[i], &op->value[idx].value[i], &ctMAC[0][i]);
			}
		}
	}
}

int
tiimat3_matrix_mac_verify(const tiimat3_Matrix *mac, const tiimat3_Matrix *op)
{
	tiimat3_Block *tmp;
	size_t row, col;
	int cmp;

	assert(mac->drow == op->drow);
	assert(mac->dcol == op->dcol);

	tmp = block_alloc(TIIMAT3_PACK);
	cmp = 0;

	for (row = 0; row < mac->drow; ++row) {
		for (col = 0; col < mac->dcol; ++col) {
			size_t idx;

			idx = (row * mac->dcol + col) * TIIMAT3_PACK;
			block_mulc(tmp, &op->value[idx], MAC[0], TIIMAT3_PACK);
			cmp += block_cmp(&mac->value[idx], tmp, 1);
		}
	}

	block_dealloc(tmp, TIIMAT3_PACK);

	return cmp;
}

int
tiimat3_matrix_mac_verifyA(const tiimat3_Matrix *mac, const tiimat3_Matrix *op)
{
	tiimat3_Message *m;
	tiimat3_Block *tmp;
	size_t row, col;
	int cmp;

	assert(mac->drow == op->drow);
	assert(mac->dcol == op->dcol);

	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	m = tiimat3_msg_alloc(1);
	tmp = block_alloc(TIIMAT3_PACK);
	cmp = 0;

	for (row = 0; row < mac->drow; ++row) {
		for (col = 0; col < mac->dcol; ++col) {
			size_t idx;

			idx = (row * mac->dcol + col) * TIIMAT3_PACK;
			block_mulc(tmp, &op->value[idx], MAC[0], TIIMAT3_PACK);

			#if TIIMAT3_HOIST || TIIMAT3_MIXED
			block_pack(m, tmp, 0, block_sigma);
			#else
			block_pack(m, tmp, TIIMAT3_D / 2 - TIIMAT3_ROTB, block_sigma);
			#endif /* TIIMAT3_HOIST */

			block_unpack(tmp, m);
			cmp += block_cmp(&mac->value[idx], tmp, 1);
		}
	}

	tiimat3_msg_dealloc(m, 1);
	block_dealloc(tmp, TIIMAT3_PACK);

	return cmp;
}

int
tiimat3_matrix_mac_verifyB(const tiimat3_Matrix *mac, const tiimat3_Matrix *op)
{
	tiimat3_Message *m;
	tiimat3_Block *tmp;
	size_t row, col;
	int cmp;

	assert(mac->drow == op->drow);
	assert(mac->dcol == op->dcol);

	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	m = tiimat3_msg_alloc(1);
	tmp = block_alloc(TIIMAT3_PACK);
	cmp = 0;

	for (row = 0; row < mac->drow; ++row) {
		for (col = 0; col < mac->dcol; ++col) {
			size_t idx;

			idx = (row * mac->dcol + col) * TIIMAT3_PACK;
			block_mulc(tmp, &op->value[idx], MAC[0], TIIMAT3_PACK);
			block_pack(m, tmp, 0, block_tau);
			block_unpack(tmp, m);
			cmp += block_cmp(&mac->value[idx], tmp, 1);
		}
	}

	tiimat3_msg_dealloc(m, 1);
	block_dealloc(tmp, TIIMAT3_PACK);

	return cmp;
}

void
tiimat3_matrix_mul(tiimat3_Matrix *rop, const tiimat3_Matrix *op1, const tiimat3_Matrix *op2)
{
	tiimat3_Block *tmp;
	size_t row, col, k;

	assert(rop->drow == op1->drow);
	assert(rop->dcol == op2->dcol);
	assert(op1->dcol == op2->drow);

	tmp = block_alloc(TIIMAT3_PACK);

	tiimat3_matrix_init_zero(rop, 1);
	for (row = 0; row < rop->drow; ++row) {
		for (col = 0; col < rop->dcol; ++col) {
			size_t idx;

			idx = (row * rop->dcol + col) * TIIMAT3_PACK;
			for (k = 0; k < op1->dcol; ++k) {
				size_t idx1, idx2;

				idx1 = (row * op1->dcol + k) * TIIMAT3_PACK;
				idx2 = (k * op2->dcol + col) * TIIMAT3_PACK;

				block_mul(tmp, &op1->value[idx1], &op2->value[idx2], TIIMAT3_PACK);
				block_add(&rop->value[idx], &rop->value[idx], tmp, TIIMAT3_PACK);
			}
		}
	}

	block_dealloc(tmp, TIIMAT3_PACK);
}

void
tiimat3_matrix_mul_enc(tiimat3_MatrixEnc *rop, tiimat3_MatrixEnc *op1, tiimat3_MatrixEnc *op2)
{
	tiimat3_MatrixEnc *rotA1, *rotA2, *rotB;
	tiimat3_BlockEnc *tmp;
	size_t row, col, k;

	assert(rop->drow == op1->drow);
	assert(rop->dcol == op2->dcol);
	assert(op1->dcol == op2->drow);

	tmp = tiimat3_util_alloc(1, sizeof *tmp);

	#if TIIMAT3_PREROT
	rotA1 = tiimat3_matrix_alloc_encdim(op1->drow, op1->dcol, 1);
	rotA2 = tiimat3_matrix_alloc_encdim(op1->drow, op1->dcol, 1);
	for (row = 0; row < op1->drow; ++row) {
		for (col = 0; col < op1->dcol; ++col) {
			size_t idx = row * op1->dcol + col;
			size_t dim = idx * TIIMAT3_DIM;

			#if TIIMAT3_HOIST
			tiimat3_block_prerotA_enc_hoisted(&rotA1->value[dim], &rotA2->value[dim], &op1->value[idx]);
			#elif TIIMAT3_MIXED
			tiimat3_block_prerotA_enc_mixed(&rotA1->value[dim], &rotA2->value[dim], &op1->value[idx]);
			#else
			tiimat3_block_prerotA_enc_reuse(&rotA1->value[dim], &rotA2->value[dim], &op1->value[idx]);
			#endif /* TIIMAT3_HOIST */
		}
	}

	rotB = tiimat3_matrix_alloc_encdim(op2->drow, op2->dcol, 1);
	for (row = 0; row < op2->drow; ++row) {
		for (col = 0; col < op2->dcol; ++col) {
			size_t idx = row * op2->dcol + col;
			size_t dim = idx * TIIMAT3_DIM;

			#if TIIMAT3_HOIST
			tiimat3_block_prerotB_enc_hoisted(&rotB->value[dim], &op2->value[idx]);
			#elif TIIMAT3_MIXED
			tiimat3_block_prerotB_enc_mixed(&rotB->value[dim], &op2->value[idx]);
			#else
			tiimat3_block_prerotB_enc_reuse(&rotB->value[dim], &op2->value[idx]);
			#endif /* TIIMAT3_HOIST */
		}
	}
	#else
	rotA1 = tiimat3_matrix_alloc_encdim(0, 0, 1);
	rotA2 = tiimat3_matrix_alloc_encdim(0, 0, 1);
	rotB = tiimat3_matrix_alloc_encdim(0, 0, 1);
	#endif /* TIIMAT3_PREROT */

	for (row = 0; row < rop->drow; ++row) {
		for (col = 0; col < rop->dcol; ++col) {
			size_t idx;

			idx = row * rop->dcol + col;
			memset(&rop->value[idx], 0, sizeof rop->value[idx]);
			for (k = 0; k < op1->dcol; ++k) {
				size_t idx1, idx2;

				idx1 = row * op1->dcol + k;
				idx2 = k * op2->dcol + col;

				#if TIIMAT3_PREROT
					idx1 *= TIIMAT3_DIM;
					idx2 *= TIIMAT3_DIM;
					tiimat3_block_mul_enc_prerot(tmp, &rotA1->value[idx1], &rotA2->value[idx1], &rotB->value[idx2]);
				#else
					#if TIIMAT3_HOIST
					tiimat3_block_mul_enc_hoisted(tmp, &op1->value[idx1], &op2->value[idx2]);
					#else
					tiimat3_block_mul_enc_reuse(tmp, &op1->value[idx1], &op2->value[idx2]);
					#endif /* TIIMAT3_HOIST */
				#endif /* TIIMAT3_PREROT */

				tiimat3_block_add_enc(&rop->value[idx], &rop->value[idx], tmp);
			}
		}
	}

	tiimat3_matrix_dealloc_enc(rotA1, 1);
	tiimat3_matrix_dealloc_enc(rotA2, 1);
	tiimat3_matrix_dealloc_enc(rotB, 1);
	tiimat3_util_dealloc(tmp);
}

void
tiimat3_matrix_print(const tiimat3_Matrix *op, size_t pack)
{
	size_t mrow, mcol, brow, bcol;

	for (mrow = 0; mrow < op->drow; ++mrow) {
		for (brow = 0; brow < TIIMAT3_DIM; ++brow) {
			for (mcol = 0; mcol < op->dcol; ++mcol) {
				for (bcol = 0; bcol < TIIMAT3_DIM; ++bcol) {
					size_t idx;

					idx = (mrow * op->dcol + mcol) * TIIMAT3_PACK + pack;
					gmp_printf("%5Zd ", op->value[idx].value[brow][bcol]);
				}
			}
			puts("");
		}
	}
}
