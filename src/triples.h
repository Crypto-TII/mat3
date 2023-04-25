#ifndef TRIPLES_H
#define TRIPLES_H

#include "libbgv.h"

#ifndef TRIPLES_HOIST
#define TRIPLES_HOIST 1
#endif /* TRIPLES_HOIST */

#ifndef TRIPLES_MIXED
#define TRIPLES_MIXED 0
#endif /* TRIPLES_MIXED */

#ifndef TRIPLES_PREROT
#define TRIPLES_PREROT 1
#endif /* TRIPLES_PREROT */

#ifndef TRIPLES_DIM
#define TRIPLES_DIM 32
#endif /* TRIPLES_DIM */

#define TRIPLES_FLOOD_BITS  80
#define TRIPLES_SHARES      5
#define TRIPLES_THREADS     16

#define TRIPLES_FLOOD_BYTES ((TRIPLES_FLOOD_BITS - 1) / CHAR_BIT + 1)
#define TRIPLES_PACK        (TIIMAT3_D / (TRIPLES_DIM * TRIPLES_DIM))
#define TRIPLES_ROTA        (TRIPLES_PACK / 2)
#define TRIPLES_ROTB        (TRIPLES_PACK / 2 * TRIPLES_DIM)

typedef struct triples_block             triples_Block;
typedef struct triples_block_ciphertext  triples_BlockEnc;
typedef struct triples_matrix_block      triples_Matrix;
typedef struct triples_matrix_ciphertext triples_MatrixEnc;

struct triples_block {
	mpz_t value[TRIPLES_DIM][TRIPLES_DIM];
};

struct triples_block_ciphertext {
	bgv_Ciphertext value[TIIMAT3_QPLEN];
};

struct triples_matrix_block {
	triples_Block *value;
	size_t drow, dcol;
};

struct triples_matrix_ciphertext {
	triples_BlockEnc *value;
	size_t drow, dcol;
};

void triples_deinit(void);
void triples_deinit_hoisted(void);
void triples_deinit_mixed(void);
void triples_init(void);
void triples_init_hoisted(void);
void triples_init_mixed(void);

void triples_matrix_add(triples_Matrix *rop, const triples_Matrix *op1, const triples_Matrix *op2);
void triples_matrix_add_enc(triples_MatrixEnc *rop, triples_MatrixEnc *op1, triples_MatrixEnc *op2);
triples_Matrix *triples_matrix_alloc(size_t drow, size_t dcol, size_t len);
triples_MatrixEnc *triples_matrix_alloc_enc(size_t drow, size_t dcol, size_t len);
int  triples_matrix_cmp(const triples_Matrix *op1, const triples_Matrix *op2);
void triples_matrix_dealloc(triples_Matrix *op, size_t len);
void triples_matrix_dealloc_enc(triples_MatrixEnc *op, size_t len);
void triples_matrix_decrypt(triples_Matrix *rop, triples_MatrixEnc *op);
void triples_matrix_encrypt(triples_MatrixEnc *rop, const triples_Matrix *op);
void triples_matrix_encryptA(triples_MatrixEnc *rop, const triples_Matrix *op);
void triples_matrix_encryptB(triples_MatrixEnc *rop, const triples_Matrix *op);
void triples_matrix_init_asc(triples_Matrix *op, size_t len);
void triples_matrix_init_rand(triples_Matrix *op, size_t len);
void triples_matrix_init_zero(triples_Matrix *op, size_t len);
void triples_matrix_init_zero_enc(triples_MatrixEnc *op, size_t len);
void triples_matrix_mac(triples_MatrixEnc *rop, triples_MatrixEnc *op);
int  triples_matrix_mac_verify(const triples_Matrix *mac, const triples_Matrix *op);
int  triples_matrix_mac_verifyA(const triples_Matrix *mac, const triples_Matrix *op);
int  triples_matrix_mac_verifyB(const triples_Matrix *mac, const triples_Matrix *op);
void triples_matrix_mul(triples_Matrix *rop, const triples_Matrix *op1, const triples_Matrix *op2);
void triples_matrix_mul_enc(triples_MatrixEnc *rop, triples_MatrixEnc *op1, triples_MatrixEnc *op2);
void triples_matrix_print(const triples_Matrix *op, size_t len);

#endif /* TRIPLES_H */


#ifdef TRIPLES_IMPL
static gmp_randstate_t randstate;
static bgv_KeySecret sk;
static bgv_KeyPublic pk[TIIMAT3_QLEN];
static bgv_KeySwitch ksw2[TIIMAT3_QPLEN];
static bgv_KeySwitch *kswA;
static bgv_KeySwitch *kswAneg;
static bgv_KeySwitch *kswB;

static mpz_t MAC[TRIPLES_SHARES + 1];
static bgv_Ciphertext ctMAC[TRIPLES_SHARES + 1][TIIMAT3_QLEN];
static tiimat3_Message Uphi[TRIPLES_DIM][2];

static void block_add(triples_Block *rop, const triples_Block *op1, const triples_Block *op2, size_t len);
static void block_add_enc(triples_BlockEnc *rop, triples_BlockEnc *op1, triples_BlockEnc *op2);
static triples_Block *block_alloc(size_t len);
static int  block_cmp(const triples_Block *op1, const triples_Block *op2, size_t len);
static void block_cpy(triples_Block *rop, const triples_Block *op);
static void block_dealloc(triples_Block *op, size_t len);
static void block_decrypt(triples_Block *rop, triples_BlockEnc *op);
static void block_encrypt(triples_BlockEnc *rop, const triples_Block *blocks, size_t rot, void (*transform)(triples_Block *, const triples_Block *));
static void block_init_rand(triples_Block *rop, size_t len);
static void block_init_zero(triples_Block *rop, size_t len);
static void block_mul(triples_Block *rop, const triples_Block *op1, const triples_Block *op2, size_t len);
static void block_mulc(triples_Block *rop, const triples_Block *op1, const mpz_t op2, size_t len);
static void block_mul_enc_hoisted(triples_BlockEnc *rop, triples_BlockEnc *A, triples_BlockEnc *B);
static void block_mul_enc_prerot(triples_BlockEnc *rop, triples_BlockEnc *A1, triples_BlockEnc *A2, triples_BlockEnc *B);
static void block_mul_enc_reuse(triples_BlockEnc *rop, triples_BlockEnc *A, triples_BlockEnc *B);
static void block_pack(tiimat3_Message *m, const triples_Block *op, size_t rot, void (*transform)(triples_Block *, const triples_Block *));
static void block_prerotA_enc_hoisted(triples_BlockEnc *rotA1, triples_BlockEnc *rotA2, triples_BlockEnc *A);
static void block_prerotA_enc_mixed(triples_BlockEnc *rotA1, triples_BlockEnc *rotA2, triples_BlockEnc *A);
static void block_prerotA_enc_reuse(triples_BlockEnc *rotA1, triples_BlockEnc *rotA2, triples_BlockEnc *A);
static void block_prerotB_enc_hoisted(triples_BlockEnc *rotB, triples_BlockEnc *B);
static void block_prerotB_enc_mixed(triples_BlockEnc *rotB, triples_BlockEnc *B);
static void block_prerotB_enc_reuse(triples_BlockEnc *rotB, triples_BlockEnc *B);
static void block_sigma(triples_Block *rop, const triples_Block *op);
static void block_tau(triples_Block *rop, const triples_Block *op);
static void block_unpack(triples_Block *vec, const tiimat3_Message *m);

static void ciphertext_cpy(bgv_Ciphertext *rop, bgv_Ciphertext *ct, size_t len);
static void ciphertext_rotA(bgv_Ciphertext *rop, bgv_Ciphertext *ct);
static void ciphertext_rotA_hoisted(bgv_Ciphertext *rop, bgv_CiphertextSwitch *csw, int k);
static void ciphertext_rotB(bgv_Ciphertext *rop, bgv_Ciphertext *ct);
static void ciphertext_rotB_hoisted(bgv_Ciphertext *rop, bgv_CiphertextSwitch *csw, int k);

static triples_MatrixEnc *matrix_alloc_encdim(size_t drow, size_t dcol, size_t len);

static void msg_cpy(tiimat3_Message *rop, const tiimat3_Message *m);
static void msg_deinit(tiimat3_Message *m);
static void msg_init(tiimat3_Message *m);
static void msg_rot(tiimat3_Message *rop, const tiimat3_Message *m, size_t steps);

static void mpz_addmod(mpz_t rop, const mpz_t op1, const mpz_t op2);
static void mpz_mulmod(mpz_t rop, const mpz_t op1, const mpz_t op2);

static void poly_flood(size_t idx, tiimat3_Poly *p, tiimat3_Seed *seed);

static void
block_add(triples_Block *rop, const triples_Block *op1, const triples_Block *op2, size_t len)
{
	size_t i, row, col;

	for (i = 0; i < len; ++i)
		for (row = 0; row < TRIPLES_DIM; ++row)
			for (col = 0; col < TRIPLES_DIM; ++col)
				mpz_addmod(rop[i].value[row][col], op1[i].value[row][col], op2[i].value[row][col]);
}

static void
block_add_enc(triples_BlockEnc *rop, triples_BlockEnc *op1, triples_BlockEnc *op2)
{
	size_t i;

	for (i = 0; i < TIIMAT3_QLEN; ++i)
		bgv_add(i, &rop->value[i], &op1->value[i], &op2->value[i]);
}

static triples_Block *
block_alloc(size_t len)
{
	triples_Block *block;
	size_t i, row, col;

	block = bgv_alloc(len, sizeof *block);
	for (i = 0; i < len; ++i)
		for (row = 0; row < TRIPLES_DIM; ++row)
			for (col = 0; col < TRIPLES_DIM; ++col)
				mpz_init(block[i].value[row][col]);

	return block;
}

static int
block_cmp(const triples_Block *op1, const triples_Block *op2, size_t len)
{
	size_t i, row, col;

	for (i = 0; i < len; ++i)
		for (row = 0; row < TRIPLES_DIM; ++row)
			for (col = 0; col < TRIPLES_DIM; ++col)
				if (mpz_cmp(op1[i].value[row][col], op2[i].value[row][col]) != 0)
					return 1;

	return 0;
}

static void
block_cpy(triples_Block *rop, const triples_Block *op)
{
	size_t row, col;

	for (row = 0; row < TRIPLES_DIM; ++row)
		for (col = 0; col < TRIPLES_DIM; ++col)
			mpz_set(rop->value[row][col], op->value[row][col]);
}

static void
block_dealloc(triples_Block *op, size_t len)
{
	size_t i, row, col;

	for (i = 0; i < len; ++i)
		for (row = 0; row < TRIPLES_DIM; ++row)
			for (col = 0; col < TRIPLES_DIM; ++col)
				mpz_clear(op[i].value[row][col]);
	bgv_dealloc(op);
}

static void
block_decrypt(triples_Block *rop, triples_BlockEnc *op)
{
	tiimat3_Seed seed;
	tiimat3_Message *m;
	tiimat3_Poly *p;
	size_t i;

	m = tiimat3_alloc_msg(1);
	p = bgv_alloc(TIIMAT3_QLEN, sizeof *p);
	tiimat3_random_seed(&seed);

	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		if (tiimat3_q[i].drop)
			continue;

		bgv_decrypt(i, &p[i], &sk, &op->value[i]);
		poly_flood(i, &p[i], &seed);
		seed.init = 1;
	}
	bgv_decode(m, p);
	block_unpack(rop, m);

	tiimat3_dealloc_msg(m, 1);
	bgv_dealloc(p);
}

static void
block_encrypt(triples_BlockEnc *rop, const triples_Block *blocks, size_t rot, void (*transform)(triples_Block *, const triples_Block *))
{
	tiimat3_Seed seed;
	tiimat3_Message *m;
	tiimat3_Poly *p;
	size_t i;

	m = tiimat3_alloc_msg(1);
	p = bgv_alloc(1, sizeof *p);
	tiimat3_random_seed(&seed);

	block_pack(m, blocks, rot, transform);
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		bgv_encode(i, p, m);
		bgv_encrypt(i, &rop->value[i], &pk[i], p, &seed);
	}

	tiimat3_dealloc_msg(m, 1);
	bgv_dealloc(p);
}

static void
block_init_rand(triples_Block *rop, size_t len)
{
	size_t i, row, col;

	for (i = 0; i < len; ++i)
		for (row = 0; row < TRIPLES_DIM; ++row)
			for (col = 0; col < TRIPLES_DIM; ++col)
				mpz_urandomm(rop[i].value[row][col], randstate, tiimat3_t.value);
}

static void
block_init_zero(triples_Block *rop, size_t len)
{
	size_t i, row, col;

	for (i = 0; i < len; ++i)
		for (row = 0; row < TRIPLES_DIM; ++row)
			for (col = 0; col < TRIPLES_DIM; ++col)
				mpz_set_ui(rop[i].value[row][col], 0);
}

void
block_mul(triples_Block *rop, const triples_Block *op1, const triples_Block *op2, size_t len)
{
	triples_Block *cpy;
	size_t i, row, col, k;
	mpz_t tmp;

	cpy = block_alloc(1);
	mpz_init(tmp);

	for (i = 0; i < len; ++i) {
		for (row = 0; row < TRIPLES_DIM; ++row) {
			for (col = 0; col < TRIPLES_DIM; ++col) {
				mpz_set_ui(cpy->value[row][col], 0);
				for (k = 0; k < TRIPLES_DIM; ++k) {
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
block_mulc(triples_Block *rop, const triples_Block *op1, const mpz_t op2, size_t len)
{
	size_t i, row, col;

	for (i = 0; i < len; ++i)
		for (row = 0; row < TRIPLES_DIM; ++row)
			for (col = 0; col < TRIPLES_DIM; ++col)
				mpz_mulmod(rop[i].value[row][col], op1[i].value[row][col], op2);
}

__attribute__((unused)) static void
block_mul_enc_hoisted(triples_BlockEnc *rop, triples_BlockEnc *A, triples_BlockEnc *B)
{
	bgv_CiphertextSwitch *cswA, *cswB;
	bgv_Ciphertext *rotA, *rotAneg, *rotB, *ct;
	bgv_Delta *delta;
	int i, k;

	cswA = bgv_alloc(TIIMAT3_QPLEN, sizeof *cswA);
	cswB = bgv_alloc(TIIMAT3_QPLEN, sizeof *cswB);
	rotA = bgv_alloc(TIIMAT3_QPLEN, sizeof *rotA);
	rotAneg = bgv_alloc(TIIMAT3_QPLEN, sizeof *rotAneg);
	rotB = bgv_alloc(TIIMAT3_QPLEN, sizeof *rotB);
	ct = bgv_alloc(2, sizeof *ct);
	delta = bgv_alloc(2, sizeof *delta);

	/* reset dropped moduli */
	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	/* ciphertext extension */
	bgv_keyswitch_ext(cswA, A->value, 1);
	bgv_keyswitch_ext(cswB, B->value, 1);

	/* first multiplication */
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		ciphertext_cpy(&ct[0], &A->value[i], 1);
		ciphertext_cpy(&ct[1], &B->value[i], 1);

		if (i == 0) {
			bgv_modswitch_delta(0, &delta[0], &ct[0]);
			bgv_modswitch_delta(0, &delta[1], &ct[1]);
			continue;
		} else {
			bgv_modswitch_ext(i, &ct[0], &delta[0]);
			bgv_modswitch_ext(i, &ct[1], &delta[1]);
		}

		bgv_mul(i, &rop->value[i], &ct[0], &ct[1]);
	}

	/* squared block multiplication */
	for (k = 1; k < TRIPLES_DIM; ++k) {
		ciphertext_rotA_hoisted(rotA, cswA, k);
		ciphertext_rotA_hoisted(rotAneg, cswA, -k);
		ciphertext_rotB_hoisted(rotB, cswB, k);

		/* multiplication */
		for (i = 0; i < TIIMAT3_QLEN; ++i) {
			tiimat3_Poly p;

			bgv_encode(i, &p, &Uphi[k][0]);
			bgv_mulc(i, &rotA[i], &rotA[i], &p);

			bgv_encode(i, &p, &Uphi[k][1]);
			bgv_mulc(i, &rotAneg[i], &rotAneg[i], &p);

			bgv_add(i, &rotA[i], &rotA[i], &rotAneg[i]);

			if (i == 0) {
				bgv_modswitch_delta(0, &delta[0], &rotA[0]);
				bgv_modswitch_delta(0, &delta[1], &rotB[0]);
				continue;
			} else {
				bgv_modswitch_ext(i, &rotA[i], &delta[0]);
				bgv_modswitch_ext(i, &rotB[i], &delta[1]);
			}

			bgv_mul(i, ct, &rotA[i], &rotB[i]);
			bgv_add(i, &rop->value[i], &rop->value[i], ct);
		}
	}
	tiimat3_mod_drop(0, 2);

	bgv_keyswitch(rop->value, ksw2, 2);

	for (i = 1; i < TIIMAT3_QLEN; ++i) {
		bgv_Delta delta;

		if (i == 1) {
			bgv_modswitch_delta(1, &delta, &rop->value[1]);
			continue;
		} else
			bgv_modswitch_ext(i, &rop->value[i], &delta);
	}
	tiimat3_mod_drop(1, 1);

	bgv_dealloc(cswA);
	bgv_dealloc(cswB);
	bgv_dealloc(rotA);
	bgv_dealloc(rotAneg);
	bgv_dealloc(rotB);
	bgv_dealloc(ct);
	bgv_dealloc(delta);
}

__attribute__((unused)) static void
block_mul_enc_prerot(triples_BlockEnc *rop, triples_BlockEnc *A1, triples_BlockEnc *A2, triples_BlockEnc *B)
{
	triples_BlockEnc *mul;
	size_t i, k;

	mul = bgv_alloc(TRIPLES_DIM, sizeof *mul);

	/* reset dropped moduli */
	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	/* squared block multiplication */
	#pragma omp parallel for
	for (k = 0; k < TRIPLES_DIM; ++k) {
		bgv_Ciphertext *ct;
		bgv_Delta *delta;
		size_t i;

		ct = bgv_alloc(2, sizeof *ct);
		delta = bgv_alloc(2, sizeof *delta);

		if (k == 0) {
			for (i = 0; i < TIIMAT3_QLEN; ++i) {
				ciphertext_cpy(&ct[0], &A2[0].value[i], 1);
				ciphertext_cpy(&ct[1], &B[0].value[i], 1);

				if (i == 0) {
					bgv_modswitch_delta(0, &delta[0], &ct[0]);
					bgv_modswitch_delta(0, &delta[1], &ct[1]);
				} else {
					bgv_modswitch_ext(i, &ct[0], &delta[0]);
					bgv_modswitch_ext(i, &ct[1], &delta[1]);
					bgv_mul(i, &mul[k].value[i], &ct[0], &ct[1]);
				}
			}
		} else {
			for (i = 0; i < TIIMAT3_QLEN; ++i) {
				tiimat3_Poly p;

				bgv_encode(i, &p, &Uphi[k][0]);
				bgv_mulc(i, &ct[0], &A2[k].value[i], &p);

				bgv_encode(i, &p, &Uphi[k][1]);
				bgv_mulc(i, &ct[1], &A1[k].value[i], &p);

				bgv_add(i, &ct[0], &ct[0], &ct[1]);
				ciphertext_cpy(&ct[1], &B[k].value[i], 1);

				if (i == 0) {
					bgv_modswitch_delta(0, &delta[0], &ct[0]);
					bgv_modswitch_delta(0, &delta[1], &ct[1]);
				} else {
					bgv_modswitch_ext(i, &ct[0], &delta[0]);
					bgv_modswitch_ext(i, &ct[1], &delta[1]);
					bgv_mul(i, &mul[k].value[i], &ct[0], &ct[1]);
				}
			}
		}

		bgv_dealloc(ct);
		bgv_dealloc(delta);
	}

	for (k = 1; k < TRIPLES_DIM; ++k) {
		for (i = 0; i < TIIMAT3_QLEN; ++i) {
			if (k == 1)
				bgv_add(i, &rop->value[i], &mul[0].value[i], &mul[1].value[i]);
			else
				bgv_add(i, &rop->value[i], &rop->value[i], &mul[k].value[i]);
		}
	}
	tiimat3_mod_drop(0, 2);
	bgv_keyswitch(rop->value, ksw2, 2);

	for (i = 1; i < TIIMAT3_QLEN; ++i) {
		bgv_Delta delta;

		if (i == 1) {
			bgv_modswitch_delta(1, &delta, &rop->value[1]);
			continue;
		} else
			bgv_modswitch_ext(i, &rop->value[i], &delta);
	}
	tiimat3_mod_drop(1, 1);

	bgv_dealloc(mul);
}

__attribute__((unused)) static void
block_mul_enc_reuse(triples_BlockEnc *rop, triples_BlockEnc *A, triples_BlockEnc *B)
{
	bgv_Ciphertext *rotA1, *rotA2, *rotB, *ct;
	bgv_Delta *delta;
	size_t i, k;

	delta = bgv_alloc(2, sizeof *delta);
	rotA1 = bgv_alloc(TIIMAT3_QPLEN, sizeof *rotA1);
	rotA2 = bgv_alloc(TIIMAT3_QPLEN, sizeof *rotA2);
	rotB = bgv_alloc(TIIMAT3_QPLEN, sizeof *rotB);
	ct = bgv_alloc(2, sizeof *ct);

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
			bgv_modswitch_delta(0, &delta[0], &ct[0]);
			bgv_modswitch_delta(0, &delta[1], &ct[1]);
			continue;
		} else {
			bgv_modswitch_ext(i, &ct[0], &delta[0]);
			bgv_modswitch_ext(i, &ct[1], &delta[1]);
		}

		bgv_mul(i, &rop->value[i], &ct[0], &ct[1]);
	}

	/* squared block multiplication */
	for (k = 1; k < TRIPLES_DIM; ++k) {
		ciphertext_rotA(rotA1, rotA1);
		ciphertext_rotA(rotA2, rotA2);
		ciphertext_rotB(rotB, rotB);

		/* multiplication */
		for (i = 0; i < TIIMAT3_QLEN; ++i) {
			tiimat3_Poly p;

			bgv_encode(i, &p, &Uphi[k][0]);
			bgv_mulc(i, &ct[0], &rotA2[i], &p);

			bgv_encode(i, &p, &Uphi[k][1]);
			bgv_mulc(i, &ct[1], &rotA1[i], &p);

			bgv_add(i, &ct[0], &ct[0], &ct[1]);
			ciphertext_cpy(&ct[1], &rotB[i], 1);

			if (i == 0) {
				bgv_modswitch_delta(0, &delta[0], &ct[0]);
				bgv_modswitch_delta(0, &delta[1], &ct[1]);
				continue;
			} else {
				bgv_modswitch_ext(i, &ct[0], &delta[0]);
				bgv_modswitch_ext(i, &ct[1], &delta[1]);
			}

			bgv_mul(i, &ct[0], &ct[0], &ct[1]);
			bgv_add(i, &rop->value[i], &rop->value[i], &ct[0]);
		}
	}
	tiimat3_mod_drop(0, 2);

	bgv_keyswitch(rop->value, ksw2, 2);

	for (i = 1; i < TIIMAT3_QLEN; ++i) {
		bgv_Delta delta;

		if (i == 1) {
			bgv_modswitch_delta(1, &delta, &rop->value[1]);
			continue;
		} else
			bgv_modswitch_ext(i, &rop->value[i], &delta);
	}
	tiimat3_mod_drop(1, 1);

	bgv_dealloc(delta);
	bgv_dealloc(rotA1);
	bgv_dealloc(rotA2);
	bgv_dealloc(rotB);
	bgv_dealloc(ct);
}

static void
block_pack(tiimat3_Message *m, const triples_Block *op, size_t rot, void (*transform)(triples_Block *, const triples_Block *))
{
	triples_Block *cpy;
	size_t row, col, i;

	cpy = block_alloc(1);

	for (i = 0; i < TRIPLES_PACK / 2; ++i) {
		transform(cpy, &op[i]);

		for (row = 0; row < TRIPLES_DIM; ++row) {
			for (col = 0; col < TRIPLES_DIM; ++col) {
				size_t j = (row * TRIPLES_DIM + col) * TRIPLES_PACK / 2 + i;
				mpz_set(m->value[j], cpy->value[row][col]);
			}
		}
	}

	for (i = 0; i < TRIPLES_PACK / 2; ++i) {
		transform(cpy, &op[TRIPLES_PACK / 2 + i]);

		for (row = 0; row < TRIPLES_DIM; ++row) {
			for (col = 0; col < TRIPLES_DIM; ++col) {
				size_t j = TIIMAT3_D / 2 + (row * TRIPLES_DIM + col) * TRIPLES_PACK / 2 + i;
				mpz_set(m->value[j], cpy->value[row][col]);
			}
		}
	}

	msg_rot(m, m, rot);
	tiimat3_msg_pack(m, m);

	block_dealloc(cpy, 1);
}

__attribute__((unused)) static void
block_prerotA_enc_hoisted(triples_BlockEnc *rotA1, triples_BlockEnc *rotA2, triples_BlockEnc *A)
{
	bgv_CiphertextSwitch *csw;
	size_t k;

	csw = bgv_alloc(TIIMAT3_QPLEN, sizeof *csw);
	bgv_keyswitch_ext(csw, A->value, 1);

	#pragma omp parallel for
	for (k = 0; k < TRIPLES_DIM; ++k) {
		if (k == 0) {
			ciphertext_cpy(rotA2[0].value, A->value, TIIMAT3_QPLEN);
		} else {
			ciphertext_rotA_hoisted(rotA1[k].value, csw, -k);
			ciphertext_rotA_hoisted(rotA2[k].value, csw, k);
		}
	}

	bgv_dealloc(csw);
}

__attribute__((unused)) static void
block_prerotA_enc_mixed(triples_BlockEnc *rotA1, triples_BlockEnc *rotA2, triples_BlockEnc *A)
{
	bgv_CiphertextSwitch *csw;
	size_t k1, k2;

	csw = bgv_alloc(TIIMAT3_QPLEN, sizeof *csw);
	ciphertext_cpy(rotA2[0].value, A->value, TIIMAT3_QPLEN);

	for (k1 = 0; k1 < TRIPLES_DIM / TRIPLES_THREADS; ++k1) {
		bgv_keyswitch_ext(csw, rotA2[k1 * TRIPLES_THREADS].value, 1);

		#pragma omp parallel for
		for (k2 = 1; k2 <= TRIPLES_THREADS; ++k2) {
			size_t k = k1 * TRIPLES_THREADS + k2;

			if (k < TRIPLES_DIM) {
				ciphertext_rotA_hoisted(rotA1[k].value, csw, -k2);
				ciphertext_rotA_hoisted(rotA2[k].value, csw, k2);
			}
		}
	}

	bgv_dealloc(csw);
}


__attribute__((unused)) static void
block_prerotA_enc_reuse(triples_BlockEnc *rotA1, triples_BlockEnc *rotA2, triples_BlockEnc *A)
{
	size_t k;

	ciphertext_rotB(rotA2[0].value, A->value);
	ciphertext_cpy(rotA1[0].value, A->value, TIIMAT3_QPLEN);

	for (k = 1; k < TRIPLES_DIM; ++k) {
		ciphertext_rotA(rotA1[k].value, rotA1[k - 1].value);
		ciphertext_rotA(rotA2[k].value, rotA2[k - 1].value);
	}
}

__attribute__((unused)) static void
block_prerotB_enc_hoisted(triples_BlockEnc *rotB, triples_BlockEnc *B)
{
	bgv_CiphertextSwitch *csw;
	size_t k;

	csw = bgv_alloc(TIIMAT3_QPLEN, sizeof *csw);
	bgv_keyswitch_ext(csw, B->value, 1);

	#pragma omp parallel for
	for (k = 0; k < TRIPLES_DIM; ++k)
		if (k == 0)
			ciphertext_cpy(rotB[0].value, B->value, TIIMAT3_QPLEN);
		else
			ciphertext_rotB_hoisted(rotB[k].value, csw, k);

	bgv_dealloc(csw);
}

__attribute__((unused)) static void
block_prerotB_enc_mixed(triples_BlockEnc *rotB, triples_BlockEnc *B)
{
	bgv_CiphertextSwitch *csw;
	size_t k1, k2;

	csw = bgv_alloc(TIIMAT3_QPLEN, sizeof *csw);
	ciphertext_cpy(rotB[0].value, B->value, TIIMAT3_QPLEN);

	for (k1 = 0; k1 < TRIPLES_DIM / TRIPLES_THREADS; ++k1) {
		bgv_keyswitch_ext(csw, rotB[k1 * TRIPLES_THREADS].value, 1);

		#pragma omp parallel for
		for (k2 = 1; k2 <= TRIPLES_THREADS; ++k2) {
			size_t k = k1 * TRIPLES_THREADS + k2;

			if (k < TRIPLES_DIM) {
				ciphertext_rotB_hoisted(rotB[k].value, csw, k2);
			}
		}
	}

	bgv_dealloc(csw);
}


__attribute__((unused)) static void
block_prerotB_enc_reuse(triples_BlockEnc *rotB, triples_BlockEnc *B)
{
	size_t k;

	ciphertext_cpy(rotB[0].value, B->value, TIIMAT3_QPLEN);

	for (k = 1; k < TRIPLES_DIM; ++k)
		ciphertext_rotB(rotB[k].value, rotB[k - 1].value);
}

static void
block_sigma(triples_Block *rop, const triples_Block *op)
{
	triples_Block *cpy;
	size_t row, col;

	cpy = block_alloc(1);

	for (row = 0; row < TRIPLES_DIM; ++row)
		for (col = 0; col < TRIPLES_DIM; ++col)
			mpz_set(cpy->value[row][col], op->value[row][(row + col) % TRIPLES_DIM]);
	block_cpy(rop, cpy);

	block_dealloc(cpy, 1);
}

static void
block_tau(triples_Block *rop, const triples_Block *op)
{
	triples_Block *cpy;
	size_t row, col;

	cpy = block_alloc(1);

	for (row = 0; row < TRIPLES_DIM; ++row)
		for (col = 0; col < TRIPLES_DIM; ++col)
			mpz_set(cpy->value[row][col], op->value[(row + col) % TRIPLES_DIM][col]);
	block_cpy(rop, cpy);

	block_dealloc(cpy, 1);
}

static void
block_unpack(triples_Block *vec, const tiimat3_Message *m)
{
	tiimat3_Message *cpy;
	size_t i, row, col;

	cpy = tiimat3_alloc_msg(1);
	tiimat3_msg_unpack(cpy, m);

	for (i = 0; i < TRIPLES_PACK / 2; ++i) {
		for (row = 0; row < TRIPLES_DIM; ++row) {
			for (col = 0; col < TRIPLES_DIM; ++col) {
				size_t j = (row * TRIPLES_DIM + col) * TRIPLES_PACK / 2 + i;
				mpz_set(vec[i].value[row][col], cpy->value[j]);
			}
		}
	}

	for (i = 0; i < TRIPLES_PACK / 2; ++i) {
		for (row = 0; row < TRIPLES_DIM; ++row) {
			for (col = 0; col < TRIPLES_DIM; ++col) {
				size_t j = TIIMAT3_D / 2 + (row * TRIPLES_DIM + col) * TRIPLES_PACK / 2 + i;
				mpz_set(vec[TRIPLES_PACK / 2 + i].value[row][col], cpy->value[j]);
			}
		}
	}

	tiimat3_dealloc_msg(cpy, 1);
}

static void
ciphertext_cpy(bgv_Ciphertext *rop, bgv_Ciphertext *ct, size_t len)
{

	memcpy(rop, ct, len * sizeof *ct);
}

static void
ciphertext_rotA(bgv_Ciphertext *rop, bgv_Ciphertext *ct)
{
	size_t i;

	if (rop == ct) {
		for (i = 0; i < TIIMAT3_QLEN; ++i)
			bgv_rot_inplace(i, &ct[i], TRIPLES_ROTA);
	} else {
		for (i = 0; i < TIIMAT3_QLEN; ++i)
			bgv_rot(i, &rop[i], &ct[i], TRIPLES_ROTA);
	}

	bgv_keyswitch(rop, &kswA[TIIMAT3_QPLEN], 1);
}

static void
ciphertext_rotA_hoisted(bgv_Ciphertext *rop, bgv_CiphertextSwitch *csw, int k)
{
	bgv_CiphertextSwitch *cswr;
	bgv_Delta *delta;
	size_t i;

	cswr = bgv_alloc(1, sizeof *cswr);
	delta = bgv_alloc(TIIMAT3_PLEN, sizeof *delta);

	if (k < 0) {
		k = -k;
		for (i = 0; i < TIIMAT3_QPLEN; ++i) {
			bgv_rot_csw(i, cswr, &csw[i], TIIMAT3_D / 2 - TRIPLES_ROTB + k * TRIPLES_ROTA);
			bgv_keyswitch_dot(i, &rop[i], cswr, &kswAneg[k * TIIMAT3_QPLEN + i]);
		}
	} else {
		for (i = 0; i < TIIMAT3_QPLEN; ++i) {
			bgv_rot_csw(i, cswr, &csw[i], k * TRIPLES_ROTA);
			bgv_keyswitch_dot(i, &rop[i], cswr, &kswA[k * TIIMAT3_QPLEN + i]);
		}
	}

	for (i = TIIMAT3_QLEN; i < TIIMAT3_QPLEN; ++i)
		bgv_keyswitch_delta(i, &delta[i - TIIMAT3_QLEN], &rop[i]);
	for (i = 0; i < TIIMAT3_QLEN; ++i)
		bgv_keyswitch_switch(i, &rop[i], delta);

	bgv_dealloc(cswr);
	bgv_dealloc(delta);
}

static void
ciphertext_rotB(bgv_Ciphertext *rop, bgv_Ciphertext *ct)
{
	size_t i;

	if (rop == ct) {
		for (i = 0; i < TIIMAT3_QLEN; ++i)
			bgv_rot_inplace(i, &ct[i], TRIPLES_ROTB);
	} else {
		for (i = 0; i < TIIMAT3_QLEN; ++i)
			bgv_rot(i, &rop[i], &ct[i], TRIPLES_ROTB);
	}

	bgv_keyswitch(rop, &kswB[TIIMAT3_QPLEN], 1);
}

static void
ciphertext_rotB_hoisted(bgv_Ciphertext *rop, bgv_CiphertextSwitch *csw, int k)
{
	bgv_CiphertextSwitch *cswr;
	bgv_Delta *delta;
	size_t i;

	cswr = bgv_alloc(1, sizeof *cswr);
	delta = bgv_alloc(TIIMAT3_PLEN, sizeof *delta);


	for (i = 0; i < TIIMAT3_QPLEN; ++i) {
		bgv_rot_csw(i, cswr, &csw[i], k * TRIPLES_ROTB);
		bgv_keyswitch_dot(i, &rop[i], cswr, &kswB[k * TIIMAT3_QPLEN + i]);
	}

	for (i = TIIMAT3_QLEN; i < TIIMAT3_QPLEN; ++i)
		bgv_keyswitch_delta(i, &delta[i - TIIMAT3_QLEN], &rop[i]);
	for (i = 0; i < TIIMAT3_QLEN; ++i)
		bgv_keyswitch_switch(i, &rop[i], delta);

	bgv_dealloc(cswr);
	bgv_dealloc(delta);
}

static triples_MatrixEnc *
matrix_alloc_encdim(size_t drow, size_t dcol, size_t len)
{
	triples_MatrixEnc *mat;
	size_t i;

	mat = bgv_alloc(len, sizeof *mat);
	for (i = 0; i < len; ++i) {
		mat[i].value = bgv_alloc(drow * dcol * TRIPLES_DIM, sizeof *mat[i].value);
		mat[i].drow = drow;
		mat[i].dcol = dcol;
	}

	return mat;
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

	cpy = tiimat3_alloc_msg(1);
	msg_cpy(cpy, m);

	for (j = 0; j < TIIMAT3_D / 2; ++j) {
		size_t idx = (j + steps) % (TIIMAT3_D / 2);
		mpz_set(rop->value[j], cpy->value[idx]);
	}
	for (; j < TIIMAT3_D; ++j) {
		size_t idx = TIIMAT3_D / 2 + ((j + steps) % (TIIMAT3_D / 2));
		mpz_set(rop->value[j], cpy->value[idx]);
	}

	tiimat3_dealloc_msg(cpy, 1);
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
	unsigned char rand[TRIPLES_FLOOD_BYTES];
	tiimat3_Poly *e;
	mpz_t tmp;
	size_t j;

	e = bgv_alloc(1, sizeof *e);
	mpz_init(tmp);

	for (j = 0; j < TIIMAT3_D; ++j) {
		tiimat3_random_bytes(rand, sizeof rand, seed);
		rand[0] >>= (CHAR_BIT - (TRIPLES_FLOOD_BITS % CHAR_BIT)) % CHAR_BIT;
		mpz_import(tmp, sizeof rand, 1, 1, 1, 0, rand);

		mpz_mod_ui(tmp, tmp, tiimat3_q[idx].value);
		e->value[j] = mpz_get_ui(tmp);
	}
	tiimat3_poly_init_error(idx, e, e);
	tiimat3_poly_add(idx, p, p, e);

	bgv_dealloc(e);
	mpz_clear(tmp);
}

void
triples_deinit(void)
{
	size_t i, k;

	bgv_deinit();
	gmp_randclear(randstate);

	bgv_dealloc(kswA);
	bgv_dealloc(kswAneg);
	bgv_dealloc(kswB);

	for (i = 0; i <= TRIPLES_SHARES; ++i)
		mpz_clear(MAC[i]);
	for (k = 1; k < TRIPLES_DIM; ++k) {
		msg_deinit(&Uphi[k][0]);
		msg_deinit(&Uphi[k][1]);
	}
}

void
triples_deinit_hoisted(void)
{

	;
}

void
triples_deinit_mixed(void)
{

	;
}

void
triples_init(void)
{
	tiimat3_Seed seed[TIIMAT3_QPLEN];
	tiimat3_Message *m;
	tiimat3_Poly p;
	size_t i, j, k;

	bgv_init();
	gmp_randinit_default(randstate);

	#if TRIPLES_HOIST
	kswA = bgv_alloc(TRIPLES_DIM * TIIMAT3_QPLEN, sizeof *kswA);
	kswAneg = bgv_alloc(TRIPLES_DIM * TIIMAT3_QPLEN, sizeof *kswAneg);
	kswB = bgv_alloc(TRIPLES_DIM * TIIMAT3_QPLEN, sizeof *kswB);
	#elif TRIPLES_MIXED
	kswA = bgv_alloc((TRIPLES_THREADS + 1) * TIIMAT3_QPLEN, sizeof *kswA);
	kswAneg = bgv_alloc((TRIPLES_THREADS + 1) * TIIMAT3_QPLEN, sizeof *kswAneg);
	kswB = bgv_alloc((TRIPLES_THREADS + 1) * TIIMAT3_QPLEN, sizeof *kswB);
	#else
	kswA = bgv_alloc(2 * TIIMAT3_QPLEN, sizeof *kswA);
	kswAneg = bgv_alloc(2 * TIIMAT3_QPLEN, sizeof *kswAneg);
	kswB = bgv_alloc(2 * TIIMAT3_QPLEN, sizeof *kswB);
	#endif /* TRIPLES_HOIST */

	for (i = 0; i <= TRIPLES_SHARES; ++i)
		mpz_init(MAC[i]);
	for (k = 1; k < TRIPLES_DIM; ++k) {
		msg_init(&Uphi[k][0]);
		msg_init(&Uphi[k][1]);
	}

	/* BGV keys */
	bgv_keygen_secret(&sk);
	tiimat3_random_seed(seed);
	for (i = 0; i < TIIMAT3_QLEN; ++i)
		bgv_keygen_public(i, &pk[i], &sk, seed);

	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_random_seed(&seed[i]);
	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		bgv_keygen_switchr(i, &kswA[TIIMAT3_QPLEN + i], &sk, TRIPLES_ROTA, seed);

	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_random_seed(&seed[i]);
	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		bgv_keygen_switchr(i, &kswB[TIIMAT3_QPLEN + i], &sk, TRIPLES_ROTB, seed);

	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_random_seed(&seed[i]);
	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		bgv_keygen_switch2(i, &ksw2[i], &sk, seed);

	/* MAC key */
	m = tiimat3_alloc_msg(1);
	mpz_set_ui(MAC[0], 0);
	memset(ctMAC[0], 0, sizeof ctMAC[0]);
	for (i = 1; i <= TRIPLES_SHARES; ++i) {
		mpz_urandomm(MAC[i], randstate, tiimat3_t.value);
		mpz_addmod(MAC[0], MAC[0], MAC[i]);

		for (j = 0; j < TIIMAT3_D; ++j)
			mpz_set(m->value[j], MAC[i]);
		tiimat3_msg_pack(m, m);

		tiimat3_random_seed(seed);
		for (j = 0; j < TIIMAT3_QLEN; ++j) {
			bgv_encode(j, &p, m);
			bgv_encrypt(j, &ctMAC[i][j], &pk[j], &p, seed);
			bgv_add(j, &ctMAC[0][j], &ctMAC[0][j], &ctMAC[i][j]);
		}
	}
	tiimat3_dealloc_msg(m, 1);

	/* Uphi */
	for (k = 1; k < TRIPLES_DIM; ++k) {
		for (j = 0; j < TIIMAT3_D; ++j)
			if (j % (TRIPLES_PACK / 2 * TRIPLES_DIM) < TRIPLES_PACK / 2 * (TRIPLES_DIM - k))
				mpz_set_ui(Uphi[k][0].value[j], 1);
			else
				mpz_set_ui(Uphi[k][0].value[j], 0);
		tiimat3_msg_pack(&Uphi[k][0], &Uphi[k][0]);

		for (j = 0; j < TIIMAT3_D; ++j)
			if (j % (TRIPLES_PACK / 2 * TRIPLES_DIM) < TRIPLES_PACK / 2 * (TRIPLES_DIM - k))
				mpz_set_ui(Uphi[k][1].value[j], 0);
			else
				mpz_set_ui(Uphi[k][1].value[j], 1);
		tiimat3_msg_pack(&Uphi[k][1], &Uphi[k][1]);
	}
}

void
triples_init_hoisted(void)
{
	tiimat3_Seed seed[TIIMAT3_QPLEN];
	size_t i, k;

	for (k = 2; k < TRIPLES_DIM; ++k) {
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_random_seed(&seed[i]);
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			bgv_keygen_switchr(i, &kswA[k * TIIMAT3_QPLEN + i], &sk, k * TRIPLES_ROTA, seed);
	}

	for (k = 1; k < TRIPLES_DIM; ++k) {
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_random_seed(&seed[i]);
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			bgv_keygen_switchr(i, &kswAneg[k * TIIMAT3_QPLEN + i], &sk, TIIMAT3_D / 2 - TRIPLES_ROTB + k * TRIPLES_ROTA, seed);
	}

	for (k = 2; k < TRIPLES_DIM; ++k) {
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_random_seed(&seed[i]);
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			bgv_keygen_switchr(i, &kswB[k * TIIMAT3_QPLEN + i], &sk, k * TRIPLES_ROTB, seed);
	}
}

void
triples_init_mixed(void)
{
	tiimat3_Seed seed[TIIMAT3_QPLEN];
	size_t i, k;

	for (k = 2; k <= TRIPLES_THREADS; ++k) {
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_random_seed(&seed[i]);
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			bgv_keygen_switchr(i, &kswA[k * TIIMAT3_QPLEN + i], &sk, k * TRIPLES_ROTA, seed);
	}

	for (k = 1; k <= TRIPLES_THREADS; ++k) {
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_random_seed(&seed[i]);
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			bgv_keygen_switchr(i, &kswAneg[k * TIIMAT3_QPLEN + i], &sk, TIIMAT3_D / 2 - TRIPLES_ROTB + k * TRIPLES_ROTA, seed);
	}

	for (k = 2; k <= TRIPLES_THREADS; ++k) {
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			tiimat3_random_seed(&seed[i]);
		for (i = 0; i < TIIMAT3_QPLEN; ++i)
			bgv_keygen_switchr(i, &kswB[k * TIIMAT3_QPLEN + i], &sk, k * TRIPLES_ROTB, seed);
	}
}

void
triples_matrix_add(triples_Matrix *rop, const triples_Matrix *op1, const triples_Matrix *op2)
{
	size_t row, col;

	assert(op1->drow == op2->drow);
	assert(op1->dcol == op2->dcol);

	for (row = 0; row < op1->drow; ++row) {
		for (col = 0; col < op1->dcol; ++col) {
			size_t idx = (row * op1->dcol + col) * TRIPLES_PACK;
			block_add(&rop->value[idx], &op1->value[idx], &op2->value[idx], TRIPLES_PACK);
		}
	}

}

void
triples_matrix_add_enc(triples_MatrixEnc *rop, triples_MatrixEnc *op1, triples_MatrixEnc *op2)
{
	size_t row, col;

	assert(op1->drow == op2->drow);
	assert(op1->dcol == op2->dcol);

	for (row = 0; row < op1->drow; ++row) {
		for (col = 0; col < op1->dcol; ++col) {
			size_t idx = row * op1->dcol + col;
			block_add_enc(&rop->value[idx], &op1->value[idx], &op2->value[idx]);
		}
	}

}

triples_Matrix *
triples_matrix_alloc(size_t drow, size_t dcol, size_t len)
{
	triples_Matrix *mat;
	size_t i;

	mat = bgv_alloc(len, sizeof *mat);
	for (i = 0; i < len; ++i) {
		mat[i].value = block_alloc(drow * dcol * TRIPLES_PACK);
		mat[i].drow = drow;
		mat[i].dcol = dcol;
	}

	return mat;
}

triples_MatrixEnc *
triples_matrix_alloc_enc(size_t drow, size_t dcol, size_t len)
{
	triples_MatrixEnc *mat;
	size_t i;

	mat = bgv_alloc(len, sizeof *mat);
	for (i = 0; i < len; ++i) {
		mat[i].value = bgv_alloc(drow * dcol, sizeof *mat[i].value);
		mat[i].drow = drow;
		mat[i].dcol = dcol;
	}

	return mat;
}

int
triples_matrix_cmp(const triples_Matrix *op1, const triples_Matrix *op2)
{
	size_t row, col;
	int cmp;

	assert(op1->drow == op2->drow);
	assert(op1->dcol == op2->dcol);

	cmp = 0;
	for (row = 0; row < op1->drow; ++row) {
		for (col = 0; col < op1->dcol; ++col) {
			size_t idx = (row * op1->dcol + col) * TRIPLES_PACK;
			cmp += block_cmp(&op1->value[idx], &op2->value[idx], TRIPLES_PACK);
		}
	}

	return cmp;
}

void
triples_matrix_dealloc(triples_Matrix *op, size_t len)
{
	size_t i;

	for (i = 0; i < len; ++i)
		block_dealloc(op[i].value, op->drow * op->dcol * TRIPLES_PACK);
	bgv_dealloc(op);
}

void
triples_matrix_dealloc_enc(triples_MatrixEnc *op, size_t len)
{
	size_t i;

	for (i = 0; i < len; ++i)
		bgv_dealloc(op[i].value);
	bgv_dealloc(op);
}

void
triples_matrix_decrypt(triples_Matrix *rop, triples_MatrixEnc *op)
{
	size_t row, col;

	for (row = 0; row < op->drow; ++row) {
		for (col = 0; col < op->dcol; ++col) {
			size_t idx = row * op->dcol + col;
			block_decrypt(&rop->value[idx * TRIPLES_PACK], &op->value[idx]);
		}
	}
}

void
triples_matrix_encrypt(triples_MatrixEnc *rop, const triples_Matrix *op)
{
	size_t row, col;

	for (row = 0; row < op->drow; ++row) {
		for (col = 0; col < op->dcol; ++col) {
			size_t idx = row * op->dcol + col;
			block_encrypt(&rop->value[idx], &op->value[idx * TRIPLES_PACK], 0, block_cpy);
		}
	}
}

void
triples_matrix_encryptA(triples_MatrixEnc *rop, const triples_Matrix *op)
{
	size_t row, col;

	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	for (row = 0; row < op->drow; ++row) {
		for (col = 0; col < op->dcol; ++col) {
			size_t idx = row * op->dcol + col;
			#if TRIPLES_HOIST || TRIPLES_MIXED
			block_encrypt(&rop->value[idx], &op->value[idx * TRIPLES_PACK], 0, block_sigma);
			#else
			block_encrypt(&rop->value[idx], &op->value[idx * TRIPLES_PACK], TIIMAT3_D / 2 - TRIPLES_ROTB, block_sigma);
			#endif /* TRIPLES_HOIST */
		}
	}
}

void
triples_matrix_encryptB(triples_MatrixEnc *rop, const triples_Matrix *op)
{
	size_t row, col;

	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	for (row = 0; row < op->drow; ++row) {
		for (col = 0; col < op->dcol; ++col) {
			size_t idx = row * op->dcol + col;
			block_encrypt(&rop->value[idx], &op->value[idx * TRIPLES_PACK], 0, block_tau);
		}
	}
}

void
triples_matrix_init_asc(triples_Matrix *op, size_t len)
{
	size_t mrow, mcol, brow, bcol, i, ii;

	for (i = 0; i < len; ++i) {
		for (mrow = 0; mrow < op[i].drow; ++mrow) {
			for (mcol = 0; mcol < op[i].dcol; ++mcol) {
				for (ii = 0; ii < TRIPLES_PACK; ++ii) {
					for (brow = 0; brow < TRIPLES_DIM; ++brow) {
						for (bcol = 0; bcol < TRIPLES_DIM; ++bcol) {
							size_t idx, val;

							idx = (mrow * op[i].dcol + mcol) * TRIPLES_PACK + ii;
							val = (mrow * op[i].dcol + brow) * op[i].dcol * TRIPLES_DIM + mcol * TRIPLES_DIM + bcol;
							mpz_set_ui(op[i].value[idx].value[brow][bcol], val);
						}
					}
				}
			}
		}
	}
}

void
triples_matrix_init_rand(triples_Matrix *op, size_t len)
{
	size_t i;

	for (i = 0; i < len; ++i)
		block_init_rand(op[i].value, op[i].drow * op[i].dcol * TRIPLES_PACK);
}

void
triples_matrix_init_zero(triples_Matrix *op, size_t len)
{
	size_t i;

	for (i = 0; i < len; ++i)
		block_init_zero(op[i].value, op[i].drow * op[i].dcol * TRIPLES_PACK);
}

void
triples_matrix_init_zero_enc(triples_MatrixEnc *op, size_t len)
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
triples_matrix_mac(triples_MatrixEnc *rop, triples_MatrixEnc *op)
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
				bgv_mul(i, &rop->value[idx].value[i], &op->value[idx].value[i], &ctMAC[0][i]);
			}
		}
	}
}

int
triples_matrix_mac_verify(const triples_Matrix *mac, const triples_Matrix *op)
{
	triples_Block *tmp;
	size_t row, col;
	int cmp;

	assert(mac->drow == op->drow);
	assert(mac->dcol == op->dcol);

	tmp = block_alloc(TRIPLES_PACK);
	cmp = 0;

	for (row = 0; row < mac->drow; ++row) {
		for (col = 0; col < mac->dcol; ++col) {
			size_t idx;

			idx = (row * mac->dcol + col) * TRIPLES_PACK;
			block_mulc(tmp, &op->value[idx], MAC[0], TRIPLES_PACK);
			cmp += block_cmp(&mac->value[idx], tmp, 1);
		}
	}

	block_dealloc(tmp, TRIPLES_PACK);

	return cmp;
}

int
triples_matrix_mac_verifyA(const triples_Matrix *mac, const triples_Matrix *op)
{
	tiimat3_Message *m;
	triples_Block *tmp;
	size_t row, col;
	int cmp;

	assert(mac->drow == op->drow);
	assert(mac->dcol == op->dcol);

	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	m = tiimat3_alloc_msg(1);
	tmp = block_alloc(TRIPLES_PACK);
	cmp = 0;

	for (row = 0; row < mac->drow; ++row) {
		for (col = 0; col < mac->dcol; ++col) {
			size_t idx;

			idx = (row * mac->dcol + col) * TRIPLES_PACK;
			block_mulc(tmp, &op->value[idx], MAC[0], TRIPLES_PACK);

			#if TRIPLES_HOIST || TRIPLES_MIXED
			block_pack(m, tmp, 0, block_sigma);
			#else
			block_pack(m, tmp, TIIMAT3_D / 2 - TRIPLES_ROTB, block_sigma);
			#endif /* TRIPLES_HOIST */

			block_unpack(tmp, m);
			cmp += block_cmp(&mac->value[idx], tmp, 1);
		}
	}

	tiimat3_dealloc_msg(m, 1);
	block_dealloc(tmp, TRIPLES_PACK);

	return cmp;
}

int
triples_matrix_mac_verifyB(const triples_Matrix *mac, const triples_Matrix *op)
{
	tiimat3_Message *m;
	triples_Block *tmp;
	size_t row, col;
	int cmp;

	assert(mac->drow == op->drow);
	assert(mac->dcol == op->dcol);

	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	m = tiimat3_alloc_msg(1);
	tmp = block_alloc(TRIPLES_PACK);
	cmp = 0;

	for (row = 0; row < mac->drow; ++row) {
		for (col = 0; col < mac->dcol; ++col) {
			size_t idx;

			idx = (row * mac->dcol + col) * TRIPLES_PACK;
			block_mulc(tmp, &op->value[idx], MAC[0], TRIPLES_PACK);
			block_pack(m, tmp, 0, block_tau);
			block_unpack(tmp, m);
			cmp += block_cmp(&mac->value[idx], tmp, 1);
		}
	}

	tiimat3_dealloc_msg(m, 1);
	block_dealloc(tmp, TRIPLES_PACK);

	return cmp;
}

void
triples_matrix_mul(triples_Matrix *rop, const triples_Matrix *op1, const triples_Matrix *op2)
{
	triples_Block *tmp;
	size_t row, col, k;

	assert(rop->drow == op1->drow);
	assert(rop->dcol == op2->dcol);
	assert(op1->dcol == op2->drow);

	tmp = block_alloc(TRIPLES_PACK);

	triples_matrix_init_zero(rop, 1);
	for (row = 0; row < rop->drow; ++row) {
		for (col = 0; col < rop->dcol; ++col) {
			size_t idx;

			idx = (row * rop->dcol + col) * TRIPLES_PACK;
			for (k = 0; k < op1->dcol; ++k) {
				size_t idx1, idx2;

				idx1 = (row * op1->dcol + k) * TRIPLES_PACK;
				idx2 = (k * op2->dcol + col) * TRIPLES_PACK;

				block_mul(tmp, &op1->value[idx1], &op2->value[idx2], TRIPLES_PACK);
				block_add(&rop->value[idx], &rop->value[idx], tmp, TRIPLES_PACK);
			}
		}
	}

	block_dealloc(tmp, TRIPLES_PACK);
}

void
triples_matrix_mul_enc(triples_MatrixEnc *rop, triples_MatrixEnc *op1, triples_MatrixEnc *op2)
{
	triples_MatrixEnc *rotA1, *rotA2, *rotB;
	triples_BlockEnc *tmp;
	size_t row, col, k;

	assert(rop->drow == op1->drow);
	assert(rop->dcol == op2->dcol);
	assert(op1->dcol == op2->drow);

	tmp = bgv_alloc(1, sizeof *tmp);

	#if TRIPLES_PREROT
	rotA1 = matrix_alloc_encdim(op1->drow, op1->dcol, 1);
	rotA2 = matrix_alloc_encdim(op1->drow, op1->dcol, 1);
	for (row = 0; row < op1->drow; ++row) {
		for (col = 0; col < op1->dcol; ++col) {
			size_t idx = row * op1->dcol + col;
			size_t dim = idx * TRIPLES_DIM;

			#if TRIPLES_HOIST
			block_prerotA_enc_hoisted(&rotA1->value[dim], &rotA2->value[dim], &op1->value[idx]);
			#elif TRIPLES_MIXED
			block_prerotA_enc_mixed(&rotA1->value[dim], &rotA2->value[dim], &op1->value[idx]);
			#else
			block_prerotA_enc_reuse(&rotA1->value[dim], &rotA2->value[dim], &op1->value[idx]);
			#endif /* TRIPLES_HOIST */
		}
	}

	rotB = matrix_alloc_encdim(op2->drow, op2->dcol, 1);
	for (row = 0; row < op2->drow; ++row) {
		for (col = 0; col < op2->dcol; ++col) {
			size_t idx = row * op2->dcol + col;
			size_t dim = idx * TRIPLES_DIM;

			#if TRIPLES_HOIST
			block_prerotB_enc_hoisted(&rotB->value[dim], &op2->value[idx]);
			#elif TRIPLES_MIXED
			block_prerotB_enc_mixed(&rotB->value[dim], &op2->value[idx]);
			#else
			block_prerotB_enc_reuse(&rotB->value[dim], &op2->value[idx]);
			#endif /* TRIPLES_HOIST */
		}
	}
	#else
	rotA1 = matrix_alloc_encdim(0, 0, 1);
	rotA2 = matrix_alloc_encdim(0, 0, 1);
	rotB = matrix_alloc_encdim(0, 0, 1);
	#endif /* TRIPLES_PREROT */

	for (row = 0; row < rop->drow; ++row) {
		for (col = 0; col < rop->dcol; ++col) {
			size_t idx;

			idx = row * rop->dcol + col;
			memset(&rop->value[idx], 0, sizeof rop->value[idx]);
			for (k = 0; k < op1->dcol; ++k) {
				size_t idx1, idx2;

				idx1 = row * op1->dcol + k;
				idx2 = k * op2->dcol + col;

				#if TRIPLES_PREROT
					idx1 *= TRIPLES_DIM;
					idx2 *= TRIPLES_DIM;
					block_mul_enc_prerot(tmp, &rotA1->value[idx1], &rotA2->value[idx1], &rotB->value[idx2]);
				#else
					#if TRIPLES_HOIST
					block_mul_enc_hoisted(tmp, &op1->value[idx1], &op2->value[idx2]);
					#else
					block_mul_enc_reuse(tmp, &op1->value[idx1], &op2->value[idx2]);
					#endif /* TRIPLES_HOIST */
				#endif /* TRIPLES_PREROT */

				block_add_enc(&rop->value[idx], &rop->value[idx], tmp);
			}
		}
	}

	triples_matrix_dealloc_enc(rotA1, 1);
	triples_matrix_dealloc_enc(rotA2, 1);
	triples_matrix_dealloc_enc(rotB, 1);
	bgv_dealloc(tmp);
}

void
triples_matrix_print(const triples_Matrix *op, size_t pack)
{
	size_t mrow, mcol, brow, bcol;

	for (mrow = 0; mrow < op->drow; ++mrow) {
		for (brow = 0; brow < TRIPLES_DIM; ++brow) {
			for (mcol = 0; mcol < op->dcol; ++mcol) {
				for (bcol = 0; bcol < TRIPLES_DIM; ++bcol) {
					size_t idx;

					idx = (mrow * op->dcol + mcol) * TRIPLES_PACK + pack;
					gmp_printf("%5Zd ", op->value[idx].value[brow][bcol]);
				}
			}
			puts("");
		}
	}
}

#endif /* TRIPLES_IMPL */
