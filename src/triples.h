#ifndef TIIMAT3_H
#define TIIMAT3_H

#include "bgv.h"

#ifndef TIIMAT3_HOIST
#define TIIMAT3_HOIST 1
#endif /* TIIMAT3_HOIST */

#ifndef TIIMAT3_MIXED
#define TIIMAT3_MIXED 0
#endif /* TIIMAT3_MIXED */

#ifndef TIIMAT3_PREROT
#define TIIMAT3_PREROT 1
#endif /* TIIMAT3_PREROT */

#ifndef TIIMAT3_DIM
#define TIIMAT3_DIM 32
#endif /* TIIMAT3_DIM */

#define TIIMAT3_FLOOD_BITS  80
#define TIIMAT3_SHARES      5
#define TIIMAT3_THREADS     16

#define TIIMAT3_FLOOD_BYTES ((TIIMAT3_FLOOD_BITS - 1) / CHAR_BIT + 1)
#define TIIMAT3_PACK        (TIIMAT3_D / (TIIMAT3_DIM * TIIMAT3_DIM))
#define TIIMAT3_ROTA        (TIIMAT3_PACK / 2)
#define TIIMAT3_ROTB        (TIIMAT3_PACK / 2 * TIIMAT3_DIM)

typedef struct tiimat3_block             tiimat3_Block;
typedef struct tiimat3_block_ciphertext  tiimat3_BlockEnc;
typedef struct tiimat3_matrix_block      tiimat3_Matrix;
typedef struct tiimat3_matrix_ciphertext tiimat3_MatrixEnc;

struct tiimat3_block {
	mpz_t value[TIIMAT3_DIM][TIIMAT3_DIM];
};

struct tiimat3_block_ciphertext {
	tiimat3_Ciphertext value[TIIMAT3_QPLEN];
};

struct tiimat3_matrix_block {
	tiimat3_Block *value;
	size_t drow, dcol;
};

struct tiimat3_matrix_ciphertext {
	tiimat3_BlockEnc *value;
	size_t drow, dcol;
};

void tiimat3_block_add_enc(tiimat3_BlockEnc *rop, tiimat3_BlockEnc *op1, tiimat3_BlockEnc *op2);
void tiimat3_block_mul_enc_hoisted(tiimat3_BlockEnc *rop, tiimat3_BlockEnc *A, tiimat3_BlockEnc *B);
void tiimat3_block_mul_enc_prerot(tiimat3_BlockEnc *rop, tiimat3_BlockEnc *A1, tiimat3_BlockEnc *A2, tiimat3_BlockEnc *B);
void tiimat3_block_mul_enc_reuse(tiimat3_BlockEnc *rop, tiimat3_BlockEnc *A, tiimat3_BlockEnc *B);
void tiimat3_block_prerotA_enc_hoisted(tiimat3_BlockEnc *rotA1, tiimat3_BlockEnc *rotA2, tiimat3_BlockEnc *A);
void tiimat3_block_prerotA_enc_mixed(tiimat3_BlockEnc *rotA1, tiimat3_BlockEnc *rotA2, tiimat3_BlockEnc *A);
void tiimat3_block_prerotA_enc_reuse(tiimat3_BlockEnc *rotA1, tiimat3_BlockEnc *rotA2, tiimat3_BlockEnc *A);
void tiimat3_block_prerotB_enc_hoisted(tiimat3_BlockEnc *rotB, tiimat3_BlockEnc *B);
void tiimat3_block_prerotB_enc_mixed(tiimat3_BlockEnc *rotB, tiimat3_BlockEnc *B);
void tiimat3_block_prerotB_enc_reuse(tiimat3_BlockEnc *rotB, tiimat3_BlockEnc *B);

void tiimat3_deinit(void);
void tiimat3_deinit_hoisted(void);
void tiimat3_deinit_mixed(void);
void tiimat3_init(void);
void tiimat3_init_hoisted(void);
void tiimat3_init_mixed(void);

void tiimat3_matrix_add(tiimat3_Matrix *rop, const tiimat3_Matrix *op1, const tiimat3_Matrix *op2);
void tiimat3_matrix_add_enc(tiimat3_MatrixEnc *rop, tiimat3_MatrixEnc *op1, tiimat3_MatrixEnc *op2);
tiimat3_Matrix *tiimat3_matrix_alloc(size_t drow, size_t dcol, size_t len);
tiimat3_MatrixEnc *tiimat3_matrix_alloc_enc(size_t drow, size_t dcol, size_t len);
tiimat3_MatrixEnc *tiimat3_matrix_alloc_encdim(size_t drow, size_t dcol, size_t len);
int  tiimat3_matrix_cmp(const tiimat3_Matrix *op1, const tiimat3_Matrix *op2);
void tiimat3_matrix_dealloc(tiimat3_Matrix *op, size_t len);
void tiimat3_matrix_dealloc_enc(tiimat3_MatrixEnc *op, size_t len);
void tiimat3_matrix_decrypt(tiimat3_Matrix *rop, tiimat3_MatrixEnc *op);
void tiimat3_matrix_encrypt(tiimat3_MatrixEnc *rop, const tiimat3_Matrix *op);
void tiimat3_matrix_encryptA(tiimat3_MatrixEnc *rop, const tiimat3_Matrix *op);
void tiimat3_matrix_encryptB(tiimat3_MatrixEnc *rop, const tiimat3_Matrix *op);
void tiimat3_matrix_init_asc(tiimat3_Matrix *op, size_t len);
void tiimat3_matrix_init_rand(tiimat3_Matrix *op, size_t len);
void tiimat3_matrix_init_zero(tiimat3_Matrix *op, size_t len);
void tiimat3_matrix_init_zero_enc(tiimat3_MatrixEnc *op, size_t len);
void tiimat3_matrix_mac(tiimat3_MatrixEnc *rop, tiimat3_MatrixEnc *op);
int  tiimat3_matrix_mac_verify(const tiimat3_Matrix *mac, const tiimat3_Matrix *op);
int  tiimat3_matrix_mac_verifyA(const tiimat3_Matrix *mac, const tiimat3_Matrix *op);
int  tiimat3_matrix_mac_verifyB(const tiimat3_Matrix *mac, const tiimat3_Matrix *op);
void tiimat3_matrix_mul(tiimat3_Matrix *rop, const tiimat3_Matrix *op1, const tiimat3_Matrix *op2);
void tiimat3_matrix_mul_enc(tiimat3_MatrixEnc *rop, tiimat3_MatrixEnc *op1, tiimat3_MatrixEnc *op2);
void tiimat3_matrix_print(const tiimat3_Matrix *op, size_t len);

#endif /* TIIMAT3_H */
