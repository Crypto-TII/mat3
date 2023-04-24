#define TRIPLES_IMPL
#include "triples.h"

#define DIM1 1
#define DIM2 2
#define DIM3 3

void
test_mat3(void)
{
	triples_Matrix *A, *B, *C, *AB, *macA, *macB, *macAB;
	triples_MatrixEnc *encA, *encB, *encAB, *encmacA, *encmacB, *encmacAB;

	A = triples_matrix_alloc(DIM1, DIM2, 1);
	B = triples_matrix_alloc(DIM2, DIM3, 1);
	C = triples_matrix_alloc(DIM1, DIM3, 1);
	AB = triples_matrix_alloc(DIM1, DIM3, 1);
	macA = triples_matrix_alloc(DIM1, DIM2, 1);
	macB = triples_matrix_alloc(DIM2, DIM3, 1);
	macAB = triples_matrix_alloc(DIM1, DIM3, 1);
	encA = triples_matrix_alloc_enc(DIM1, DIM2, 1);
	encB = triples_matrix_alloc_enc(DIM2, DIM3, 1);
	encAB = triples_matrix_alloc_enc(DIM1, DIM3, 1);
	encmacA = triples_matrix_alloc_enc(DIM1, DIM2, 1);
	encmacB = triples_matrix_alloc_enc(DIM2, DIM3, 1);
	encmacAB = triples_matrix_alloc_enc(DIM1, DIM3, 1);

	triples_matrix_init_rand(A, 1);
	triples_matrix_init_rand(B, 1);
	triples_matrix_mul(C, A, B);

	triples_matrix_encryptA(encA, A);
	triples_matrix_mac(encmacA, encA);
	triples_matrix_decrypt(macA, encmacA);
	assert(triples_matrix_mac_verifyA(macA, A) == 0);

	triples_matrix_encryptB(encB, B);
	triples_matrix_mac(encmacB, encB);
	triples_matrix_decrypt(macB, encmacB);
	assert(triples_matrix_mac_verifyB(macB, B) == 0);

	triples_matrix_mul_enc(encAB, encA, encB);
	triples_matrix_decrypt(AB, encAB);
	assert(triples_matrix_cmp(C, AB) == 0);

	triples_matrix_mac(encmacAB, encAB);
	triples_matrix_decrypt(macAB, encmacAB);
	assert(triples_matrix_mac_verify(macAB, AB) == 0);

	triples_matrix_dealloc(A, 1);
	triples_matrix_dealloc(B, 1);
	triples_matrix_dealloc(C, 1);
	triples_matrix_dealloc(AB, 1);
	triples_matrix_dealloc(macA, 1);
	triples_matrix_dealloc(macB, 1);
	triples_matrix_dealloc(macAB, 1);
	triples_matrix_dealloc_enc(encA, 1);
	triples_matrix_dealloc_enc(encB, 1);
	triples_matrix_dealloc_enc(encAB, 1);
	triples_matrix_dealloc_enc(encmacA, 1);
	triples_matrix_dealloc_enc(encmacB, 1);
	triples_matrix_dealloc_enc(encmacAB, 1);
}

void
test_mat3_shared(void)
{
	triples_Matrix *A, *B, *C, *AB, *macA, *macB, *macAB;
	triples_MatrixEnc *encA, *encB, *encAB, *encmacA, *encmacB, *encmacAB;
	size_t i;

	A = triples_matrix_alloc(DIM1, DIM2, TRIPLES_SHARES + 1);
	B = triples_matrix_alloc(DIM2, DIM3, TRIPLES_SHARES + 1);
	C = triples_matrix_alloc(DIM1, DIM3, 1);
	AB = triples_matrix_alloc(DIM1, DIM3, 1);
	macA = triples_matrix_alloc(DIM1, DIM2, 1);
	macB = triples_matrix_alloc(DIM2, DIM3, 1);
	macAB = triples_matrix_alloc(DIM1, DIM3, 1);
	encA = triples_matrix_alloc_enc(DIM1, DIM2, TRIPLES_SHARES + 1);
	encB = triples_matrix_alloc_enc(DIM2, DIM3, TRIPLES_SHARES + 1);
	encAB = triples_matrix_alloc_enc(DIM1, DIM3, 1);
	encmacA = triples_matrix_alloc_enc(DIM1, DIM2, 1);
	encmacB = triples_matrix_alloc_enc(DIM2, DIM3, 1);
	encmacAB = triples_matrix_alloc_enc(DIM1, DIM3, 1);

	triples_matrix_init_zero(&A[0], 1);
	triples_matrix_init_rand(&A[1], TRIPLES_SHARES);
	for (i = 1; i <= TRIPLES_SHARES; ++i)
		triples_matrix_add(&A[0], &A[0], &A[i]);
	triples_matrix_init_zero(&B[0], 1);
	triples_matrix_init_rand(&B[1], TRIPLES_SHARES);
	for (i = 1; i <= TRIPLES_SHARES; ++i)
		triples_matrix_add(&B[0], &B[0], &B[i]);
	triples_matrix_mul(C, A, B);

	triples_matrix_init_zero_enc(&encA[0], 1);
	for (i = 1; i <= TRIPLES_SHARES; ++i) {
		triples_matrix_encryptA(&encA[i], &A[i]);
		triples_matrix_add_enc(&encA[0], &encA[0], &encA[i]);
	}
	triples_matrix_mac(encmacA, encA);
	triples_matrix_decrypt(macA, encmacA);
	assert(triples_matrix_mac_verifyA(macA, A) == 0);

	triples_matrix_init_zero_enc(&encB[0], 1);
	for (i = 1; i <= TRIPLES_SHARES; ++i) {
		triples_matrix_encryptB(&encB[i], &B[i]);
		triples_matrix_add_enc(&encB[0], &encB[0], &encB[i]);
	}
	triples_matrix_mac(encmacB, encB);
	triples_matrix_decrypt(macB, encmacB);
	assert(triples_matrix_mac_verifyB(macB, B) == 0);

	triples_matrix_mul_enc(encAB, encA, encB);
	triples_matrix_decrypt(AB, encAB);
	assert(triples_matrix_cmp(C, AB) == 0);

	triples_matrix_mac(encmacAB, encAB);
	triples_matrix_decrypt(macAB, encmacAB);
	assert(triples_matrix_mac_verify(macAB, AB) == 0);

	triples_matrix_dealloc(A, TRIPLES_SHARES + 1);
	triples_matrix_dealloc(B, TRIPLES_SHARES + 1);
	triples_matrix_dealloc(C, 1);
	triples_matrix_dealloc(AB, 1);
	triples_matrix_dealloc(macA, 1);
	triples_matrix_dealloc(macB, 1);
	triples_matrix_dealloc(macAB, 1);
	triples_matrix_dealloc_enc(encA, TRIPLES_SHARES + 1);
	triples_matrix_dealloc_enc(encB, TRIPLES_SHARES + 1);
	triples_matrix_dealloc_enc(encAB, 1);
	triples_matrix_dealloc_enc(encmacA, 1);
	triples_matrix_dealloc_enc(encmacB, 1);
	triples_matrix_dealloc_enc(encmacAB, 1);
}

int
main(void)
{

	triples_init();
	triples_deinit();
	triples_init();

	#if TRIPLES_HOIST
	triples_init_hoisted();
	triples_deinit_hoisted();
	triples_init_hoisted();
	#elif TRIPLES_MIXED
	triples_init_mixed();
	triples_deinit_mixed();
	triples_init_mixed();
	#endif /* TRIPLES_HOIST */

	fputs("[+] Testing mat3:        ", stderr);
	test_mat3(), fputs("4/4.\n", stderr);

	fputs("[+] Testing mat3 shared: ", stderr);
	test_mat3_shared(), fputs("4/4.\n", stderr);

	return 0;
}
