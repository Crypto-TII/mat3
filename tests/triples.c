#include "../src/triples.c"

#define DIM1 1
#define DIM2 2
#define DIM3 3

void
test_mat3(void)
{
	tiimat3_Matrix *A, *B, *C, *AB, *macA, *macB, *macAB;
	tiimat3_MatrixEnc *encA, *encB, *encAB, *encmacA, *encmacB, *encmacAB;

	A = tiimat3_matrix_alloc(DIM1, DIM2, 1);
	B = tiimat3_matrix_alloc(DIM2, DIM3, 1);
	C = tiimat3_matrix_alloc(DIM1, DIM3, 1);
	AB = tiimat3_matrix_alloc(DIM1, DIM3, 1);
	macA = tiimat3_matrix_alloc(DIM1, DIM2, 1);
	macB = tiimat3_matrix_alloc(DIM2, DIM3, 1);
	macAB = tiimat3_matrix_alloc(DIM1, DIM3, 1);
	encA = tiimat3_matrix_alloc_enc(DIM1, DIM2, 1);
	encB = tiimat3_matrix_alloc_enc(DIM2, DIM3, 1);
	encAB = tiimat3_matrix_alloc_enc(DIM1, DIM3, 1);
	encmacA = tiimat3_matrix_alloc_enc(DIM1, DIM2, 1);
	encmacB = tiimat3_matrix_alloc_enc(DIM2, DIM3, 1);
	encmacAB = tiimat3_matrix_alloc_enc(DIM1, DIM3, 1);

	tiimat3_matrix_init_rand(A, 1);
	tiimat3_matrix_init_rand(B, 1);
	tiimat3_matrix_mul(C, A, B);

	tiimat3_matrix_encryptA(encA, A);
	tiimat3_matrix_mac(encmacA, encA);
	tiimat3_matrix_decrypt(macA, encmacA);
	assert(tiimat3_matrix_mac_verifyA(macA, A) == 0);

	tiimat3_matrix_encryptB(encB, B);
	tiimat3_matrix_mac(encmacB, encB);
	tiimat3_matrix_decrypt(macB, encmacB);
	assert(tiimat3_matrix_mac_verifyB(macB, B) == 0);

	tiimat3_matrix_mul_enc(encAB, encA, encB);
	tiimat3_matrix_decrypt(AB, encAB);
	assert(tiimat3_matrix_cmp(C, AB) == 0);

	tiimat3_matrix_mac(encmacAB, encAB);
	tiimat3_matrix_decrypt(macAB, encmacAB);
	assert(tiimat3_matrix_mac_verify(macAB, AB) == 0);

	tiimat3_matrix_dealloc(A, 1);
	tiimat3_matrix_dealloc(B, 1);
	tiimat3_matrix_dealloc(C, 1);
	tiimat3_matrix_dealloc(AB, 1);
	tiimat3_matrix_dealloc(macA, 1);
	tiimat3_matrix_dealloc(macB, 1);
	tiimat3_matrix_dealloc(macAB, 1);
	tiimat3_matrix_dealloc_enc(encA, 1);
	tiimat3_matrix_dealloc_enc(encB, 1);
	tiimat3_matrix_dealloc_enc(encAB, 1);
	tiimat3_matrix_dealloc_enc(encmacA, 1);
	tiimat3_matrix_dealloc_enc(encmacB, 1);
	tiimat3_matrix_dealloc_enc(encmacAB, 1);
}

void
test_mat3_shared(void)
{
	tiimat3_Matrix *A, *B, *C, *AB, *macA, *macB, *macAB;
	tiimat3_MatrixEnc *encA, *encB, *encAB, *encmacA, *encmacB, *encmacAB;
	size_t i;

	A = tiimat3_matrix_alloc(DIM1, DIM2, TIIMAT3_SHARES + 1);
	B = tiimat3_matrix_alloc(DIM2, DIM3, TIIMAT3_SHARES + 1);
	C = tiimat3_matrix_alloc(DIM1, DIM3, 1);
	AB = tiimat3_matrix_alloc(DIM1, DIM3, 1);
	macA = tiimat3_matrix_alloc(DIM1, DIM2, 1);
	macB = tiimat3_matrix_alloc(DIM2, DIM3, 1);
	macAB = tiimat3_matrix_alloc(DIM1, DIM3, 1);
	encA = tiimat3_matrix_alloc_enc(DIM1, DIM2, TIIMAT3_SHARES + 1);
	encB = tiimat3_matrix_alloc_enc(DIM2, DIM3, TIIMAT3_SHARES + 1);
	encAB = tiimat3_matrix_alloc_enc(DIM1, DIM3, 1);
	encmacA = tiimat3_matrix_alloc_enc(DIM1, DIM2, 1);
	encmacB = tiimat3_matrix_alloc_enc(DIM2, DIM3, 1);
	encmacAB = tiimat3_matrix_alloc_enc(DIM1, DIM3, 1);

	tiimat3_matrix_init_zero(&A[0], 1);
	tiimat3_matrix_init_rand(&A[1], TIIMAT3_SHARES);
	for (i = 1; i <= TIIMAT3_SHARES; ++i)
		tiimat3_matrix_add(&A[0], &A[0], &A[i]);
	tiimat3_matrix_init_zero(&B[0], 1);
	tiimat3_matrix_init_rand(&B[1], TIIMAT3_SHARES);
	for (i = 1; i <= TIIMAT3_SHARES; ++i)
		tiimat3_matrix_add(&B[0], &B[0], &B[i]);
	tiimat3_matrix_mul(C, A, B);

	tiimat3_matrix_init_zero_enc(&encA[0], 1);
	for (i = 1; i <= TIIMAT3_SHARES; ++i) {
		tiimat3_matrix_encryptA(&encA[i], &A[i]);
		tiimat3_matrix_add_enc(&encA[0], &encA[0], &encA[i]);
	}
	tiimat3_matrix_mac(encmacA, encA);
	tiimat3_matrix_decrypt(macA, encmacA);
	assert(tiimat3_matrix_mac_verifyA(macA, A) == 0);

	tiimat3_matrix_init_zero_enc(&encB[0], 1);
	for (i = 1; i <= TIIMAT3_SHARES; ++i) {
		tiimat3_matrix_encryptB(&encB[i], &B[i]);
		tiimat3_matrix_add_enc(&encB[0], &encB[0], &encB[i]);
	}
	tiimat3_matrix_mac(encmacB, encB);
	tiimat3_matrix_decrypt(macB, encmacB);
	assert(tiimat3_matrix_mac_verifyB(macB, B) == 0);

	tiimat3_matrix_mul_enc(encAB, encA, encB);
	tiimat3_matrix_decrypt(AB, encAB);
	assert(tiimat3_matrix_cmp(C, AB) == 0);

	tiimat3_matrix_mac(encmacAB, encAB);
	tiimat3_matrix_decrypt(macAB, encmacAB);
	assert(tiimat3_matrix_mac_verify(macAB, AB) == 0);

	tiimat3_matrix_dealloc(A, TIIMAT3_SHARES + 1);
	tiimat3_matrix_dealloc(B, TIIMAT3_SHARES + 1);
	tiimat3_matrix_dealloc(C, 1);
	tiimat3_matrix_dealloc(AB, 1);
	tiimat3_matrix_dealloc(macA, 1);
	tiimat3_matrix_dealloc(macB, 1);
	tiimat3_matrix_dealloc(macAB, 1);
	tiimat3_matrix_dealloc_enc(encA, TIIMAT3_SHARES + 1);
	tiimat3_matrix_dealloc_enc(encB, TIIMAT3_SHARES + 1);
	tiimat3_matrix_dealloc_enc(encAB, 1);
	tiimat3_matrix_dealloc_enc(encmacA, 1);
	tiimat3_matrix_dealloc_enc(encmacB, 1);
	tiimat3_matrix_dealloc_enc(encmacAB, 1);
}

int
main(void)
{

	tiimat3_init();
	tiimat3_deinit();
	tiimat3_init();

	#if TIIMAT3_HOIST
	tiimat3_init_hoisted();
	tiimat3_deinit_hoisted();
	tiimat3_init_hoisted();
	#elif TIIMAT3_MIXED
	tiimat3_init_mixed();
	tiimat3_deinit_mixed();
	tiimat3_init_mixed();
	#endif /* TIIMAT3_HOIST */

	fputs("[+] Testing mat3:        ", stderr);
	test_mat3(), fputs("4/4.\n", stderr);

	fputs("[+] Testing mat3 shared: ", stderr);
	test_mat3_shared(), fputs("4/4.\n", stderr);

	return 0;
}
