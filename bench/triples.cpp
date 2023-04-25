#define TIIMAT3_DIM 128

#include <benchmark/benchmark.h>
#include <gmp.h>

extern "C" {
	#include "triples.h"
}


static void
BM_triples(benchmark::State& state)
{
	const size_t dim = TIIMAT3_BLOCKS;

	auto A = tiimat3_matrix_alloc(dim, dim, TIIMAT3_SHARES);
	auto B = tiimat3_matrix_alloc(dim, dim, TIIMAT3_SHARES);

	tiimat3_init();
	#if TIIMAT3_HOIST
	tiimat3_init_hoisted();
	#elif TIIMAT3_MIXED
	tiimat3_init_mixed();
	#endif /* TIIMAT3_HOIST */

	auto encA = tiimat3_matrix_alloc_enc(dim, dim, TIIMAT3_SHARES);
	tiimat3_matrix_init_rand(&A[1], TIIMAT3_SHARES - 1);
	for (size_t i = 1; i < TIIMAT3_SHARES; ++i)
		tiimat3_matrix_encryptA(&encA[i], &A[i]);

	auto encB = tiimat3_matrix_alloc_enc(dim, dim, TIIMAT3_SHARES);
	tiimat3_matrix_init_rand(&B[1], TIIMAT3_SHARES - 1);
	for (size_t i = 1; i < TIIMAT3_SHARES; ++i)
		tiimat3_matrix_encryptB(&encB[i], &B[i]);

	auto AB = tiimat3_matrix_alloc(dim, dim, 1);
	auto macA = tiimat3_matrix_alloc(dim, dim, 1);
	auto macB = tiimat3_matrix_alloc(dim, dim, 1);
	auto macAB = tiimat3_matrix_alloc(dim, dim, 1);
	auto encAB = tiimat3_matrix_alloc_enc(dim, dim, 1);
	auto encmacA = tiimat3_matrix_alloc_enc(dim, dim, 1);
	auto encmacB = tiimat3_matrix_alloc_enc(dim, dim, 1);
	auto encmacAB = tiimat3_matrix_alloc_enc(dim, dim, 1);

	for (auto _ : state) {
		tiimat3_matrix_init_rand(&A[0], 1);
		tiimat3_matrix_encryptA(&encA[0], &A[0]);
		for (size_t i = 1; i < TIIMAT3_SHARES; ++i)
			tiimat3_matrix_add_enc(&encA[0], &encA[0], &encA[i]);
		tiimat3_matrix_mac(encmacA, encA);
		tiimat3_matrix_decrypt(macA, encmacA);

		tiimat3_matrix_init_rand(&B[0], 1);
		tiimat3_matrix_encryptA(&encB[0], &B[0]);
		for (size_t i = 1; i < TIIMAT3_SHARES; ++i)
			tiimat3_matrix_add_enc(&encB[0], &encB[0], &encB[i]);
		tiimat3_matrix_mac(encmacB, encB);
		tiimat3_matrix_decrypt(macB, encmacB);

		tiimat3_matrix_mul_enc(encAB, encA, encB);
		tiimat3_matrix_decrypt(AB, encAB);

		tiimat3_matrix_mac(encmacAB, encAB);
		tiimat3_matrix_decrypt(macAB, encmacAB);
	}

	tiimat3_matrix_dealloc(A, TIIMAT3_SHARES);
	tiimat3_matrix_dealloc(B, TIIMAT3_SHARES);
	tiimat3_matrix_dealloc(AB, 1);
	tiimat3_matrix_dealloc(macA, 1);
	tiimat3_matrix_dealloc(macB, 1);
	tiimat3_matrix_dealloc(macAB, 1);

	tiimat3_matrix_dealloc_enc(encA, TIIMAT3_SHARES);
	tiimat3_matrix_dealloc_enc(encB, TIIMAT3_SHARES);
	tiimat3_matrix_dealloc_enc(encAB, 1);
	tiimat3_matrix_dealloc_enc(encmacA, 1);
	tiimat3_matrix_dealloc_enc(encmacB, 1);
	tiimat3_matrix_dealloc_enc(encmacAB, 1);
}

BENCHMARK(BM_triples)->Unit(benchmark::kSecond);
BENCHMARK_MAIN();
