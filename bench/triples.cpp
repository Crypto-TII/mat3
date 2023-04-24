#include <benchmark/benchmark.h>
#include <gmp.h>

#define TRIPLES_DIM 128

extern "C" {
#define TRIPLES_IMPL
#include "triples.h"
}


static void
BM_triples(benchmark::State& state)
{
	const size_t dim = TRIPLES_BLOCKS;

	auto A = triples_matrix_alloc(dim, dim, TRIPLES_SHARES);
	auto B = triples_matrix_alloc(dim, dim, TRIPLES_SHARES);

	triples_init();
	#if TRIPLES_HOIST
	triples_init_hoisted();
	#elif TRIPLES_MIXED
	triples_init_mixed();
	#endif /* TRIPLES_HOIST */

	auto encA = triples_matrix_alloc_enc(dim, dim, TRIPLES_SHARES);
	triples_matrix_init_rand(&A[1], TRIPLES_SHARES - 1);
	for (size_t i = 1; i < TRIPLES_SHARES; ++i)
		triples_matrix_encryptA(&encA[i], &A[i]);

	auto encB = triples_matrix_alloc_enc(dim, dim, TRIPLES_SHARES);
	triples_matrix_init_rand(&B[1], TRIPLES_SHARES - 1);
	for (size_t i = 1; i < TRIPLES_SHARES; ++i)
		triples_matrix_encryptB(&encB[i], &B[i]);

	auto AB = triples_matrix_alloc(dim, dim, 1);
	auto macA = triples_matrix_alloc(dim, dim, 1);
	auto macB = triples_matrix_alloc(dim, dim, 1);
	auto macAB = triples_matrix_alloc(dim, dim, 1);
	auto encAB = triples_matrix_alloc_enc(dim, dim, 1);
	auto encmacA = triples_matrix_alloc_enc(dim, dim, 1);
	auto encmacB = triples_matrix_alloc_enc(dim, dim, 1);
	auto encmacAB = triples_matrix_alloc_enc(dim, dim, 1);

	for (auto _ : state) {
		triples_matrix_init_rand(&A[0], 1);
		triples_matrix_encryptA(&encA[0], &A[0]);
		for (size_t i = 1; i < TRIPLES_SHARES; ++i)
			triples_matrix_add_enc(&encA[0], &encA[0], &encA[i]);
		triples_matrix_mac(encmacA, encA);
		triples_matrix_decrypt(macA, encmacA);

		triples_matrix_init_rand(&B[0], 1);
		triples_matrix_encryptA(&encB[0], &B[0]);
		for (size_t i = 1; i < TRIPLES_SHARES; ++i)
			triples_matrix_add_enc(&encB[0], &encB[0], &encB[i]);
		triples_matrix_mac(encmacB, encB);
		triples_matrix_decrypt(macB, encmacB);

		triples_matrix_mul_enc(encAB, encA, encB);
		triples_matrix_decrypt(AB, encAB);

		triples_matrix_mac(encmacAB, encAB);
		triples_matrix_decrypt(macAB, encmacAB);
	}

	triples_matrix_dealloc(A, TRIPLES_SHARES);
	triples_matrix_dealloc(B, TRIPLES_SHARES);
	triples_matrix_dealloc(AB, 1);
	triples_matrix_dealloc(macA, 1);
	triples_matrix_dealloc(macB, 1);
	triples_matrix_dealloc(macAB, 1);

	triples_matrix_dealloc_enc(encA, TRIPLES_SHARES);
	triples_matrix_dealloc_enc(encB, TRIPLES_SHARES);
	triples_matrix_dealloc_enc(encAB, 1);
	triples_matrix_dealloc_enc(encmacA, 1);
	triples_matrix_dealloc_enc(encmacB, 1);
	triples_matrix_dealloc_enc(encmacAB, 1);
}

BENCHMARK(BM_triples)->Unit(benchmark::kSecond);
BENCHMARK_MAIN();
