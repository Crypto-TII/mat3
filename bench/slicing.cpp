#define TRIPLES_HOIST  1
#define TRIPLES_PREROT 1

#include <benchmark/benchmark.h>
#include <gmp.h>

extern "C" {
#define TRIPLES_IMPL
#include "triples.h"
}

#define CEIL(a, b) (((a) - 1) / (b) + 1)


static void
BM_slicing(benchmark::State& state)
{
	const size_t dim = CEIL(192, TRIPLES_DIM);

	auto A = triples_matrix_alloc(dim, dim, TRIPLES_SHARES);
	auto B = triples_matrix_alloc(dim, dim, TRIPLES_SHARES);

	triples_init();
	triples_init_hoisted();

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

BENCHMARK(BM_slicing)->Unit(benchmark::kSecond);
BENCHMARK_MAIN();
