#define TRIPLES_HOIST  1
#define TRIPLES_PREROT 1
#define TRIPLES_DIM    128

#include <benchmark/benchmark.h>
#include <gmp.h>

extern "C" {
#define TRIPLES_IMPL
#include "triples.h"
}


static void
BM_init(benchmark::State& state)
{
	for (auto _ : state) {
		state.ResumeTiming();
		triples_init();

		state.PauseTiming();
		triples_deinit();
	}
}

static void
BM_init_hoisted(benchmark::State& state)
{

	triples_init();
	triples_init_hoisted();

	for (auto _ : state) {
		state.PauseTiming();
		triples_deinit();
		triples_deinit_hoisted();

		state.ResumeTiming();
		triples_init();
		triples_init_hoisted();
	}
}

static void
BM_encryptA(benchmark::State& state)
{
	const size_t dim = state.range(0);
	auto mat = triples_matrix_alloc(dim, dim, 1);
	auto enc = triples_matrix_alloc_enc(dim, dim, 1);

	for (auto _ : state) {
		state.PauseTiming();
		triples_matrix_init_rand(mat, 1);

		state.ResumeTiming();
		triples_matrix_encryptA(enc, mat);
	}

	triples_matrix_dealloc(mat, 1);
	triples_matrix_dealloc_enc(enc, 1);
}

static void
BM_encryptB(benchmark::State& state)
{
	const size_t dim = state.range(0);
	auto mat = triples_matrix_alloc(dim, dim, 1);
	auto enc = triples_matrix_alloc_enc(dim, dim, 1);

	for (auto _ : state) {
		state.PauseTiming();
		triples_matrix_init_rand(mat, 1);

		state.ResumeTiming();
		triples_matrix_encryptB(enc, mat);
	}

	triples_matrix_dealloc(mat, 1);
	triples_matrix_dealloc_enc(enc, 1);
}

static void
BM_permuteA(benchmark::State& state)
{
	const size_t dim = state.range(0);
	auto mat = triples_matrix_alloc(dim, dim, 1);
	auto enc = triples_matrix_alloc_enc(dim, dim, 1);
	auto rot1 = matrix_alloc_encdim(dim, dim, 1);
	auto rot2 = matrix_alloc_encdim(dim, dim, 1);

	for (auto _ : state) {
		state.PauseTiming();
		triples_matrix_init_rand(mat, 1);
		triples_matrix_encryptA(enc, mat);

		state.ResumeTiming();
		for (size_t row = 0; row < dim; ++row) {
			for (size_t col = 0; col < dim; ++col) {
				size_t idx = row * dim + col;
				size_t idxdim = idx * TRIPLES_DIM;

				#if TRIPLES_HOIST
				block_prerotA_enc_hoisted(&rot1->value[idxdim], &rot2->value[idxdim], &enc->value[idx]);
				#else
				block_prerotA_enc(&rot1->value[idxdim], &rot2->value[idxdim], &enc->value[idx]);
				#endif /* TRIPLES_HOIST */
			}
		}
	}

	triples_matrix_dealloc(mat, 1);
	triples_matrix_dealloc_enc(enc, 1);
	triples_matrix_dealloc_enc(rot1, 1);
	triples_matrix_dealloc_enc(rot2, 1);
}

static void
BM_permuteB(benchmark::State& state)
{
	const size_t dim = state.range(0);
	auto mat = triples_matrix_alloc(dim, dim, 1);
	auto enc = triples_matrix_alloc_enc(dim, dim, 1);
	auto rot = matrix_alloc_encdim(dim, dim, 1);

	for (auto _ : state) {
		state.PauseTiming();
		triples_matrix_init_rand(mat, 1);
		triples_matrix_encryptB(enc, mat);

		state.ResumeTiming();
		for (size_t row = 0; row < dim; ++row) {
			for (size_t col = 0; col < dim; ++col) {
				size_t idx = row * dim + col;
				size_t idxdim = idx * TRIPLES_DIM;

				#if TRIPLES_HOIST
				block_prerotB_enc_hoisted(&rot->value[idxdim], &enc->value[idx]);
				#else
				block_prerotB_enc(&rot->value[idxdim], &enc->value[idx]);
				#endif /* TRIPLES_HOIST */
			}
		}
	}

	triples_matrix_dealloc(mat, 1);
	triples_matrix_dealloc_enc(enc, 1);
	triples_matrix_dealloc_enc(rot, 1);
}

static void
BM_blockmul(benchmark::State& state)
{
	const size_t dim = state.range(0);

	auto A = triples_matrix_alloc(dim, dim, 1);
	auto B = triples_matrix_alloc(dim, dim, 1);
	auto encA = triples_matrix_alloc_enc(dim, dim, 1);
	auto encB = triples_matrix_alloc_enc(dim, dim, 1);
	auto encAB = triples_matrix_alloc_enc(dim, dim, 1);
	auto tmp = (triples_BlockEnc *)bgv_alloc(1, sizeof(triples_BlockEnc));

	#if TRIPLES_PREROT
	auto rotA1 = matrix_alloc_encdim(dim, dim, 1);
	auto rotA2 = matrix_alloc_encdim(dim, dim, 1);
	auto rotB = matrix_alloc_encdim(dim, dim, 1);
	#else
	auto rotA1 = matrix_alloc_encdim(0, 0, 1);
	auto rotA2 = matrix_alloc_encdim(0, 0, 1);
	auto rotB = matrix_alloc_encdim(0, 0, 1);
	#endif /* TRIPLES_PREROT */

	for (auto _ : state) {
		state.PauseTiming();
		triples_matrix_init_rand(A, 1);
		triples_matrix_init_rand(B, 1);
		triples_matrix_encryptA(encA, A);
		triples_matrix_encryptB(encB, B);

		#if TRIPLES_PREROT
		for (size_t row = 0; row < dim; ++row) {
			for (size_t col = 0; col < dim; ++col) {
				size_t idx = row * dim + col;
				size_t idxdim = idx * TRIPLES_DIM;

				#if TRIPLES_HOIST
				block_prerotA_enc_hoisted(&rotA1->value[idxdim], &rotA2->value[idxdim], &encA->value[idx]);
				block_prerotB_enc_hoisted(&rotB->value[idxdim], &encB->value[idx]);
				#else
				block_prerotA_enc(&rotA1->value[idxdim], &rotA2->value[idxdim], &encA->value[idx]);
				block_prerotB_enc(&rot->value[idxdim], &enc->value[idx]);
				#endif /* TRIPLES_HOIST */
			}
		}
		#endif /* TRIPLES_PREROT */

		state.ResumeTiming();
		for (size_t row = 0; row < dim; ++row) {
			for (size_t col = 0; col < dim; ++col) {
				size_t idx = row * dim + col;

				memset(&encAB->value[idx], 0, sizeof encAB->value[idx]);
				for (size_t k = 0; k < dim; ++k) {
					size_t idxA = row * dim + k;
					size_t idxB = k * dim + col;

					#if TRIPLES_PREROT
					idxA *= TRIPLES_DIM;
					idxB *= TRIPLES_DIM;
					block_mul_enc_prerot(tmp, &rotA1->value[idxA], &rotA2->value[idxA], &rotB->value[idxB]);
					#else
					#if TRIPLES_HOIST
					block_mul_enc_hoisted(tmp, &encA->value[idxA], &encB->value[idxB]);
					#else
					block_mul_enc_reuse(tmp, &encA->value[idxA], &encB->value[idxB]);
					#endif /* TRIPLES_HOIST */
					#endif /* TRIPLES_PREROT */

					block_add_enc(&encAB->value[idx], &encAB->value[idx], tmp);
				}
			}
		}
	}

	bgv_dealloc(tmp);
	triples_matrix_dealloc(A, 1);
	triples_matrix_dealloc(B, 1);
	triples_matrix_dealloc_enc(encA, 1);
	triples_matrix_dealloc_enc(encB, 1);
	triples_matrix_dealloc_enc(encAB, 1);
	triples_matrix_dealloc_enc(rotA1, 1);
	triples_matrix_dealloc_enc(rotA2, 1);
	triples_matrix_dealloc_enc(rotB, 1);
}

static void
BM_addmac(benchmark::State& state)
{
	const size_t dim = state.range(0);
	auto mat = triples_matrix_alloc(dim, dim, 1);
	auto enc = triples_matrix_alloc_enc(dim, dim, 1);
	auto mac = triples_matrix_alloc_enc(dim, dim, 1);

	for (auto _ : state) {
		state.PauseTiming();
		triples_matrix_init_rand(mat, 1);
		triples_matrix_encryptA(enc, mat);

		state.ResumeTiming();
		triples_matrix_mac(mac, enc);
	}

	triples_matrix_dealloc(mat, 1);
	triples_matrix_dealloc_enc(enc, 1);
	triples_matrix_dealloc_enc(mac, 1);
}

static void
BM_decrypt(benchmark::State& state)
{
	const size_t dim = state.range(0);
	auto mat = triples_matrix_alloc(dim, dim, 1);
	auto enc = triples_matrix_alloc_enc(dim, dim, 1);

	for (auto _ : state) {
		state.PauseTiming();
		triples_matrix_init_rand(mat, 1);
		triples_matrix_encryptA(enc, mat);

		state.ResumeTiming();
		triples_matrix_decrypt(mat, enc);
	}

	triples_matrix_dealloc(mat, 1);
	triples_matrix_dealloc_enc(enc, 1);
}

BENCHMARK(BM_init)->Unit(benchmark::kSecond);
BENCHMARK(BM_init_hoisted)->Unit(benchmark::kSecond);
BENCHMARK(BM_encryptA)->Unit(benchmark::kSecond)
	->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(8);
BENCHMARK(BM_encryptB)->Unit(benchmark::kSecond)
	->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(8);
BENCHMARK(BM_permuteA)->Unit(benchmark::kSecond)
	->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(8);
BENCHMARK(BM_permuteB)->Unit(benchmark::kSecond)
	->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(8);
BENCHMARK(BM_blockmul)->Unit(benchmark::kSecond)
	->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(8);
BENCHMARK(BM_addmac)->Unit(benchmark::kSecond)
	->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(8);
BENCHMARK(BM_decrypt)->Unit(benchmark::kSecond)
	->Arg(1)->Arg(2)->Arg(3)->Arg(4)->Arg(8);
BENCHMARK_MAIN();
