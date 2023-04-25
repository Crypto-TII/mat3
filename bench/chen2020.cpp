#define TIIMAT3_HOIST  1
#define TIIMAT3_PREROT 1
#define TIIMAT3_DIM    128

#include <benchmark/benchmark.h>
#include <gmp.h>

extern "C" {
	#include "triples.h"
}


static void
BM_init(benchmark::State& state)
{
	for (auto _ : state) {
		state.ResumeTiming();
		tiimat3_init();

		state.PauseTiming();
		tiimat3_deinit();
	}
}

static void
BM_init_hoisted(benchmark::State& state)
{

	tiimat3_init();
	tiimat3_init_hoisted();

	for (auto _ : state) {
		state.PauseTiming();
		tiimat3_deinit();
		tiimat3_deinit_hoisted();

		state.ResumeTiming();
		tiimat3_init();
		tiimat3_init_hoisted();
	}
}

static void
BM_encryptA(benchmark::State& state)
{
	const size_t dim = state.range(0);
	auto mat = tiimat3_matrix_alloc(dim, dim, 1);
	auto enc = tiimat3_matrix_alloc_enc(dim, dim, 1);

	for (auto _ : state) {
		state.PauseTiming();
		tiimat3_matrix_init_rand(mat, 1);

		state.ResumeTiming();
		tiimat3_matrix_encryptA(enc, mat);
	}

	tiimat3_matrix_dealloc(mat, 1);
	tiimat3_matrix_dealloc_enc(enc, 1);
}

static void
BM_encryptB(benchmark::State& state)
{
	const size_t dim = state.range(0);
	auto mat = tiimat3_matrix_alloc(dim, dim, 1);
	auto enc = tiimat3_matrix_alloc_enc(dim, dim, 1);

	for (auto _ : state) {
		state.PauseTiming();
		tiimat3_matrix_init_rand(mat, 1);

		state.ResumeTiming();
		tiimat3_matrix_encryptB(enc, mat);
	}

	tiimat3_matrix_dealloc(mat, 1);
	tiimat3_matrix_dealloc_enc(enc, 1);
}

static void
BM_permuteA(benchmark::State& state)
{
	const size_t dim = state.range(0);
	auto mat = tiimat3_matrix_alloc(dim, dim, 1);
	auto enc = tiimat3_matrix_alloc_enc(dim, dim, 1);
	auto rot1 = tiimat3_matrix_alloc_encdim(dim, dim, 1);
	auto rot2 = tiimat3_matrix_alloc_encdim(dim, dim, 1);

	for (auto _ : state) {
		state.PauseTiming();
		tiimat3_matrix_init_rand(mat, 1);
		tiimat3_matrix_encryptA(enc, mat);

		state.ResumeTiming();
		for (size_t row = 0; row < dim; ++row) {
			for (size_t col = 0; col < dim; ++col) {
				size_t idx = row * dim + col;
				size_t idxdim = idx * TIIMAT3_DIM;

				#if TIIMAT3_HOIST
				tiimat3_block_prerotA_enc_hoisted(&rot1->value[idxdim], &rot2->value[idxdim], &enc->value[idx]);
				#else
				tiimat3_block_prerotA_enc(&rot1->value[idxdim], &rot2->value[idxdim], &enc->value[idx]);
				#endif /* TIIMAT3_HOIST */
			}
		}
	}

	tiimat3_matrix_dealloc(mat, 1);
	tiimat3_matrix_dealloc_enc(enc, 1);
	tiimat3_matrix_dealloc_enc(rot1, 1);
	tiimat3_matrix_dealloc_enc(rot2, 1);
}

static void
BM_permuteB(benchmark::State& state)
{
	const size_t dim = state.range(0);
	auto mat = tiimat3_matrix_alloc(dim, dim, 1);
	auto enc = tiimat3_matrix_alloc_enc(dim, dim, 1);
	auto rot = tiimat3_matrix_alloc_encdim(dim, dim, 1);

	for (auto _ : state) {
		state.PauseTiming();
		tiimat3_matrix_init_rand(mat, 1);
		tiimat3_matrix_encryptB(enc, mat);

		state.ResumeTiming();
		for (size_t row = 0; row < dim; ++row) {
			for (size_t col = 0; col < dim; ++col) {
				size_t idx = row * dim + col;
				size_t idxdim = idx * TIIMAT3_DIM;

				#if TIIMAT3_HOIST
				tiimat3_block_prerotB_enc_hoisted(&rot->value[idxdim], &enc->value[idx]);
				#else
				tiimat3_block_prerotB_enc(&rot->value[idxdim], &enc->value[idx]);
				#endif /* TIIMAT3_HOIST */
			}
		}
	}

	tiimat3_matrix_dealloc(mat, 1);
	tiimat3_matrix_dealloc_enc(enc, 1);
	tiimat3_matrix_dealloc_enc(rot, 1);
}

static void
BM_blockmul(benchmark::State& state)
{
	const size_t dim = state.range(0);

	auto A = tiimat3_matrix_alloc(dim, dim, 1);
	auto B = tiimat3_matrix_alloc(dim, dim, 1);
	auto encA = tiimat3_matrix_alloc_enc(dim, dim, 1);
	auto encB = tiimat3_matrix_alloc_enc(dim, dim, 1);
	auto encAB = tiimat3_matrix_alloc_enc(dim, dim, 1);
	auto tmp = (tiimat3_BlockEnc *)tiimat3_util_alloc(1, sizeof(tiimat3_BlockEnc));

	#if TIIMAT3_PREROT
	auto rotA1 = tiimat3_matrix_alloc_encdim(dim, dim, 1);
	auto rotA2 = tiimat3_matrix_alloc_encdim(dim, dim, 1);
	auto rotB = tiimat3_matrix_alloc_encdim(dim, dim, 1);
	#else
	auto rotA1 = matrix_alloc_encdim(0, 0, 1);
	auto rotA2 = matrix_alloc_encdim(0, 0, 1);
	auto rotB = matrix_alloc_encdim(0, 0, 1);
	#endif /* TIIMAT3_PREROT */

	for (auto _ : state) {
		state.PauseTiming();
		tiimat3_matrix_init_rand(A, 1);
		tiimat3_matrix_init_rand(B, 1);
		tiimat3_matrix_encryptA(encA, A);
		tiimat3_matrix_encryptB(encB, B);

		#if TIIMAT3_PREROT
		for (size_t row = 0; row < dim; ++row) {
			for (size_t col = 0; col < dim; ++col) {
				size_t idx = row * dim + col;
				size_t idxdim = idx * TIIMAT3_DIM;

				#if TIIMAT3_HOIST
				tiimat3_block_prerotA_enc_hoisted(&rotA1->value[idxdim], &rotA2->value[idxdim], &encA->value[idx]);
				tiimat3_block_prerotB_enc_hoisted(&rotB->value[idxdim], &encB->value[idx]);
				#else
				tiimat3_block_prerotA_enc(&rotA1->value[idxdim], &rotA2->value[idxdim], &encA->value[idx]);
				tiimat3_block_prerotB_enc(&rot->value[idxdim], &enc->value[idx]);
				#endif /* TIIMAT3_HOIST */
			}
		}
		#endif /* TIIMAT3_PREROT */

		state.ResumeTiming();
		for (size_t row = 0; row < dim; ++row) {
			for (size_t col = 0; col < dim; ++col) {
				size_t idx = row * dim + col;

				memset(&encAB->value[idx], 0, sizeof encAB->value[idx]);
				for (size_t k = 0; k < dim; ++k) {
					size_t idxA = row * dim + k;
					size_t idxB = k * dim + col;

					#if TIIMAT3_PREROT
					idxA *= TIIMAT3_DIM;
					idxB *= TIIMAT3_DIM;
					tiimat3_block_mul_enc_prerot(tmp, &rotA1->value[idxA], &rotA2->value[idxA], &rotB->value[idxB]);
					#else
					#if TIIMAT3_HOIST
					tiimat3_block_mul_enc_hoisted(tmp, &encA->value[idxA], &encB->value[idxB]);
					#else
					tiimat3_block_mul_enc_reuse(tmp, &encA->value[idxA], &encB->value[idxB]);
					#endif /* TIIMAT3_HOIST */
					#endif /* TIIMAT3_PREROT */

					tiimat3_block_add_enc(&encAB->value[idx], &encAB->value[idx], tmp);
				}
			}
		}
	}

	tiimat3_util_dealloc(tmp);
	tiimat3_matrix_dealloc(A, 1);
	tiimat3_matrix_dealloc(B, 1);
	tiimat3_matrix_dealloc_enc(encA, 1);
	tiimat3_matrix_dealloc_enc(encB, 1);
	tiimat3_matrix_dealloc_enc(encAB, 1);
	tiimat3_matrix_dealloc_enc(rotA1, 1);
	tiimat3_matrix_dealloc_enc(rotA2, 1);
	tiimat3_matrix_dealloc_enc(rotB, 1);
}

static void
BM_addmac(benchmark::State& state)
{
	const size_t dim = state.range(0);
	auto mat = tiimat3_matrix_alloc(dim, dim, 1);
	auto enc = tiimat3_matrix_alloc_enc(dim, dim, 1);
	auto mac = tiimat3_matrix_alloc_enc(dim, dim, 1);

	for (auto _ : state) {
		state.PauseTiming();
		tiimat3_matrix_init_rand(mat, 1);
		tiimat3_matrix_encryptA(enc, mat);

		state.ResumeTiming();
		tiimat3_matrix_mac(mac, enc);
	}

	tiimat3_matrix_dealloc(mat, 1);
	tiimat3_matrix_dealloc_enc(enc, 1);
	tiimat3_matrix_dealloc_enc(mac, 1);
}

static void
BM_decrypt(benchmark::State& state)
{
	const size_t dim = state.range(0);
	auto mat = tiimat3_matrix_alloc(dim, dim, 1);
	auto enc = tiimat3_matrix_alloc_enc(dim, dim, 1);

	for (auto _ : state) {
		state.PauseTiming();
		tiimat3_matrix_init_rand(mat, 1);
		tiimat3_matrix_encryptA(enc, mat);

		state.ResumeTiming();
		tiimat3_matrix_decrypt(mat, enc);
	}

	tiimat3_matrix_dealloc(mat, 1);
	tiimat3_matrix_dealloc_enc(enc, 1);
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
