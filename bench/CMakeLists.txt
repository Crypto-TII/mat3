set(SRC ${PROJECT_SOURCE_DIR}/src)
include_directories(${SRC})

add_executable(bench_chen2020 chen2020.cpp
	${SRC}/util.h ${SRC}/util.c
	${SRC}/random.h ${SRC}/random.c
	${SRC}/msg_mpz.h ${SRC}/msg_mpz.c
	${SRC}/bgv.h ${SRC}/bgv.c
	${SRC}/triples.h ${SRC}/triples.c)
if(TIIMAT3_USE_HEXL)
	target_sources(bench_chen2020 PRIVATE
		${SRC}/ring_hexl.h ${SRC}/ring_hexl.cpp)
	target_link_libraries(bench_chen2020 hexl)
else()
	target_sources(bench_chen2020 PRIVATE
		${SRC}/ring_i64.h ${SRC}/ring_i64.cpp)
endif()
target_compile_options(bench_chen2020 PRIVATE -fopenmp)
target_link_options(bench_chen2020 PRIVATE -fopenmp)
target_link_libraries(bench_chen2020 benchmark crypto gmp ssl)

foreach(DIM IN ITEMS 32 64 128)
	add_executable(bench_slicing_${DIM} slicing.cpp
		${SRC}/util.h ${SRC}/util.c
		${SRC}/random.h ${SRC}/random.c
		${SRC}/msg_mpz.h ${SRC}/msg_mpz.c
		${SRC}/bgv.h ${SRC}/bgv.c
		${SRC}/triples.h ${SRC}/triples.c)
	if(TIIMAT3_USE_HEXL)
		target_sources(bench_slicing_${DIM} PRIVATE
			${SRC}/ring_hexl.h ${SRC}/ring_hexl.cpp)
		target_link_libraries(bench_slicing_${DIM} hexl)
	else()
		target_sources(bench_slicing_${DIM} PRIVATE
			${SRC}/ring_i64.h ${SRC}/ring_i64.cpp)
	endif()
	target_compile_definitions(bench_slicing_${DIM} PRIVATE
		-DTIIMAT3_DIM=${DIM})
	target_compile_options(bench_slicing_${DIM} PRIVATE -fopenmp)
	target_link_options(bench_slicing_${DIM} PRIVATE -fopenmp)
	target_link_libraries(bench_slicing_${DIM} benchmark crypto gmp ssl)
endforeach()

foreach(BLOCKS IN ITEMS 1 2 3 4 8)
	add_executable(bench_triples_hoist_${BLOCKS} triples.cpp
		${SRC}/util.h ${SRC}/util.c
		${SRC}/random.h ${SRC}/random.c
		${SRC}/msg_mpz.h ${SRC}/msg_mpz.c
		${SRC}/bgv.h ${SRC}/bgv.c
		${SRC}/triples.h ${SRC}/triples.c)
	if(TIIMAT3_USE_HEXL)
		target_sources(bench_triples_hoist_${BLOCKS} PRIVATE
			${SRC}/ring_hexl.h ${SRC}/ring_hexl.cpp)
		target_link_libraries(bench_triples_hoist_${BLOCKS} hexl)
	else()
		target_sources(bench_triples_hoist_${BLOCKS} PRIVATE
			${SRC}/ring_i64.h ${SRC}/ring_i64.cpp)
	endif()
	target_compile_definitions(bench_triples_hoist_${BLOCKS} PRIVATE
		-DTIIMAT3_BLOCKS=${BLOCKS}
		-DTIIMAT3_HOIST=1
		-DTIIMAT3_PREROT=1)
	target_compile_options(bench_triples_hoist_${BLOCKS} PRIVATE -fopenmp)
	target_link_options(bench_triples_hoist_${BLOCKS} PRIVATE -fopenmp)
	target_link_libraries(bench_triples_hoist_${BLOCKS} benchmark crypto gmp ssl)

	add_executable(bench_triples_mixed_${BLOCKS} triples.cpp
		${SRC}/util.h ${SRC}/util.c
		${SRC}/random.h ${SRC}/random.c
		${SRC}/msg_mpz.h ${SRC}/msg_mpz.c
		${SRC}/bgv.h ${SRC}/bgv.c
		${SRC}/triples.h ${SRC}/triples.c)
	if(TIIMAT3_USE_HEXL)
		target_sources(bench_triples_mixed_${BLOCKS} PRIVATE
			${SRC}/ring_hexl.h ${SRC}/ring_hexl.cpp)
		target_link_libraries(bench_triples_mixed_${BLOCKS} hexl)
	else()
		target_sources(bench_triples_mixed_${BLOCKS} PRIVATE
			${SRC}/ring_i64.h ${SRC}/ring_i64.cpp)
	endif()
	target_compile_definitions(bench_triples_mixed_${BLOCKS} PRIVATE
		-DTIIMAT3_BLOCKS=${BLOCKS}
		-DTIIMAT3_HOIST=0
		-DTIIMAT3_MIXED=1
		-DTIIMAT3_PREROT=1)
	target_compile_options(bench_triples_mixed_${BLOCKS} PRIVATE -fopenmp)
	target_link_options(bench_triples_mixed_${BLOCKS} PRIVATE -fopenmp)
	target_link_libraries(bench_triples_mixed_${BLOCKS} benchmark crypto gmp ssl)

	add_executable(bench_triples_recomp_${BLOCKS} triples.cpp
		${SRC}/util.h ${SRC}/util.c
		${SRC}/random.h ${SRC}/random.c
		${SRC}/msg_mpz.h ${SRC}/msg_mpz.c
		${SRC}/bgv.h ${SRC}/bgv.c
		${SRC}/triples.h ${SRC}/triples.c)
	if(TIIMAT3_USE_HEXL)
		target_sources(bench_triples_recomp_${BLOCKS} PRIVATE
			${SRC}/ring_hexl.h ${SRC}/ring_hexl.cpp)
		target_link_libraries(bench_triples_recomp_${BLOCKS} hexl)
	else()
		target_sources(bench_triples_recomp_${BLOCKS} PRIVATE
			${SRC}/ring_i64.h ${SRC}/ring_i64.cpp)
	endif()
	target_compile_definitions(bench_triples_recomp_${BLOCKS} PRIVATE
		-DTIIMAT3_BLOCKS=${BLOCKS}
		-DTIIMAT3_HOIST=1
		-DTIIMAT3_PREROT=0)
	target_compile_options(bench_triples_recomp_${BLOCKS} PRIVATE -fopenmp)
	target_link_options(bench_triples_recomp_${BLOCKS} PRIVATE -fopenmp)
	target_link_libraries(bench_triples_recomp_${BLOCKS} benchmark crypto gmp ssl)

	add_executable(bench_triples_reuse_${BLOCKS} triples.cpp
		${SRC}/util.h ${SRC}/util.c
		${SRC}/random.h ${SRC}/random.c
		${SRC}/msg_mpz.h ${SRC}/msg_mpz.c
		${SRC}/bgv.h ${SRC}/bgv.c
		${SRC}/triples.h ${SRC}/triples.c)
	if(TIIMAT3_USE_HEXL)
		target_sources(bench_triples_reuse_${BLOCKS} PRIVATE
			${SRC}/ring_hexl.h ${SRC}/ring_hexl.cpp)
		target_link_libraries(bench_triples_reuse_${BLOCKS} hexl)
	else()
		target_sources(bench_triples_reuse_${BLOCKS} PRIVATE
			${SRC}/ring_i64.h ${SRC}/ring_i64.cpp)
	endif()
	target_compile_definitions(bench_triples_reuse_${BLOCKS} PRIVATE
		-DTIIMAT3_BLOCKS=${BLOCKS}
		-DTIIMAT3_HOIST=0
		-DTIIMAT3_PREROT=1)
	target_compile_options(bench_triples_reuse_${BLOCKS} PRIVATE -fopenmp)
	target_link_options(bench_triples_reuse_${BLOCKS} PRIVATE -fopenmp)
	target_link_libraries(bench_triples_reuse_${BLOCKS} benchmark crypto gmp ssl)
endforeach()
