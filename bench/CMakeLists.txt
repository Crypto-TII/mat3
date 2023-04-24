set(SRC ${PROJECT_SOURCE_DIR}/src)
include_directories(${SRC})

add_executable(bench_chen2020 chen2020.cpp
	${SRC}/triples.h)
target_compile_options(bench_chen2020 PRIVATE -fopenmp -fpermissive)
target_link_options(bench_chen2020 PRIVATE -fopenmp)
target_link_libraries(bench_chen2020 libbgv benchmark)

foreach(DIM IN ITEMS 32 64 128)
	add_executable(bench_slicing_${DIM} slicing.cpp
		${SRC}/triples.h)
	target_compile_definitions(bench_slicing_${DIM} PRIVATE
		-DTRIPLES_DIM=${DIM})
	target_compile_options(bench_slicing_${DIM} PRIVATE -fopenmp -fpermissive)
	target_link_options(bench_slicing_${DIM} PRIVATE -fopenmp)
	target_link_libraries(bench_slicing_${DIM} libbgv benchmark)
endforeach()

foreach(BLOCKS IN ITEMS 1 2 3 4 8)
	add_executable(bench_triples_hoist_${BLOCKS} triples.cpp
		${SRC}/triples.h)
	target_compile_definitions(bench_triples_hoist_${BLOCKS} PRIVATE
		-DTRIPLES_BLOCKS=${BLOCKS}
		-DTRIPLES_HOIST=1
		-DTRIPLES_PREROT=1)
	target_compile_options(bench_triples_hoist_${BLOCKS} PRIVATE -fopenmp -fpermissive)
	target_link_options(bench_triples_hoist_${BLOCKS} PRIVATE -fopenmp)
	target_link_libraries(bench_triples_hoist_${BLOCKS} libbgv benchmark)

	add_executable(bench_triples_mixed_${BLOCKS} triples.cpp
		${SRC}/triples.h)
	target_compile_definitions(bench_triples_mixed_${BLOCKS} PRIVATE
		-DTRIPLES_BLOCKS=${BLOCKS}
		-DTRIPLES_HOIST=0
		-DTRIPLES_MIXED=1
		-DTRIPLES_PREROT=1)
	target_compile_options(bench_triples_mixed_${BLOCKS} PRIVATE -fopenmp -fpermissive)
	target_link_options(bench_triples_mixed_${BLOCKS} PRIVATE -fopenmp)
	target_link_libraries(bench_triples_mixed_${BLOCKS} libbgv benchmark)

	add_executable(bench_triples_recomp_${BLOCKS} triples.cpp
		${SRC}/triples.h)
	target_compile_definitions(bench_triples_recomp_${BLOCKS} PRIVATE
		-DTRIPLES_BLOCKS=${BLOCKS}
		-DTRIPLES_HOIST=1
		-DTRIPLES_PREROT=0)
	target_compile_options(bench_triples_recomp_${BLOCKS} PRIVATE -fopenmp -fpermissive)
	target_link_options(bench_triples_recomp_${BLOCKS} PRIVATE -fopenmp)
	target_link_libraries(bench_triples_recomp_${BLOCKS} libbgv benchmark)

	add_executable(bench_triples_reuse_${BLOCKS} triples.cpp
		${SRC}/triples.h)
	target_compile_definitions(bench_triples_reuse_${BLOCKS} PRIVATE
		-DTRIPLES_BLOCKS=${BLOCKS}
		-DTRIPLES_HOIST=0
		-DTRIPLES_PREROT=1)
	target_compile_options(bench_triples_reuse_${BLOCKS} PRIVATE -fopenmp -fpermissive)
	target_link_options(bench_triples_reuse_${BLOCKS} PRIVATE -fopenmp)
	target_link_libraries(bench_triples_reuse_${BLOCKS} libbgv benchmark)
endforeach()
