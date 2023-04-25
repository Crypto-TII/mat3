set(SRC ${PROJECT_SOURCE_DIR}/src)
include_directories(${SRC})

add_executable(test_random random.c)
target_compile_options(test_random PRIVATE -fsanitize=address)
target_link_options(test_random PRIVATE -fsanitize=address)

add_executable(test_ring_i64 ring_i64.c
	${SRC}/ring_i64.h
	${SRC}/random.h ${SRC}/random.c
	${SRC}/util.h ${SRC}/util.c)
target_compile_options(test_ring_i64 PRIVATE -fsanitize=address)
target_link_libraries(test_ring_i64 gmp)
target_link_options(test_ring_i64 PRIVATE -fsanitize=address)

add_executable(test_msg_mpz msg_mpz.c
	${SRC}/msg_mpz.h
	${SRC}/util.h ${SRC}/util.c)
target_compile_options(test_msg_mpz PRIVATE -fsanitize=address)
target_link_libraries(test_msg_mpz gmp)
target_link_options(test_msg_mpz PRIVATE -fsanitize=address)

add_executable(test_bgv bgv.c
	${SRC}/bgv.h
	${SRC}/msg_mpz.h
	${SRC}/random.h ${SRC}/random.c
	${SRC}/util.h ${SRC}/util.c)
target_compile_options(test_bgv PRIVATE -fsanitize=address)
target_link_libraries(test_bgv gmp)
target_link_options(test_bgv PRIVATE -fsanitize=address)

if(BGV_USE_HEXL)
	target_sources(test_bgv PRIVATE ${SRC}/ring_hexl.h ${SRC}/ring_hexl.cpp)
	target_link_libraries(test_bgv hexl)
else()
	target_sources(test_bgv PRIVATE ${SRC}/ring_i64.h)
endif()

foreach(HOIST IN ITEMS 0 1)
	foreach(PREROT IN ITEMS 0 1)
		add_executable(test_triples_${HOIST}${PREROT} triples.c
			${SRC}/triples.h)
		target_compile_definitions(test_triples_${HOIST}${PREROT} PRIVATE
			-DTRIPLES_DIM=4
			-DTRIPLES_HOIST=${HOIST}
			-DTRIPLES_PREROT=${PREROT})
		target_compile_options(test_triples_${HOIST}${PREROT} PRIVATE
			-fopenmp -fsanitize=address -Wno-type-limits)
		target_link_options(test_triples_${HOIST}${PREROT} PRIVATE -fopenmp -fsanitize=address)
		target_link_libraries(test_triples_${HOIST}${PREROT} libbgv)
	endforeach()
endforeach()
