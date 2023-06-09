set(SRC ${PROJECT_SOURCE_DIR}/src)
include_directories(${SRC})

add_executable(test_random random.c
	${SRC}/util.h ${SRC}/util.c)
target_compile_options(test_random PRIVATE -fsanitize=address)
target_link_libraries(test_random crypto ssl)
target_link_options(test_random PRIVATE -fsanitize=address)

add_executable(test_ring_i64 ring_i64.c
	${SRC}/util.h ${SRC}/util.c
	${SRC}/random.h ${SRC}/random.c)
target_compile_options(test_ring_i64 PRIVATE -fsanitize=address)
target_link_libraries(test_ring_i64 crypto gmp ssl)
target_link_options(test_ring_i64 PRIVATE -fsanitize=address)

add_executable(test_msg_mpz msg_mpz.c
	${SRC}/util.h ${SRC}/util.c)
target_compile_options(test_msg_mpz PRIVATE -fsanitize=address)
target_link_libraries(test_msg_mpz crypto gmp ssl)
target_link_options(test_msg_mpz PRIVATE -fsanitize=address)

add_executable(test_bgv bgv.c
	${SRC}/util.h ${SRC}/util.c
	${SRC}/random.h ${SRC}/random.c
	${SRC}/msg_mpz.h ${SRC}/msg_mpz.c)
target_compile_options(test_bgv PRIVATE -fsanitize=address)
target_link_libraries(test_bgv crypto gmp ssl)
target_link_options(test_bgv PRIVATE -fsanitize=address)

if(TIIMAT3_USE_HEXL)
	target_sources(test_bgv PRIVATE
		${SRC}/ring_hexl.h ${SRC}/ring_hexl.cpp)
	target_link_libraries(test_bgv hexl)
else()
	target_sources(test_bgv PRIVATE
		${SRC}/ring_i64.h ${SRC}/ring_i64.c)
endif()

foreach(HOIST IN ITEMS 0 1)
	foreach(PREROT IN ITEMS 0 1)
		add_executable(test_triples_${HOIST}${PREROT} triples.c
			${SRC}/util.h ${SRC}/util.c
			${SRC}/random.h ${SRC}/random.c
			${SRC}/msg_mpz.h ${SRC}/msg_mpz.c
			${SRC}/bgv.h ${SRC}/bgv.c)
		if(TIIMAT3_USE_HEXL)
			target_sources(test_triples_${HOIST}${PREROT} PRIVATE
				${SRC}/ring_hexl.h ${SRC}/ring_hexl.cpp)
			target_link_libraries(test_triples_${HOIST}${PREROT}
				crypto gmp hexl ssl)
		else()
			target_sources(test_triples_${HOIST}${PREROT} PRIVATE
				${SRC}/ring_i64.h ${SRC}/ring_i64.c)
			target_link_libraries(test_triples_${HOIST}${PREROT}
				crypto gmp ssl)
		endif()

		target_compile_definitions(test_triples_${HOIST}${PREROT} PRIVATE
			-DTIIMAT3_DIM=4
			-DTIIMAT3_HOIST=${HOIST}
			-DTIIMAT3_PREROT=${PREROT})
		target_compile_options(test_triples_${HOIST}${PREROT} PRIVATE
			-fopenmp -fsanitize=address -Wno-type-limits)
		target_link_options(test_triples_${HOIST}${PREROT} PRIVATE -fopenmp -fsanitize=address)
	endforeach()
endforeach()
