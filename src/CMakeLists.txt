add_library(libbgv STATIC libbgv.h libbgv.c
	bgv.h
	msg_mpz.h
	random.h
	util.h util.c)
target_link_libraries(libbgv gmp)

if(BGV_USE_HEXL)
	target_sources(libbgv PRIVATE ring_hexl.h ring_hexl.cpp)
	target_link_libraries(libbgv hexl)
else()
	target_sources(libbgv PRIVATE ring_i64.h)
endif()
