add_library(bgv STATIC bgv.h bgv.c
	msg_mpz.h msg_mpz.c
	random.h random.c
	util.h util.c)
target_link_libraries(bgv gmp)

if(BGV_USE_HEXL)
	target_sources(bgv PRIVATE ring_hexl.h ring_hexl.cpp)
	target_link_libraries(bgv hexl)
else()
	target_sources(bgv PRIVATE ring_i64.h)
endif()
