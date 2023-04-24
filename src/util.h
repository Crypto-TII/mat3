#ifndef BGV_UTIL_H
#define BGV_UTIL_H

#include <stdint.h>
#include <stdlib.h>

void *bgv_util_alloc(size_t len, size_t size);
void  bgv_util_dealloc(void *ptr);

/* https://en.wikipedia.org/wiki/Hacker%27s_Delight */
uint64_t bgv_util_bitrev(uint64_t a, size_t log2);
uint64_t bgv_util_invmod(uint64_t a, uint64_t mod);

#endif /* BGV_UTIL_H */


#ifdef BGV_UTIL_IMPL

void *
bgv_util_alloc(size_t len, size_t size)
{
	void *p;

	p = calloc(len, size);
	if (p == 0) {
		fputs("bgv: calloc\n", stderr);
		exit(1);
	}

	return p;
}

void
bgv_util_dealloc(void *ptr)
{

	free(ptr);
}

uint64_t
bgv_util_bitrev(uint64_t a, size_t log2)
{

	a = ((a & 0x5555555555555555) <<  1) | ((a & 0xaaaaaaaaaaaaaaaa) >>  1);
	a = ((a & 0x3333333333333333) <<  2) | ((a & 0xcccccccccccccccc) >>  2);
	a = ((a & 0x0f0f0f0f0f0f0f0f) <<  4) | ((a & 0xf0f0f0f0f0f0f0f0) >>  4);
	a = ((a & 0x00ff00ff00ff00ff) <<  8) | ((a & 0xff00ff00ff00ff00) >>  8);
	a = ((a & 0x0000ffff0000ffff) << 16) | ((a & 0xffff0000ffff0000) >> 16);

	return ((a << 32) | (a >> 32)) >> (64 - log2);
}

uint64_t
bgv_util_invmod(uint64_t a, uint64_t mod)
{
	uint64_t q, r, r0, r1, t, t0, t1;

	r0 = mod, r1 = a;
	t0 = 0, t1 = 1;
	do {
		q = r0 / r1;
		r = r0 % r1;
		t = t0 - q * t1;

		r0 = r1, r1 = r;
		t0 = t1, t1 = t;
	} while (r != 0);

	return t0;
}

#endif /* BGV_UTIL_IMPL */
