#include "util.h"

void *
tiimat3_util_alloc(size_t len, size_t size)
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
tiimat3_util_dealloc(void *ptr)
{

	free(ptr);
}

uint64_t
tiimat3_util_bitrev(uint64_t a, size_t log2)
{

	a = ((a & 0x5555555555555555) <<  1) | ((a & 0xaaaaaaaaaaaaaaaa) >>  1);
	a = ((a & 0x3333333333333333) <<  2) | ((a & 0xcccccccccccccccc) >>  2);
	a = ((a & 0x0f0f0f0f0f0f0f0f) <<  4) | ((a & 0xf0f0f0f0f0f0f0f0) >>  4);
	a = ((a & 0x00ff00ff00ff00ff) <<  8) | ((a & 0xff00ff00ff00ff00) >>  8);
	a = ((a & 0x0000ffff0000ffff) << 16) | ((a & 0xffff0000ffff0000) >> 16);

	return ((a << 32) | (a >> 32)) >> (64 - log2);
}

uint64_t
tiimat3_util_invmod(uint64_t a, uint64_t mod)
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

size_t
tiimat3_util_msb64(uint64_t a)
{
	size_t i;

	i = 0;
	while (a >>= 1)
		++i;

	return i;
}

uint32_t
tiimat3_util_popcount32(uint32_t a)
{

	a = (a & 0x55555555) + ((a >>  1) & 0x55555555);
	a = (a & 0x33333333) + ((a >>  2) & 0x33333333);
	a = (a & 0x0f0f0f0f) + ((a >>  4) & 0x0f0f0f0f);
	a = (a & 0x00ff00ff) + ((a >>  8) & 0x00ff00ff);
	a = (a & 0x0000ffff) + ((a >> 16) & 0x0000ffff);

	return a;
}

uint64_t
tiimat3_util_rol64(uint64_t a, size_t r)
{

	return (a << r) | (a >> (64 - r));
}
