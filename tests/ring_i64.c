#include "../src/ring_i64.c"

void
test_barrett(void)
{
	int64_t mod, barrett, i;

	mod = 0xfff1;
	barrett = barrett_precomp(mod);

	for (i = -mod; i < 0; ++i)
		assert(barrett_reduce(i, mod, barrett) == i + mod);
	for (i = 0; i < mod; ++i)
		assert(barrett_reduce(i, mod, barrett) == i);
	for (i = 0; i < mod; ++i)
		assert(barrett_reduce(0x1234 * mod + i, mod, barrett) == i);
	for (i = 0; i < mod; ++i)
		assert(barrett_reduce(-0xabcd * mod + i, mod, barrett) == i);
}

void
test_montgomery(void)
{
	int64_t mod, inv, pow, tmp, a, b, c, i;

	mod = 0xfff1;
	inv = montgomery_precomp_inv(mod);
	pow = montgomery_precomp_pow(mod);

	assert((uint64_t)mod * inv == 1);

	for (i = 0; i < mod; ++i) {
		tmp = montgomery_init(i, mod, inv, pow);
		assert(montgomery_deinit(tmp, mod, inv) == i);
	}

	a = 0x1234;
	b = -0xabcd;
	tmp = (a * b) % mod + mod;

	a = montgomery_init(a, mod, inv, pow);
	b = montgomery_init(b, mod, inv, pow);
	c = montgomery_mulred(a, b, mod, inv);
	assert(montgomery_deinit(c, mod, inv) == tmp);

	c = montgomery_powred(a, 0, mod, inv, pow);
	assert(montgomery_deinit(c, mod, inv) == 1);

	c = montgomery_powred(a, 1, mod, inv, pow);
	assert(c == a);

	c = montgomery_powred(a, mod - 1, mod, inv, pow);
	assert(montgomery_deinit(c, mod, inv) == 1);
}

void
test_seiler(void)
{
	int64_t poly1[32768], poly2[32768] = { 0 };
	int64_t mod, root, barrett, pow, inv;
	size_t d, i;

	d = 32768;
	mod = 0x1230001;
	barrett = barrett_precomp(mod);
	inv = montgomery_precomp_inv(mod);
	pow = montgomery_precomp_pow(mod);
	root = seiler_precomp(d, mod, barrett, inv, pow);

	for (i = 0; i < d; ++i)
		poly1[i] = montgomery_init(i, mod, inv, pow);
	seiler_ntt(poly1, d, mod, root, barrett, inv, pow);
	seiler_intt(poly1, d, mod, root, barrett, inv, pow);
	for (i = 0; i < d; ++i)
		assert(montgomery_deinit(poly1[i], mod, inv) == (int64_t)i);

	for (i = 0; i < d; ++i)
		poly1[i] = montgomery_init(1, mod, inv, pow);
	poly2[1] = montgomery_init(2, mod, inv, pow);
	seiler_ntt(poly1, d, mod, root, barrett, inv, pow);
	seiler_ntt(poly2, d, mod, root, barrett, inv, pow);
	for (i = 0; i < d; ++i)
		poly1[i] = montgomery_mulred(poly1[i], poly2[i], mod, inv);
	seiler_intt(poly1, d, mod, root, barrett, inv, pow);
	assert(montgomery_deinit(poly1[0], mod, inv) == mod - 2);
	for (i = 1; i < d; ++i)
		assert(montgomery_deinit(poly1[i], mod, inv) == 2);
}

int
main(void)
{

	fputs("[+] Testing Barrett (i64):    ", stderr);
	test_barrett(), fputs("4/4.\n", stderr);

	fputs("[+] Testing Montgomery (i64): ", stderr);
	test_montgomery(), fputs("6/6.\n", stderr);

	fputs("[+] Testing Seiler (i64):     ", stderr);
	test_seiler(), fputs("2/2.\n", stderr);

	return 0;
}
