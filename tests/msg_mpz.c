#include "../src/msg_mpz.c"

void
test_cmod(void)
{
	mpz_t a, mod;

	mpz_inits(a, mod, 0);

	mpz_set_si(mod, 60);
	mpz_set_si(a, -59);
	mpz_cmod(a, a, mod);
	assert(mpz_cmp_si(a, 1) == 0);

	mpz_set_si(a, 29);
	mpz_cmod(a, a, mod);
	assert(mpz_cmp_si(a, 29) == 0);

	mpz_set_si(a, 30);
	mpz_cmod(a, a, mod);
	assert(mpz_cmp_si(a, -30) == 0);

	mpz_set_si(a, 59);
	mpz_cmod(a, a, mod);
	assert(mpz_cmp_si(a, -1) == 0);

	mpz_clears(a, mod, 0);
}

void
test_crt(void)
{
	mpz_t poly[1024];
	const uint64_t *crns[4];
	uint64_t *rns[4], mods[4] = { 3, 5, 7, 11 };
	size_t d, i;


	d = 1024;
	for (i = 0; i < 4; ++i)
		rns[i] = calloc(d, sizeof *rns);
	for (i = 0; i < d; ++i)
		mpz_init(poly[i]);

	for (i = 0; i < d; ++i) {
		rns[0][i] = i % mods[0];
		rns[1][i] = i % mods[1];
		rns[2][i] = i % mods[2];
		rns[3][i] = i % mods[3];
	}

	crns[0] = rns[0];
	crns[1] = rns[1];
	crns[2] = rns[2];
	crns[3] = rns[3];
	mpz_crt(poly, d, crns, mods, 4);

	for (i = 0; i < d; ++i)
		assert(mpz_cmp_ui(poly[i], i) == 0);

	for (i = 0; i < 4; ++i)
		free(rns[i]);
	for (i = 0; i < d; ++i)
		mpz_clear(poly[i]);
}

void
test_seiler(void)
{
	mpz_t poly[32768];
	mpz_t mod, root;
	size_t d, i;

	mpz_inits(mod, root, 0);

	d = 32768;
	mpz_set_ui(mod, 0x1230001);
	for (i = 0; i < d; ++i)
		mpz_init_set_ui(poly[i], i);
	mpz_precomp(root, d, mod);
	mpz_ntt(poly, d, mod, root);
	mpz_intt(poly, d, mod, root);
	for (i = 0; i < d; ++i)
		assert(mpz_cmp_ui(poly[i], i) == 0);

	for (i = 0; i < d; ++i)
		mpz_clear(poly[i]);
	mpz_clears(mod, root, 0);
}

int
main(void)
{

	fputs("[+] Testing cmod (mpz):       ", stderr);
	test_cmod(), fputs("4/4.\n", stderr);

	fputs("[+] Testing CRT (mpz):        ", stderr);
	test_crt(), fputs("1/1.\n", stderr);

	fputs("[+] Testing Seiler (mpz):     ", stderr);
	test_seiler(), fputs("1/1.\n", stderr);

	return 0;
}
