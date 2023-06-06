#include "../src/random.c"

#define SAMPLES 10000000
int8_t samples_i8[SAMPLES];
uint64_t samples_u64[SAMPLES];

static double
binomial(int n, int k)
{
	long bin;
	int i;

	if (k > n - k)
		k = n - k;

	bin = 1;
	for (i = 1; i <= k; ++i, --n) {
		if (n % i == 0)
			bin = n / i * bin;
		else if (bin % i == 0)
			bin = bin / i * n;
		else
			bin = bin * n / i;
	}

	return bin;
}

static double
chi2(double *exp, double *obs, size_t len)
{
	double chi;
	size_t i;

	chi = 0;
	for (i = 0; i < len; ++i)
		if (exp[i] > 0)
			chi += (obs[i] - exp[i]) / exp[i] * (obs[i] - exp[i]);

	return chi;
}

static void
test_keccak(void)
{
	const char *in = "fhe";
	const char *exp = "\xe3\xe5\x0b\x8b\x71\xc7\xd3\x9b\xbf\x2b\xf2\xea"
	                  "\x6c\x9a\x1f\x3b\x73\x83\x8a\xcb\xeb\x73\xd7\x22"
	                  "\x31\x44\x09\xfa\x5d\xac\x1e\x63\xdd\xec\xd8\xa2"
	                  "\xde\x22\xf4\xb5\x80\x4f\x56\x70\x45\xf8\xc5\x1e";

	uint64_t state[KECCAK_STATE_LEN];
	unsigned char hash[SHA3_384_RATE];

	sha3_384_absorb(state, (const unsigned char *)in, strlen(in));
	sha3_384_squeeze(state, hash);

	assert(memcmp(hash, exp, SHA3_384_BYTES) == 0);
}

static void
test_random(void)
{
	/* https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm */
	tiimat3_Seed seed;
	double obs[256], exp[256];
	signed char r;
	size_t sum, i;

	tiimat3_random_seed(&seed);
	random_bytes(samples_i8, SAMPLES, &seed);

	memset(obs, 0, sizeof obs);
	for (i = 0; i < SAMPLES; ++i)
		++obs[128 + samples_i8[i]];

	sum = 0;
	for (i = 0; i < 256; ++i)
		sum += obs[i];
	assert(sum == SAMPLES);

	for (i = 0; i < 256; ++i)
		exp[i] = SAMPLES / 256.0;
	assert(chi2(exp, obs, 8) < 26.125);

	seed.init = 1;
	for (i = 0; i < SAMPLES; ++i) {
		random_bytes(&r, 1, &seed);
		assert(r == samples_i8[i]);
	}
}

static void
test_sample(void)
{
	/* https://www.itl.nist.gov/div898/handbook/eda/section3/eda3674.htm */
	tiimat3_Seed seed;
	double obs[101], exp[101];
	size_t sum, i;

	tiimat3_random_seed(&seed);

	tiimat3_random_cbd1_i8(samples_i8, SAMPLES, &seed);
	memset(obs, 0, sizeof obs);
	for (i = 0; i < SAMPLES; ++i)
		++obs[1 + samples_i8[i]];
	sum = 0;
	for (i = 0; i <= 2; ++i)
		sum += obs[i];
	assert(sum == SAMPLES);
	for (i = 0; i <= 2; ++i)
		exp[i] = binomial(2, i) * SAMPLES / (1 << 2);
	assert(chi2(obs, exp, 3) < 13.816);

	tiimat3_random_cbd21_i8(samples_i8, SAMPLES, &seed);
	memset(obs, 0, sizeof obs);
	for (i = 0; i < SAMPLES; ++i)
		++obs[21 + samples_i8[i]];
	sum = 0;
	for (i = 0; i <= 42; ++i)
		sum += obs[i];
	assert(sum == SAMPLES);
	for (i = 0; i <= 42; ++i)
		exp[i] = (double)binomial(42, i) * SAMPLES / (1l << 42);
	assert(chi2(obs, exp, 43) < 76.084);

	tiimat3_random_uniform_u64(samples_u64, SAMPLES, 101, &seed);
	memset(obs, 0, sizeof obs);
	for (i = 0; i < SAMPLES; ++i)
		++obs[samples_u64[i]];
	sum = 0;
	for (i = 0; i < 101; ++i)
		sum += obs[i];
	assert(sum == SAMPLES);
	for (i = 0; i < 101; ++i)
		exp[i] = SAMPLES / 101.0;
	assert(chi2(obs, exp, 101) < 149.449);
}

int
main(void)
{

	fputs("[+] Testing Keccak:           ", stderr);
	test_keccak(), fputs("1/1.\n", stderr);

	fputs("[+] Testing random:           ", stderr);
	test_random(), fputs("3/3.\n", stderr);

	fputs("[+] Testing sample:           ", stderr);
	test_sample(), fputs("6/6.\n", stderr);

	return 0;
}
