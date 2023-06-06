#include "../src/bgv.c"

#define LEN(a) (sizeof a / sizeof *a)

void
msg_add(tiimat3_Message *rop, const tiimat3_Message *op1, const tiimat3_Message *op2)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j) {
		mpz_add(rop->value[j], op1->value[j], op2->value[j]);
		mpz_mod(rop->value[j], rop->value[j], tiimat3_t.value);
	}
}

tiimat3_Message *
msg_alloc_rand(void)
{
	gmp_randstate_t state;
	tiimat3_Message *m;
	size_t j;

	gmp_randinit_default(state);

	m = tiimat3_msg_alloc(1);
	for (j = 0; j < TIIMAT3_D; ++j)
		mpz_urandomm(m->value[j], state, tiimat3_t.value);

	gmp_randclear(state);

	return m;
}

int
msg_equal(const tiimat3_Message *m1, const tiimat3_Message *m2)
{
	size_t j;
	int cmp;

	cmp = 0;
	for (j = 0; j < TIIMAT3_D; ++j)
		cmp += mpz_cmp(m1->value[j], m2->value[j]);

	return cmp == 0;
}

void
msg_mul(tiimat3_Message *rop, const tiimat3_Message *op1, const tiimat3_Message *op2)
{
	size_t j;

	for (j = 0; j < TIIMAT3_D; ++j) {
		mpz_mul(rop->value[j], op1->value[j], op2->value[j]);
		mpz_mod(rop->value[j], rop->value[j], tiimat3_t.value);
	}
}

void
msg_rot(tiimat3_Message *rop, const tiimat3_Message *m, size_t steps)
{
	mpz_t tmp;
	size_t i, j;

	mpz_init(tmp);

	for (i = 0; i < steps; ++i) {
		mpz_set(tmp, m->value[0]);
		for (j = 0; j < TIIMAT3_D / 2 - 1; ++j)
			mpz_set(rop->value[j], m->value[j + 1]);
		mpz_set(rop->value[j], tmp);

		mpz_set(tmp, m->value[TIIMAT3_D / 2]);
		for (j = TIIMAT3_D / 2; j < TIIMAT3_D - 1; ++j)
			mpz_set(rop->value[j], m->value[j + 1]);
		mpz_set(rop->value[j], tmp);
	}

	mpz_clear(tmp);
}

void
test_pack(void)
{
	tiimat3_Message *cmp, *m;
	tiimat3_Poly *p;
	size_t i;

	cmp = msg_alloc_rand();
	m = tiimat3_msg_alloc(1);
	p = tiimat3_util_alloc(TIIMAT3_QLEN, sizeof *p);

	tiimat3_msg_pack(m, cmp);
	for (i = 0; i < TIIMAT3_QLEN; ++i)
		tiimat3_bgv_encode(i, &p[i], m);
	tiimat3_bgv_decode(m, p);
	tiimat3_msg_unpack(m, m);

	assert(msg_equal(cmp, m));

	tiimat3_util_dealloc(p);
	tiimat3_msg_dealloc(cmp, 1);
	tiimat3_msg_dealloc(m, 1);
}

void
test_encrypt(void)
{
	tiimat3_Seed seed[2];
	tiimat3_KeySecret *sk;
	tiimat3_KeyPublic *pk;

	tiimat3_Message *cmp, *m;
	tiimat3_Ciphertext *ct;
	tiimat3_Poly *p;
	size_t i;

	sk = tiimat3_util_alloc(1, sizeof *sk);
	pk = tiimat3_util_alloc(1, sizeof *pk);
	cmp = msg_alloc_rand();
	m = tiimat3_msg_alloc(1);
	ct = tiimat3_util_alloc(1, sizeof *ct);
	p = tiimat3_util_alloc(TIIMAT3_QLEN, sizeof *p);

	tiimat3_bgv_keygen_secret(sk);
	for (i = 0; i < LEN(seed); ++i)
		tiimat3_random_seed(&seed[i]);

	tiimat3_msg_pack(m, cmp);
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		tiimat3_bgv_encode(i, &p[i], m);
		tiimat3_bgv_keygen_public(i, pk, sk, &seed[0]);
		tiimat3_bgv_encrypt(i, ct, pk, &p[i], &seed[1]);
		tiimat3_bgv_decrypt(i, &p[i], sk, ct);
	}
	tiimat3_bgv_decode(m, p);
	tiimat3_msg_unpack(m, m);
	assert(msg_equal(cmp, m));

	tiimat3_util_dealloc(sk);
	tiimat3_util_dealloc(pk);
	tiimat3_msg_dealloc(cmp, 1);
	tiimat3_msg_dealloc(m, 1);
	tiimat3_util_dealloc(ct);
	tiimat3_util_dealloc(p);
}

void
test_arithmetic(void)
{
	tiimat3_Seed seed[4];
	tiimat3_KeySecret *sk;
	tiimat3_KeyPublic *pk;

	tiimat3_Message *cmp, *m[4];
	tiimat3_Ciphertext *ct;
	tiimat3_Poly *p;
	size_t i;

	sk = tiimat3_util_alloc(1, sizeof *sk);
	pk = tiimat3_util_alloc(1, sizeof *pk);
	cmp = tiimat3_msg_alloc(1);
	m[0] = msg_alloc_rand();
	m[1] = msg_alloc_rand();
	m[2] = msg_alloc_rand();
	m[3] = msg_alloc_rand();
	ct = tiimat3_util_alloc(3, sizeof *ct);
	p = tiimat3_util_alloc(TIIMAT3_QLEN, sizeof *p);

	tiimat3_bgv_keygen_secret(sk);
	for (i = 0; i < LEN(seed); ++i)
		tiimat3_random_seed(&seed[i]);

	msg_mul(cmp, m[0], m[3]);
	msg_add(cmp, cmp, m[1]);
	msg_mul(cmp, cmp, m[2]);

	tiimat3_msg_pack(m[0], m[0]);
	tiimat3_msg_pack(m[1], m[1]);
	tiimat3_msg_pack(m[2], m[2]);
	tiimat3_msg_pack(m[3], m[3]);

	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		tiimat3_bgv_keygen_public(i, pk, sk, &seed[0]);

		tiimat3_bgv_encode(i, &p[i], m[0]);
		tiimat3_bgv_encrypt(i, &ct[0], pk, &p[i], &seed[1]);
		tiimat3_bgv_encode(i, &p[i], m[1]);
		tiimat3_bgv_encrypt(i, &ct[1], pk, &p[i], &seed[2]);
		tiimat3_bgv_encode(i, &p[i], m[2]);
		tiimat3_bgv_encrypt(i, &ct[2], pk, &p[i], &seed[3]);
		tiimat3_bgv_encode(i, &p[i], m[3]);

		tiimat3_bgv_mulc(i, &ct[0], &ct[0], &p[i]);
		tiimat3_bgv_add(i, &ct[0], &ct[0], &ct[1]);
		tiimat3_bgv_mul(i, &ct[0], &ct[0], &ct[2]);

		tiimat3_bgv_decrypt(i, &p[i], sk, &ct[0]);
	}

	tiimat3_bgv_decode(m[0], p);
	tiimat3_msg_unpack(m[0], m[0]);
	assert(msg_equal(cmp, m[0]));

	tiimat3_util_dealloc(sk);
	tiimat3_util_dealloc(pk);
	tiimat3_msg_dealloc(cmp, 1);
	for (i = 0; i < 4; ++i)
		tiimat3_msg_dealloc(m[i], 1);
	tiimat3_util_dealloc(ct);
	tiimat3_util_dealloc(p);
}

void
test_modswitch(void)
{
	tiimat3_Seed seed[2];
	tiimat3_KeySecret *sk;
	tiimat3_KeyPublic *pk;

	tiimat3_Message *cmp, *m;
	tiimat3_Ciphertext *ct;
	tiimat3_Delta *delta;
	tiimat3_Poly *p;
	size_t i;

	sk = tiimat3_util_alloc(1, sizeof *sk);
	pk = tiimat3_util_alloc(1, sizeof *pk);
	cmp = msg_alloc_rand();
	m = tiimat3_msg_alloc(1);
	ct = tiimat3_util_alloc(TIIMAT3_QLEN, sizeof *ct);
	delta = tiimat3_util_alloc(2, sizeof *delta);
	p = tiimat3_util_alloc(TIIMAT3_QLEN, sizeof *p);

	tiimat3_bgv_keygen_secret(sk);
	for (i = 0; i < LEN(seed); ++i)
		tiimat3_random_seed(&seed[i]);

	tiimat3_msg_pack(m, cmp);
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		tiimat3_bgv_keygen_public(i, pk, sk, &seed[0]);

		tiimat3_bgv_encode(i, &p[i], m);
		tiimat3_bgv_encrypt(i, &ct[i], pk, &p[i], &seed[1]);

		if (i == 0) {
			tiimat3_bgv_modswitch_delta(0, &delta[0], &ct[0]);
			continue;
		} else
			tiimat3_bgv_modswitch_ext(i, &ct[i], &delta[0]);

		if (i == 1) {
			tiimat3_bgv_modswitch_delta(1, &delta[1], &ct[1]);
			continue;
		} else
			tiimat3_bgv_modswitch_ext(i, &ct[i], &delta[1]);

		tiimat3_bgv_decrypt(i, &p[i], sk, &ct[i]);
	}

	tiimat3_mod_drop(0, 1);
	tiimat3_mod_drop(1, 1);

	tiimat3_bgv_decode(m, p);
	tiimat3_msg_unpack(m, m);
	assert(msg_equal(cmp, m));

	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	tiimat3_util_dealloc(delta);
	tiimat3_util_dealloc(ct);
	tiimat3_util_dealloc(p);
	tiimat3_msg_dealloc(cmp, 1);
	tiimat3_msg_dealloc(m, 1);
	tiimat3_util_dealloc(pk);
	tiimat3_util_dealloc(sk);
}

void
test_keyswitch(void)
{
	tiimat3_Seed seed[TIIMAT3_OMEGA + 3];
	tiimat3_KeySecret *sk;
	tiimat3_KeyPublic *pk;
	tiimat3_KeySwitch *ksw, *ksw2, *kswr;

	tiimat3_Message *cmp, *m;
	tiimat3_Ciphertext *ct;
	tiimat3_CiphertextSwitch *csw, *cswr;
	tiimat3_Delta *delta;
	tiimat3_Poly *p;
	size_t i;

	sk = tiimat3_util_alloc(1, sizeof *sk);
	pk = tiimat3_util_alloc(1, sizeof *pk);
	ksw = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *ksw);
	ksw2 = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *ksw2);
	kswr = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *kswr);
	cmp = msg_alloc_rand();
	m = tiimat3_msg_alloc(1);
	ct = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *ct);
	csw = tiimat3_util_alloc(TIIMAT3_QPLEN, sizeof *csw);
	cswr = tiimat3_util_alloc(1, sizeof *cswr);
	delta = tiimat3_util_alloc(TIIMAT3_PLEN, sizeof *delta);
	p = tiimat3_util_alloc(TIIMAT3_QLEN, sizeof *p);

	tiimat3_bgv_keygen_secret(sk);

	for (i = 0; i < LEN(seed); ++i)
		tiimat3_random_seed(&seed[i]);
	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_bgv_keygen_switch(i, &ksw[i], sk, sk, &seed[3]);

	for (i = 3; i < LEN(seed); ++i)
		tiimat3_random_seed(&seed[i]);
	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_bgv_keygen_switchr(i, &kswr[i], sk, 1, &seed[3]);

	for (i = 3; i < LEN(seed); ++i)
		tiimat3_random_seed(&seed[i]);
	for (i = 0; i < TIIMAT3_QPLEN; ++i)
		tiimat3_bgv_keygen_switch2(i, &ksw2[i], sk, &seed[3]);

	tiimat3_msg_pack(m, cmp);
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		tiimat3_bgv_keygen_public(i, pk, sk, &seed[0]);
		tiimat3_bgv_encode(i, &p[i], m);
		tiimat3_bgv_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
	}
	tiimat3_bgv_keyswitch(ct, ksw, 1);
	for (i = 0; i < TIIMAT3_QLEN; ++i)
		tiimat3_bgv_decrypt(i, &p[i], sk, &ct[i]);
	tiimat3_bgv_decode(m, p);
	tiimat3_msg_unpack(m, m);
	assert(msg_equal(cmp, m));

	tiimat3_msg_pack(m, cmp);
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		tiimat3_bgv_keygen_public(i, pk, sk, &seed[0]);
		tiimat3_bgv_encode(i, &p[i], m);
		tiimat3_bgv_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
		tiimat3_bgv_rot_inplace(i, &ct[i], 1);
	}
	tiimat3_bgv_keyswitch(ct, kswr, 1);
	for (i = 0; i < TIIMAT3_QLEN; ++i)
		tiimat3_bgv_decrypt(i, &p[i], sk, &ct[i]);
	tiimat3_bgv_decode(m, p);
	tiimat3_msg_unpack(m, m);
	msg_rot(cmp, cmp, 1);
	assert(msg_equal(cmp, m));

	tiimat3_msg_pack(m, cmp);
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		tiimat3_bgv_keygen_public(i, pk, sk, &seed[0]);
		tiimat3_bgv_encode(i, &p[i], m);
		tiimat3_bgv_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
		tiimat3_bgv_mul(i, &ct[i], &ct[i], &ct[i]);
	}
	tiimat3_bgv_keyswitch(ct, ksw2, 2);
	for (i = 0; i < TIIMAT3_QLEN; ++i)
		tiimat3_bgv_decrypt(i, &p[i], sk, &ct[i]);
	tiimat3_bgv_decode(m, p);
	tiimat3_msg_unpack(m, m);
	msg_mul(cmp, cmp, cmp);
	assert(msg_equal(cmp, m));

	tiimat3_msg_pack(m, cmp);
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		tiimat3_bgv_keygen_public(i, pk, sk, &seed[0]);
		tiimat3_bgv_encode(i, &p[i], m);
		tiimat3_bgv_encrypt(i, &ct[i], pk, &p[i], &seed[1]);

		if (i == 0) {
			tiimat3_bgv_modswitch_delta(0, &delta[0], &ct[0]);
			continue;
		} else
			tiimat3_bgv_modswitch_ext(i, &ct[i], &delta[0]);

		if (i == 1) {
			tiimat3_bgv_modswitch_delta(1, &delta[1], &ct[1]);
			continue;
		} else
			tiimat3_bgv_modswitch_ext(i, &ct[i], &delta[1]);
	}
	tiimat3_mod_drop(0, 1);
	tiimat3_mod_drop(1, 1);
	tiimat3_bgv_keyswitch(ct, ksw, 1);
	for (i = 2; i < TIIMAT3_QLEN; ++i)
		tiimat3_bgv_decrypt(i, &p[i], sk, &ct[i]);
	tiimat3_bgv_decode(m, p);
	tiimat3_msg_unpack(m, m);
	assert(msg_equal(cmp, m));

	tiimat3_mod_drop(0, 0);
	tiimat3_mod_drop(1, 0);

	tiimat3_msg_pack(m, cmp);
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		tiimat3_bgv_keygen_public(i, pk, sk, &seed[0]);
		tiimat3_bgv_encode(i, &p[i], m);
		tiimat3_bgv_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
	}
	tiimat3_bgv_keyswitch_ext(csw, ct, 1);
	for (i = 0; i < TIIMAT3_QPLEN; ++i) {
		tiimat3_bgv_rot_csw(i, cswr, &csw[i], 1);
		tiimat3_bgv_keyswitch_dot(i, &ct[i], cswr, &kswr[i]);
	}
	for (i = TIIMAT3_QLEN; i < TIIMAT3_QPLEN; ++i)
		tiimat3_bgv_keyswitch_delta(i, &delta[i - TIIMAT3_QLEN], &ct[i]);
	for (i = 0; i < TIIMAT3_QLEN; ++i) {
		tiimat3_bgv_keyswitch_switch(i, &ct[i], delta);
		tiimat3_bgv_decrypt(i, &p[i], sk, &ct[i]);
	}
	tiimat3_bgv_decode(m, p);
	tiimat3_msg_unpack(m, m);
	msg_rot(cmp, cmp, 1);
	assert(msg_equal(cmp, m));

	tiimat3_util_dealloc(sk);
	tiimat3_util_dealloc(pk);
	tiimat3_util_dealloc(ksw);
	tiimat3_util_dealloc(ksw2);
	tiimat3_util_dealloc(kswr);

	tiimat3_msg_dealloc(cmp, 1);
	tiimat3_msg_dealloc(m, 1);
	tiimat3_util_dealloc(ct);
	tiimat3_util_dealloc(csw);
	tiimat3_util_dealloc(cswr);
	tiimat3_util_dealloc(delta);
	tiimat3_util_dealloc(p);
}

int
main(void)
{

	tiimat3_bgv_init();

	fputs("[+] Testing BGV packing:      ", stderr);
	test_pack(), fputs("1/1.\n", stderr);

	fputs("[+] Testing BGV encrypt:      ", stderr);
	test_encrypt(), fputs("1/1.\n", stderr);

	fputs("[+] Testing BGV arithmetic:   ", stderr);
	test_arithmetic(), fputs("1/1.\n", stderr);

	fputs("[+] Testing BGV modswitch:    ", stderr);
	test_modswitch(), fputs("1/1.\n", stderr);

	fputs("[+] Testing BGV keyswitch:    ", stderr);
	test_keyswitch(), fputs("5/5.\n", stderr);

	tiimat3_bgv_deinit();

	return 0;
}
