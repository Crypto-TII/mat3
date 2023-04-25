#define BGV_BGV_IMPL
#define BGV_MSG_IMPL
#define BGV_RANDOM_IMPL
#define BGV_RING_IMPL
#include "bgv.h"

#define BGV_UTIL_IMPL
#include "util.h"

#include "ring_hexl.h"

#define LEN(a) (sizeof a / sizeof *a)

void
msg_add(bgv_Message *rop, const bgv_Message *op1, const bgv_Message *op2)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j) {
		mpz_add(rop->value[j], op1->value[j], op2->value[j]);
		mpz_mod(rop->value[j], rop->value[j], bgv_t.value);
	}
}

bgv_Message *
msg_alloc_rand(void)
{
	gmp_randstate_t state;
	bgv_Message *m;
	size_t j;

	gmp_randinit_default(state);

	m = bgv_alloc_msg(1);
	for (j = 0; j < BGV_D; ++j)
		mpz_urandomm(m->value[j], state, bgv_t.value);

	gmp_randclear(state);

	return m;
}

int
msg_equal(const bgv_Message *m1, const bgv_Message *m2)
{
	size_t j;
	int cmp;

	cmp = 0;
	for (j = 0; j < BGV_D; ++j)
		cmp += mpz_cmp(m1->value[j], m2->value[j]);

	return cmp == 0;
}

void
msg_mul(bgv_Message *rop, const bgv_Message *op1, const bgv_Message *op2)
{
	size_t j;

	for (j = 0; j < BGV_D; ++j) {
		mpz_mul(rop->value[j], op1->value[j], op2->value[j]);
		mpz_mod(rop->value[j], rop->value[j], bgv_t.value);
	}
}

void
msg_rot(bgv_Message *rop, const bgv_Message *m, size_t steps)
{
	mpz_t tmp;
	size_t i, j;

	mpz_init(tmp);

	for (i = 0; i < steps; ++i) {
		mpz_set(tmp, m->value[0]);
		for (j = 0; j < BGV_D / 2 - 1; ++j)
			mpz_set(rop->value[j], m->value[j + 1]);
		mpz_set(rop->value[j], tmp);

		mpz_set(tmp, m->value[BGV_D / 2]);
		for (j = BGV_D / 2; j < BGV_D - 1; ++j)
			mpz_set(rop->value[j], m->value[j + 1]);
		mpz_set(rop->value[j], tmp);
	}

	mpz_clear(tmp);
}

void
test_pack(void)
{
	bgv_Message *cmp, *m;
	bgv_Poly *p;
	size_t i;

	cmp = msg_alloc_rand();
	m = bgv_alloc_msg(1);
	p = bgv_alloc(BGV_QLEN, sizeof *p);

	bgv_pack(m, cmp);
	for (i = 0; i < BGV_QLEN; ++i)
		bgv_encode(i, &p[i], m);
	bgv_decode(m, p);
	bgv_unpack(m, m);

	assert(msg_equal(cmp, m));

	bgv_dealloc(p);
	bgv_dealloc_msg(cmp, 1);
	bgv_dealloc_msg(m, 1);
}

void
test_encrypt(void)
{
	tiimat3_Seed seed[2];
	bgv_KeySecret *sk;
	bgv_KeyPublic *pk;

	bgv_Message *cmp, *m;
	bgv_Ciphertext *ct;
	bgv_Poly *p;
	size_t i;

	sk = bgv_alloc(1, sizeof *sk);
	pk = bgv_alloc(1, sizeof *pk);
	cmp = msg_alloc_rand();
	m = bgv_alloc_msg(1);
	ct = bgv_alloc(1, sizeof *ct);
	p = bgv_alloc(BGV_QLEN, sizeof *p);

	bgv_keygen_secret(sk);
	for (i = 0; i < LEN(seed); ++i)
		tiimat3_random_seed(&seed[i]);

	bgv_pack(m, cmp);
	for (i = 0; i < BGV_QLEN; ++i) {
		bgv_encode(i, &p[i], m);
		bgv_keygen_public(i, pk, sk, &seed[0]);
		bgv_encrypt(i, ct, pk, &p[i], &seed[1]);
		bgv_decrypt(i, &p[i], sk, ct);
	}
	bgv_decode(m, p);
	bgv_unpack(m, m);
	assert(msg_equal(cmp, m));

	bgv_dealloc(sk);
	bgv_dealloc(pk);
	bgv_dealloc_msg(cmp, 1);
	bgv_dealloc_msg(m, 1);
	bgv_dealloc(ct);
	bgv_dealloc(p);
}

void
test_arithmetic(void)
{
	tiimat3_Seed seed[4];
	bgv_KeySecret *sk;
	bgv_KeyPublic *pk;

	bgv_Message *cmp, *m[4];
	bgv_Ciphertext *ct;
	bgv_Poly *p;
	size_t i;

	sk = bgv_alloc(1, sizeof *sk);
	pk = bgv_alloc(1, sizeof *pk);
	cmp = bgv_alloc_msg(1);
	m[0] = msg_alloc_rand();
	m[1] = msg_alloc_rand();
	m[2] = msg_alloc_rand();
	m[3] = msg_alloc_rand();
	ct = bgv_alloc(3, sizeof *ct);
	p = bgv_alloc(BGV_QLEN, sizeof *p);

	bgv_keygen_secret(sk);
	for (i = 0; i < LEN(seed); ++i)
		tiimat3_random_seed(&seed[i]);

	msg_mul(cmp, m[0], m[3]);
	msg_add(cmp, cmp, m[1]);
	msg_mul(cmp, cmp, m[2]);

	bgv_pack(m[0], m[0]);
	bgv_pack(m[1], m[1]);
	bgv_pack(m[2], m[2]);
	bgv_pack(m[3], m[3]);

	for (i = 0; i < BGV_QLEN; ++i) {
		bgv_keygen_public(i, pk, sk, &seed[0]);

		bgv_encode(i, &p[i], m[0]);
		bgv_encrypt(i, &ct[0], pk, &p[i], &seed[1]);
		bgv_encode(i, &p[i], m[1]);
		bgv_encrypt(i, &ct[1], pk, &p[i], &seed[2]);
		bgv_encode(i, &p[i], m[2]);
		bgv_encrypt(i, &ct[2], pk, &p[i], &seed[3]);
		bgv_encode(i, &p[i], m[3]);

		bgv_mulc(i, &ct[0], &ct[0], &p[i]);
		bgv_add(i, &ct[0], &ct[0], &ct[1]);
		bgv_mul(i, &ct[0], &ct[0], &ct[2]);

		bgv_decrypt(i, &p[i], sk, &ct[0]);
	}

	bgv_decode(m[0], p);
	bgv_unpack(m[0], m[0]);
	assert(msg_equal(cmp, m[0]));

	bgv_dealloc(sk);
	bgv_dealloc(pk);
	bgv_dealloc_msg(cmp, 1);
	for (i = 0; i < LEN(m); ++i)
		bgv_dealloc_msg(m[i], 1);
	bgv_dealloc(ct);
	bgv_dealloc(p);
}

void
test_modswitch(void)
{
	tiimat3_Seed seed[2];
	bgv_KeySecret *sk;
	bgv_KeyPublic *pk;

	bgv_Message *cmp, *m;
	bgv_Ciphertext *ct;
	bgv_Delta *delta;
	bgv_Poly *p;
	size_t i;

	sk = bgv_alloc(1, sizeof *sk);
	pk = bgv_alloc(1, sizeof *pk);
	cmp = msg_alloc_rand();
	m = bgv_alloc_msg(1);
	ct = bgv_alloc(BGV_QLEN, sizeof *ct);
	delta = bgv_alloc(2, sizeof *delta);
	p = bgv_alloc(BGV_QLEN, sizeof *p);

	bgv_keygen_secret(sk);
	for (i = 0; i < LEN(seed); ++i)
		tiimat3_random_seed(&seed[i]);

	bgv_pack(m, cmp);
	for (i = 0; i < BGV_QLEN; ++i) {
		bgv_keygen_public(i, pk, sk, &seed[0]);

		bgv_encode(i, &p[i], m);
		bgv_encrypt(i, &ct[i], pk, &p[i], &seed[1]);

		if (i == 0) {
			bgv_modswitch_delta(0, &delta[0], &ct[0]);
			continue;
		} else
			bgv_modswitch_ext(i, &ct[i], &delta[0]);

		if (i == 1) {
			bgv_modswitch_delta(1, &delta[1], &ct[1]);
			continue;
		} else
			bgv_modswitch_ext(i, &ct[i], &delta[1]);

		bgv_decrypt(i, &p[i], sk, &ct[i]);
	}

	bgv_mod_drop(0, 1);
	bgv_mod_drop(1, 1);

	bgv_decode(m, p);
	bgv_unpack(m, m);
	assert(msg_equal(cmp, m));

	bgv_mod_drop(0, 0);
	bgv_mod_drop(1, 0);

	bgv_dealloc(delta);
	bgv_dealloc(ct);
	bgv_dealloc(p);
	bgv_dealloc_msg(cmp, 1);
	bgv_dealloc_msg(m, 1);
	bgv_dealloc(pk);
	bgv_dealloc(sk);
}

void
test_keyswitch(void)
{
	tiimat3_Seed seed[BGV_OMEGA + 3];
	bgv_KeySecret *sk;
	bgv_KeyPublic *pk;
	bgv_KeySwitch *ksw, *ksw2, *kswr;

	bgv_Message *cmp, *m;
	bgv_Ciphertext *ct;
	bgv_CiphertextSwitch *csw, *cswr;
	bgv_Delta *delta;
	bgv_Poly *p;
	size_t i;

	sk = bgv_alloc(1, sizeof *sk);
	pk = bgv_alloc(1, sizeof *pk);
	ksw = bgv_alloc(BGV_QPLEN, sizeof *ksw);
	ksw2 = bgv_alloc(BGV_QPLEN, sizeof *ksw2);
	kswr = bgv_alloc(BGV_QPLEN, sizeof *kswr);
	cmp = msg_alloc_rand();
	m = bgv_alloc_msg(1);
	ct = bgv_alloc(BGV_QPLEN, sizeof *ct);
	csw = bgv_alloc(BGV_QPLEN, sizeof *csw);
	cswr = bgv_alloc(1, sizeof *cswr);
	delta = bgv_alloc(BGV_PLEN, sizeof *delta);
	p = bgv_alloc(BGV_QLEN, sizeof *p);

	bgv_keygen_secret(sk);

	for (i = 0; i < LEN(seed); ++i)
		tiimat3_random_seed(&seed[i]);
	for (i = 0; i < BGV_QPLEN; ++i)
		bgv_keygen_switch(i, &ksw[i], sk, sk, &seed[3]);

	for (i = 3; i < LEN(seed); ++i)
		tiimat3_random_seed(&seed[i]);
	for (i = 0; i < BGV_QPLEN; ++i)
		bgv_keygen_switchr(i, &kswr[i], sk, 1, &seed[3]);

	for (i = 3; i < LEN(seed); ++i)
		tiimat3_random_seed(&seed[i]);
	for (i = 0; i < BGV_QPLEN; ++i)
		bgv_keygen_switch2(i, &ksw2[i], sk, &seed[3]);

	bgv_pack(m, cmp);
	for (i = 0; i < BGV_QLEN; ++i) {
		bgv_keygen_public(i, pk, sk, &seed[0]);
		bgv_encode(i, &p[i], m);
		bgv_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
	}
	bgv_keyswitch(ct, ksw, 1);
	for (i = 0; i < BGV_QLEN; ++i)
		bgv_decrypt(i, &p[i], sk, &ct[i]);
	bgv_decode(m, p);
	bgv_unpack(m, m);
	assert(msg_equal(cmp, m));

	bgv_pack(m, cmp);
	for (i = 0; i < BGV_QLEN; ++i) {
		bgv_keygen_public(i, pk, sk, &seed[0]);
		bgv_encode(i, &p[i], m);
		bgv_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
		bgv_rot_inplace(i, &ct[i], 1);
	}
	bgv_keyswitch(ct, kswr, 1);
	for (i = 0; i < BGV_QLEN; ++i)
		bgv_decrypt(i, &p[i], sk, &ct[i]);
	bgv_decode(m, p);
	bgv_unpack(m, m);
	msg_rot(cmp, cmp, 1);
	assert(msg_equal(cmp, m));

	bgv_pack(m, cmp);
	for (i = 0; i < BGV_QLEN; ++i) {
		bgv_keygen_public(i, pk, sk, &seed[0]);
		bgv_encode(i, &p[i], m);
		bgv_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
		bgv_mul(i, &ct[i], &ct[i], &ct[i]);
	}
	bgv_keyswitch(ct, ksw2, 2);
	for (i = 0; i < BGV_QLEN; ++i)
		bgv_decrypt(i, &p[i], sk, &ct[i]);
	bgv_decode(m, p);
	bgv_unpack(m, m);
	msg_mul(cmp, cmp, cmp);
	assert(msg_equal(cmp, m));

	bgv_pack(m, cmp);
	for (i = 0; i < BGV_QLEN; ++i) {
		bgv_keygen_public(i, pk, sk, &seed[0]);
		bgv_encode(i, &p[i], m);
		bgv_encrypt(i, &ct[i], pk, &p[i], &seed[1]);

		if (i == 0) {
			bgv_modswitch_delta(0, &delta[0], &ct[0]);
			continue;
		} else
			bgv_modswitch_ext(i, &ct[i], &delta[0]);

		if (i == 1) {
			bgv_modswitch_delta(1, &delta[1], &ct[1]);
			continue;
		} else
			bgv_modswitch_ext(i, &ct[i], &delta[1]);
	}
	bgv_mod_drop(0, 1);
	bgv_mod_drop(1, 1);
	bgv_keyswitch(ct, ksw, 1);
	for (i = 2; i < BGV_QLEN; ++i)
		bgv_decrypt(i, &p[i], sk, &ct[i]);
	bgv_decode(m, p);
	bgv_unpack(m, m);
	assert(msg_equal(cmp, m));

	bgv_mod_drop(0, 0);
	bgv_mod_drop(1, 0);

	bgv_pack(m, cmp);
	for (i = 0; i < BGV_QLEN; ++i) {
		bgv_keygen_public(i, pk, sk, &seed[0]);
		bgv_encode(i, &p[i], m);
		bgv_encrypt(i, &ct[i], pk, &p[i], &seed[1]);
	}
	bgv_keyswitch_ext(csw, ct, 1);
	for (i = 0; i < BGV_QPLEN; ++i) {
		bgv_rot_csw(i, cswr, &csw[i], 1);
		bgv_keyswitch_dot(i, &ct[i], cswr, &kswr[i]);
	}
	for (i = BGV_QLEN; i < BGV_QPLEN; ++i)
		bgv_keyswitch_delta(i, &delta[i - BGV_QLEN], &ct[i]);
	for (i = 0; i < BGV_QLEN; ++i) {
		bgv_keyswitch_switch(i, &ct[i], delta);
		bgv_decrypt(i, &p[i], sk, &ct[i]);
	}
	bgv_decode(m, p);
	bgv_unpack(m, m);
	msg_rot(cmp, cmp, 1);
	assert(msg_equal(cmp, m));

	bgv_dealloc(sk);
	bgv_dealloc(pk);
	bgv_dealloc(ksw);
	bgv_dealloc(ksw2);
	bgv_dealloc(kswr);

	bgv_dealloc_msg(cmp, 1);
	bgv_dealloc_msg(m, 1);
	bgv_dealloc(ct);
	bgv_dealloc(csw);
	bgv_dealloc(cswr);
	bgv_dealloc(delta);
	bgv_dealloc(p);
}

int
main(void)
{

	bgv_init();

	fputs("[+] Testing BGV packing:    ", stderr);
	test_pack(), fputs("1/1.\n", stderr);

	fputs("[+] Testing BGV encrypt:    ", stderr);
	test_encrypt(), fputs("1/1.\n", stderr);

	fputs("[+] Testing BGV arithmetic: ", stderr);
	test_arithmetic(), fputs("1/1.\n", stderr);

	fputs("[+] Testing BGV modswitch:  ", stderr);
	test_modswitch(), fputs("1/1.\n", stderr);

	fputs("[+] Testing BGV keyswitch:  ", stderr);
	test_keyswitch(), fputs("5/5.\n", stderr);

	bgv_deinit();

	return 0;
}
