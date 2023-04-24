#ifndef BGV_RANDOM_H
#define BGV_RANDOM_H

#include <sys/random.h>
#include <assert.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <x86intrin.h>

typedef struct bgv_seed bgv_Seed;
struct bgv_seed {
	uint64_t state[25];
	unsigned char value[16];
	size_t pos;
	int init;
};

void bgv_sample_bytes(unsigned char *out, size_t len, bgv_Seed *seed);
void bgv_sample_cbd1_i8(int8_t *out, size_t len, bgv_Seed *seed);
void bgv_sample_cbd1_i64(int64_t *out, size_t len, bgv_Seed *seed);
void bgv_sample_cbd21_i8(int8_t *out, size_t len, bgv_Seed *seed);
void bgv_sample_cbd21_i64(int64_t *out, size_t len, bgv_Seed *seed);
void bgv_sample_seed(bgv_Seed *seed);
void bgv_sample_uniform_i64(int64_t *out, size_t len, uint64_t mod, bgv_Seed *seed);
void bgv_sample_uniform_u64(uint64_t *out, size_t len, uint64_t mod, bgv_Seed *seed);

#endif /* BGV_RANDOM_H */


#ifdef BGV_RANDOM_IMPL

/* https://github.com/pq-crystals/dilithium/blob/master/ref/fips202.c */
#define KECCAK_SEED_LEN  16
#define KECCAK_STATE_LEN 25
#define SHA3_384_BYTES   48
#define SHA3_384_RATE    104
#define SHAKE128_RATE    168

#define BSR64(x)    __bsrq(x)
#define LEN(a)      (sizeof a / sizeof *a)
#define POPCNT32(x) __popcntd(x)
#define ROL64(x, r) __rolq(x, r)

#define sha3_384_absorb(s, in, len) keccak_absorb(s, in, len, 104, 0x06)
#define sha3_384_squeeze(s, out)    keccak_squeeze(s, out, 104)
#define shake128_absorb(s, in, len) keccak_absorb(s, in, len, 168, 0x1f)
#define shake128_squeeze(s, out)    keccak_squeeze(s, out, 168)

static void keccak_absorb(uint64_t *state, const unsigned char *in, size_t len, size_t rate, unsigned char suffix);
static void keccak_permute(uint64_t *state);
static void keccak_squeeze(uint64_t *state, unsigned char *out, size_t rate);
static void random_bytes(void *out, size_t len, bgv_Seed *seed);

static void
keccak_absorb(uint64_t *state, const unsigned char *in, size_t len, size_t rate, unsigned char suffix)
{
	unsigned char *s;
	size_t i;

	for (i = 0; i < KECCAK_STATE_LEN; ++i)
		state[i] = 0;

	s = (unsigned char *)state;
	while (len > rate) {
		for (i = 0; i < rate; ++i)
			s[i] ^= in[i];
		in += rate;
		len -= rate;

		keccak_permute(state);
	}
	for (i = 0; i < len; ++i)
		s[i] ^= in[i];

	if (len != 0)
		s[len] ^= suffix;
	s[rate - 1] ^= 0x80;
}

static void
keccak_permute(uint64_t *state)
{
	const uint64_t constants[] = {
		0x0000000000000001,
		0x0000000000008082,
		0x800000000000808a,
		0x8000000080008000,
		0x000000000000808b,
		0x0000000080000001,
		0x8000000080008081,
		0x8000000000008009,
		0x000000000000008a,
		0x0000000000000088,
		0x0000000080008009,
		0x000000008000000a,
		0x000000008000808b,
		0x800000000000008b,
		0x8000000000008089,
		0x8000000000008003,
		0x8000000000008002,
		0x8000000000000080,
		0x000000000000800a,
		0x800000008000000a,
		0x8000000080008081,
		0x8000000000008080,
		0x0000000080000001,
		0x8000000080008008,
	};

	uint64_t Aba, Abe, Abi, Abo, Abu;
	uint64_t Aga, Age, Agi, Ago, Agu;
	uint64_t Aka, Ake, Aki, Ako, Aku;
	uint64_t Ama, Ame, Ami, Amo, Amu;
	uint64_t Asa, Ase, Asi, Aso, Asu;
	uint64_t BCa, BCe, BCi, BCo, BCu;
	uint64_t Da, De, Di, Do, Du;
	uint64_t Eba, Ebe, Ebi, Ebo, Ebu;
	uint64_t Ega, Ege, Egi, Ego, Egu;
	uint64_t Eka, Eke, Eki, Eko, Eku;
	uint64_t Ema, Eme, Emi, Emo, Emu;
	uint64_t Esa, Ese, Esi, Eso, Esu;
	size_t round;

	Aba = state[0];
	Abe = state[1];
	Abi = state[2];
	Abo = state[3];
	Abu = state[4];
	Aga = state[5];
	Age = state[6];
	Agi = state[7];
	Ago = state[8];
	Agu = state[9];
	Aka = state[10];
	Ake = state[11];
	Aki = state[12];
	Ako = state[13];
	Aku = state[14];
	Ama = state[15];
	Ame = state[16];
	Ami = state[17];
	Amo = state[18];
	Amu = state[19];
	Asa = state[20];
	Ase = state[21];
	Asi = state[22];
	Aso = state[23];
	Asu = state[24];

	for (round = 0; round < 24; round += 2) {
		BCa = Aba ^ Aga ^ Aka ^ Ama ^ Asa;
		BCe = Abe ^ Age ^ Ake ^ Ame ^ Ase;
		BCi = Abi ^ Agi ^ Aki ^ Ami ^ Asi;
		BCo = Abo ^ Ago ^ Ako ^ Amo ^ Aso;
		BCu = Abu ^ Agu ^ Aku ^ Amu ^ Asu;

		Da = BCu ^ ROL64(BCe, 1);
		De = BCa ^ ROL64(BCi, 1);
		Di = BCe ^ ROL64(BCo, 1);
		Do = BCi ^ ROL64(BCu, 1);
		Du = BCo ^ ROL64(BCa, 1);

		Aba ^= Da;
		BCa = Aba;
		Age ^= De;
		BCe = ROL64(Age, 44);
		Aki ^= Di;
		BCi = ROL64(Aki, 43);
		Amo ^= Do;
		BCo = ROL64(Amo, 21);
		Asu ^= Du;
		BCu = ROL64(Asu, 14);
		Eba = BCa ^((~BCe) & BCi);
		Eba ^= constants[round];
		Ebe = BCe ^ ((~BCi) & BCo);
		Ebi = BCi ^ ((~BCo) & BCu);
		Ebo = BCo ^ ((~BCu) & BCa);
		Ebu = BCu ^ ((~BCa) & BCe);

		Abo ^= Do;
		BCa = ROL64(Abo, 28);
		Agu ^= Du;
		BCe = ROL64(Agu, 20);
		Aka ^= Da;
		BCi = ROL64(Aka, 3);
		Ame ^= De;
		BCo = ROL64(Ame, 45);
		Asi ^= Di;
		BCu = ROL64(Asi, 61);
		Ega = BCa ^ ((~BCe) & BCi);
		Ege = BCe ^ ((~BCi) & BCo);
		Egi = BCi ^ ((~BCo) & BCu);
		Ego = BCo ^ ((~BCu) & BCa);
		Egu = BCu ^ ((~BCa) & BCe);

		Abe ^= De;
		BCa = ROL64(Abe, 1);
		Agi ^= Di;
		BCe = ROL64(Agi, 6);
		Ako ^= Do;
		BCi = ROL64(Ako, 25);
		Amu ^= Du;
		BCo = ROL64(Amu, 8);
		Asa ^= Da;
		BCu = ROL64(Asa, 18);
		Eka = BCa ^ ((~BCe) & BCi);
		Eke = BCe ^ ((~BCi) & BCo);
		Eki = BCi ^ ((~BCo) & BCu);
		Eko = BCo ^ ((~BCu) & BCa);
		Eku = BCu ^ ((~BCa) & BCe);

		Abu ^= Du;
		BCa = ROL64(Abu, 27);
		Aga ^= Da;
		BCe = ROL64(Aga, 36);
		Ake ^= De;
		BCi = ROL64(Ake, 10);
		Ami ^= Di;
		BCo = ROL64(Ami, 15);
		Aso ^= Do;
		BCu = ROL64(Aso, 56);
		Ema = BCa ^ ((~BCe) & BCi);
		Eme = BCe ^ ((~BCi) & BCo);
		Emi = BCi ^ ((~BCo) & BCu);
		Emo = BCo ^ ((~BCu) & BCa);
		Emu = BCu ^ ((~BCa) & BCe);

		Abi ^= Di;
		BCa = ROL64(Abi, 62);
		Ago ^= Do;
		BCe = ROL64(Ago, 55);
		Aku ^= Du;
		BCi = ROL64(Aku, 39);
		Ama ^= Da;
		BCo = ROL64(Ama, 41);
		Ase ^= De;
		BCu = ROL64(Ase, 2);
		Esa = BCa ^ ((~BCe) & BCi);
		Ese = BCe ^ ((~BCi) & BCo);
		Esi = BCi ^ ((~BCo) & BCu);
		Eso = BCo ^ ((~BCu) & BCa);
		Esu = BCu ^ ((~BCa) & BCe);

		BCa = Eba ^ Ega ^ Eka ^ Ema ^ Esa;
		BCe = Ebe ^ Ege ^ Eke ^ Eme ^ Ese;
		BCi = Ebi ^ Egi ^ Eki ^ Emi ^ Esi;
		BCo = Ebo ^ Ego ^ Eko ^ Emo ^ Eso;
		BCu = Ebu ^ Egu ^ Eku ^ Emu ^ Esu;

		Da = BCu ^ ROL64(BCe, 1);
		De = BCa ^ ROL64(BCi, 1);
		Di = BCe ^ ROL64(BCo, 1);
		Do = BCi ^ ROL64(BCu, 1);
		Du = BCo ^ ROL64(BCa, 1);

		Eba ^= Da;
		BCa = Eba;
		Ege ^= De;
		BCe = ROL64(Ege, 44);
		Eki ^= Di;
		BCi = ROL64(Eki, 43);
		Emo ^= Do;
		BCo = ROL64(Emo, 21);
		Esu ^= Du;
		BCu = ROL64(Esu, 14);
		Aba = BCa ^ ((~BCe) & BCi);
		Aba ^= constants[round + 1];
		Abe = BCe ^ ((~BCi) & BCo);
		Abi = BCi ^ ((~BCo) & BCu);
		Abo = BCo ^ ((~BCu) & BCa);
		Abu = BCu ^ ((~BCa) & BCe);

		Ebo ^= Do;
		BCa = ROL64(Ebo, 28);
		Egu ^= Du;
		BCe = ROL64(Egu, 20);
		Eka ^= Da;
		BCi = ROL64(Eka, 3);
		Eme ^= De;
		BCo = ROL64(Eme, 45);
		Esi ^= Di;
		BCu = ROL64(Esi, 61);
		Aga = BCa ^ ((~BCe) & BCi);
		Age = BCe ^ ((~BCi) & BCo);
		Agi = BCi ^ ((~BCo) & BCu);
		Ago = BCo ^ ((~BCu) & BCa);
		Agu = BCu ^ ((~BCa) & BCe);

		Ebe ^= De;
		BCa = ROL64(Ebe, 1);
		Egi ^= Di;
		BCe = ROL64(Egi, 6);
		Eko ^= Do;
		BCi = ROL64(Eko, 25);
		Emu ^= Du;
		BCo = ROL64(Emu, 8);
		Esa ^= Da;
		BCu = ROL64(Esa, 18);
		Aka = BCa ^ ((~BCe) & BCi);
		Ake = BCe ^ ((~BCi) & BCo);
		Aki = BCi ^ ((~BCo) & BCu);
		Ako = BCo ^ ((~BCu) & BCa);
		Aku = BCu ^ ((~BCa) & BCe);

		Ebu ^= Du;
		BCa = ROL64(Ebu, 27);
		Ega ^= Da;
		BCe = ROL64(Ega, 36);
		Eke ^= De;
		BCi = ROL64(Eke, 10);
		Emi ^= Di;
		BCo = ROL64(Emi, 15);
		Eso ^= Do;
		BCu = ROL64(Eso, 56);
		Ama = BCa ^ ((~BCe) & BCi);
		Ame = BCe ^ ((~BCi) & BCo);
		Ami = BCi ^ ((~BCo) & BCu);
		Amo = BCo ^ ((~BCu) & BCa);
		Amu = BCu ^ ((~BCa) & BCe);

		Ebi ^= Di;
		BCa = ROL64(Ebi, 62);
		Ego ^= Do;
		BCe = ROL64(Ego, 55);
		Eku ^= Du;
		BCi = ROL64(Eku, 39);
		Ema ^= Da;
		BCo = ROL64(Ema, 41);
		Ese ^= De;
		BCu = ROL64(Ese, 2);
		Asa = BCa ^ ((~BCe) & BCi);
		Ase = BCe ^ ((~BCi) & BCo);
		Asi = BCi ^ ((~BCo) & BCu);
		Aso = BCo ^ ((~BCu) & BCa);
		Asu = BCu ^ ((~BCa) & BCe);
	}

	state[0] = Aba;
	state[1] = Abe;
	state[2] = Abi;
	state[3] = Abo;
	state[4] = Abu;
	state[5] = Aga;
	state[6] = Age;
	state[7] = Agi;
	state[8] = Ago;
	state[9] = Agu;
	state[10] = Aka;
	state[11] = Ake;
	state[12] = Aki;
	state[13] = Ako;
	state[14] = Aku;
	state[15] = Ama;
	state[16] = Ame;
	state[17] = Ami;
	state[18] = Amo;
	state[19] = Amu;
	state[20] = Asa;
	state[21] = Ase;
	state[22] = Asi;
	state[23] = Aso;
	state[24] = Asu;
}

static void
keccak_squeeze(uint64_t *state, unsigned char *out, size_t rate)
{
	const unsigned char *s;
	size_t i;

	keccak_permute(state);

	s = (unsigned char *)state;
	for (i = 0; i < rate; ++i)
		out[i] = s[i];
}

static void
random_bytes(void *out, size_t len, bgv_Seed *seed)
{
	unsigned char *ptr, *s;
	size_t unused, i;

	s = (unsigned char *)seed->state;
	if (seed->init) {
		shake128_absorb(seed->state, seed->value, KECCAK_SEED_LEN);
		shake128_squeeze(seed->state, s);
		seed->init = 0;
		seed->pos = 0;
	}

	ptr = out;
	unused = SHAKE128_RATE - seed->pos;

	if (len < unused) {
		for (i = 0; i < len; ++i)
			ptr[i] = s[seed->pos++];
		return;
	}

	for (i = 0; i < unused; ++i)
		ptr[i] = s[seed->pos++];
	ptr += unused;
	len -= unused;

	while (len > SHAKE128_RATE) {
		shake128_squeeze(seed->state, ptr);
		ptr += SHAKE128_RATE;
		len -= SHAKE128_RATE;
	}

	shake128_squeeze(seed->state, s);
	for (seed->pos = 0; seed->pos < len; ++seed->pos)
		ptr[seed->pos] = s[seed->pos];
}

void
bgv_sample_bytes(unsigned char *out, size_t len, bgv_Seed *seed)
{

	random_bytes(out, len, seed);
}

void
bgv_sample_cbd1_i8(int8_t *out, size_t len, bgv_Seed *seed)
{
	size_t i;

	assert(len % 4 == 0);

	for (i = 0; i < len; /* noop */) {
		int8_t cbd[] = { -1, 0, 0, 1 };
		unsigned char r;

		random_bytes(&r, 1, seed);
		out[i++] = cbd[r & 3];
		out[i++] = cbd[(r >>= 2) & 3];
		out[i++] = cbd[(r >>= 2) & 3];
		out[i++] = cbd[r >> 2];
	}
}

void
bgv_sample_cbd1_i64(int64_t *out, size_t len, bgv_Seed *seed)
{
	size_t i;

	assert(len % 4 == 0);

	for (i = 0; i < len; /* noop */) {
		int8_t cbd[] = { -1, 0, 0, 1 };
		unsigned char r;

		random_bytes(&r, 1, seed);
		out[i++] = cbd[r & 3];
		out[i++] = cbd[(r >>= 2) & 3];
		out[i++] = cbd[(r >>= 2) & 3];
		out[i++] = cbd[r >> 2];
	}
}

void
bgv_sample_cbd21_i8(int8_t *out, size_t len, bgv_Seed *seed)
{
	size_t i;

	for (i = 0; i < len; ++i) {
		uint32_t r[2];
		random_bytes(r, 2 * sizeof *r, seed);

		r[0] &= 0x001fffff;
		r[1] &= 0x001fffff;

		out[i] = POPCNT32(r[0]) - POPCNT32(r[1]);
	}
}

void
bgv_sample_cbd21_i64(int64_t *out, size_t len, bgv_Seed *seed)
{
	size_t i;

	for (i = 0; i < len; ++i) {
		uint32_t r[2];
		random_bytes(r, 2 * sizeof *r, seed);

		r[0] &= 0x001fffff;
		r[1] &= 0x001fffff;

		out[i] = POPCNT32(r[0]) - POPCNT32(r[1]);
	}
}

void
bgv_sample_seed(bgv_Seed *seed)
{

	if (getrandom(seed->value, KECCAK_SEED_LEN, 0) != KECCAK_SEED_LEN) {
		perror("getrandom");
		exit(1);
	}
	seed->init = 1;
}

void
bgv_sample_uniform_i64(int64_t *out, size_t len, uint64_t mod, bgv_Seed *seed)
{

	bgv_sample_uniform_u64((uint64_t *)out, len, mod, seed);
}

void
bgv_sample_uniform_u64(uint64_t *out, size_t len, uint64_t mod, bgv_Seed *seed)
{
	uint64_t mask;
	size_t i;

	mask = ((uint64_t)1 << (BSR64(mod) + 1)) - 1;
	for (i = 0; i < len; ++i) {
		uint64_t r;

		do {
			random_bytes(&r, sizeof r, seed);
			r &= mask;
		} while (r >= mod);

		out[i] = r;
	}
}

#endif /* BGV_RANDOM_IMPL */
