#ifndef TIIMAT3_H
#define TIIMAT3_H

#define TIIMAT3_OMEGA 3
#define TIIMAT3_KMAX  4

#include "util.h"
#include "random.h"
#include "msg_mpz.h"
#if TIIMAT3_USE_HEXL
	#include "ring_hexl.h"
	#define tiimat3_mod_init tiimat3_mod_init_mpz
	#define tiimat3_msg_crt  tiimat3_msg_crt_u64
	#define tiimat3_poly_mod tiimat3_poly_mod_mpz
#else
	#include "ring_i64.h"
	#define tiimat3_mod_init tiimat3_mod_init_mpz
	#define tiimat3_msg_crt  tiimat3_msg_crt_i64
	#define tiimat3_poly_mod tiimat3_poly_mod_mpz
#endif /* TIIMAT3_USE_HEXL */

typedef struct tiimat3_ciphertext        tiimat3_Ciphertext;
typedef struct tiimat3_ciphertext_switch tiimat3_CiphertextSwitch;
typedef struct tiimat3_delta             tiimat3_Delta;
typedef struct tiimat3_key_public        tiimat3_KeyPublic;
typedef struct tiimat3_key_secret        tiimat3_KeySecret;
typedef struct tiimat3_key_switch        tiimat3_KeySwitch;

struct tiimat3_ciphertext {
	tiimat3_Poly poly[3];
	size_t degree;
	int ntt;
};

struct tiimat3_ciphertext_switch {
	tiimat3_Ciphertext ct;
	tiimat3_Poly ext[TIIMAT3_OMEGA];
	size_t poly;
};

struct tiimat3_delta {
	tiimat3_Ciphertext ct;
	tiimat3_Digit inv;
	size_t idx;
};

struct tiimat3_key_public {
	tiimat3_Seed seed;
	tiimat3_Poly poly;
};

struct tiimat3_key_secret {
	tiimat3_Poly poly;
};

struct tiimat3_key_switch {
	tiimat3_Seed seed;
	tiimat3_Poly poly[TIIMAT3_OMEGA];
};

extern size_t tiimat3_chunk_idx[TIIMAT3_OMEGA][TIIMAT3_KMAX];
extern size_t tiimat3_chunk_len[TIIMAT3_OMEGA];

void tiimat3_bgv_add(size_t idx, tiimat3_Ciphertext *rop, tiimat3_Ciphertext *op1, tiimat3_Ciphertext *op2);
void tiimat3_bgv_decode(tiimat3_Message *m, const tiimat3_Poly *p);
void tiimat3_bgv_decrypt(size_t idx, tiimat3_Poly *p, const tiimat3_KeySecret *sk, tiimat3_Ciphertext *ct);
void tiimat3_bgv_deinit(void);
void tiimat3_bgv_encode(size_t idx, tiimat3_Poly *p, const tiimat3_Message *m);
void tiimat3_bgv_encrypt(size_t idx, tiimat3_Ciphertext *ct, tiimat3_KeyPublic *pk, const tiimat3_Poly *p, tiimat3_Seed *seed);
void tiimat3_bgv_init(void);
void tiimat3_bgv_keygen_public(size_t idx, tiimat3_KeyPublic *pk, const tiimat3_KeySecret *sk, tiimat3_Seed *seed);
void tiimat3_bgv_keygen_secret(tiimat3_KeySecret *sk);
void tiimat3_bgv_keygen_switch(size_t idx, tiimat3_KeySwitch *ksw, const tiimat3_KeySecret *sk, const tiimat3_KeySecret *sk2, tiimat3_Seed *seed);
void tiimat3_bgv_keygen_switch2(size_t idx, tiimat3_KeySwitch *ksw, const tiimat3_KeySecret *sk, tiimat3_Seed *seed);
void tiimat3_bgv_keygen_switchr(size_t idx, tiimat3_KeySwitch *ksw, const tiimat3_KeySecret *sk, size_t steps, tiimat3_Seed *seed);
void tiimat3_bgv_keyswitch(tiimat3_Ciphertext *ct, tiimat3_KeySwitch *ksw, size_t dsw);
void tiimat3_bgv_keyswitch_delta(size_t idx, tiimat3_Delta *delta, tiimat3_Ciphertext *ct);
void tiimat3_bgv_keyswitch_dot(size_t idx, tiimat3_Ciphertext *ct, tiimat3_CiphertextSwitch *csw, tiimat3_KeySwitch *ksw);
void tiimat3_bgv_keyswitch_ext(tiimat3_CiphertextSwitch *csw, tiimat3_Ciphertext *ct, size_t poly);
void tiimat3_bgv_keyswitch_switch(size_t idx, tiimat3_Ciphertext *ct, tiimat3_Delta *delta);
void tiimat3_bgv_modswitch(size_t idx, tiimat3_Ciphertext *ct, tiimat3_Delta *delta);
void tiimat3_bgv_modswitch_delta(size_t idx, tiimat3_Delta *delta, tiimat3_Ciphertext *ct);
void tiimat3_bgv_modswitch_ext(size_t idx, tiimat3_Ciphertext *ct, tiimat3_Delta *delta);
void tiimat3_bgv_mulc(size_t idx, tiimat3_Ciphertext *rop, tiimat3_Ciphertext *ct, const tiimat3_Poly *p);
void tiimat3_bgv_mul(size_t idx, tiimat3_Ciphertext *rop, tiimat3_Ciphertext *op1, tiimat3_Ciphertext *op2);
void tiimat3_bgv_rot(size_t idx, tiimat3_Ciphertext *rop, tiimat3_Ciphertext *ct, size_t steps);
void tiimat3_bgv_rot_csw(size_t idx, tiimat3_CiphertextSwitch *rop, tiimat3_CiphertextSwitch *csw, size_t steps);
void tiimat3_bgv_rot_csw_inplace(size_t idx, tiimat3_CiphertextSwitch *csw, size_t steps);
void tiimat3_bgv_rot_inplace(size_t idx, tiimat3_Ciphertext *ct, size_t steps);

#endif /* TIIMAT3_H */
