#ifndef TIIMAT3_MSG_H
#define TIIMAT3_MSG_H

#include <assert.h>
#include <gmp.h>
#include <stdio.h>
#include <x86intrin.h>

#include "util.h"

#define TIIMAT3_D 32768
#define TIIMAT3_T "340282366920938463463374607431773454337"

typedef struct tiimat3_message     tiimat3_Message;
typedef struct tiimat3_message_mod tiimat3_MessageMod;

struct tiimat3_message {
	mpz_t value[TIIMAT3_D];
};

struct tiimat3_message_mod {
	mpz_t value;
	mpz_t root;
};

extern tiimat3_MessageMod tiimat3_t;

tiimat3_Message *tiimat3_alloc_msg(size_t len);
void tiimat3_dealloc_msg(tiimat3_Message *m, size_t len);
void tiimat3_msg_cmod(tiimat3_Message *rop, const tiimat3_Message *m, const mpz_t mod);
void tiimat3_msg_crt_i64(tiimat3_Message *m, const int64_t **rns, const int64_t *mods, size_t len);
void tiimat3_msg_crt_u64(tiimat3_Message *m, const uint64_t **rns, const uint64_t *mods, size_t len);
void tiimat3_msg_mulc(tiimat3_Message *rop, const tiimat3_Message *m, uint64_t c);
void tiimat3_msg_pack(tiimat3_Message *rop, const tiimat3_Message *m);
void tiimat3_msg_unpack(tiimat3_Message *rop, const tiimat3_Message *m);

void tiimat3_msgmod_deinit(void);
void tiimat3_msgmod_init(void);

#endif /* TIIMAT3_MSG_H */
