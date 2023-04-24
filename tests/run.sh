#!/bin/sh

export OMP_NUM_THREADS=16

echo "[*] Running BGV tests"
build/bin/test_random
build/bin/test_ring_i64
build/bin/test_msg_mpz
build/bin/test_bgv

for config in 00 01 10 11; do
	echo "[*] Running triples test (${config})"
	build/bin/test_triples_${config}
done
