# Matrix Triple Generation

To build this code, use the following commands:
```
cmake -S . -B build
cmake --build build
```

There also is a helper script running all tests:
```
tests/run.sh
```

Finally, if you want to run benchmarks, please specify
`-DTIIMAT3_BUILD_BENCH=1` when running the first `cmake` command.


## Tests

We generate multiple tests for the BGV library and the triples implementation.
For the BGV code, we individually test the random number implementation, the
signed ring arithmetic and the mpz plaintext implementation. We also provide
integration tests for the BGV scheme itself that show how to use the library.

For the triples, we provide multiple test files for the different combinations
of hoisting and pre-rotation:
- `test_triples_00`: no hoisting, no pre-rotation
- `test_triples_01`: no hoisting, pre-rotation
- `test_triples_10`: hoisting, no pre-rotation
- `test_triples_11`: hoisting, pre-rotation

## Citations

The accompanying paper is published at [AsiaCCS 2023] and available on [ePrint].

[AsiaCCS 2023]: https://doi.org/10.1145/3579856.3590344
[ePrint]: https://eprint.iacr.org/2023/593
