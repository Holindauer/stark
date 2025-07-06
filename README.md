# Scalable Transparent Arguments of Knowledge

Rust implementation of STARKs (Scalable Transparent ARguments of Knowledge) with Fast Reed-Solomon Interactive Oracle Proofs (FRI). Based on [Alan Szepieniec's excellent STARK Anatomy tutorial](https://aszepieniec.github.io/stark-anatomy/).

## Overview

A fully functional STARK proof system that generates and verifies proofs for computations. The implementation uses the Rescue-Prime hash function as an example computation and includes all necessary components for a complete STARK proof system.

## Quick Start

```bash
# Run all tests
cargo test

# Run benchmarks
cargo bench

# Run specific benchmark suite
cargo bench stark_benchmarks
cargo bench polynomial_benchmarks
```
### Core Components

1. **Field Module** (`field.rs`)
   - Prime field arithmetic over 2^128 - 45 * 2^40 + 1
   - Optimized modular operations
   - Primitive root computation

2. **Polynomial Module** (`univariate_poly.rs`)
   - Polynomial arithmetic and evaluation
   - NTT-based fast multiplication and evaluation
   - Polynomial composition for transition constraints

3. **NTT Module** (`ntt.rs`)
   - Number Theoretic Transform implementation
   - O(n log n) interpolation and evaluation
   - Cooley-Tukey algorithm with bit-reversal

4. **FRI Module** (`fri.rs`)
   - Fast Reed-Solomon IOP implementation
   - Folding and query phases
   - Merkle tree commitments

5. **STARK Module** (`stark.rs`)
   - Complete STARK prover and verifier
   - Rescue-Prime hash function as example computation
   - Fiat-Shamir transform for non-interactivity

## Testing

```bash
# Run all tests
cargo test

# Run specific test suite
cargo test field
cargo test polynomial
cargo test fri
cargo test stark

# Run with output
cargo test -- --nocapture
```

## Benchmarking

The project includes comprehensive Criterion benchmarks:

```bash
# Run all benchmarks
cargo bench

# Generate HTML reports
cargo bench -- --output-format html

# Compare with baseline
cargo bench -- --save-baseline main
cargo bench -- --baseline main
```

Benchmark reports are generated in `target/criterion/`.

## Development

### Building
```bash
cargo build --release
```

### Documentation
```bash
cargo doc --open
```

### Code Coverage
```bash
cargo tarpaulin -o html
```

## Learning Resources

- [STARK Anatomy Tutorial](https://aszepieniec.github.io/stark-anatomy/) - Original tutorial this implementation is based on
- [STARK Paper](https://eprint.iacr.org/2018/046.pdf) - Original STARK paper by Ben-Sasson et al.
- [FRI Paper](https://eccc.weizmann.ac.il/report/2017/134/) - Fast Reed-Solomon IOP paper

## Acknowledgments

- Alan Szepieniec for the excellent STARK Anatomy tutorial
- The StarkWare team for pioneering STARK technology
- Claude (Anthropic) for assistance in debugging and optimization
