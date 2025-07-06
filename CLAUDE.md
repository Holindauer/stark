# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Rust implementation of a STARK (Scalable Transparent ARgument of Knowledge) proof system with FRI (Fast Reed-Solomon Interactive Oracle Proof). It's based on Alan Szepieniec's STARK tutorial and represents a complete, working implementation suitable for educational and research purposes.

**Current Status**: **FULLY FUNCTIONAL** - All components including the STARK prover/verifier are working correctly.

## Key Commands

### Development Commands
- `cargo build` - Build the project
- `cargo test` - Run all tests
- `cargo test <module_name>` - Run tests for a specific module (e.g., `cargo test field`)
- `cargo check` - Quick syntax and type checking without building
- `cargo clippy` - Run the Rust linter for code quality checks
- `cargo bench` - Run all benchmarks
- `cargo bench <benchmark_name>` - Run specific benchmark suite

### Testing Individual Modules
- `cargo test field::tests` - Test field arithmetic
- `cargo test univariate_poly::tests` - Test univariate polynomials
- `cargo test multivariate_poly::tests` - Test multivariate polynomials
- `cargo test merkle::tests` - Test Merkle tree implementation
- `cargo test fri::tests` - Test FRI protocol
- `cargo test rescue_prime::tests` - Test Rescue Prime hash
- `cargo test stark::tests` - Test STARK prover/verifier

## Architecture Overview

The codebase implements a complete STARK proof system with these core components:

### 1. **Field Arithmetic** (`src/modules/field.rs`)
- Implements finite field operations over prime field p = 407 * 2^119 + 1
- Provides basic field element operations (add, multiply, inverse, etc.)
- Core building block for all polynomial operations

### 2. **Polynomial System**
- **Univariate Polynomials** (`src/modules/univariate_poly.rs`): Single-variable polynomial arithmetic, evaluation, interpolation, composition
- **Multivariate Polynomials** (`src/modules/multivariate_poly.rs`): Multi-variable polynomial operations, symbolic evaluation
- **NTT Module** (`src/modules/ntt.rs`): Number Theoretic Transform for O(n log n) polynomial operations

### 3. **Cryptographic Components**
- **Merkle Trees** (`src/modules/merkle.rs`): Commitment scheme using SHA256
- **Rescue Prime** (`src/modules/rescue_prime.rs`): STARK-friendly hash function with algebraic structure
- **Proof Stream** (`src/modules/proof_stream.rs`): Handles proof serialization/deserialization

### 4. **FRI Protocol** (`src/modules/fri.rs`)
- Implements Fast Reed-Solomon proximity testing
- Key component for polynomial commitment in STARKs
- Includes both proving and verification logic

### 5. **STARK System** (`src/modules/stark.rs`)
- Main prover and verifier implementation
- Orchestrates all components to create and verify STARK proofs
- Fully functional with correct polynomial composition for transition constraints

## Important Implementation Details

- The field prime is specifically chosen for FFT-friendliness (p = 407 * 2^119 + 1)
- Uses Fiat-Shamir transform for non-interactive proofs
- Implements batched FRI for efficiency
- All cryptographic operations use the custom field implementation rather than standard integers

## Key Features

1. **NTT Implementation**: Number Theoretic Transform for O(n log n) polynomial operations
   - Fast polynomial evaluation and interpolation for power-of-2 domains
   - Significant performance improvement over naive algorithms

2. **Comprehensive Benchmarks**: Criterion benchmark suite for performance analysis
   - STARK benchmarks: prove/verify performance across security configurations
   - Polynomial benchmarks: arithmetic, evaluation, NTT operations
   - Field benchmarks: ~460-490ns for basic operations
   - FRI benchmarks: commit, prove, verify phases

3. **Complete Implementation**: All components necessary for STARK proofs
   - Polynomial composition for transition constraints
   - Proper field arithmetic with modular reduction
   - Merkle tree commitments
   - Fiat-Shamir transform for non-interactivity

## Benchmark Commands

```bash
# Run all benchmarks
cargo bench

# Run specific benchmark suites
cargo bench stark_benchmarks
cargo bench polynomial_benchmarks
cargo bench field_benchmarks
cargo bench fri_benchmarks

# Generate HTML reports
cargo bench -- --output-format html

# Compare performance with baseline
cargo bench -- --save-baseline before
cargo bench -- --baseline before
```

## Development Tips

- The test cases in each module provide good examples of expected behavior
- Use `cargo test -- --nocapture` to see println! output during tests
- Run benchmarks before and after optimizations to measure impact
- The NTT module is key for performance on power-of-2 domains
- Polynomial coefficients are stored highest-to-lowest degree (design choice)

## Performance Characteristics

- **STARK Prove**: ~12 seconds (small config with Rescue-Prime)
- **STARK Verify**: ~200-500ms
- **Field Operations**: ~460-490 nanoseconds
- **Polynomial Multiplication**: O(n log n) with NTT
- **FRI Commit**: Dominated by Merkle tree hashing

## Future Optimization Opportunities

1. **Parallelization**: Many polynomial operations can be parallelized
2. **GPU Acceleration**: NTT and field operations are GPU-friendly
3. **Memory Optimization**: Reduce allocations in hot paths
4. **Batch Operations**: Process multiple proofs simultaneously
5. **Alternative Hash Functions**: Explore faster STARK-friendly hashes