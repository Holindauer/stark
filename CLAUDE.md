# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Rust implementation of a STARK (Scalable Transparent ARgument of Knowledge) proof system with FRI (Fast Reed-Solomon Interactive Oracle Proof). It's based on Alan Szepieniec's STARK tutorial and is intended for educational/research purposes rather than production use.

**Current Status**: There's a known bug in the STARK prover or verifier causing proof verification to fail, though all component tests pass.

## Key Commands

### Development Commands
- `cargo build` - Build the project
- `cargo test` - Run all tests
- `cargo test <module_name>` - Run tests for a specific module (e.g., `cargo test field`)
- `cargo check` - Quick syntax and type checking without building
- `cargo clippy` - Run the Rust linter for code quality checks

### Testing Individual Modules
- `cargo test field::tests` - Test field arithmetic
- `cargo test univariate_poly::tests` - Test univariate polynomials
- `cargo test multivariate_poly::tests` - Test multivariate polynomials
- `cargo test merkle::tests` - Test Merkle tree implementation
- `cargo test fri::tests` - Test FRI protocol
- `cargo test rescue_prime::tests` - Test Rescue Prime hash
- `cargo test stark::tests` - Test STARK prover/verifier (currently failing)

## Architecture Overview

The codebase implements a complete STARK proof system with these core components:

### 1. **Field Arithmetic** (`src/modules/field.rs`)
- Implements finite field operations over prime field p = 407 * 2^119 + 1
- Provides basic field element operations (add, multiply, inverse, etc.)
- Core building block for all polynomial operations

### 2. **Polynomial System**
- **Univariate Polynomials** (`src/modules/univariate_poly.rs`): Single-variable polynomial arithmetic, evaluation, interpolation
- **Multivariate Polynomials** (`src/modules/multivariate_poly.rs`): Multi-variable polynomial operations, symbolic evaluation

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
- Contains the bug that needs fixing

## Important Implementation Details

- The field prime is specifically chosen for FFT-friendliness (p = 407 * 2^119 + 1)
- Uses Fiat-Shamir transform for non-interactive proofs
- Implements batched FRI for efficiency
- All cryptographic operations use the custom field implementation rather than standard integers

## Known Issues

1. **STARK Verification Bug**: The complete STARK proof verification fails despite all component tests passing
2. **Empty main.rs**: The binary entry point is not implemented (library-only usage)
3. **Unused imports**: Some modules have warnings about unused imports that could be cleaned up

## Development Tips

- When debugging the STARK verification issue, focus on the interaction between `stark.rs` and the FRI protocol
- The test cases in each module provide good examples of expected behavior
- Use `cargo test -- --nocapture` to see println! output during tests