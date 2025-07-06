# STARK Proof System Benchmarks

This directory contains comprehensive Criterion benchmarks for the STARK proof system implementation.

## Benchmark Suites

### 1. STARK Benchmarks (`stark_benchmarks.rs`)
- **STARK Prove**: Measures proof generation time across different security configurations
- **STARK Verify**: Measures verification time for pre-generated proofs  
- **End-to-End**: Complete prove + verify cycle benchmarks
- **Rescue Prime**: Hash function and trace generation performance

**Configurations tested:**
- Small: expansion_factor=4, colinearity_tests=2, security_level=2
- Medium: expansion_factor=4, colinearity_tests=4, security_level=16  
- Large: expansion_factor=8, colinearity_tests=8, security_level=32

### 2. Polynomial Benchmarks (`polynomial_benchmarks.rs`)
- **Basic Operations**: Addition, multiplication, division across polynomial degrees
- **Single Evaluation**: Horner's method evaluation performance
- **Domain Evaluation**: Traditional vs NTT-based bulk evaluation
- **Lagrange Interpolation**: O(n²) interpolation performance
- **NTT Operations**: Fast interpolation and evaluation using Number Theoretic Transform
- **Composition**: Polynomial composition p(q(x)) performance
- **Zeroifier**: Vanishing polynomial construction

**Degree ranges tested:** 16, 64, 256, 1024 (varies by operation)

### 3. Field Element Benchmarks (`field_benchmarks.rs`)
- **Arithmetic**: Basic field operations (+, -, *, /, -)
- **Inverse**: Modular inverse computation using extended Euclidean algorithm
- **Power**: Exponentiation with various exponent sizes
- **Sampling**: Field element creation from random bytes
- **Primitive Roots**: nth root of unity computation
- **Serialization**: Bincode and string conversion performance
- **Batch Operations**: Bulk arithmetic operations

### 4. FRI Benchmarks (`fri_benchmarks.rs`)
- **Commit**: Merkle tree commitment phase
- **Prove**: Complete FRI proving including folding rounds
- **Verify**: FRI verification with authentication paths
- **Index Sampling**: Fiat-Shamir challenge generation
- **Domain Evaluation**: Domain point generation
- **Query**: Opening proofs at challenged indices
- **Colinearity Test**: Three-point colinearity verification

## Running Benchmarks

### Run All Benchmarks
```bash
cargo bench
```

### Run Specific Benchmark Suite
```bash
cargo bench stark_benchmarks
cargo bench polynomial_benchmarks
cargo bench field_benchmarks
cargo bench fri_benchmarks
```

### Run Specific Benchmark
```bash
cargo bench stark_prove
cargo bench ntt_operations
cargo bench field_arithmetic
```

### Generate HTML Reports
```bash
cargo bench -- --output-format html
```

Reports will be generated in `target/criterion/` with interactive HTML visualizations.

## Benchmark Configuration

The benchmarks use Criterion's default configuration with some customizations:

- **Sample Size**: Reduced for expensive operations (STARK prove/verify, Lagrange interpolation)
- **Batch Size**: `LargeInput` for operations requiring setup
- **Measurement**: Wall-clock time with statistical analysis
- **Iterations**: Automatically determined by Criterion for statistical significance

## Performance Expectations

### Typical Performance (Release Mode)
- **STARK Prove (Small)**: ~10-15 seconds
- **STARK Verify (Small)**: ~100-500ms  
- **Field Arithmetic**: ~1-10 microseconds
- **Polynomial Operations**: 
  - Addition/Subtraction: ~1-100 microseconds
  - Multiplication: ~10 microseconds - 10ms (degree dependent)
  - NTT Evaluation: ~100 microseconds - 10ms
- **FRI Operations**:
  - Prove: ~1-10 seconds
  - Verify: ~10-100ms

### Performance Factors
- **Polynomial Degree**: Operations scale with degree (linear for evaluation, quadratic for naive multiplication)
- **Domain Size**: NTT operations scale as O(n log n)
- **Security Parameters**: Higher security requires more colinearity tests and larger domains
- **CPU Architecture**: Performance varies significantly between architectures

## Optimization Notes

### Performance Bottlenecks Identified
1. **Lagrange Interpolation**: O(n²) complexity, replaced with NTT for power-of-2 domains
2. **Polynomial Division**: Expensive for high-degree polynomials
3. **Field Inverse**: Uses extended Euclidean algorithm, inherently expensive
4. **Merkle Tree Operations**: Hashing dominates for large trees

### Optimizations Implemented
1. **NTT-Based Operations**: O(n log n) polynomial operations for power-of-2 domains
2. **Coefficient Caching**: Avoid recomputation where possible
3. **Batch Operations**: Reduce per-operation overhead
4. **Memory Management**: Minimize allocations in hot paths

## Interpreting Results

### Key Metrics
- **Mean**: Average execution time
- **Std Dev**: Performance consistency (lower is better)
- **Median**: Typical performance (less affected by outliers)
- **MAD**: Median Absolute Deviation, robust consistency measure

### Regression Detection
Criterion automatically detects performance regressions between runs:
- **Green**: Performance improved or stable
- **Yellow**: Minor performance change
- **Red**: Significant performance regression

### Comparison
Use `cargo bench -- --save-baseline <name>` to save baselines and compare:
```bash
cargo bench -- --save-baseline before_optimization
# Make changes
cargo bench -- --baseline before_optimization
```

## Continuous Integration

For CI environments, consider:
```bash
# Shorter benchmarks for CI
cargo bench -- --quick
# Save baseline on main branch
cargo bench -- --save-baseline main
# Compare on feature branches  
cargo bench -- --baseline main
```

## Hardware Requirements

### Minimum Requirements
- **CPU**: 2+ cores, 2GHz+
- **Memory**: 4GB RAM
- **Storage**: 1GB for benchmark data and reports

### Recommended
- **CPU**: 4+ cores, 3GHz+ with consistent boost clocks
- **Memory**: 8GB+ RAM for large polynomial operations
- **Storage**: SSD for faster I/O during benchmarking

### Notes
- Disable CPU frequency scaling for consistent results
- Close other applications to reduce noise
- Run benchmarks multiple times for statistical significance
- Consider thermal throttling on mobile/laptop hardware