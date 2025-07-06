use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use zk_stark::modules::fri::Fri;
use zk_stark::modules::field::FieldElement;
use zk_stark::modules::univariate_poly::Polynomial;
use zk_stark::modules::proof_stream::ProofStream;
use num_bigint::BigInt;
use rand::Rng;

fn generate_test_polynomial(degree: usize) -> Polynomial {
    let mut rng = rand::thread_rng();
    let coeffs: Vec<BigInt> = (0..=degree)
        .map(|_| BigInt::from(rng.gen::<u32>()))
        .collect();
    Polynomial::new(coeffs)
}

fn setup_fri_instance(domain_length: usize, expansion_factor: usize, num_colinearity_tests: usize) -> (Fri, Polynomial, Vec<FieldElement>) {
    // Setup FRI parameters
    let generator = FieldElement::generator();
    let omega = FieldElement::primitive_nth_root(domain_length as u128);
    
    let fri = Fri::new(
        generator,
        omega.clone(),
        domain_length,
        expansion_factor,
        num_colinearity_tests
    );
    
    // Create a test polynomial
    let degree = (domain_length / expansion_factor) - 1;
    let polynomial = generate_test_polynomial(degree);
    
    // Generate domain
    let mut domain = Vec::new();
    for i in 0..domain_length {
        domain.push(omega.pow(i as u128));
    }
    
    (fri, polynomial, domain)
}

fn bench_fri_commit(c: &mut Criterion) {
    let mut group = c.benchmark_group("fri_commit");
    
    let configs = vec![
        ("small", 64, 4, 2),
        ("medium", 256, 4, 4),
        ("large", 1024, 8, 8),
    ];
    
    for (name, domain_length, expansion_factor, num_colinearity_tests) in configs {
        group.bench_with_input(
            BenchmarkId::new("commit", name),
            &(domain_length, expansion_factor, num_colinearity_tests),
            |b, &(domain_len, exp_factor, col_tests)| {
                b.iter_batched(
                    || {
                        let (fri, polynomial, domain) = setup_fri_instance(domain_len, exp_factor, col_tests);
                        let codeword = polynomial.eval_domain(domain);
                        let mut proof_stream = ProofStream::new();
                        (fri, codeword, proof_stream)
                    },
                    |(fri, codeword, mut proof_stream)| {
                        black_box(fri.commit(codeword, &mut proof_stream))
                    },
                    criterion::BatchSize::LargeInput,
                )
            },
        );
    }
    
    group.finish();
}

fn bench_fri_prove(c: &mut Criterion) {
    let mut group = c.benchmark_group("fri_prove");
    group.sample_size(10); // Fewer samples for expensive operations
    
    let configs = vec![
        ("small", 64, 4, 2),
        ("medium", 256, 4, 4),
        ("large", 1024, 8, 8),
    ];
    
    for (name, domain_length, expansion_factor, num_colinearity_tests) in configs {
        group.bench_with_input(
            BenchmarkId::new("prove", name),
            &(domain_length, expansion_factor, num_colinearity_tests),
            |b, &(domain_len, exp_factor, col_tests)| {
                b.iter_batched(
                    || {
                        let (fri, polynomial, domain) = setup_fri_instance(domain_len, exp_factor, col_tests);
                        let codeword = polynomial.eval_domain(domain);
                        let mut proof_stream = ProofStream::new();
                        (fri, codeword, proof_stream)
                    },
                    |(fri, codeword, mut proof_stream)| {
                        black_box(fri.prove(codeword, &mut proof_stream))
                    },
                    criterion::BatchSize::LargeInput,
                )
            },
        );
    }
    
    group.finish();
}

fn bench_fri_verify(c: &mut Criterion) {
    let mut group = c.benchmark_group("fri_verify");
    
    let configs = vec![
        ("small", 64, 4, 2),
        ("medium", 256, 4, 4),
    ];
    
    for (name, domain_length, expansion_factor, num_colinearity_tests) in configs {
        group.bench_with_input(
            BenchmarkId::new("verify", name),
            &(domain_length, expansion_factor, num_colinearity_tests),
            |b, &(domain_len, exp_factor, col_tests)| {
                // Pre-generate proof for verification
                let (fri, polynomial, domain) = setup_fri_instance(domain_len, exp_factor, col_tests);
                let codeword = polynomial.eval_domain(domain);
                let mut proof_stream = ProofStream::new();
                fri.prove(codeword, &mut proof_stream);
                
                b.iter(|| {
                    let mut proof_stream_copy = proof_stream.clone();
                    let mut polynomial_values = Vec::new();
                    black_box(fri.verify(&mut proof_stream_copy, &mut polynomial_values))
                });
            },
        );
    }
    
    group.finish();
}

fn bench_fri_sample_indices(c: &mut Criterion) {
    let mut group = c.benchmark_group("fri_sample_indices");
    
    let fri = {
        let generator = FieldElement::generator();
        let omega = FieldElement::primitive_nth_root(256);
        Fri::new(generator, omega, 256, 4, 8)
    };
    
    let seed = vec![1, 2, 3, 4, 5, 6, 7, 8, 9, 10];
    let sizes = vec![128, 256, 512, 1024];
    let reduced_sizes = vec![32, 64, 128, 256];
    let numbers = vec![2, 4, 8, 16];
    
    for ((size, reduced_size), number) in sizes.into_iter().zip(reduced_sizes).zip(numbers) {
        group.bench_with_input(
            BenchmarkId::new("sample_indices", format!("size_{}_reduced_{}_num_{}", size, reduced_size, number)),
            &(size, reduced_size, number),
            |b, &(s, rs, n)| {
                b.iter(|| {
                    black_box(fri.sample_indices(&seed, s, rs, n))
                })
            },
        );
    }
    
    group.finish();
}

fn bench_fri_eval_domain(c: &mut Criterion) {
    let mut group = c.benchmark_group("fri_eval_domain");
    
    let domain_lengths = vec![64, 256, 1024, 4096];
    
    for domain_length in domain_lengths {
        let generator = FieldElement::generator();
        let omega = FieldElement::primitive_nth_root(domain_length as u128);
        let fri = Fri::new(generator, omega, domain_length, 4, 2);
        
        group.bench_with_input(
            BenchmarkId::new("eval_domain", domain_length),
            &domain_length,
            |b, _| {
                b.iter(|| {
                    black_box(fri.eval_domain())
                })
            },
        );
    }
    
    group.finish();
}

fn bench_fri_query(c: &mut Criterion) {
    let mut group = c.benchmark_group("fri_query");
    
    let configs = vec![
        ("small", 64, 4, 2),
        ("medium", 256, 4, 4),
    ];
    
    for (name, domain_length, expansion_factor, num_colinearity_tests) in configs {
        group.bench_with_input(
            BenchmarkId::new("query", name),
            &(domain_length, expansion_factor, num_colinearity_tests),
            |b, &(domain_len, exp_factor, col_tests)| {
                // Setup for query benchmark
                let (fri, polynomial, domain) = setup_fri_instance(domain_len, exp_factor, col_tests);
                let current_codeword = polynomial.eval_domain(domain);
                
                // Create a "next" codeword (simplified for benchmarking)
                let next_codeword: Vec<FieldElement> = current_codeword.iter()
                    .take(current_codeword.len() / 2)
                    .cloned()
                    .collect();
                
                let c_indices = vec![0, 1]; // Simple indices for testing
                
                b.iter(|| {
                    let mut proof_stream = ProofStream::new();
                    black_box(fri.query(
                        current_codeword.clone(),
                        next_codeword.clone(),
                        c_indices.clone(),
                        &mut proof_stream
                    ))
                });
            },
        );
    }
    
    group.finish();
}

fn bench_fri_colinearity_test(c: &mut Criterion) {
    let mut group = c.benchmark_group("fri_colinearity_test");
    
    // Test different point configurations
    let point_sets = vec![
        vec![
            (BigInt::from(1), BigInt::from(1)),
            (BigInt::from(2), BigInt::from(2)),
            (BigInt::from(3), BigInt::from(3))
        ],
        vec![
            (BigInt::from(0), BigInt::from(5)),
            (BigInt::from(1), BigInt::from(7)),
            (BigInt::from(2), BigInt::from(9))
        ],
        vec![
            (BigInt::from(10), BigInt::from(100)),
            (BigInt::from(20), BigInt::from(200)),
            (BigInt::from(30), BigInt::from(300))
        ],
    ];
    
    for (i, points) in point_sets.iter().enumerate() {
        group.bench_with_input(
            BenchmarkId::new("colinearity_test", format!("set_{}", i)),
            points,
            |b, points| {
                b.iter(|| {
                    black_box(Polynomial::test_colinearity(black_box(points.clone())))
                })
            },
        );
    }
    
    group.finish();
}

criterion_group!(
    benches,
    bench_fri_commit,
    bench_fri_prove,
    bench_fri_verify,
    bench_fri_sample_indices,
    bench_fri_eval_domain,
    bench_fri_query,
    bench_fri_colinearity_test
);
criterion_main!(benches);