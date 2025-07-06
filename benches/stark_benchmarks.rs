use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use zk_stark::modules::stark::Stark;
use zk_stark::modules::rescue_prime::RescuePrime;
use zk_stark::modules::field::FieldElement;

fn bench_stark_prove(c: &mut Criterion) {
    let mut group = c.benchmark_group("stark_prove");
    
    // Different configurations to benchmark
    let configs = vec![
        ("small", 4, 2, 2),      // expansion_factor, num_colinearity_tests, security_level
        ("medium", 4, 4, 16),
        ("large", 8, 8, 32),
    ];
    
    for (name, expansion_factor, num_colinearity_tests, security_level) in configs {
        group.bench_with_input(
            BenchmarkId::new("prove", name),
            &(expansion_factor, num_colinearity_tests, security_level),
            |b, &(exp_factor, col_tests, sec_level)| {
                b.iter_batched(
                    || {
                        // Setup - not included in timing
                        let rp = RescuePrime::new();
                        let input_element = FieldElement::sample(vec![48, 120, 100, 101, 97, 100, 98, 101, 101, 102]);
                        let num_cycles = rp.N + 1;
                        let state_width = rp.m;
                        
                        let stark = Stark::new(
                            exp_factor,
                            col_tests,
                            sec_level,
                            state_width,
                            num_cycles
                        );
                        
                        let trace = rp.trace(input_element.clone());
                        let air = rp.transition_constraints(stark.omicron.clone());
                        let output_element = rp.hash(input_element.clone());
                        let boundary = rp.boundary_constraints(output_element);
                        
                        (stark, trace, air, boundary)
                    },
                    |(stark, trace, air, boundary)| {
                        // Actual benchmark - prove
                        black_box(stark.prove(trace, air, boundary))
                    },
                    criterion::BatchSize::LargeInput,
                )
            },
        );
    }
    
    group.finish();
}

fn bench_stark_verify(c: &mut Criterion) {
    let mut group = c.benchmark_group("stark_verify");
    
    // Different configurations to benchmark
    let configs = vec![
        ("small", 4, 2, 2),
        ("medium", 4, 4, 16),
    ];
    
    for (name, expansion_factor, num_colinearity_tests, security_level) in configs {
        group.bench_with_input(
            BenchmarkId::new("verify", name),
            &(expansion_factor, num_colinearity_tests, security_level),
            |b, &(exp_factor, col_tests, sec_level)| {
                // Pre-generate proof for verification benchmark
                let rp = RescuePrime::new();
                let input_element = FieldElement::sample(vec![48, 120, 100, 101, 97, 100, 98, 101, 101, 102]);
                let num_cycles = rp.N + 1;
                let state_width = rp.m;
                
                let stark = Stark::new(
                    exp_factor,
                    col_tests,
                    sec_level,
                    state_width,
                    num_cycles
                );
                
                let trace = rp.trace(input_element.clone());
                let air = rp.transition_constraints(stark.omicron.clone());
                let output_element = rp.hash(input_element.clone());
                let boundary = rp.boundary_constraints(output_element);
                
                let proof = stark.prove(trace, air.clone(), boundary.clone());
                
                b.iter(|| {
                    black_box(stark.verify(proof.clone(), air.clone(), boundary.clone()))
                });
            },
        );
    }
    
    group.finish();
}

fn bench_stark_end_to_end(c: &mut Criterion) {
    let mut group = c.benchmark_group("stark_end_to_end");
    group.sample_size(10); // Fewer samples for long-running benchmarks
    
    let configs = vec![
        ("small", 4, 2, 2),
        ("medium", 4, 4, 16),
    ];
    
    for (name, expansion_factor, num_colinearity_tests, security_level) in configs {
        group.bench_with_input(
            BenchmarkId::new("prove_and_verify", name),
            &(expansion_factor, num_colinearity_tests, security_level),
            |b, &(exp_factor, col_tests, sec_level)| {
                b.iter_batched(
                    || {
                        let rp = RescuePrime::new();
                        let input_element = FieldElement::sample(vec![48, 120, 100, 101, 97, 100, 98, 101, 101, 102]);
                        let num_cycles = rp.N + 1;
                        let state_width = rp.m;
                        
                        let stark = Stark::new(
                            exp_factor,
                            col_tests,
                            sec_level,
                            state_width,
                            num_cycles
                        );
                        
                        let trace = rp.trace(input_element.clone());
                        let air = rp.transition_constraints(stark.omicron.clone());
                        let output_element = rp.hash(input_element.clone());
                        let boundary = rp.boundary_constraints(output_element);
                        
                        (stark, trace, air, boundary)
                    },
                    |(stark, trace, air, boundary)| {
                        let proof = stark.prove(trace, air.clone(), boundary.clone());
                        let verified = stark.verify(proof, air, boundary);
                        black_box(verified)
                    },
                    criterion::BatchSize::LargeInput,
                )
            },
        );
    }
    
    group.finish();
}

fn bench_rescue_prime_hash(c: &mut Criterion) {
    let mut group = c.benchmark_group("rescue_prime");
    
    let rp = RescuePrime::new();
    let input_element = FieldElement::sample(vec![48, 120, 100, 101, 97, 100, 98, 101, 101, 102]);
    
    group.bench_function("hash", |b| {
        b.iter(|| {
            black_box(rp.hash(black_box(input_element.clone())))
        })
    });
    
    group.bench_function("trace_generation", |b| {
        b.iter(|| {
            black_box(rp.trace(black_box(input_element.clone())))
        })
    });
    
    group.finish();
}

criterion_group!(
    benches,
    bench_stark_prove,
    bench_stark_verify,
    bench_stark_end_to_end,
    bench_rescue_prime_hash
);
criterion_main!(benches);