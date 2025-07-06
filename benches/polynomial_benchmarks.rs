use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use zk_stark::modules::univariate_poly::Polynomial;
use zk_stark::modules::field::FieldElement;
use zk_stark::modules::ntt::NTT;
use num_bigint::BigInt;
use rand::Rng;

fn generate_random_polynomial(degree: usize) -> Polynomial {
    let mut rng = rand::thread_rng();
    let coeffs: Vec<BigInt> = (0..=degree)
        .map(|_| BigInt::from(rng.gen::<u64>()))
        .collect();
    Polynomial::new(coeffs)
}

fn generate_random_field_elements(n: usize) -> Vec<FieldElement> {
    let mut rng = rand::thread_rng();
    (0..n)
        .map(|_| FieldElement::new(BigInt::from(rng.gen::<u64>())))
        .collect()
}

fn bench_polynomial_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("polynomial_operations");
    
    let degrees = vec![16, 64, 256, 1024];
    
    for degree in degrees {
        let poly1 = generate_random_polynomial(degree);
        let poly2 = generate_random_polynomial(degree);
        
        group.bench_with_input(
            BenchmarkId::new("addition", degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    black_box(poly1.clone() + poly2.clone())
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("multiplication", degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    black_box(poly1.clone() * poly2.clone())
                })
            },
        );
        
        if degree <= 256 { // Division is expensive, limit to smaller degrees
            group.bench_with_input(
                BenchmarkId::new("division", degree),
                &degree,
                |b, _| {
                    b.iter(|| {
                        black_box(poly1.clone() / poly2.clone())
                    })
                },
            );
        }
    }
    
    group.finish();
}

fn bench_polynomial_evaluation(c: &mut Criterion) {
    let mut group = c.benchmark_group("polynomial_evaluation");
    
    let degrees = vec![16, 64, 256, 1024];
    
    for degree in degrees {
        let poly = generate_random_polynomial(degree);
        let point = FieldElement::new(BigInt::from(42));
        
        group.bench_with_input(
            BenchmarkId::new("single_evaluation", degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    black_box(poly.eval(black_box(point.clone())))
                })
            },
        );
    }
    
    group.finish();
}

fn bench_polynomial_domain_evaluation(c: &mut Criterion) {
    let mut group = c.benchmark_group("polynomial_domain_evaluation");
    
    let sizes = vec![16, 64, 256, 1024];
    
    for size in sizes {
        let poly = generate_random_polynomial(size);
        let domain = generate_random_field_elements(size);
        
        group.bench_with_input(
            BenchmarkId::new("traditional_evaluation", size),
            &size,
            |b, _| {
                b.iter(|| {
                    let domain_clone = domain.clone();
                    black_box(domain_clone.iter().map(|x| poly.eval(x.clone())).collect::<Vec<_>>())
                })
            },
        );
        
        // For power-of-2 sizes, also benchmark NTT evaluation
        if (size as u32).is_power_of_two() && size >= 8 {
            group.bench_with_input(
                BenchmarkId::new("ntt_evaluation", size),
                &size,
                |b, _| {
                    b.iter(|| {
                        black_box(poly.eval_domain_ntt(size))
                    })
                },
            );
        }
    }
    
    group.finish();
}

fn bench_lagrange_interpolation(c: &mut Criterion) {
    let mut group = c.benchmark_group("lagrange_interpolation");
    group.sample_size(10); // Fewer samples as this can be slow
    
    let sizes = vec![8, 16, 32, 64]; // Smaller sizes as O(n²) is expensive
    
    for size in sizes {
        let domain = generate_random_field_elements(size);
        let values = generate_random_field_elements(size);
        
        group.bench_with_input(
            BenchmarkId::new("lagrange", size),
            &size,
            |b, _| {
                b.iter(|| {
                    black_box(Polynomial::lagrange(
                        black_box(domain.clone()),
                        black_box(values.clone())
                    ))
                })
            },
        );
    }
    
    group.finish();
}

fn bench_ntt_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("ntt_operations");
    
    let sizes = vec![8, 16, 32, 64, 128, 256, 512, 1024];
    
    for size in sizes {
        if !(size as u32).is_power_of_two() {
            continue;
        }
        
        let coeffs = generate_random_field_elements(size);
        let evaluations = generate_random_field_elements(size);
        
        group.bench_with_input(
            BenchmarkId::new("ntt_interpolate", size),
            &size,
            |b, _| {
                let ntt = NTT::new(size);
                b.iter(|| {
                    black_box(ntt.interpolate(black_box(evaluations.clone())))
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("ntt_evaluate", size),
            &size,
            |b, _| {
                let ntt = NTT::new(size);
                b.iter(|| {
                    black_box(ntt.evaluate(black_box(coeffs.clone())))
                })
            },
        );
    }
    
    group.finish();
}

fn bench_polynomial_composition(c: &mut Criterion) {
    let mut group = c.benchmark_group("polynomial_composition");
    
    let degrees = vec![8, 16, 32, 64]; // Composition can be expensive
    
    for degree in degrees {
        let poly1 = generate_random_polynomial(degree);
        let poly2 = generate_random_polynomial(degree / 2); // Smaller inner polynomial
        
        group.bench_with_input(
            BenchmarkId::new("compose", degree),
            &degree,
            |b, _| {
                b.iter(|| {
                    black_box(poly1.compose(&poly2))
                })
            },
        );
    }
    
    group.finish();
}

fn bench_polynomial_zeroifier(c: &mut Criterion) {
    let mut group = c.benchmark_group("polynomial_zeroifier");
    
    let domain_sizes = vec![8, 16, 32, 64, 128];
    
    for size in domain_sizes {
        let domain = generate_random_field_elements(size);
        
        group.bench_with_input(
            BenchmarkId::new("zeroifier_domain", size),
            &size,
            |b, _| {
                b.iter(|| {
                    black_box(Polynomial::zeroifier_domain(black_box(domain.clone())))
                })
            },
        );
    }
    
    group.finish();
}

criterion_group!(
    benches,
    bench_polynomial_operations,
    bench_polynomial_evaluation,
    bench_polynomial_domain_evaluation,
    bench_lagrange_interpolation,
    bench_ntt_operations,
    bench_polynomial_composition,
    bench_polynomial_zeroifier
);
criterion_main!(benches);