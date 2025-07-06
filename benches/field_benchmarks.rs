use criterion::{black_box, criterion_group, criterion_main, Criterion, BenchmarkId};
use zk_stark::modules::field::FieldElement;
use num_bigint::BigInt;
use rand::Rng;

fn generate_random_field_element() -> FieldElement {
    let mut rng = rand::thread_rng();
    FieldElement::new(BigInt::from(rng.gen::<u64>()))
}

fn bench_field_arithmetic(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_arithmetic");
    
    let a = generate_random_field_element();
    let b = generate_random_field_element();
    
    group.bench_function("addition", |bench| {
        bench.iter(|| {
            black_box(a.clone() + b.clone())
        })
    });
    
    group.bench_function("subtraction", |bench| {
        bench.iter(|| {
            black_box(a.clone() - b.clone())
        })
    });
    
    group.bench_function("multiplication", |bench| {
        bench.iter(|| {
            black_box(a.clone() * b.clone())
        })
    });
    
    group.bench_function("division", |bench| {
        bench.iter(|| {
            black_box(a.clone() / b.clone())
        })
    });
    
    group.bench_function("negation", |bench| {
        bench.iter(|| {
            black_box(-a.clone())
        })
    });
    
    group.finish();
}

fn bench_field_inverse(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_inverse");
    
    let elements = (0..100).map(|_| generate_random_field_element()).collect::<Vec<_>>();
    
    group.bench_function("inverse", |b| {
        let mut idx = 0;
        b.iter(|| {
            let result = black_box(elements[idx % elements.len()].inverse());
            idx += 1;
            result
        })
    });
    
    group.finish();
}

fn bench_field_power(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_power");
    
    let base = generate_random_field_element();
    let exponents = vec![2u128, 4, 8, 16, 32, 64, 128, 256, 1024];
    
    for exp in exponents {
        group.bench_with_input(
            BenchmarkId::new("pow", exp),
            &exp,
            |b, &exponent| {
                b.iter(|| {
                    black_box(base.pow(black_box(exponent)))
                })
            },
        );
    }
    
    group.finish();
}

fn bench_field_sample(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_sample");
    
    let random_bytes = vec![
        vec![1, 2, 3, 4, 5, 6, 7, 8],
        vec![9, 10, 11, 12, 13, 14, 15, 16, 17],
        vec![48, 120, 100, 101, 97, 100, 98, 101, 101, 102],
    ];
    
    for (_i, bytes) in random_bytes.iter().enumerate() {
        group.bench_with_input(
            BenchmarkId::new("sample", format!("{}_bytes", bytes.len())),
            bytes,
            |b, bytes| {
                b.iter(|| {
                    black_box(FieldElement::sample(black_box(bytes.clone())))
                })
            },
        );
    }
    
    group.finish();
}

fn bench_field_primitive_roots(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_primitive_roots");
    
    let orders = vec![4u128, 8, 16, 32, 64, 128, 256, 512, 1024];
    
    for order in orders {
        group.bench_with_input(
            BenchmarkId::new("primitive_nth_root", order),
            &order,
            |b, &n| {
                b.iter(|| {
                    black_box(FieldElement::primitive_nth_root(black_box(n)))
                })
            },
        );
    }
    
    group.finish();
}

fn bench_field_creation(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_creation");
    
    let big_ints = vec![
        BigInt::from(42),
        BigInt::from(1234567890u64),
        BigInt::parse_bytes(b"123456789012345678901234567890", 10).unwrap(),
    ];
    
    for (i, big_int) in big_ints.iter().enumerate() {
        group.bench_with_input(
            BenchmarkId::new("from_bigint", format!("size_{}", i)),
            big_int,
            |b, big_int| {
                b.iter(|| {
                    black_box(FieldElement::new(black_box(big_int.clone())))
                })
            },
        );
    }
    
    group.finish();
}

fn bench_field_comparison(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_comparison");
    
    let elements = (0..100).map(|_| generate_random_field_element()).collect::<Vec<_>>();
    
    group.bench_function("equality", |b| {
        let mut idx = 0;
        b.iter(|| {
            let a = &elements[idx % elements.len()];
            let b = &elements[(idx + 1) % elements.len()];
            idx += 1;
            black_box(a == b)
        })
    });
    
    group.finish();
}

fn bench_field_serialization(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_serialization");
    
    let element = generate_random_field_element();
    
    group.bench_function("to_string", |b| {
        b.iter(|| {
            black_box(element.to_string())
        })
    });
    
    // Test bincode serialization if available
    group.bench_function("bincode_serialize", |b| {
        b.iter(|| {
            black_box(bincode::serialize(&element).unwrap())
        })
    });
    
    let serialized = bincode::serialize(&element).unwrap();
    group.bench_function("bincode_deserialize", |b| {
        b.iter(|| {
            black_box(bincode::deserialize::<FieldElement>(&serialized).unwrap())
        })
    });
    
    group.finish();
}

fn bench_field_batch_operations(c: &mut Criterion) {
    let mut group = c.benchmark_group("field_batch_operations");
    
    let sizes = vec![10, 100, 1000];
    
    for size in sizes {
        let elements = (0..size).map(|_| generate_random_field_element()).collect::<Vec<_>>();
        
        group.bench_with_input(
            BenchmarkId::new("batch_addition", size),
            &size,
            |b, _| {
                b.iter(|| {
                    let mut sum = FieldElement::zero();
                    for elem in &elements {
                        sum = sum + elem.clone();
                    }
                    black_box(sum)
                })
            },
        );
        
        group.bench_with_input(
            BenchmarkId::new("batch_multiplication", size),
            &size,
            |b, _| {
                b.iter(|| {
                    let mut product = FieldElement::one();
                    for elem in &elements {
                        product = product * elem.clone();
                    }
                    black_box(product)
                })
            },
        );
    }
    
    group.finish();
}

criterion_group!(
    benches,
    bench_field_arithmetic,
    bench_field_inverse,
    bench_field_power,
    bench_field_sample,
    bench_field_primitive_roots,
    bench_field_creation,
    bench_field_comparison,
    bench_field_serialization,
    bench_field_batch_operations
);
criterion_main!(benches);