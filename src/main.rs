use zk_stark::modules::merkle::{*};
use zk_stark::modules::field::{*};
use zk_stark::modules::univariate_poly::{*};
use zk_stark::modules::fri::{*};
use zk_stark::modules::proof_stream::{*};

use num_bigint::BigInt;

fn main() {

        // setup
        let degree: usize = 63;
        let expansion_factor: usize = 4;
        let num_colinearity_tests: usize = 17;

        let initial_codeword_length: usize = (degree + 1) * expansion_factor;
        let mut log_codeword_length: usize = 0;
        let mut codeword_length: usize = initial_codeword_length;
        while codeword_length > 1 {
            codeword_length /= 2;
            log_codeword_length += 1;
        }

        assert!(1 << log_codeword_length == initial_codeword_length, "log incorrectly computed");

        let omega = FieldElement::primitive_nth_root(initial_codeword_length as i128);
        let generator = FieldElement::generator();

        // ensure valid properties for field eleemnts
        assert_eq!(omega.pow(1 << log_codeword_length), FieldElement::one(), "omega not the nth root of unity");
        assert!(omega.pow(1 << (log_codeword_length - 1)) != generator, "omega not primitive");

        // setup fri object
        let fri = Fri::new(generator, omega.clone(), initial_codeword_length, expansion_factor, num_colinearity_tests);

        // setup poly
        let mut coeffs: Vec<FieldElement> = vec![];
        for i in 0..degree+1 {
            coeffs.push(FieldElement::new(BigInt::from(i)));
        }
        let polynomial = Polynomial{coeffs};
        println!("polynomial: {}", polynomial);

        // setup domain for evaluation
        let mut domain: Vec<FieldElement> = vec![];
        for i in 0..initial_codeword_length {
            domain.push(omega.pow(i as i128).clone());
        }

        // get codeword
        let codeword = polynomial.eval_domain(domain);

        // setup proof stream
        let mut proof_stream = ProofStream::new();

        // query test
        fri.prove(codeword, &mut proof_stream);

        let mut points: Vec<(usize, FieldElement)> = vec![];
        let verdict: bool = fri.verify(&mut proof_stream, &mut points);
        if verdict {
            println!("FRI test passed");
        } else {
            println!("FRI test failed");
        }
    }
