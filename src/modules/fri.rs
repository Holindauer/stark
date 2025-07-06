use crate::modules::field::{*};
use crate::modules::merkle::{*};
use crate::modules::proof_stream::{*};
use crate::modules::univariate_poly::{*};
use blake2::{Blake2b, Digest};
use typenum::U64;
use num_bigint::BigInt;
extern crate num_traits;
extern crate num_bigint;
use hex;
use num_traits::Num; 
use generic_array::typenum::U32;
use generic_array::GenericArray;

// Merkle
type OutputSize = U32; // 32 bytes (256 bits) output
type HashOutput = GenericArray<u8, OutputSize>; // Fixed-size hash output


pub struct Fri {
    pub offset: FieldElement,
    pub omega: FieldElement,
    pub domain_length: usize,
    pub expansion_factor: usize,
    pub num_colinearity_tests: usize,
}

impl Fri {
    
    // constructor 
    pub fn new( offset: FieldElement, omega: FieldElement, initial_domain_length: usize, expansion_factor: usize, num_colinearity_tests: usize ) -> Fri {        

        let fri_obj = Fri { offset, omega, domain_length: initial_domain_length, expansion_factor, num_colinearity_tests };
        assert!(fri_obj.num_rounds() > 0, "Cannot do FRI with less than 1 round");
        fri_obj
    }

    // computes num rounds of the fri protocol
    pub fn num_rounds(&self) -> usize {

        let mut codeword_length = self.domain_length.clone();
        let mut num_rounds = 0; 

        // while codeword length is greater than the expansion factor(reed-soloman) and the num
        // colinearity checks does not exceed a quarter of the length of the working codeword
        while codeword_length > self.expansion_factor && 4 * self.num_colinearity_tests.clone() < codeword_length {
            codeword_length /= 2;
            num_rounds += 1;
        }
        num_rounds
    }

    // computes the indices to sample for the fri protocol
    pub fn sample_indices(&self, seed: &Vec<u8>, size: usize, reduced_size: usize, number: usize) -> Vec<usize> {
        let mut indices = Vec::new();
        let mut reduced_indices = Vec::new();
        let mut counter: u64 = 0;
        
        while indices.len() < number {
            let mut hasher = Blake2b::new();
            hasher.update(seed);
            hasher.update(&counter.to_le_bytes());
            
            // Explicitly typing the hash_result with the expected output size
            let hash_result: GenericArray<u8, U64> = hasher.finalize();
    
            let index = Fri::sample_index(hash_result.to_vec(), size);
            let reduced_index = index % reduced_size;
            counter += 1;
    
            if !reduced_indices.contains(&reduced_index) {
                indices.push(index);
                reduced_indices.push(reduced_index);
            }
        }
    
        indices
    }

    // helper function to sample index from byte array of random bytes
    pub fn sample_index(byte_array: Vec<u8>, size: usize) -> usize {
        let mut acc: usize = 0;
        for b in byte_array {
            acc = (acc << 8) ^ (b as usize);
        }
        acc % size
    }

    // computes domain for polynomial evaluation
    pub fn eval_domain(&self) -> Vec<FieldElement> {

        let mut domain: Vec<FieldElement> = Vec::new();
        for i in 0..self.domain_length {
            domain.push( self.offset.clone() * self.omega.pow(i as u128));
        }
        domain
    }

    // commit 
    pub fn commit(&self, mut codeword: Vec<FieldElement>, proof_stream: &mut ProofStream) -> Vec<Vec<FieldElement>> {
        
        // shorthand
        let one = FieldElement::one();
        let two = FieldElement::new(BigInt::from(2));

        // working values
        let mut omega = self.omega.clone();
        let mut offset = self.offset.clone();
        let mut codewords = Vec::new();

        // split and fold rounds
        for r in 0..self.num_rounds() {
            let N: usize = codeword.len();

            // ensure n has the right order
            assert!(omega.pow(N as u128 - 1) == omega.inverse(), "error in commit: omega does not have the right order!");

            // Commit merkle root to proof stream
            let root: HashOutput = Merkle::commit(&codeword.iter().map(|x| bincode::serialize(x).unwrap()).collect::<Vec<Vec<u8>>>());            
            proof_stream.push(hex::encode(root));   

            // prepare next round only if necessary
            if r == self.num_rounds() - 1 { break; }

            // get challenge
            let fiat_shamir: Vec<u8> = proof_stream.prover_fiat_shamir(32); // 32 bytes
            let alpha: FieldElement = FieldElement::sample(fiat_shamir);            

            // collect codewords
            codewords.push(codeword.clone());

            // split and fold
            let mut new_codeword: Vec<FieldElement> = vec![];
            for i in 0..(N/2) {
                new_codeword.push(
                    two.inverse() * ( (one.clone() + alpha.clone() / (offset.clone() * omega.pow(i as u128))) * codeword.get(i).unwrap().clone() + (one.clone() - alpha.clone() / (offset.clone() * omega.pow(i as u128)))  * codeword.get(N/2 + i).unwrap().clone()) // long boi
                );
            }
            codeword = new_codeword;

            omega = omega.pow(2);
            offset = offset.pow(2);
        }

        // send last codeword as json
        proof_stream.push(serde_json::to_string(&codeword).unwrap());

        // collect last codeword too
        codewords.push(codeword.clone());

        codewords
    }

    // query
    pub fn query(&self, current_codeword: Vec<FieldElement>, next_codeword: Vec<FieldElement>, c_indices: Vec<usize>, proof_stream: &mut ProofStream) -> Vec<usize> {

        // infer a and b indices
        let mut a_indices: Vec<usize> = Vec::new();
        let mut b_indices: Vec<usize> = Vec::new();
        for idx in c_indices.clone() {
            a_indices.push(idx);
            b_indices.push(idx + current_codeword.len()/2);
        }

        // reveal leaves
        for s in 0..self.num_colinearity_tests {

            // push points to proof stream
            proof_stream.push(
                serde_json::to_string(
                    // tuple of three colinear points to check
                    &(
                        current_codeword.get(a_indices[s]).unwrap().to_string(),
                        current_codeword.get(b_indices[s]).unwrap().to_string(),
                        next_codeword.get(c_indices[s]).unwrap().to_string()
                    )
                ).unwrap()
            );
        }

        // reveal authentication path
        for s in 0..self.num_colinearity_tests {

            // a indices
            let serialized_codeword: Vec<Vec<u8>> = current_codeword.iter() // convert to bytes
            .map(|element| bincode::serialize(element).unwrap()).collect();
            let open_a: Vec<HashOutput> = Merkle::open(a_indices[s], &serialized_codeword);
            let serialized_open_a = serde_json::to_string(&open_a).unwrap();
            proof_stream.push(serialized_open_a);

            // b indices
            let serialized_codeword: Vec<Vec<u8>> = current_codeword.iter() // convert to bytes
            .map(|element| bincode::serialize(element).unwrap()).collect();
            let open_b: Vec<HashOutput> = Merkle::open(b_indices[s], &serialized_codeword);
            let serialized_open_b = serde_json::to_string(&open_b).unwrap();
            proof_stream.push(serialized_open_b);

            // c indices
            let serialized_codeword: Vec<Vec<u8>> = next_codeword.iter() // convert to bytes
            .map(|element| bincode::serialize(element).unwrap()).collect();
            let open_c: Vec<HashOutput> = Merkle::open(c_indices[s], &serialized_codeword);
            let serialized_open_c = serde_json::to_string(&open_c).unwrap();
            proof_stream.push(serialized_open_c);
        }

        // return a + b indices
        a_indices.extend(b_indices);
        a_indices
    }

    // prove
    pub fn prove(&self, codeword: Vec<FieldElement>, proof_stream: &mut ProofStream) -> Vec<usize> {

        // commit
        let codewords: Vec<Vec<FieldElement>> = self.commit(codeword, proof_stream);

        // get indices
        let top_level_indices = self.sample_indices(
            &proof_stream.prover_fiat_shamir(32),            // seed
            codewords.get(0).unwrap().len() / 2,             // size  
            codewords.get(codewords.len()-1).unwrap().len(), // reduced size
            self.num_colinearity_tests                       // number of tests
        );
        let mut indices = top_level_indices.clone(); // working indices


        // query phase
        for i in 0..codewords.len()-1 {
            
            // fold
            let mut folded_indices = vec![];
            for idx in indices.clone() { folded_indices.push( idx % (codewords.get(i).unwrap().len()/2) );} 
            indices = folded_indices;
            
            // query
            indices = self.query(
                codewords.get(i).unwrap().clone(),
                codewords.get(i+1).unwrap().clone(),
                indices,
                proof_stream
            );
        }

        // return the indices from the first round of queries
        // We need to reconstruct them from the top level since indices were folded
        let mut a_indices: Vec<usize> = Vec::new();
        let mut b_indices: Vec<usize> = Vec::new();
        for idx in top_level_indices {
            a_indices.push(idx);
            b_indices.push(idx + codewords.get(0).unwrap().len()/2);
        }
        a_indices.extend(b_indices);
        a_indices
    }

    pub fn verify(&self, proof_stream: &mut ProofStream, polynomial_values: &mut Vec<(usize, FieldElement)>) -> bool {

        // short hand
        let mut omega: FieldElement = self.omega.clone();
        let mut offset: FieldElement = self.offset.clone();

        // extract merkle roots and alphas
        let mut roots: Vec<String> = Vec::new();
        let mut alphas: Vec<FieldElement> = Vec::new();
        for r in 0..self.num_rounds(){
            roots.push(proof_stream.pull());
            alphas.push(FieldElement::sample(proof_stream.verifier_fiat_shamir(32)))
        }

        // extract last codeword
        let last_codeword: Vec<FieldElement> = serde_json::from_str(&proof_stream.pull()).unwrap();

        // compute last codeword's merkle root
        let serialized_codeword_commit: Vec<Vec<u8>> = last_codeword.iter() // to bytes
        .map(|x| bincode::serialize(x).unwrap()).collect::<Vec<Vec<u8>>>();
        let root: String = hex::encode(Merkle::commit(&serialized_codeword_commit));

        // verify the last codeword's root is the same as the last root
        if *roots.get(roots.len()-1).unwrap() != root {
            println!("last codeword is not well formed");
            return false;
        }
        
        // check that it is low degree
        let degree = (last_codeword.len() / self.expansion_factor) - 1;
        let mut last_omega: FieldElement = omega.clone();
        let mut last_offset: FieldElement = offset.clone();
        for r in 0..self.num_rounds()-1 {
            last_omega = last_omega.pow(2);
            last_offset = last_offset.pow(2);
        }

        // assert that last_omega has the right order
        assert!(last_omega.inverse() == last_omega.pow((last_codeword.len()-1) as u128));

        // compute interpolant
        let mut last_domain: Vec<FieldElement> = vec![];
        for i in 0..last_codeword.len() {
            last_domain.push(last_offset.clone() * last_omega.pow(i as u128));
        }
        let poly = Polynomial::lagrange(last_domain.clone(), last_codeword.clone());

        // verify by evaluating
        let poly_eval: Vec<FieldElement> = poly.eval_domain(last_domain);
        assert_eq!(poly_eval, last_codeword, "last codeword is not low degree");
        if poly.degree() > degree {
            println!("last codeword does not correspond to polynomial of low enough degree");
            println!("degree: {}, expected: {}", poly.degree(), degree);
            return false;
        }

        // get indices
        let top_level_indices = self.sample_indices(
            &proof_stream.verifier_fiat_shamir(32),      // seed
            self.domain_length >> 1,                     // size
            self.domain_length >> (self.num_rounds()-1), // reduced size
            self.num_colinearity_tests                   // num tests
        );

        // for every round, check consistency of subsequent layers
        for r in 0..self.num_rounds()-1{

            // fold c indices
            let mut c_indices = vec![];
            for idx in top_level_indices.clone() {
                c_indices.push(idx % (self.domain_length >> (r+1)));
            }

            // infer a and b indices
            let mut a_indices: Vec<usize> = Vec::new();
            let mut b_indices: Vec<usize> = Vec::new();
            for idx in c_indices.clone() { 
                a_indices.push(idx); 
                b_indices.push(idx + (self.domain_length >> (r+1)));    
            }

            // read values and check colinearity
            let mut aa: Vec<FieldElement> = vec![];
            let mut bb: Vec<FieldElement> = vec![];
            let mut cc: Vec<FieldElement> = vec![];
            for s in 0..self.num_colinearity_tests {

                // read points from proof stream
                let (ay, by, cy): (String, String, String) = serde_json::from_str(&proof_stream.pull()).unwrap();

                // convert to field elements
                let ay = FieldElement::new(BigInt::from_str_radix(&ay, 10).expect("BigInt parse err"));
                let by = FieldElement::new(BigInt::from_str_radix(&by, 10).expect("BigInt parse err"));
                let cy = FieldElement::new(BigInt::from_str_radix(&cy, 10).expect("BigInt parse err"));

                // push to respective vectors
                aa.push(ay.clone()); bb.push(by.clone()); cc.push(cy.clone());

                // record top-layer values for layer verification
                if r == 0 {
                    polynomial_values.push( (a_indices[s], ay.clone()) );
                    polynomial_values.push( (b_indices[s], by.clone()) );
                }

                // colinearity check 
                let ax = offset.clone() * omega.clone().pow(a_indices[s] as u128);
                let bx = offset.clone() * omega.pow(b_indices[s] as u128);
                let cx = alphas.get(r).unwrap();
                let points = vec![(ax.value, ay.value), (bx.value, by.value), (cx.value.clone(), cy.value)];
                if false == Polynomial::test_colinearity  (points) {
                    println!("colinearity test failed");
                    return false;
                }
            }   

            // verify authentication paths
            for i in 0..self.num_colinearity_tests {

                // pull authentication path from proof stream for a indices
                let path: String  = proof_stream.pull();
                let deserialized_path: Vec<HashOutput> = serde_json::from_str(&path).unwrap();

                // verify path for a indices
                if false == Merkle::verify(
                    &GenericArray::clone_from_slice(&hex::decode(roots.get(r).unwrap()).unwrap()), // decode --> GenericArray
                    *a_indices.get(i).unwrap(),
                    &deserialized_path,
                    &bincode::serialize(aa.get(i).unwrap()).unwrap()
                ) { return false; }    

                // pull authentication path from proof stream for b indices
                let path: String  = proof_stream.pull();
                let deserialized_path: Vec<HashOutput> = serde_json::from_str(&path).unwrap();

                // verify path for b indices
                if false == Merkle::verify(
                    &GenericArray::clone_from_slice(&hex::decode(roots.get(r).unwrap()).unwrap()), // decode --> GenericArray
                    *b_indices.get(i).unwrap(),
                    &deserialized_path,
                    &bincode::serialize(bb.get(i).unwrap()).unwrap()
                ) { return false; }

                // pull authentication path for c indices
                let path: String  = proof_stream.pull();
                let deserialized_path: Vec<HashOutput> = serde_json::from_str(&path).unwrap();

                // verify path for c indices
                if false == Merkle::verify(
                    &GenericArray::clone_from_slice(&hex::decode(roots.get(r+1).unwrap()).unwrap()), // decode --> GenericArray
                    *c_indices.get(i).unwrap(),
                    &deserialized_path,
                    &bincode::serialize(cc.get(i).unwrap()).unwrap()
                ) { return false; }
            }

            // square omega and offset for next round
            omega = omega.pow(2);
            offset = offset.pow(2);
        }

        true
    }
}


#[cfg(test)]
mod tests {
    use core::panic;

    use super::*;

    #[test]
    fn test_fri() {

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

        let omega = FieldElement::primitive_nth_root(initial_codeword_length as u128);
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

        // setup domain for evaluation
        let mut domain: Vec<FieldElement> = vec![];
        for i in 0..initial_codeword_length {
            domain.push(omega.pow(i as u128).clone());
        }

        // get codeword
        let mut codeword = polynomial.eval_domain(domain);

        // setup proof stream
        println!("Testing Valid Codeword");
        let mut proof_stream = ProofStream::new();

        // prove
        fri.prove(codeword.clone(), &mut proof_stream);

        // verify 
        let mut points: Vec<(usize, FieldElement)> = vec![];
        let verdict: bool = fri.verify(&mut proof_stream, &mut points);

        // ensure correct verdict
        if verdict == false{ panic!("Fri proof should be valid but is not"); } 


        for (idx, val) in points {
            if polynomial.eval(omega.pow(idx as u128)) != val {
                panic!("Polynomial evals to wrong value");
            }
        }
        println!("Success!");

        println!("Testing Invalid Codeword...");

        // disturb codeword and test again
        let mut proof_stream = ProofStream::new();
        for i in 0..degree/3{
            codeword[i] = FieldElement::zero();
        }

        fri.prove(codeword.clone(), &mut proof_stream);
        points = vec![];
        assert_eq!(false, fri.verify(&mut proof_stream, &mut points));
        println!("Success!!");
    }
}