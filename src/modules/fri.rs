use crate::modules::field::{*};
use crate::modules::merkle::{*};
use crate::modules::proof_stream::{*};
use crate::modules::univariate_poly::{*};
use blake2::{Blake2b, Digest};
use typenum::U64;
use num_bigint::BigInt;
use num_traits::{Zero, One};
use hex;
use serde::{Serialize, Deserialize};

use generic_array::typenum::U32;
use generic_array::GenericArray;

type OutputSize = U32; // 32 bytes (256 bits) output
type HashOutput = GenericArray<u8, OutputSize>; // Fixed-size hash output


struct Fri {
    offset: FieldElement,
    omega: FieldElement,
    domain_length: usize,
    expansion_factor: usize,
    num_colinearity_tests: usize,
}

impl Fri {
    
    // constructor 
    fn new( offset: FieldElement, omega: FieldElement, initial_domain_length: usize, expansion_factor: usize, num_colinearity_tests: usize ) -> Fri {        

        let fri_obj = Fri { offset, omega, domain_length: initial_domain_length, expansion_factor, num_colinearity_tests };
        assert!(fri_obj.num_rounds() > 0, "Cannot do FRI with less than 1 round");
        fri_obj
    }

    // computes num rounds of the fri protocol
    fn num_rounds(&self) -> usize {

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
    fn sample_indices(&self, seed: &[u8], size: usize, reduced_size: usize, number: usize) -> Vec<usize> {
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
    fn sample_index(byte_array: Vec<u8>, size: usize) -> usize {
        let mut acc: usize = 0;
        for b in byte_array {
            acc = (acc << 8) ^ (b as usize);
        }
        acc % size
    }

    // computes domain for polynomial evaluation
    fn eval_domain(&self) -> Vec<FieldElement> {

        let mut domain: Vec<FieldElement> = Vec::new();
        for i in 0..self.domain_length {
            domain.push( self.offset.clone() * self.omega.pow(i as i128));
        }
        domain
    }

    // commit phase
    pub fn commit(&self, mut codeword: Vec<FieldElement>, proof_stream: &mut ProofStream) -> Vec<Vec<FieldElement>> {
        
        // shorthand
        let one = FieldElement::one();
        let two = FieldElement::new(BigInt::from(2));

        // working values
        let mut omega = self.omega.clone();
        let mut offset = self.offset.clone();
        let mut codewords = Vec::new();

        println!("rounds: {}", self.num_rounds());

        // split and fold rounds
        for r in 0..self.num_rounds() {
            let N: usize = codeword.len();

            // ensure n has the right order
            assert!(omega.pow(N as i128 - 1) == omega.inverse(), "error in commit: omega does not have the right order!");

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
                    two.inverse() * ( (one.clone() + alpha.clone() / (offset.clone() * omega.pow(i as i128))) * codeword.get(i).unwrap().clone() + (one.clone() - alpha.clone() / (offset.clone() * omega.pow(i as i128)))  * codeword.get(N/2 + i).unwrap().clone()) // long boi
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
}


#[cfg(test)]
mod tests {
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

        // setup domain for evaluation
        let mut domain: Vec<FieldElement> = vec![];
        for i in 0..initial_codeword_length {
            domain.push(omega.pow(i as i128).clone());
        }

        // get codeword
        let codeword = polynomial.eval_domain(domain);

        // setup proof stream
        let mut proof_stream = ProofStream::new();





        // commit test
        let codewords = fri.commit(codeword, &mut proof_stream);


        // print len of each codword
        for c in codewords {
            println!("Codeword length: {}", c.len());
        }

    }
}