use crate::modules::field::{*};
use crate::modules::merkle::{*};
use crate::modules::proof_stream::{*};
use blake2::{Blake2b, Digest};
use generic_array::GenericArray;
use typenum::U64;

use num_bigint::RandBigInt;
use num_bigint::BigInt;
use num_traits::{Zero, One};
use num_traits::ToPrimitive;

struct Fri {
    offset: FieldElement,
    omega: FieldElement,
    domain_length: i32,
    expansion_factor: i32,
    num_colinearity_tests: i32,
}

impl Fri {
    
    // constructor 
    fn new( 
        offset: FieldElement, 
        omega: FieldElement, 
        initial_domain_length: i32,
        expansion_factor: i32, 
        num_colinearity_tests: i32 
        ) -> Fri {        


        let fri_obj = Fri { offset, omega, domain_length: initial_domain_length, expansion_factor, num_colinearity_tests };
        assert!(fri_obj.num_rounds() > 0, "Cannot do FRI with less than 1 round");
        fri_obj
    }

    // computes num rounds of the fri protocol
    fn num_rounds(&self) -> i32 {

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
    fn eval_domain(&self){

        let mut domain: Vec<FieldElement> = Vec::new();
        for i in 0..self.domain_length {
            domain.push( self.offset.clone() * self.omega.pow(i as i128));
        }
        domain
    }


}