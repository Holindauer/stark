use crate::modules::field::{*};
use crate::modules::fri::{*};
use crate::modules::rescue_prime::{*};
use crate::modules::univariate_poly::{*};
use crate::modules::multivariate_poly::{*};
use crate::modules::proof_stream::{*};
use crate::modules::merkle::{*};

use rand::{random, RngCore};
use num_bigint::BigInt;
use blake2::{Blake2b, Digest};
use std::collections::HashMap;
use generic_array::typenum::U32;
use generic_array::GenericArray;


pub struct Stark {
    pub expansion_factor: usize,
    pub num_colinearity_tests: usize,
    pub security_level: usize,
    pub num_randomizers : usize,
    pub num_registers: usize,
    pub original_trace_length: usize,   
    pub generator: FieldElement,
    pub omega: FieldElement,
    pub omicron: FieldElement,
    pub omicron_domain: Vec<FieldElement>,
    pub fri: Fri,
}   


impl Stark {

    pub fn new(
        expansion_factor: usize, 
        num_colinearity_tests: usize, 
        security_level: usize, 
        num_registers: usize, 
        original_trace_length: usize
    ) -> Stark {
        // TODO impl asserts here

        // transition constraints degree
        let transition_constraints_degree = 2;

        // randomizers
        let num_randomizers = 4 * num_colinearity_tests;

        // compute domain lengths
        let randomized_trace_length = original_trace_length + num_randomizers;

        // omicron domain length
        let product = randomized_trace_length * transition_constraints_degree;
        let bit_length = (std::mem::size_of_val(&product) * 8) - product.leading_zeros() as usize; 
        let omicron_domain_length = 1 << bit_length;   

        // fri domain length
        let fri_domain_length = omicron_domain_length * expansion_factor;

        // field elements
        let omega: FieldElement = FieldElement::primitive_nth_root(fri_domain_length as u128);
        let omicron: FieldElement = FieldElement::primitive_nth_root(omicron_domain_length as u128);
        let omicron_domain: Vec<FieldElement> = (0..omicron_domain_length).map(|i| omicron.pow(i as u128)).collect();

        assert_eq!(omicron.value, BigInt::from(65907963977709178563567092354521124432 as u128));

        // generator
        let generator = FieldElement::generator();

        // init fri
        let fri = Fri::new(generator.clone(), omega.clone(), fri_domain_length, expansion_factor, num_colinearity_tests);

        Stark {
            expansion_factor,
            num_colinearity_tests,
            security_level,
            num_randomizers,
            num_registers,
            original_trace_length,
            generator,
            omega,
            omicron,
            omicron_domain,
            fri,  
        }
    }
    
    // boundary zeroifiers
    fn boundary_zeroifiers(&self, boundary: Vec<(usize, usize, FieldElement)>) -> Vec<Polynomial> {
        
        let mut zeroifiers: Vec<Polynomial> = vec![];
        for s in 0..self.num_registers {

            // get points
            let mut points: Vec<FieldElement>= vec![];
            for (c, r, v) in boundary.clone() {
                if r == s {
                    points.push(self.omicron.pow(c as u128));
                }
            }

            // create zeroifier
            let zeroifier = Polynomial::zeroifier_domain(points.clone());
            zeroifiers.push(zeroifier.clone());
        }   
        zeroifiers
    } 

    // boundary interpolants
    fn boundary_interpolants(&self, boundary: Vec<(usize, usize, FieldElement)>) -> Vec<Polynomial> {

        let mut interpolants: Vec<Polynomial> = vec![];
        for s in 0..self.num_registers {

            // get interpolation domain and values
            let mut domain: Vec<FieldElement> = vec![];
            let mut values: Vec<FieldElement> = vec![];
            for (c, r, v) in boundary.clone() {
                if r == s {
                    domain.push(self.omicron.pow(c as u128));
                    values.push(v.clone());
                }
            }

            // create interpolant
            let interpolatnt = Polynomial::lagrange(domain.clone(), values.clone());
            interpolants.push(interpolatnt.clone());
        }   

        interpolants
    }  

    // zeroifier over omicron domain
    fn transition_zeroifier(&self) -> Polynomial {
        let domain: Vec<FieldElement> = self.omicron_domain[0..(self.original_trace_length-1)].to_vec();
        Polynomial::zeroifier_domain(domain)
    }

    // transition poly degree bounds
    fn transition_degree_bounds(&self, transition_constraints: Vec<MPolynomial>) -> Vec<usize> {

        // reconstruct point degrees
        let mut point_degrees: Vec<usize> = vec![];
        point_degrees.push(1);
        for i in 0..2*self.num_randomizers {
            point_degrees.push(self.original_trace_length + self.num_randomizers-1);
        }

        // get max from sums
        let mut maxes: Vec<usize> = vec![];
        for a in transition_constraints {

            let mut sums: Vec<usize> = vec![];
            for (k, v) in a.dict.iter() {
                
                let mut sum: usize = 0;
                let zipped: Vec<(&usize, &u128)>= point_degrees.iter().zip(k.iter()).collect();
                for (r, l) in zipped {
                    sum += r * (*l as usize);
                }
                sums.push(sum);
            }

            // push max
            maxes.push(sums.iter().max().unwrap().clone());
        }
        maxes
    }

    // transition quotient poly degree bounds
    fn transition_quotient_degree_bounds(&self, transition_constraints: Vec<MPolynomial>) -> Vec<usize> {
        
        let mut bounds: Vec<usize> = vec![];
        for d in self.transition_degree_bounds(transition_constraints){
            bounds.push(d - (self.original_trace_length - 1))
        }
        bounds
    }

    // boundary quotient poly degree bounds
    fn boundary_quotient_degree_bounds(&self, randomized_trace_length: usize, boundary: Vec<(usize, usize, FieldElement)>) -> Vec<usize> {
        let randomized_trace_degree = randomized_trace_length - 1;
        let mut bounds: Vec<usize> = vec![];
        for bz in self.boundary_zeroifiers(boundary) {
            bounds.push(randomized_trace_degree - bz.degree());
        }
        bounds
    }

    // max degree
    fn max_degree(&self, transition_constraints: Vec<MPolynomial>) -> usize {
       let md: usize = *self.transition_quotient_degree_bounds(transition_constraints).iter().max().unwrap(); 
        // Calculate the next power of 2 minus 1
        // This is equivalent to Python's (1 << (len(bin(md)[2:]))) - 1
        if md == 0 {
            0
        } else {
            // Number of bits needed to represent md
            let bits_needed = 64 - md.leading_zeros() as usize;
            (1 << bits_needed) - 1
        }
    }

    // smaple weights 
    fn sample_weights(&self, number: usize, randomness: Vec<u8>) -> Vec<FieldElement> {

        let mut weights: Vec<FieldElement> = Vec::with_capacity(number);
        for i in 0..number {
            // hash w/ blake2b 
            let mut hasher = Blake2b::new();
            hasher.update(&randomness);
            hasher.update(i.to_le_bytes()); // Convert index to bytes in little-endian format
            let hash_result: HashOutput = hasher.finalize(); 

            // sample field element
            weights.push(FieldElement::sample(hash_result.to_vec()));
        }

        weights
    }

    // stark prover
    pub fn prove(
        &self, 
        trace: Vec<Vec<FieldElement>>, 
        transition_constraints: Vec<MPolynomial>, 
        boundary: Vec<(usize, usize, FieldElement)>
    ) -> Vec<u8> {

        // make trace mutable
        let mut trace = trace;

        // create proof stream object 
        let mut proof_stream: ProofStream = ProofStream::new();

        // concatenate randomizers
        for k in 0..self.num_randomizers {

            // gen randomizer
            let mut randomizer: Vec<FieldElement> = vec![];
            for s in 0..self.num_registers {

                // create random bytes
                let mut rng = rand::thread_rng();
                let mut random_bytes: Vec<u8> = vec![0; 17];
                rng.fill_bytes(&mut random_bytes);       

                // sample random value         
                let value: FieldElement = FieldElement::sample(random_bytes);
                randomizer.push(value);
            }
            trace.push(randomizer);
        }

        // create trace domain
        let mut trace_domain = vec![];
        for i in 0..trace.len() { trace_domain.push(self.omicron.pow(i as u128)); }

        // interpolate trace polynomials over domain
        let mut trace_polynomials: Vec<Polynomial> = vec![];
        for s in 0..self.num_registers {

            // seperate trace for single register
            let mut single_trace: Vec<FieldElement> = vec![];
            for c in 0..trace.len() { 
                single_trace.push(trace[c][s].clone());
            }

            // interpolate trace polynomial
            let trace_poly = Polynomial::lagrange(trace_domain.clone(), single_trace.clone());
            assert!(trace_domain.len() == single_trace.len());

            trace_polynomials.push(trace_poly.clone());

        }

        // get boundary interpolants and zeroifiers
        let interpolants = self.boundary_interpolants(boundary.clone());
        let zeroifiers = self.boundary_zeroifiers(boundary.clone());

        // subtract boundary inteprolants and divide out zeroifiers
        let mut boundary_quotients: Vec<Polynomial> = vec![];
        for s in 0..self.num_registers {
            let interpolant = interpolants[s].clone();
            let zeroifier = zeroifiers[s].clone();
            let quotient = (trace_polynomials[s].clone() - interpolant.clone()) / zeroifier.clone();
            boundary_quotients.push(quotient.clone());
        }

        // commit to each boundary quotient
        let fri_domain = self.fri.eval_domain();
        let mut boundary_quotient_codewords: Vec<Vec<FieldElement>> = vec![];
        for s in 0..self.num_registers {


            // evaluate boundary quotient
            let codeword: Vec<FieldElement> = boundary_quotients[s].eval_domain(fri_domain.clone());
            boundary_quotient_codewords.push(codeword.clone());


            // merkleize boundary quotient
            let merkle_root = Merkle::commit(&codeword.iter().map(|x| bincode::serialize(x).unwrap()).collect::<Vec<Vec<u8>>>());            
            
            // push to proof stream
            proof_stream.push(hex::encode(merkle_root));   
        }

        // setup point for symbolic evaluation
        let mut point: Vec<Polynomial> = vec![];
        // The polynomial X has coefficients [0, 1] following Python reference (lowest to highest)
        // But our Rust implementation uses highest to lowest, so we need [1, 0]
        // This represents 1*x + 0 = x
        point.push(Polynomial::new(vec![BigInt::from(1), BigInt::from(0)]));


        point.extend(trace_polynomials.clone());
        for tp in trace_polynomials.clone() {
            point.push(tp.scale(&self.omicron.clone()));
        }   

        // symbolically evaluate transition constraints
        let mut transition_polynomials: Vec<Polynomial> = vec![];
        for a in transition_constraints.clone() {
            transition_polynomials.push(
                a.eval_symbolic(&point.clone())
            )
        }    

        // divide out zeroifier
        let mut transition_quotients: Vec<Polynomial> = vec![];
        for tp in transition_polynomials {
            transition_quotients.push(tp / self.transition_zeroifier());
        }

        // create randomizer polynomial
        let mut rand_coeffs: Vec<FieldElement> = vec![];
        let mut rng = rand::thread_rng();
        for i in 0..self.max_degree(transition_constraints.clone())+1{

            // random field element coefficient          
            let mut random_bytes: Vec<u8> = vec![0; 17];
            rng.fill_bytes(&mut random_bytes);  
            rand_coeffs.push(FieldElement::sample(random_bytes));
        }
        let randomizer_poly = Polynomial{coeffs: rand_coeffs};

        // commit to randomizer polynomial
        let randomizer_codeword: Vec<FieldElement> = randomizer_poly.eval_domain(fri_domain.clone());
        let randomizer_root: HashOutput = Merkle::commit(
            &randomizer_codeword.iter().map(|x| bincode::serialize(x).unwrap()).collect::<Vec<Vec<u8>>>()
        );
        proof_stream.push(hex::encode(randomizer_root));

        // println!("after randomizer polynomial proof stream len {}", proof_stream.objects.len());



        /*
            get weights for nonlinear combination
              - 1 randomizer
              - 2 for every transition quotient
              - 3 for every boundary quotient
         */
        let fiat_shamir_bytes = proof_stream.prover_fiat_shamir(32);
        let weights = self.sample_weights(
            1 + 2*transition_quotients.len() + 2*boundary_quotients.len(),
            fiat_shamir_bytes.clone()
        );
        println!("DEBUG Prover: Number of weights sampled: {}", weights.len());
        println!("DEBUG Prover: Num transition quotients: {}", transition_quotients.len());
        println!("DEBUG Prover: Num boundary quotients: {}", boundary_quotients.len());
        println!("DEBUG Prover: First 4 bytes of Fiat-Shamir: {:?}", &fiat_shamir_bytes[..4]);
        println!("DEBUG Prover: First weight: {}", weights[0].value);

        // ensure transition quotient degrees match degree bounds
        let tq_degrees: Vec<usize> = transition_quotients.iter().map(|tq| tq.degree()).collect();
        assert!(tq_degrees == self.transition_quotient_degree_bounds(transition_constraints.clone()));
    
        // compute terms of nonlinear combination polynomial
        // The polynomial x has coefficients [1, 0] (highest to lowest in our Rust impl)
        // This represents 1*x + 0 = x
        let x = Polynomial{coeffs: vec![FieldElement::one(), FieldElement::zero()]};
        let max_degree = self.max_degree(transition_constraints.clone()); 
        let mut terms: Vec<Polynomial> = vec![];
        terms.push(randomizer_poly);
        for i in 0..transition_quotients.len() {
            terms.push(transition_quotients[i].clone());
            let shift = max_degree - self.transition_quotient_degree_bounds(transition_constraints.clone())[i];
            terms.push(x.pow(shift as u128) * transition_quotients[i].clone());
        }
        for i in 0..self.num_registers {
            terms.push(boundary_quotients[i].clone());
            let shift = max_degree - self.boundary_quotient_degree_bounds(trace.len(), boundary.clone())[i];
            terms.push(x.pow(shift as u128) * boundary_quotients[i].clone())
        }

        // take weighted sum
        let mut combination = Polynomial{coeffs: vec![FieldElement::zero()]};
        for i in 0..terms.len() {
            combination = combination + (Polynomial{coeffs: vec![weights[i].clone()]} * terms[i].clone());
        }
        println!("DEBUG: Combination polynomial degree: {}", combination.degree());
        println!("DEBUG: Max degree: {}", max_degree);
        // compute matching codeword
        let combined_codeword = combination.eval_domain(fri_domain.clone());
        println!("DEBUG: FRI domain length: {}", fri_domain.len());
        println!("DEBUG: First few FRI domain points: {:?}", 
                 fri_domain.iter().take(3).map(|x| x.value.to_string()).collect::<Vec<_>>());
        println!("DEBUG: First few values of combined codeword: {:?}", 
                 combined_codeword.iter().take(5).map(|x| x.value.to_string()).collect::<Vec<_>>());

        // prove low degree of combination polynomial w/ fri, and collect indices
        let mut indices = self.fri.prove(combined_codeword.clone(), &mut proof_stream);
        println!("DEBUG: FRI returned indices: {:?}", indices);
        indices.sort(); // Sort indices to match verifier's expectation
        println!("DEBUG: Sorted indices: {:?}", indices);
        
        // Debug: Check values at the indices FRI is querying
        for (i, idx) in indices.iter().enumerate().take(2) {
            println!("DEBUG Prover: combined_codeword[{}] = {}", idx, combined_codeword[*idx].value);
            // Also check the evaluation domain point
            let domain_point = self.fri.offset.clone() * self.fri.omega.pow(*idx as u128);
            println!("DEBUG Prover: domain point at index {} = {}", idx, domain_point.value);
            println!("DEBUG Prover: combination.eval(domain_point) = {}", combination.eval(domain_point));
        }

        println!("after fri proof stream len {}", proof_stream.objects.len());

        // process indices
        let mut duplicated_indices: Vec<usize> = indices.clone();
        for i in indices.clone() {
            duplicated_indices.push( (i + self.expansion_factor) % self.fri.domain_length );
        }
        duplicated_indices.sort(); // Sort to match verifier's expectation
        
        let mut quadrupled_indices: Vec<usize> = duplicated_indices.clone();
        for i in duplicated_indices.clone() {
            quadrupled_indices.push( (i + (self.fri.domain_length / 2)) % self.fri.domain_length )
        }
        quadrupled_indices.sort();


        // open indicated positions in the boudnary quotient codeword
        for bqc in boundary_quotient_codewords {
            for i in duplicated_indices.clone() {

                // commit to element of codeword at index
                proof_stream.push(serde_json::to_string(&bqc[i]).unwrap());

                // get authentication path for i'th word in codeword
                let serialized_bqc: Vec<Vec<u8>> = bqc.iter() 
                .map(|element| bincode::serialize(element).unwrap()).collect();// convert to bytes
                let auth_path: String = serde_json::to_string(&Merkle::open(i, &serialized_bqc)).unwrap();

                // commit to auth path
                proof_stream.push(auth_path);
            }
        }

        println!("after pushing bqc stream len {}", proof_stream.objects.len());


        // ... as well as in the randomizer
        for (idx_pos, i) in indices.clone().into_iter().enumerate() {

            // commit to i'th word in codeword
            let randomizer_value = &randomizer_codeword[i];
            if idx_pos == 0 {
                println!("DEBUG Prover: Sending randomizer[{}] = {}", i, randomizer_value.value);
            }
            proof_stream.push(serde_json::to_string(randomizer_value).unwrap());
            
            // get authentication path for i'th word in codeword
            let serialized_rc: Vec<Vec<u8>> = randomizer_codeword.iter() 
            .map(|element| bincode::serialize(element).unwrap()).collect();// convert to bytes
            let auth_path: String = serde_json::to_string(&Merkle::open(i, &serialized_rc)).unwrap();

            // commit to auth path
            proof_stream.push(auth_path);
        }

        println!("after pushing randomizer proof stream len {}", proof_stream.objects.len());


        // final proof is the serialized proof stream
        proof_stream.serialize()
    }


    pub fn verify(
        &self,
        proof: Vec<u8>,
        transition_constraints: Vec<MPolynomial>,
        boundary: Vec<(usize, usize, FieldElement)>
    ) -> bool{

        // infer trace length from boundary conditions
        let mut cycles: Vec<usize> = vec![];
        for (c, r, v) in boundary.clone() { cycles.push(c); }        
        let original_trace_length = 1 + cycles.iter().max().unwrap();
        let randomized_trace_length = original_trace_length + self.num_randomizers;

        // deserialize with right proof stream
        let mut proof_stream = ProofStream::deserialize(&proof);

        // get merkle roots of boundary quotient codewords
        let mut boundary_quotient_roots: Vec<String> = vec![];
        for s in 0..self.num_registers {
            boundary_quotient_roots.push(proof_stream.pull());
        }

        // get merkle root of randomizer polynomial
        let randomizer_root: String = proof_stream.pull();

        println!("read idx after randomizer root {}", proof_stream.read_idx);

        // get weights for nonlinear combination
        let fiat_shamir_bytes = proof_stream.verifier_fiat_shamir(32);
        let weights = self.sample_weights(
            1 + 2*transition_constraints.len() + 2*self.boundary_interpolants(boundary.clone()).len(),
            fiat_shamir_bytes.clone()
        );
        println!("DEBUG Verifier: Number of weights sampled: {}", weights.len());
        println!("DEBUG Verifier: First 4 bytes of Fiat-Shamir: {:?}", &fiat_shamir_bytes[..4]);
        println!("DEBUG Verifier: First weight: {}", weights[0].value);

        // verify low degree of combination polynomial
        let mut polynomial_values: Vec<(usize, FieldElement)> = vec![];
        let mut verifier_accepts: bool = self.fri.verify( &mut proof_stream, &mut polynomial_values);
        polynomial_values.sort_by_key(|iv| iv.0); // sort
        println!("DEBUG Verifier: FRI returned {} polynomial values", polynomial_values.len());
        if polynomial_values.len() > 0 {
            println!("DEBUG Verifier: First polynomial value: ({}, {})", 
                     polynomial_values[0].0, polynomial_values[0].1.value);
        }


        println!("read idx after fri {}", proof_stream.read_idx);

        if verifier_accepts == false { 
            println!("DEBUG: FRI verification failed");
            return false; 
        }
        
        // isolate indices and values
        let indices: Vec<usize> = polynomial_values.iter().map(|iv| iv.0).collect();
        let values: Vec<FieldElement> = polynomial_values.iter().map(|iv| iv.1.clone()).collect();
        
        println!("DEBUG Verifier: FRI returned indices: {:?}", indices);


        // read and verify leafs which are elements of boundary quotient codewords
        let mut duplicated_indices: Vec<usize> = indices.clone();
        for i in indices.clone() { duplicated_indices.push((i + self.expansion_factor) % self.fri.domain_length)};
        duplicated_indices.sort();

        // collect leafs/auth paths into vec, verify auth paths
        let mut leafs: Vec<HashMap<usize, FieldElement>> = vec![];
        for r in 0..boundary_quotient_roots.len() {

            let mut hash_map = HashMap::new();
            for i in duplicated_indices.clone() {
                
                // get leaf value and insert into map
                let leaf_value: FieldElement = serde_json::from_str(&proof_stream.pull()).unwrap();
                hash_map.insert(i.clone(), leaf_value.clone());

                // get auth_path for leaf
                let auth_path: Vec<HashOutput> = serde_json::from_str(&proof_stream.pull()).unwrap();

                // verify auth path
                let merkle_verification: bool = Merkle::verify(
                    &GenericArray::clone_from_slice(&hex::decode(boundary_quotient_roots[r].clone()).unwrap()), // decode --> GenericArray
                    i,
                    &auth_path,
                    &bincode::serialize(&leaf_value).unwrap()
                );
                
                // update verifier acceptance
                let verifier_accepts = verifier_accepts && merkle_verification;
                if !verifier_accepts {
                    return false;
                }
            }
            leafs.push(hash_map);
        }

        // read and verify randomizer leafs
        let mut randomizer: HashMap<usize, FieldElement> = HashMap::new();
        for i in indices.clone() {

            // get leaf value and insert into map
            let leaf_value: FieldElement = serde_json::from_str(&proof_stream.pull()).unwrap();
            randomizer.insert(i.clone(), leaf_value.clone());
            
            // get auth_path for leaf
            let auth_path: Vec<HashOutput> = serde_json::from_str(&proof_stream.pull()).unwrap();
        
            // verify auth path
            let merkle_verification: bool = Merkle::verify(
                &GenericArray::clone_from_slice(&hex::decode(randomizer_root.clone()).unwrap()), // decode --> GenericArray
                i,
                &auth_path,
                &bincode::serialize(&leaf_value).unwrap()
            );
            
            // update verifier acceptance
            let verifier_accepts = verifier_accepts && merkle_verification;
            if !verifier_accepts {
                return false;
            }
            
        }

        // verify leafs of combination polynomial
        for i in 0..indices.len() {
            let current_index = indices[i];

            // get trace values by applying a correction ot the boundary quotient values (which are the leafs)
            let domain_current_index = self.generator.clone() * (self.omega.pow(current_index as u128));
            let next_index = (current_index + self.expansion_factor) % self.fri.domain_length;
            let domain_next_index = self.generator.clone() * (self.omega.pow(next_index as u128));
            
            if i == 0 {
                println!("DEBUG Verifier: Evaluating combination at index {} (domain point {})", 
                         current_index, domain_current_index.value);
                println!("DEBUG Verifier: Generator: {}, Omega: {}", self.generator.value, self.omega.value);
                println!("DEBUG Verifier: FRI offset: {}, FRI omega: {}", self.fri.offset.value, self.fri.omega.value);
                println!("DEBUG Verifier: Expected FRI value at this index: {}", values[i].value);
            }

            // current trace
            let mut current_trace: Vec<FieldElement> = vec![];
            for i in 0..self.num_registers { current_trace.push(FieldElement::zero());}

            // next trace
            let mut next_trace: Vec<FieldElement> = vec![];
            for i in 0..self.num_registers { next_trace.push(FieldElement::zero());}

            for s in 0..self.num_registers {
                // retrieve zeroifier and interpolant
                let zeroifier = self.boundary_zeroifiers(boundary.clone())[s].clone();
                let interpolant = self.boundary_interpolants(boundary.clone())[s].clone();
                
                
                let current_trace_value = leafs[s].get(&current_index).unwrap().clone() * zeroifier.eval(domain_current_index.clone()) + interpolant.eval(domain_current_index.clone());
                let next_trace_value = leafs[s].get(&next_index).unwrap().clone() * zeroifier.eval(domain_next_index.clone()) + interpolant.eval(domain_next_index.clone());
                current_trace[s] = current_trace_value;
                next_trace[s] = next_trace_value;
            }
            
            // get eval point
            let mut point: Vec<FieldElement> = vec![];
            point.push(domain_current_index.clone());
            point.extend(current_trace.clone());
            point.extend(next_trace.clone());

            // eval transition constraints
            let mut transition_constraints_values: Vec<FieldElement> = vec![];
            for s in 0..transition_constraints.len() {
                let tc_value = transition_constraints[s].clone().eval(&point);
                transition_constraints_values.push(tc_value.clone());
                if i == 0 && s == 0 {
                    println!("DEBUG Verifier: Evaluating transition constraint {} at point", s);
                    println!("  Point has {} elements", point.len());
                    println!("  Point[0] (x): {}", point[0].value);
                    println!("  Point[1] (current state[0]): {}", point[1].value);
                    println!("  Point[2] (current state[1]): {}", point[2].value);
                    println!("  Point[3] (next state[0]): {}", point[3].value);
                    println!("  Point[4] (next state[1]): {}", point[4].value);
                    println!("  Result: {}", tc_value.value);
                }
            }

            let counter = 0;
            let mut terms: Vec<FieldElement> = vec![];
            let randomizer_term = randomizer.get(&current_index).unwrap().clone();
            terms.push(randomizer_term.clone());
            if i == 0 {
                println!("DEBUG Verifier: Term 0 (randomizer): {}", randomizer_term.value);
            }
            
            for s in 0..transition_constraints_values.len(){
                let tcv = transition_constraints_values[s].clone();
                let tz_eval = self.transition_zeroifier().eval(domain_current_index.clone());
                if i == 0 && s == 0 {
                    println!("  Transition constraint value {}: {}", s, tcv.value);
                    println!("  Transition zeroifier eval: {}", tz_eval.value);
                }
                let quotient = tcv / tz_eval;
                terms.push(quotient.clone());
                let shift = self.max_degree(transition_constraints.clone()) - self.transition_quotient_degree_bounds(transition_constraints.clone())[s];
                let shifted_term = quotient.clone() * (domain_current_index.pow(shift as u128));
                terms.push(shifted_term.clone());
                if i == 0 {
                    println!("  Term {} (tq {}): {}", 1 + s*2, s, quotient.value);
                    println!("  Term {} (tq {} shifted by {}): {}", 2 + s*2, s, shift, shifted_term.value);
                }
            }
            for s in 0..self.num_registers {
                let bqv = leafs[s].get(&current_index).unwrap().clone();
                terms.push(bqv.clone());
                let shift = self.max_degree(transition_constraints.clone()) - self.boundary_quotient_degree_bounds(randomized_trace_length, boundary.clone())[s];
                let shifted_term = bqv.clone() * (domain_current_index.pow(shift as u128));
                terms.push(shifted_term.clone());
                if i == 0 {
                    println!("  Term {} (bq {} from Merkle): {}", 5 + s*2, s, bqv.value);
                    println!("  Term {} (bq {} shifted by {}): {}", 6 + s*2, s, shift, shifted_term.value);
                    
                    // Debug: What should this boundary quotient value be?
                    let trace_value = current_trace[s].clone();
                    let interpolant = self.boundary_interpolants(boundary.clone())[s].eval(domain_current_index.clone());
                    let zeroifier = self.boundary_zeroifiers(boundary.clone())[s].eval(domain_current_index.clone());
                    let expected_bq = (trace_value - interpolant) / zeroifier;
                    println!("  DEBUG: Reconstructed bq {}: {}", s, expected_bq.value);
                }
            }
            
            // construct combination
            let mut combination: FieldElement = FieldElement::zero();
            for j in 0..terms.len() {
                combination = combination + (terms[j].clone() * weights[j].clone())
            }
            
            // verify against combination polynomail value
            if i == 0 {
                println!("DEBUG Verifier: Computed combination: {}", combination.value);
                println!("DEBUG Verifier: Expected (from FRI): {}", values[i].value);
            }
            verifier_accepts = verifier_accepts && (combination == values[i]);
            if !verifier_accepts {
                println!("Verification failed at index {}: combination {} != values[i] {}", i, combination.value, values[i].value);
                println!("  Current index: {}", current_index);
                println!("  Domain current index: {}", domain_current_index.value);
                println!("  Number of terms: {}", terms.len());
                println!("  Weights length: {}", weights.len());
                return false;
            }

        }

        println!("DEBUG: Final verifier_accepts = {}", verifier_accepts);
        verifier_accepts
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stark() {

        // test setup
        let expansion_factor: usize = 4;
        let num_colinearity_tests: usize = 2;
        let security_level: usize = 2;

        // init rescue prime
        let rp = RescuePrime::new();
        let mut output_element = FieldElement::sample(vec![48, 120, 100, 101, 97, 100, 98, 101, 101, 102]); // 0xdeadbeef

        for trial in 0..20 {

            // setup
            let input_element = output_element.clone();
            println!("running trial with input: {}", input_element.value);
            output_element = rp.hash(input_element.clone());
            let num_cycles = rp.N+1;
            let state_width = rp.m;

            // init stark
            let stark = Stark::new(
                expansion_factor,
                num_colinearity_tests,
                security_level,
                state_width,
                num_cycles
            );

            // prove honestly
            println!("honest proof generation...");

            // prove
            let trace: Vec<Vec<FieldElement>> = rp.trace(input_element.clone());
            let air: Vec<MPolynomial> = rp.transition_constraints(stark.omicron.clone());
            let boundary: Vec<(usize, usize, FieldElement)> = rp.boundary_constraints(output_element.clone());
            let proof = stark.prove(trace, air.clone(), boundary.clone());
            println!("num bytes in proof: {:?}", proof.len());

            // verify
            let verdict = stark.verify(proof, air.clone(), boundary.clone());

            assert!(verdict);

            break;
        }
    }
}