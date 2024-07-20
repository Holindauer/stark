// use sha3::{Digest, Shake256};
// use sha3::digest::{Update, ExtendableOutput};

// use serde::{Serialize, Deserialize};
// use std::string;


// use bincode;

use sha3::{Digest, Shake256};
use sha3::digest::{Update, ExtendableOutput, XofReader};
use serde::{Serialize, Deserialize};
use bincode;
use std::vec;

#[derive(Serialize, Deserialize)]
struct ProofStream {
    objects: Vec<String>,
    read_idx: usize,
}

impl ProofStream {
    // Constructor
    fn new() -> Self {
        ProofStream {
            objects: vec![],
            read_idx: 0,
        }
    }

    // Method to add an object to the stream
    fn push(&mut self, obj: String) {
        self.objects.push(obj);
    }

    // Method to retrieve and remove the next object from the stream
    fn pull(&mut self) -> String {
        assert!(self.read_idx < self.objects.len(), "ProofStream: cannot pull object; queue empty.");
        let obj = self.objects[self.read_idx].clone();
        self.read_idx += 1;
        obj
    }

    // Method to serialize the proof stream
    fn serialize(&self) -> Vec<u8> {
        bincode::serialize(&self.objects).unwrap()
    }

    // Static method to deserialize the proof stream
    fn deserialize(bb: &[u8]) -> Self {
        let objects: Vec<String> = bincode::deserialize(bb).unwrap();
        ProofStream {
            objects,
            read_idx: 0,
        }
    }

    // Method for the prover to generate a hash
    fn prover_fiat_shamir(&self, num_bytes: usize) -> Vec<u8> {
        let data = self.serialize();
        let mut hasher = Shake256::default();
        hasher.update(&data);
        let mut output = vec![0u8; num_bytes];
        let mut xof = hasher.finalize_xof();
        xof.read(&mut output);
        output
    }

    // Method for the verifier to generate a hash
    fn verifier_fiat_shamir(&self, num_bytes: usize) -> Vec<u8> {
        let data = bincode::serialize(&self.objects[..self.read_idx]).unwrap();
        let mut hasher = Shake256::default();
        hasher.update(&data);
        let mut output = vec![0u8; num_bytes];
        let mut xof = hasher.finalize_xof();
        xof.read(&mut output);
        output
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_serialize() {
        let mut proof1 = ProofStream::new();
        proof1.push("1".to_string());
        proof1.push("2".to_string());
        proof1.push("3".to_string());
        proof1.push("4".to_string());

        let serialized = proof1.serialize();
        let mut proof2 = ProofStream::deserialize(&serialized);

        assert_eq!(proof1.pull(), proof2.pull(), "pulled object 0 don't match");
        assert_eq!(proof1.pull(), proof2.pull(), "pulled object 1 don't match");
        assert_eq!(proof1.pull(), proof2.pull(), "pulled object 2 don't match");
        assert_eq!(proof1.pull(), "4", "object 3 pulled from proof1 is not 4");
        assert_eq!(proof2.pull(), "4", "object 3 pulled from proof2 is not 4");
        assert_eq!(proof1.prover_fiat_shamir(32), proof2.prover_fiat_shamir(32), "fiat shamir is not the same");
        
        println!("All tests passed");
    }
}
