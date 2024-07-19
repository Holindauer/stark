// use crate::modules::field::{*};
// use std::vec;
// use std::collections::HashSet;
// use sha2::{Sha256, Digest};

// #[derive(Debug, Clone)]
// struct MerkleTree {
//     root: FieldElement,
//     tree: Vec<FieldElement>,
//     height: usize,
//     facts: HashSet<FieldElement>,
// }


// // impl MerkleTree {

// //     pub fn new(data: Vec<FieldElement>) -> MerkleTree {
// //         assert!(data.len() > 0);

// //         // calculate number of leaves (next power of 2 from length of data)
// //         let num_leaves: usize = ((data.len() as f64).log2().ceil() as f64).exp2() as usize;

// //         // extend data with zeros to make it a power of 2
// //         let mut data = data.clone();  
// //         data.extend((0..num_leaves - data.len()).map(|_| FieldElement::zero()));

// //         // get height
// //         let height = (num_leaves as f64).log2().ceil() as usize;

// //         // facts set for leaf data
// //         let facts = HashSet::new();

// //         // build tree
// //         let root = MerkleTree::build_tree();

// //         MerkleTree {
// //             root, height, facts,
// //         }
// //     }


// //     pub fn build_tree() -> Vec<FieldElement> {

        
// //     }
// // }





// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::modules::field::{*};

//     #[test]
//     fn test_merkle_tree() {


//         let data = vec![
//             FieldElement::new(1), FieldElement::new(2), FieldElement::new(3), FieldElement::new(4),
//             FieldElement::new(1), FieldElement::new(2), FieldElement::new(3)];
        
        
        
//         let tree = MerkleTree::new(data);

//         // println!("{:?}", tree);
//     }
// }
