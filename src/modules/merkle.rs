extern crate blake2;
extern crate generic_array;
extern crate typenum;
use blake2::{Blake2b, Digest};
use generic_array::typenum::U32;
use generic_array::GenericArray;

// define the hash output size
type OutputSize = U32; // 32 bytes (256 bits) output
pub type HashOutput = GenericArray<u8, OutputSize>; // Fixed-size hash output

pub struct Merkle;

impl Merkle {

    // Commit to a list of leafs recursively
    fn commit_leafs(leafs: &[HashOutput]) -> HashOutput {
        assert!(leafs.len().is_power_of_two(), "length must be power of two");

        // Return root
        if leafs.len() == 1 { return leafs[0].clone(); } 
        else {
            //  split leaves in two, commit to each halfS
            let mid = leafs.len() / 2;  
            let left = Self::commit_leafs(&leafs[..mid]); //
            let right = Self::commit_leafs(&leafs[mid..]);

            // Hash the two roots
            let mut hasher = Blake2b::<OutputSize>::new();
            hasher.update(left);
            hasher.update(right);
            hasher.finalize()
        }
    }


    // Commit to a list of data elements represented as byte vectors
    pub fn commit(data_array: &Vec<Vec<u8>>) -> HashOutput {

        // Hash each data element
        let leafs: Vec<HashOutput> = data_array.iter()
            .map(|da| {
                let mut hasher = Blake2b::<OutputSize>::new();
                hasher.update(da);
                hasher.finalize()
            })
            .collect();

        // return the computed root hash
        Self::commit_leafs(&leafs)
    }

    // Open a leaf at a given index, recursively collects authentication path
    fn open_(index: usize, leafs: &[HashOutput]) -> Vec<HashOutput> {
        assert!(leafs.len().is_power_of_two(), "length must be power of two");
        assert!(index < leafs.len(), "cannot open invalid index");  

        let mid = leafs.len() / 2;

        // Return the sibling of the leaf
        if leafs.len() == 2 {  
            vec![leafs[1 - index]] 
        } 
        else if index < mid {
                // open the left half, push the right sibling
                let mut result = Self::open_(index, &leafs[..mid]);
                result.push(Self::commit_leafs(&leafs[mid..]));
                result
        }
        else {
                // open the right half, push the left sibling
                let mut result = Self::open_(index - mid, &leafs[mid..]);
                result.push(Self::commit_leafs(&leafs[..mid]));
                result
        }
    }

    // Open a leaf at a given index, returns authentication path
    pub fn open(index: usize, data_array: &Vec<Vec<u8>>)  -> Vec<HashOutput> {

        // Hash each data element
        let leafs: Vec<HashOutput> = data_array.iter()
            .map(|da| {
                let mut hasher = Blake2b::<OutputSize>::new();
                hasher.update(da);
                hasher.finalize()
            })
            .collect();


        // return the computed authentication path
        Self::open_(index, &leafs)
    }

    // Recurisvely rebuilds the root hash from authentication path
    fn verify_(root: &HashOutput, index: usize, path: &[HashOutput], leaf: &HashOutput) -> bool {
        assert!(index < (1 << path.len()), "cannot verify invalid index");

        // verify root hash
        if path.len() == 1 {
            let mut hasher = Blake2b::<OutputSize>::new();
            if index == 0 {
                hasher.update(leaf);
                hasher.update(&path[0]);
            } else {
                hasher.update(&path[0]);
                hasher.update(leaf);
            }   
            // ensure provided root hash is correct compared to computed hash
            root == &hasher.finalize()
        } else { 
            // recursively rebuild the root hash from authentication path
            let mut hasher = Blake2b::<OutputSize>::new();
            if index % 2 == 0 {
                hasher.update(leaf);
                hasher.update(&path[0]);
            } else {
                hasher.update(&path[0]);
                hasher.update(leaf);
            }
            // verify the next level of the tree
            Self::verify_(root, index >> 1, &path[1..], &hasher.finalize())
        }
    }

    // Verify a leaf at a given index, via the authentication path
    pub fn verify(root: &HashOutput, index: usize, path: &[HashOutput], data_element: &Vec<u8>) -> bool {
        let mut hasher = Blake2b::<OutputSize>::new();
        hasher.update(data_element);
        let digest = hasher.finalize();
        Self::verify_(root, index, path, &digest)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, thread_rng};
    use rand::distributions::Standard;

    // Helper function to generate random Vec<u8> of random length
    fn random_data_vec_u8() -> Vec<u8> {
        let mut rng = thread_rng();
        let size = rng.gen_range(1..100); // Random size between 1 and 99
        rng.sample_iter(&Standard).take(size).collect()
    }

    // Helper function to generate random GenericArray<u8, OutputSize>
    fn random_data_gen_arr_u8() -> GenericArray<u8, OutputSize> {
        let mut rng = thread_rng();
        let data: Vec<u8> = rng.sample_iter(&Standard).take(32).collect();
        GenericArray::clone_from_slice(&data)
    }

    #[test]
    fn test_merkle() {
        let n = 64;
        let leafs: Vec<Vec<u8>> = (0..n).map(|_| random_data_vec_u8()).collect();
        let root = Merkle::commit(&leafs);

    //     // Test that opening and verifying each leaf works correctly
        for i in 0..n {
            let path = Merkle::open(i, &leafs);
            assert!(Merkle::verify(&root, i, &path, &leafs[i]), "Verification failed for correct leaf");
        }

        // Test that opening with wrong data fails
        for i in 0..n {
            let path = Merkle::open(i, &leafs);
            let wrong_data = random_data_vec_u8();
            assert!(!Merkle::verify(&root, i, &path, &wrong_data), "Verification succeeded for wrong data");
        }

        // Test that using a wrong index fails
        for i in 0..n {
            let path = Merkle::open(i, &leafs);
            let j = (i + 1) % n; // Use a different index
            assert!(!Merkle::verify(&root, j, &path, &leafs[i]), "Verification succeeded for wrong index");
        }

        // Test with a false root
        for i in 0..n {
            let path = Merkle::open(i, &leafs);
            let fake_root = random_data_gen_arr_u8();
            assert!(!Merkle::verify(&fake_root, i, &path, &leafs[i]), "Verification succeeded with a false root");
        }

        // Test each path element with a fake value
        for i in 0..n {
            let mut path = Merkle::open(i, &leafs);
            for j in 0..path.len() {
                let original = path[j].clone();
                path[j] = random_data_gen_arr_u8();
                assert!(!Merkle::verify(&root, i, &path, &leafs[i]), "Verification succeeded with a tampered path");
                path[j] = original; // Restore original path for next iteration
            }
        }

        // Test using a different root
        let fake_leafs: Vec<Vec<u8>> = (0..n).map(|_| random_data_vec_u8()).collect();
        let fake_root = Merkle::commit(&fake_leafs);
        for i in 0..n {
            let path = Merkle::open(i, &leafs);
            assert!(!Merkle::verify(&fake_root, i, &path, &leafs[i]), "Verification succeeded with a different root");
        }
    }
}
