use crate::modules::field::FieldElement;
use num_bigint::BigInt;

/// Number Theoretic Transform (NTT) implementation for efficient polynomial operations
/// This provides O(n log n) interpolation and evaluation instead of O(n²) Lagrange
pub struct NTT {
    pub n: usize,
    pub omega: FieldElement,  // primitive nth root of unity
    pub omega_inv: FieldElement,  // inverse of omega
}

impl NTT {
    /// Create a new NTT instance for given size n (must be power of 2)
    pub fn new(n: usize) -> Self {
        assert!(n.is_power_of_two(), "NTT size must be power of 2");
        assert!(n <= (1_usize << 30), "NTT size too large for practical use"); // Reasonable limit
        
        let omega = FieldElement::primitive_nth_root(n as u128);
        let omega_inv = omega.inverse();
        
        NTT { n, omega, omega_inv }
    }
    
    /// Forward NTT: polynomial coefficients -> evaluations at roots of unity
    pub fn forward(&self, coeffs: &mut Vec<FieldElement>) {
        assert_eq!(coeffs.len(), self.n, "Input size must match NTT size");
        self.ntt_recursive(coeffs, false);
    }
    
    /// Inverse NTT: evaluations at roots of unity -> polynomial coefficients  
    pub fn inverse(&self, evals: &mut Vec<FieldElement>) {
        assert_eq!(evals.len(), self.n, "Input size must match NTT size");
        self.ntt_recursive(evals, true);
        
        // Scale by 1/n
        let n_inv = FieldElement::new(BigInt::from(self.n)).inverse();
        for coeff in evals.iter_mut() {
            *coeff = coeff.clone() * n_inv.clone();
        }
    }
    
    /// Recursive NTT implementation using Cooley-Tukey algorithm
    fn ntt_recursive(&self, a: &mut Vec<FieldElement>, inverse: bool) {
        let n = a.len();
        if n <= 1 { return; }
        
        // Bit-reverse permutation
        self.bit_reverse(a);
        
        // Bottom-up iterative NTT
        let mut len = 2;
        while len <= n {
            let w = if inverse {
                self.omega_inv.pow((self.n / len) as u128)
            } else {
                self.omega.pow((self.n / len) as u128)
            };
            
            for i in (0..n).step_by(len) {
                let mut wn = FieldElement::one();
                for j in 0..(len / 2) {
                    let u = a[i + j].clone();
                    let v = a[i + j + len / 2].clone() * wn.clone();
                    a[i + j] = u.clone() + v.clone();
                    a[i + j + len / 2] = u - v;
                    wn = wn * w.clone();
                }
            }
            len *= 2;
        }
    }
    
    /// Bit-reverse permutation for NTT
    fn bit_reverse(&self, a: &mut Vec<FieldElement>) {
        let n = a.len();
        let mut j = 0;
        for i in 1..n {
            let mut bit = n >> 1;
            while j & bit != 0 {
                j ^= bit;
                bit >>= 1;
            }
            j ^= bit;
            if i < j {
                a.swap(i, j);
            }
        }
    }
    
    /// Fast polynomial interpolation using NTT
    /// Given evaluations at consecutive powers of omega, returns polynomial coefficients
    pub fn interpolate(&self, evaluations: Vec<FieldElement>) -> Vec<FieldElement> {
        assert_eq!(evaluations.len(), self.n, "Evaluations size must match NTT size");
        let mut coeffs = evaluations;
        self.inverse(&mut coeffs);
        coeffs
    }
    
    /// Fast polynomial evaluation using NTT  
    /// Given polynomial coefficients, returns evaluations at consecutive powers of omega
    pub fn evaluate(&self, coefficients: Vec<FieldElement>) -> Vec<FieldElement> {
        let mut padded_coeffs = coefficients;
        // Pad with zeros if needed
        padded_coeffs.resize(self.n, FieldElement::zero());
        self.forward(&mut padded_coeffs);
        padded_coeffs
    }
    
    /// Fast polynomial multiplication using NTT
    pub fn multiply(&self, a: &Vec<FieldElement>, b: &Vec<FieldElement>) -> Vec<FieldElement> {
        let result_size = a.len() + b.len() - 1;
        let ntt_size = result_size.next_power_of_two();
        let ntt = NTT::new(ntt_size);
        
        let mut a_padded = a.clone();
        let mut b_padded = b.clone();
        a_padded.resize(ntt_size, FieldElement::zero());
        b_padded.resize(ntt_size, FieldElement::zero());
        
        // Transform to evaluation form
        ntt.forward(&mut a_padded);
        ntt.forward(&mut b_padded);
        
        // Pointwise multiplication
        for i in 0..ntt_size {
            a_padded[i] = a_padded[i].clone() * b_padded[i].clone();
        }
        
        // Transform back to coefficient form
        ntt.inverse(&mut a_padded);
        
        // Trim to actual result size
        a_padded.truncate(result_size);
        a_padded
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::BigInt;

    #[test]
    fn test_ntt_basic() {
        let ntt = NTT::new(4);
        
        // Test polynomial 1 + 2x + 3x^2 + 4x^3
        let coeffs = vec![
            FieldElement::new(BigInt::from(1)),
            FieldElement::new(BigInt::from(2)), 
            FieldElement::new(BigInt::from(3)),
            FieldElement::new(BigInt::from(4))
        ];
        
        // Forward transform then inverse should give back original
        let evals = ntt.evaluate(coeffs.clone());
        let recovered = ntt.interpolate(evals);
        
        for i in 0..4 {
            assert_eq!(coeffs[i], recovered[i], "NTT round-trip failed at index {}", i);
        }
    }
    
    #[test]
    fn test_ntt_interpolation() {
        let ntt = NTT::new(8);
        
        // Create some test evaluations
        let domain: Vec<FieldElement> = (0..8)
            .map(|i| ntt.omega.pow(i as u128))
            .collect();
            
        let evaluations: Vec<FieldElement> = (0..8)
            .map(|i| FieldElement::new(BigInt::from(i * i + 1)))  // f(x) = x^2 + 1
            .collect();
        
        let coeffs = ntt.interpolate(evaluations.clone());
        let recovered_evals = ntt.evaluate(coeffs);
        
        for i in 0..8 {
            assert_eq!(evaluations[i], recovered_evals[i], 
                      "Interpolation test failed at index {}", i);
        }
    }
    
    #[test]
    fn test_ntt_multiplication() {
        let ntt = NTT::new(8);
        
        // Test (1 + x) * (1 + 2x) = 1 + 3x + 2x^2
        let a = vec![
            FieldElement::new(BigInt::from(1)),
            FieldElement::new(BigInt::from(1))
        ];
        let b = vec![
            FieldElement::new(BigInt::from(1)),
            FieldElement::new(BigInt::from(2))
        ];
        
        let result = ntt.multiply(&a, &b);
        
        assert_eq!(result.len(), 3);
        assert_eq!(result[0], FieldElement::new(BigInt::from(1))); // constant term
        assert_eq!(result[1], FieldElement::new(BigInt::from(3))); // x term
        assert_eq!(result[2], FieldElement::new(BigInt::from(2))); // x^2 term
    }
}