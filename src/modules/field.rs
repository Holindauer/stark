use std::ops::{Add, Sub, Mul, Div, Neg};
use std::fmt;
use rand::Rng;
use serde::{Serialize, Deserialize};

/**
 * @notice field.rs contains a finite field implementation for the modulus 3 * 2^30 + 1
 * It contains a generate with value 5. 
 */


#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct FieldElement{
    pub value: i128,
}

impl FieldElement {

    // constructor for field element at value w/ prime mod 3 * 2^30 + 1
    pub fn new(input_value: i128) -> FieldElement{
        let value: i128 = input_value % FieldElement::modulus();
        FieldElement{value}
    }

    // generator constructor
    pub fn generator() -> FieldElement{ FieldElement::new(5) }

    // modulus 
    pub fn modulus() -> i128 { 3 * 2_i128.pow(30) + 1 }

    // zero and one constructor 
    pub fn zero() -> FieldElement{ FieldElement::new(0) }
    pub fn one() -> FieldElement{ FieldElement::new(1) }

    // exponentiation
    pub fn pow(&self, exp: i128) -> FieldElement {

        let mut base = self.clone();
        let mut exp = exp.clone();
        let mut result = FieldElement::one(); 

        // exponentiation by repeated squaring
        while exp > 0 {
            if exp % 2 == 1 {
                result = result * base;
            }
            base = base * base;
            exp /= 2;
        }
        result 
    }

    // random field element
    pub fn random(&self) -> FieldElement {
        let mut rng = rand::thread_rng();
        FieldElement::new(rng.gen_range(0..3221225472))
    }

    // multiplicative inverse
    pub fn inverse(&self) -> FieldElement {
        let mut t = 0;
        let mut new_t = 1;
        let mut r = FieldElement::modulus();
        let mut new_r = self.value;
        
        while new_r != 0 {
            let quotient = r / new_r;
            (t, new_t) = (new_t, t - quotient * new_t);
            (r, new_r) = (new_r, r - quotient * new_r);
        }

        if r > 1 { // This means `element.value` and `element.modulus` are not coprime
            panic!("Element does not have an inverse");
        }

        // Ensure the result is positive
        if t < 0 {
            t += FieldElement::modulus();
        }

        FieldElement { value: t }
    }
}

// addition
impl Add for FieldElement {
    type Output = FieldElement;

    // add then mod by prime
    fn add(self, other: FieldElement) -> FieldElement {
        FieldElement{
            value: (self.value + other.value) % FieldElement::modulus(),}
    }
}

// subtraction
impl Sub for FieldElement {
    type Output = FieldElement;

    // subtract then mod by prime
    fn sub(self, other: FieldElement) -> FieldElement {
        let mut result = (self.value - other.value) % FieldElement::modulus();
        if result < 0 {
            result += FieldElement::modulus();
        }
        FieldElement { value: result  }
    }
}

// negation
impl Neg for FieldElement {
    type Output = FieldElement;

    fn neg(self) -> FieldElement {
        let mut result = FieldElement::modulus() - self.value;
        if result >= FieldElement::modulus() {
            result -= FieldElement::modulus();
        }
        FieldElement { value: result }
    }
}

// multiplication
impl Mul for FieldElement {
    type Output = FieldElement;

    // multiply then mod by prime
    fn mul(self, other: FieldElement) -> FieldElement {
        FieldElement{ value: (self.value * other.value) % FieldElement::modulus() }
    }
}

// division
impl Div for FieldElement {
    type Output = FieldElement;

    // multiply by inverse then mod by prime
    fn div(self, other: FieldElement) -> FieldElement {
        self * other.inverse()
    }
}

// prints field element value
impl fmt::Display for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Simple modulo to ensure it's within range and non-negative
        let value = (self.value % FieldElement::modulus() + FieldElement::modulus()) % FieldElement::modulus();
        write!(f, "{}", value)
    }
}

// tests
#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn addition() {
        let a = FieldElement::new(3221225472);
        let b = FieldElement::new(10);
        assert_eq!(a + b, FieldElement::new(9));
    }

    #[test]
    fn subtraction() {

        let a = FieldElement::new(3221225472);
        let b = FieldElement::new(10);
        let c = a - b;
        println!("{:?}", c);

        assert_eq!(3221225472 - 10, c.value);
    }

    #[test]
    fn inverse_test_fuzz() { 

        // get rand int within mod
        let mut rng = rand::thread_rng();
        let random_number = rng.gen_range(1..3221225472); 
    
        // ensure multiplicative inverse works
        let a = -FieldElement::new(random_number);
        let c  = a.inverse();
        assert_eq!(c * a, FieldElement::one());
    }

    #[test]
    fn invers_test_2() {
        // Test with various non-zero elements
        let test_values = vec![1, 2, 3, 5, 1234567, 3221225470]; // Various values, including near the modulus
        for val in test_values {
            let elem = FieldElement::new(val);
            let inv = elem.inverse();
            let product = elem * inv;

            // Check if the product is 1 (multiplicative identity)
            assert_eq!(product, FieldElement::one(), "Failed inverse test for value: {}", val);
        }

        // Optional: Test edge case for zero, which should not have an inverse
        // Depending on how you want to handle this, this part might be different
        let zero_elem = FieldElement::zero();
        let result = std::panic::catch_unwind(|| zero_elem.inverse());
        assert!(result.is_err(), "Inverse of zero did not panic as expected");
    }

    #[test]
    fn test_pow() {
        let a = FieldElement::new(2);
        assert_eq!(a.pow(32).value, 2_i128.pow(32) % 3221225473);
    }


    #[test]
    fn test_negative_handling() {
        let elem1 = FieldElement::new(2);
        let elem2 = FieldElement::new(3);
        assert_eq!((elem1 - elem2).value, FieldElement::modulus() - 1);
    }

}