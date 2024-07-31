
use core::panic;
use std::ops::{Add, Sub, Mul, Div, Neg};
use std::fmt;
use num_bigint::RandBigInt;
use num_bigint::BigInt;
use num_traits::{Zero, One};
use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct FieldElement {
    pub value: BigInt,
}

impl FieldElement {
    // constructor for field element with value within prime modulus
    pub fn new(value: BigInt) -> FieldElement {
        FieldElement { value: value % FieldElement::modulus() }
    }

    // generator constructor
    pub fn generator() -> FieldElement { FieldElement::new(BigInt::from(85408008396924667383611388730472331217 as i128)) }

    // prime modulus of 1 + 407 * 2^119
    pub fn modulus() -> BigInt { BigInt::from(407) * (BigInt::from(2).pow(119)) + 1 }

    // zero and one constructors
    pub fn zero() -> FieldElement { FieldElement::new(BigInt::zero()) }
    pub fn one() -> FieldElement { FieldElement::new(BigInt::one()) }

    // check if field element is zero
    pub fn is_zero(&self) -> bool { self.value.is_zero() }

    // exponentiation
    pub fn pow(&self, exp: u128) -> FieldElement {
        let mut base = self.clone();
        let mut exp = BigInt::from(exp);
        let mut result = FieldElement::one();

        //  modular exponentiation 
        while !exp.is_zero() {
            if &exp % 2 == BigInt::from(1) {
                result = result * base.clone();
            }
            base = base.clone() * base;
            exp /= 2;
        }
        result
    }

    // random field element
    pub fn random(&self) -> FieldElement {
        let mut rng = rand::thread_rng();
        let upper_bound = FieldElement::modulus();
        let random_value = rng.gen_bigint_range(&BigInt::zero(), &upper_bound);
        FieldElement::new(random_value)
    }

    // multiplicative inverse using the extended Euclidean algorithm
    pub fn inverse(&self) -> FieldElement {
        let zero = BigInt::zero();
        let one = BigInt::one();
        let mut t = zero.clone();
        let mut new_t = one.clone();
        let mut r = FieldElement::modulus();
        let mut new_r = self.value.clone();

        while new_r != zero {
            let quotient = &r / &new_r;

            let temp_t = t.clone();
            t = new_t.clone();
            new_t = temp_t - &quotient * &new_t;

            let temp_r = r.clone();
            r = new_r.clone();
            new_r = temp_r - quotient * new_r;
        }

        if r > one {
            panic!("Element does not have an inverse");
        }

        // Adjust negative result to be positive
        if t < zero {
            t = t + FieldElement::modulus();
        }

        FieldElement::new(t)
    }

    // primitive n'th root of unity c^n = 1
    pub fn primitive_nth_root(n: i128) -> FieldElement {
        assert!(n <= 1 << 119 && (n & (n-1)) == 0, "Field does not have nth root of unity where n > 2^119 or not power of two.");

        // accumulate root of unity by squaring generator
        let mut root = FieldElement::generator();
        let mut order: i128 = 1 << 119;
        while order != n {
            root = root.pow(2);
            order = order / 2;
        }
        root
    }

    // sample random field element
    pub fn sample(byte_array: Vec<u8>) -> FieldElement {
        let mut acc = BigInt::zero();
        for b in byte_array {
            acc = (acc << 8) ^ BigInt::from(b);
        }   
        FieldElement::new(acc % FieldElement::modulus())
    }
}

// Note: modulo applied in FieldElement constructor
impl Add for FieldElement {
    type Output = Self;
    fn add(self, other: Self) -> Self { FieldElement::new(self.value + other.value) }
}

impl Sub for FieldElement {
    type Output = Self;

    fn sub(self, other: Self) -> Self {
        FieldElement::new(  
            // Ensure the result is positive by adding the modulus before taking the modulo
            (self.value - other.value + FieldElement::modulus()) % FieldElement::modulus() 
        )
    }
}

impl Mul for FieldElement {
    type Output = Self;
    fn mul(self, other: Self) -> Self { FieldElement::new(self.value * other.value) }
}

impl Div for FieldElement {
    type Output = Self;
    fn div(self, other: Self) -> Self { self * other.inverse() }
}

impl Neg for FieldElement {
    type Output = Self;
    fn neg(self) -> Self { FieldElement::new(-self.value) }
}

impl fmt::Display for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use num_bigint::ToBigInt;
    use rand::Rng;

    #[test]
    fn test_addition() {
        let a = FieldElement::new(FieldElement::modulus() - 1.to_bigint().unwrap());
        let b = FieldElement::new(10.to_bigint().unwrap());
        assert_eq!(a + b, FieldElement::new(9.to_bigint().unwrap()));
    }

    #[test]
    fn test_subtraction() {
        let a = FieldElement::new(FieldElement::modulus() - 1.to_bigint().unwrap());
        let b = FieldElement::new(10.to_bigint().unwrap());
        let c = a - b;

        println!("{:?}", c);

        assert_eq!(c, FieldElement::new((FieldElement::modulus() - 11.to_bigint().unwrap()) % FieldElement::modulus()));
    }

    #[test]
    fn test_inverse_1() {
        let elem = FieldElement::new(3.to_bigint().unwrap());
        let inv_elem = elem.inverse();
        let one_elem = FieldElement::new(1.to_bigint().unwrap());
        assert_eq!(elem * inv_elem, one_elem, "Inverse test failed");
    }

    #[test]
    fn test_inverse_2() {
        let test_values: Vec<i128> = vec![1, 2, 3, 5, 1234567, 3221225470];
        for val in test_values {
            let elem = FieldElement::new(val.to_bigint().unwrap());
            let inv = elem.inverse();
            let product = elem * inv;

            assert_eq!(product, FieldElement::one(), "Failed inverse test for value: {}", val);
        }

        let zero_elem = FieldElement::zero();
        let result = std::panic::catch_unwind(|| zero_elem.inverse());
        assert!(result.is_err(), "Inverse of zero did not panic as expected");
    }

    #[test]
    fn test_pow() {
        let a = FieldElement::new(2.to_bigint().unwrap());
        let expected_value = 2.to_bigint().unwrap().pow(32) % FieldElement::modulus();
        assert_eq!(a.pow(32).value, expected_value);
    }

    #[test]
    fn test_negative_handling() {
        let elem1 = FieldElement::new(2.to_bigint().unwrap());
        let elem2 = FieldElement::new(3.to_bigint().unwrap());
        let result = elem1 - elem2;
        let expected_value = (FieldElement::modulus() - 1.to_bigint().unwrap()) % FieldElement::modulus();
        assert_eq!(result.value, expected_value);
    }

    #[test]
    fn test_primitive_nth_root() {
        let n = 16;
        let root = FieldElement::primitive_nth_root(n);
        assert_eq!(root.pow(n as u128), FieldElement::one());
    }

    #[test]
    fn test_sample() {

        // create vec of 32 random bytes
        let mut rng = rand::thread_rng();
        let byte_array: Vec<u8> = (0..32).map(|_| rng.gen()).collect();
        let elem = FieldElement::sample(byte_array);

        // check if element is within modulus
        assert!(elem.value < FieldElement::modulus());
    }
}

