use std::ops::{Add, Sub, Mul, Div, Neg};
use std::fmt;
use rand::Rng;


#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FieldElement{
    pub value: i128,
    modulus: i128,
}

// creates new field element at value w/ prime mod 3 * 2^30 + 1
pub fn new_field_element(input_value: i128) -> FieldElement{
    let modulus: i128 = 3 * 2_i128.pow(30) + 1; 
    let value: i128 = input_value % modulus;
    FieldElement{value, modulus}
}
    
// generator constructor
pub fn generator() -> FieldElement{
    new_field_element(5)
}

// zero constructor
pub fn zero() -> FieldElement{
    new_field_element(0)
}

// one constructor
pub fn one() -> FieldElement{
    new_field_element(1)
}

// addition
impl Add for FieldElement {
    type Output = FieldElement;

    // add then mod by prime
    fn add(self, other: FieldElement) -> FieldElement {
        FieldElement{
            value: (self.value + other.value) % self.modulus, 
            modulus: self.modulus
        }
    }
}

// subtraction
impl Sub for FieldElement {
    type Output = FieldElement;

    // subtract then mod by prime
    fn sub(self, other: FieldElement) -> FieldElement {
        FieldElement{
            value: (self.value - other.value) % self.modulus, 
            modulus: self.modulus
        }
    }
}

// negation
impl Neg for FieldElement {
    type Output = FieldElement;

    fn neg(self) -> FieldElement {
        let mut result = self.modulus - self.value;
        if result >= self.modulus {
            result -= self.modulus;
        }
        FieldElement { value: result, modulus: self.modulus }
    }
}

// multiplication
impl Mul for FieldElement {
    type Output = FieldElement;

    // multiply then mod by prime
    fn mul(self, other: FieldElement) -> FieldElement {
        FieldElement{
            value: (self.value * other.value) % self.modulus, 
            modulus: self.modulus
        }
    }
}

// division
impl Div for FieldElement {
    type Output = FieldElement;

    // multiply by inverse then mod by prime
    fn div(self, other: FieldElement) -> FieldElement {
        self * inverse(other)
    }
}

// exponentiation
pub fn power(base: FieldElement, exp: i128) -> FieldElement {

    let mut base = base.clone();
    let mut exp = exp.clone();
    let mut result = one(); 

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

// multiplicative inverse
pub fn inverse(element: FieldElement) -> FieldElement {
    let mut t = 0;
    let mut new_t = 1;
    let mut r = element.modulus;
    let mut new_r = element.value;

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
        t += element.modulus;
    }

    FieldElement {
        value: t,
        modulus: element.modulus
    }
}


// random element
pub fn random() -> FieldElement {
    let mut rng = rand::thread_rng();
    new_field_element(rng.gen_range(0..3221225472))
}


// prints field element value
impl fmt::Display for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Simple modulo to ensure it's within range and non-negative
        let value = (self.value % self.modulus + self.modulus) % self.modulus;
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
        let a = new_field_element(3221225472);
        let b = new_field_element(10);
        assert_eq!(a + b, new_field_element(9));
    }

    #[test]
    fn subtraction() {

        let a = new_field_element(3221225472);
        let b = new_field_element(10);
        let c = a - b;
        println!("{:?}", c);

        assert_eq!(3221225472 - 10, c.value);
    }

    #[test]
    fn inverse_test() { 

        // get rand int within mod
        let mut rng = rand::thread_rng();
        let random_number = rng.gen_range(1..3221225472); 
    
        // ensure multiplicative inverse works
        let a = -new_field_element(random_number);
        let c  = inverse(a);
        assert_eq!(c * a, one());
    }

    #[test]
    fn test_pow() {
        let a = new_field_element(2);
        assert_eq!(power(a, 32).value, 2_i128.pow(32) % 3221225473);
    }

}