use std::ops::{Add, Sub, Mul, Div, Neg};
use std::fmt;


#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FieldElement{

    value: i64,
    modulus: i64,
}

// creates new field element at value w/ prime mod 3 * 2^30 + 1
pub fn new_field_element(input_value: i64) -> FieldElement{
    let modulus: i64 = 3 * 2_i64.pow(30) + 1; 
    let value: i64 = input_value % modulus;
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

        let mut value: i64 = self.value - other.value;

        if value < 0 {
            value = self.modulus - value; // if negative, add modulus
        }

        FieldElement{
            value: self.value % self.modulus, 
            modulus: self.modulus
        }
    }
}

// negation
impl Neg for FieldElement {
    type Output = FieldElement;

    // negate then mod by prime
    fn neg(self) -> FieldElement {
        FieldElement{
            value: zero().value - self.value, 
            modulus: self.modulus
        }
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

// field element inverse func
fn inverse(element: FieldElement) -> FieldElement {
    let mut t: i64 = 0;
    let mut new_t: i64 = 1;

    let mut r: i64 = element.modulus;
    let mut new_r: i64 = element.value;

    while new_r != 0 {

        let quotient: i64 = r / new_r;

        t = new_t;
        new_t = t - (quotient * new_t);

        r = new_r;
        new_r = r - (quotient * new_r);
    }

    if r != 1 {
        panic!("{} is not invertible", element);
    }

    FieldElement{
        value: t % element.modulus,
        modulus: element.modulus
    }
}


// prints field element value
impl fmt::Display for FieldElement {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{}", self.value)
    }
}

