use crate::modules::field::{*};
use crate::modules::univariate_poly::{*};
use std::collections::HashMap;
use num_bigint::RandBigInt;
use num_bigint::BigInt;
use num_traits::{Zero, One};
use std::ops::{Add, Sub, Mul, Neg};

#[derive(Debug, Clone, PartialEq, Eq)]
struct MPolynomial {
    /**
        Multivariate polynomials are represented as dictionaries with exponents
        vectors as keys and coefficients as values.

        f(x,y,z) = 17 + 2xy + 42z - 19x^6*y^3*z^12 is represented as:
        {
            (0,0,0) => 17,
            (1,1,0) => 2,
            (0,0,1) => 42,
            (6,3,12) => -19,
        }
     */

    pub dict: HashMap<Vec<i128>, FieldElement>,
}


impl MPolynomial {

    pub fn zero() -> MPolynomial {
        MPolynomial { dict: HashMap::new() }
    }

    pub fn is_zero(&self) -> bool {
        self.dict.values().all(FieldElement::is_zero)
    }

    // exponentiation by squaring
    pub fn pow(self, exponent: usize) -> Self {
        if self.is_zero() { return MPolynomial::zero(); } // 0^n = 0

        // FieldElement has a method `one()` that returns the multiplicative identity element
        let field_one = FieldElement::one();
        let num_variables = self.dict.keys().next().unwrap().len(); // Assuming dictionary is not empty
        let mut acc = MPolynomial { dict: HashMap::from([(vec![0; num_variables], field_one)]) };

        for b in format!("{:b}", exponent).chars() {
            acc = acc.clone() * acc; // Square the accumulator
            if b == '1' {
                acc = acc.clone() * self.clone(); // Multiply by self if the bit is 1
            }
        }

        acc
    }


}


impl Add for MPolynomial {
    type Output = MPolynomial;


    fn add(self, rhs: Self) -> Self::Output {

        // result dict fomr lhs
        let mut result_dict = self.dict.clone();

        // Iterate over the terms in the rhs side poly
        for (exponents, coeff) in rhs.dict {

            // if the term exists in the lhs dict, add the coefficients
            if let Some(current_coeff) = result_dict.get_mut(&exponents) {
                *current_coeff = current_coeff.clone() + coeff;

            } else { // otherwise, insert the new term
                result_dict.insert(exponents, coeff);
            }
        }

        MPolynomial { dict: result_dict }
    }
}

impl Sub for MPolynomial {
    type Output = MPolynomial;

    fn sub(self, other: MPolynomial) -> Self::Output { 
        self + (-other) // add negated rhs
    }
}


impl Neg for MPolynomial {
    type Output = MPolynomial;

    fn neg(self) -> Self::Output {

        // Cnegate all coeffs in clone of dict
        let mut result_dict = self.dict.clone();        
        for (_exponents, coeff) in result_dict.iter_mut() {
            *coeff = -coeff.clone();
        }

        MPolynomial { dict: result_dict }
    }
}


impl Mul for MPolynomial {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        let mut result_dict: HashMap<Vec<i128>, FieldElement> = HashMap::new();

        // Iterate over each term in the left polynomial
        for (exponents_left, coeff_left) in &self.dict {
            // Iterate over each term in the right polynomial
            for (exponents_right, coeff_right) in &rhs.dict {

                // Calculate new exponents by adding corresponding exponents
                let mut new_exponents = exponents_left.clone();
                for (index, exponent) in exponents_right.iter().enumerate() {
                
                    if index < new_exponents.len() {
                        new_exponents[index] += exponent;
                    } else {
                        // Handle the case where the number of variables in both polynomials differs
                        new_exponents.push(*exponent);
                    }
                }

                // Calculate the new coefficient
                let new_coeff = coeff_left.clone() * coeff_right.clone();

                // Add or update the term in the result dictionary
                if let Some(existing_coeff) = result_dict.get_mut(&new_exponents) {
                    *existing_coeff = existing_coeff.clone() + new_coeff;
                } else {
                    result_dict.insert(new_exponents, new_coeff);
                }
            }
        }

        MPolynomial { dict: result_dict }
    }
}