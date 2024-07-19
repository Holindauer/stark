use crate::modules::field::{*};
use crate::modules::univariate_poly::{*};
use std::collections::HashMap;
use num_bigint::RandBigInt;
use num_bigint::BigInt;
use num_traits::{Zero, One};
use num_traits::ToPrimitive;
use std::ops::{Add, Sub, Mul, Neg};

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
#[derive(Debug, Clone, PartialEq, Eq)]
struct MPolynomial {
    pub dict: HashMap<Vec<i128>, FieldElement>,
}

impl MPolynomial {

    // constructors
    pub fn new(dict: HashMap<Vec<i128>, FieldElement>) -> MPolynomial { MPolynomial { dict } }
    pub fn zero() -> MPolynomial { MPolynomial { dict: HashMap::new() } }
    pub fn constant(element: i128) -> MPolynomial { 
        MPolynomial { dict: HashMap::from([(vec![0], FieldElement::new(BigInt::from(element)))]) } 
    }

    // checks
    pub fn is_zero(&self) -> bool { self.dict.values().all(FieldElement::is_zero) }

    /// Evaluates the polynomial at a given point (a vector of field elements)
    /// by accumulating the product of each term's coefficient and the point
    pub fn eval(&self, point: &[FieldElement]) -> FieldElement {
        let mut acc = FieldElement::zero();  
        for (exponents, coeff) in &self.dict {
            let mut prod = coeff.clone();


            for (i, &exponent) in exponents.iter().enumerate() {
                // handle cases where point does not cover all variables
                if i >= point.len() { continue; }
                
                let mut term = if exponent == 0 {
                    FieldElement::one() // c^0 = 1
                } else {   

                    // calculate the term by raising the point to the exponent
                    let base = &point[i];
                    let mut term = base.clone();
                    for _ in 1..exponent {
                        term = term.clone() * base.clone();
                    }
                    term
                };

                // multiply the term by the accumulated product
                prod = prod.clone() * term;
            }
            // add the product to the accumulator
            acc = acc.clone() + prod;
        }
        acc
}

    // efficient exponentiation by squaring
    pub fn pow(self, exponent: usize) -> Self {
        if self.is_zero() { return MPolynomial::zero(); } // 0^n = 0
        if exponent == 0 { return MPolynomial::constant(1); } // c^0 = 1

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

    /// Returns a vector of multivariate polynomials representing each indeterminate's linear 
    /// function with a leading coefficient of one. For example, for three indeterminates, it 
    /// returns: [f(x,y,z) = x, f(x,y,z) = y, f(x,y,z) = z]
    pub fn variables(num_variables: usize) -> Vec<Self> {
        let mut variables = Vec::new();
        for i in 0..num_variables {

            // Create a dictionary with a single term for each variable
            let mut exponent = vec![0; num_variables];

            // package the exponent 1 in to the dict at the index i 
            exponent[i] = 1;
            let coefficient = FieldElement::one();
            let mut dict = HashMap::new();
            dict.insert(exponent, coefficient);

            // push the polynomial to the vector
            variables.push(MPolynomial { dict });
        }
        variables
    }

    /// Lifts a univariate polynomial into a multivariate polynomial context.
    pub fn lift(poly: &Polynomial, variable_index: usize) -> Self {
        if poly.is_zero() { return MPolynomial::zero(); }

        let variables = MPolynomial::variables(variable_index + 1);
        let x = variables[variable_index].clone();  // Assuming variables returns a vector of MPolynomial
        let mut acc = MPolynomial::zero();

        for (i, coeff) in poly.coeffs.iter().rev().enumerate() {
            let coeff_val: i128 = coeff.value.clone().to_i128().unwrap();
            let term = MPolynomial::constant(coeff_val) * x.clone().pow(i as usize);
            acc = acc + term;
        }
        acc
    }
}

// multivariate poly addition
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

// multivariate poly subtraction
impl Sub for MPolynomial {
    type Output = MPolynomial;

    fn sub(self, other: MPolynomial) -> Self::Output { 
        self + (-other) // add negated rhs
    }
}

// multivariate poly negation
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

// multivariate poly multiplication
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



#[cfg(test)]
mod tests {
    use std::vec;

    use super::*;

    #[test]
    fn test_evaluate(){
        // polynomials
        let variables = MPolynomial::variables(4);
        let mpoly1 = MPolynomial::constant(1) * variables[0].clone() + MPolynomial::constant(2) * variables[1].clone() + MPolynomial::constant(5) * variables[2].clone().pow(3);
        let mpoly2 = MPolynomial::constant(1) * variables[0].clone() * variables[3].clone() + MPolynomial::constant(5) * (variables[3].clone().pow(3)) + MPolynomial::constant(5);  
        let mpoly3 = mpoly1.clone() * mpoly2.clone();

        // point to evaluate
        let point = vec![0, 5, 5, 2].iter().map(|&x| FieldElement::new(BigInt::from(x))).collect::<Vec<FieldElement>>();

        // evaluateS
        let eval1 = mpoly1.eval(&point);
        let eval2 = mpoly2.eval(&point);
        let eval3 = mpoly3.eval(&point);

        // check correct values
        assert_eq!(eval1.value, BigInt::from(635));
        assert_eq!(eval2.value, BigInt::from(45));
        assert_eq!(eval3.value, BigInt::from(28575));

        // check eval commutativity
        assert_eq!(eval1.clone() * eval2.clone(), eval3);
        assert_eq!(eval1 + eval2, (mpoly1 + mpoly2).eval(&point));
    }

    #[test]
    fn test_lift(){

        // interpolate a univariate polynomial
        let upoly = Polynomial::lagrange(
            vec![FieldElement::new(BigInt::from(0)), FieldElement::new(BigInt::from(1)), FieldElement::new(BigInt::from(2))],
            vec![FieldElement::new(BigInt::from(2)), FieldElement::new(BigInt::from(5)), FieldElement::new(BigInt::from(5))],
        );

        // lift the univariate polynomial to a multivariate polynomial
        let mpoly = MPolynomial::lift(&upoly, 3);

        // ensure correct coefficients in mpoly
        let key_1: Vec<i128> = vec![0];
        let key_2: Vec<i128> = vec ![0, 0, 0, 1];
        let key_3: Vec<i128> = vec![0, 0, 0, 2];

        let coeff_1: FieldElement = mpoly.dict.get(&key_1).unwrap().clone();
        let coeff_2: FieldElement = mpoly.dict.get(&key_2).unwrap().clone();
        let coeff_3: FieldElement = mpoly.dict.get(&key_3).unwrap().clone();
        
        assert_eq!(coeff_1.value, BigInt::from(2));
        assert_eq!(coeff_2.value, BigInt::from(135248948571115190067962368383525060613 as i128));
        assert_eq!(coeff_3.value, BigInt::from(135248948571115190067962368383525060607 as i128));

        // 270497897142230380135924736767050121204
        let upoly_eval = upoly.eval(FieldElement::new(BigInt::from(5)));
        let mpoly_eval = mpoly.eval(&vec![FieldElement::zero(), FieldElement::zero(), FieldElement::zero(), FieldElement::new(BigInt::from(5))]);

        assert_eq!(upoly_eval, mpoly_eval);
    }
}



