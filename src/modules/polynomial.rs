use crate::modules::field::{*};
use std::vec;
use std::ops::{Add, Sub, Mul, Div, Neg};
use std::fmt;

/**
 * polynomial.rs implements polynomial operations, evaluation, and lagrange 
 * interpolation over the specific finite field implemented in field.rs
 */

 #[derive(Debug, Clone)]
pub struct Polynomial {
    pub coeffs: Vec<FieldElement>,
}

// polynomial constructor
impl Polynomial {

    // coefficients stored highest to lowest degree
    pub fn new(coeffs: Vec<FieldElement>) -> Self {
        Polynomial { coeffs }
    }
    
    // evaluate polynomial at field element x
    pub fn eval(&self, x: FieldElement) -> FieldElement {

        // evaluate polynomial at x using Horner's method
        let mut val: i128 = 0;
        for coef in self.coeffs.iter() {
            val = (val * x.value + coef.value) % x.modulus;
        }
        new_field_element(val)
    }
}

// polynomial addition
impl Add for Polynomial {
    type Output = Polynomial;

    fn add(self, other: Polynomial) -> Polynomial {

        // len of max poly is len of result
        let max_poly: Polynomial; let min_poly: Polynomial;
        
        // len correspons to highest degree - 1
        if self.coeffs.len() > other.coeffs.len() { max_poly = self; min_poly = other; } 
        else { max_poly = other; min_poly = self; }

        // difference of length
        let diff_len: usize = max_poly.coeffs.len() - min_poly.coeffs.len();

        // push higher degree terms to result from max poly
        let mut res: Vec<FieldElement> = vec![];

        // if max degree is different, add higher degrees first
        if diff_len > 0{
            for i in 0..diff_len {
                res.push(*max_poly.coeffs.get(i).unwrap());
            }
        }
        
        // then add together the remaining shared degrees elementwise
        for i in 0..min_poly.coeffs.len() {
            let coeff_1: FieldElement = *max_poly.coeffs.get(diff_len+i).unwrap();
            let coeff_2: FieldElement = *min_poly.coeffs.get(i).unwrap();
            res.push(coeff_1 + coeff_2);
        }

        Polynomial::new(res)
    }
}

// polynomial subtraction
impl Sub for Polynomial {
    type Output = Polynomial;

    fn sub(self, rhs: Polynomial) -> Polynomial {

        // init result coeffs vec
        let mut res: Vec<FieldElement> = vec![];
        let diff_len: usize;

        let rhs_offset: usize;
        let lhs_offset: usize;

        let max_poly: Polynomial;
        let min_poly: Polynomial;

        // if subtrahend is longer, push negated rhs terms until poly degree union 
        if rhs.coeffs.len() > self.coeffs.len() {

            // set max and min polynomials, and difference in length
            max_poly = rhs.clone(); min_poly = self.clone();
            diff_len = max_poly.coeffs.len() - min_poly.coeffs.len();
            for i in 0..diff_len {res.push(zero() - *rhs.coeffs.get(i).unwrap()); }

            // set offsets for degree union subtraction
            rhs_offset = diff_len; lhs_offset = 0;
        }
        // if minuend is longer, push lhs terms until poly degree union
        else if self.coeffs.len() > rhs.coeffs.len() {

            // set max and min polynomials, and difference in length
            max_poly = self.clone(); min_poly = rhs.clone();
            diff_len = max_poly.coeffs.len() - min_poly.coeffs.len();

            // push lhs terms until poly degree union
            for i in 0..diff_len {res.push(*max_poly.coeffs.get(i).unwrap()); }

            // set offsets for degree union subtraction
            rhs_offset = 0; lhs_offset = diff_len;
        }
        else { 
            // max and min polynomials are equal in this case
            min_poly = rhs.clone();
            rhs_offset = 0; lhs_offset = 0;
        }

        // sub remaining terms elementwise
        let union_len = min_poly.coeffs.len();
        for i in 0..union_len {
            let rhs_coeff: FieldElement = *rhs.coeffs.get(rhs_offset+i).unwrap();
            let lhs_coeff: FieldElement = *self.coeffs.get(lhs_offset+i).unwrap();
            res.push(lhs_coeff - rhs_coeff);
        }

        println!("{:?}", res);

        Polynomial::new(res)
    }
}

// negation
impl Neg for Polynomial {
    type Output = Polynomial;

    // negate by subtracting from zero
    fn neg(self) -> Polynomial {   
        let mut res: Vec<FieldElement> = vec![];
        for coef in self.coeffs.iter() {
            res.push(zero() - *coef);
        }
        Polynomial::new(res)
    }
}



// prints polynomial formatted polynomial
impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let len = self.coeffs.len();

        // display coeefs from highest term to lowest
        for (i, coeff) in self.coeffs.iter().enumerate() {
            if i > 0 { write!(f, " + ")?; }
            write!(f, "{}x^{}", coeff, len - i - 1)?;
        }
        Ok(())
    }
}


// tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn add_same_degree() {

        // add polynomials
        let poly_1: Polynomial = Polynomial::new(vec![new_field_element(10), new_field_element(3), new_field_element(1)]);
        let poly_2: Polynomial = Polynomial::new(vec![new_field_element(90), new_field_element(3), new_field_element(1)]);
        let poly = poly_1.clone() + poly_2.clone();

        // ensure correct
        assert_eq!(poly.coeffs.get(0).unwrap().value, 100);
        assert_eq!(poly.coeffs.get(1).unwrap().value, 6);
        assert_eq!(poly.coeffs.get(2).unwrap().value, 2);
    }


    #[test]
    fn sub_same_degree() {

        // add polynomials
        let poly_1: Polynomial = Polynomial::new(vec![new_field_element(80), new_field_element(6), new_field_element(1)]);
        let poly_2: Polynomial = Polynomial::new(vec![new_field_element(40), new_field_element(3), new_field_element(6)]);
        let poly = poly_1.clone() - poly_2.clone();

        println!("coeffs len: {}", poly.coeffs.len());

        // ensure correct
        assert_eq!(poly.coeffs.get(0).unwrap().value, 40);
        assert_eq!(poly.coeffs.get(1).unwrap().value, 3);
        assert_eq!(poly.coeffs.get(2).unwrap().value, -5);
    }

    #[test]
    fn sub_larger_lhs() {

        // add polynomials
        let lhs: Polynomial = Polynomial::new(vec![
            new_field_element(80), new_field_element(6), new_field_element(1), 
            new_field_element(80), new_field_element(6), new_field_element(1)]);
        let rhs: Polynomial = Polynomial::new(vec![new_field_element(40), new_field_element(3), new_field_element(6)]);
        let poly = lhs.clone() - rhs.clone();

        // ensure correct
        assert_eq!(poly.coeffs.len(), 6);
        assert_eq!(poly.coeffs.get(0).unwrap().value, 80);
        assert_eq!(poly.coeffs.get(1).unwrap().value, 6);
        assert_eq!(poly.coeffs.get(2).unwrap().value, 1);
        assert_eq!(poly.coeffs.get(3).unwrap().value, 40);
        assert_eq!(poly.coeffs.get(4).unwrap().value, 3);
        assert_eq!(poly.coeffs.get(5).unwrap().value, -5);
    }

    #[test]
    fn sub_larger_rhs() {

        // add polynomials
        let lhs: Polynomial = Polynomial::new(vec![new_field_element(40), new_field_element(3), new_field_element(6)]);
        let rhs: Polynomial = Polynomial::new(vec![
            new_field_element(80), new_field_element(6), new_field_element(1), 
            new_field_element(80), new_field_element(6), new_field_element(1)]);
        let poly = lhs.clone() - rhs.clone();

        // ensure correct
        assert_eq!(poly.coeffs.len(), 6);
        assert_eq!(poly.coeffs.get(0).unwrap().value, -80);
        assert_eq!(poly.coeffs.get(1).unwrap().value, -6);
        assert_eq!(poly.coeffs.get(2).unwrap().value, -1);
        assert_eq!(poly.coeffs.get(3).unwrap().value, -40);
        assert_eq!(poly.coeffs.get(4).unwrap().value, -3);
        assert_eq!(poly.coeffs.get(5).unwrap().value, 5);
    }

    #[test]
    fn add_diff_degree() {

        // field element coefficients
        let coeffs_1: Vec<FieldElement> = vec![new_field_element(10), new_field_element(3), new_field_element(1)];
        let coeffs_2: Vec<FieldElement> = vec![
            new_field_element(1), new_field_element(2), new_field_element(3),
            new_field_element(90), new_field_element(3), new_field_element(1)
            ];

        // polynomials
        let poly_1: Polynomial = Polynomial::new(coeffs_1);
        let poly_2: Polynomial = Polynomial::new(coeffs_2);

        // add them
        let poly = poly_1.clone() + poly_2.clone();

        // ensure correct
        assert_eq!(poly.coeffs.get(0).unwrap().value, 1);
        assert_eq!(poly.coeffs.get(1).unwrap().value, 2);
        assert_eq!(poly.coeffs.get(2).unwrap().value, 3);
        assert_eq!(poly.coeffs.get(3).unwrap().value, 100);
        assert_eq!(poly.coeffs.get(4).unwrap().value, 6);
        assert_eq!(poly.coeffs.get(5).unwrap().value, 2);
    }

    #[test]
    fn negate_test() {
        // negate, evaluate, and ensure correct
        let poly: Polynomial = Polynomial::new(vec![new_field_element(10), new_field_element(3), new_field_element(1)]);
        let neg_poly = -poly.clone();
        let eval_neg_poly = neg_poly.eval(new_field_element(9));
        assert_eq!(eval_neg_poly.value,  -838);
    }


    #[test]
    fn eval_test(){
        let poly = Polynomial::new(vec![new_field_element(10), new_field_element(3), new_field_element(1)]);
        let out = poly.eval(new_field_element(2));
        assert_eq!(out.value,  47);
    }
}


