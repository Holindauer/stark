use crate::modules::field::{*};
use std::vec;
use std::ops::{Add, Sub, Mul, Div};
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


/**

def eval(self, point):
    """
    Evaluates the polynomial at the given point using Horner evaluation.
    """
    point = FieldElement.typecast(point).val
    val = 0
    for coef in self.poly[::-1]:
        val = (val * point + coef.val) % FieldElement.k_modulus
    return FieldElement(val)

 */




// tests
#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn add_same_degree() {


        // field element coefficients
        let coeffs_1: Vec<FieldElement> = vec![new_field_element(10), new_field_element(3), new_field_element(1)];
        let coeffs_2: Vec<FieldElement> = vec![new_field_element(90), new_field_element(3), new_field_element(1)];

        // polynomials
        let poly_1: Polynomial = Polynomial::new(coeffs_1);
        let poly_2: Polynomial = Polynomial::new(coeffs_2);

        // add them
        let poly = poly_1.clone() + poly_2.clone();

        // ensure correct
        assert_eq!(poly.coeffs.get(0).unwrap().value, 100);
        assert_eq!(poly.coeffs.get(1).unwrap().value, 6);
        assert_eq!(poly.coeffs.get(2).unwrap().value, 2);
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
}


