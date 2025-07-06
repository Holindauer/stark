use crate::modules::field::{*};
use std::vec;
use std::ops::{Add, Sub, Mul, Div, Neg};
use std::fmt;
use serde::{Serialize, Deserialize};
use std::fs::File;
use std::io::{self, Read, Write};
use num_bigint::RandBigInt;
use num_bigint::BigInt;
use num_traits::{Zero, One};


/**
 * polynomial.rs implements polynomial operations, evaluation, and lagrange 
 * interpolation over the specific finite field implemented in field.rs
 */

 #[derive(PartialEq, Debug, Clone, Serialize, Deserialize)]
pub struct Polynomial {
    pub coeffs: Vec<FieldElement>,
}

// polynomial constructor
impl Polynomial {

    // coefficients stored highest to lowest degree
    pub fn new(coeffs: Vec<BigInt>) -> Self {
        Polynomial { coeffs: coeffs.iter().map(|x| FieldElement::new(x.clone())).collect() }
    }
    
    // evaluate polynomial at field element x
    pub fn eval(&self, x: FieldElement) -> FieldElement {

        // evaluate polynomial at x using Horner's method
        let mut val: BigInt = BigInt::zero();
        for coef in self.coeffs.iter() {
            val = (val * x.value.clone() + coef.value.clone()) % FieldElement::modulus();
        }
        FieldElement::new(val)
    }

    // evaluate over a domain of field elements
    pub fn eval_domain(&self, domain: Vec<FieldElement>) -> Vec<FieldElement> {
        domain.iter().map(|x| self.eval(x.clone())).collect()
    }

    // get degree of polynomial
    pub fn degree(&self) -> usize {
        let mut degree = self.coeffs.len() - 1;
        let zero = FieldElement::zero();

        // remove trailing zeros
        for c in self.coeffs.iter() {
            if c != &zero { break; }
            degree -= 1;
        }
        degree
    }

    // checks for zero polynomial
    pub fn is_zero(&self) -> bool {
        let mut is_zero: bool = true;
        for i in 0..self.coeffs.len() {
            if self.coeffs.get(i).unwrap().value != BigInt::zero() { is_zero = false; }
        }
        is_zero
    }

    // Constructs a monomial x^degree with a specific coefficient
    pub fn monomial(degree: usize, coefficient: FieldElement) -> Polynomial {
        let mut coeffs = vec![FieldElement::zero(); degree + 1];
        coeffs[0] = coefficient;  // Highest degree term
        Polynomial { coeffs }
    }

    // Constructs the Lagrange basis polynomial for a given index i
    pub fn lagrange_basis(x: &Vec<FieldElement>, i: usize) -> Polynomial {
        let mut l = Polynomial::monomial(0, FieldElement::one());  // Start with 1

        for j in 0..x.len() {
            if i != j {
                // (x - x_j)
                let monomial = Polynomial::monomial(1, FieldElement::one());
                let constant = Polynomial::new(vec![-x[j].value.clone()]);  // Use negation method from FieldElement
                let num = monomial + constant;

                // Divide by (x_i - x_j), ensure subtraction handles negatives correctly
                let denom_diff = x[i].clone() - x[j].clone();
                let denom = denom_diff.inverse();  // Should handle negative by wrapping around the modulus correctly

                l = l * num;
                l = l.scale(&denom);
            }
        }

        l
    }

    // Performs Lagrange interpolation based on given x and y values
    pub fn lagrange(x: Vec<FieldElement>, y: Vec<FieldElement>) -> Polynomial {
        assert_eq!(x.len(), y.len(), "x and y must be the same length");

        let n = x.len();
        let mut p = Polynomial::new(vec![]); // Start with the zero polynomial

        for i in 0..n {
            let li = Polynomial::lagrange_basis(&x, i);
            p = p + (li.scale(&y[i]));
        }

        // addd modulus to negative elements
        for i in 0..p.coeffs.len() {
            if p.coeffs[i].value < BigInt::zero(){
                p.coeffs[i].value += FieldElement::modulus();
            }
        }

        p
    }

    // Scales all coefficients of the polynomial by a FieldElement
    pub fn scale(&self, factor: &FieldElement) -> Polynomial {
        let scaled_coeffs: Vec<FieldElement> = self.coeffs.iter().map(|coeff| coeff.clone() * factor.clone()).collect();
        Polynomial::new(scaled_coeffs.iter().map(|x| x.value.clone()).collect())
    }

    // Save the polynomial to a JSON file
    pub fn save(&self, filename: &str) -> io::Result<()> {
        let json = serde_json::to_string(&self).unwrap();
        let mut file = File::create(filename)?;
        file.write_all(json.as_bytes())?;
        Ok(())
    }

    // Load a polynomial from a JSON file
    pub fn load(filename: &str) -> io::Result<Polynomial> {
        let mut file = File::open(filename)?;
        let mut json = String::new();
        file.read_to_string(&mut json)?;
        let poly: Polynomial = serde_json::from_str(&json).unwrap();
        Ok(poly)
    }


    // vanishing polynomial (x - c_1)(x - c_2) . . . (c - c_n)
    pub fn zeroifier(domain: Vec<BigInt>) -> Polynomial {
        
        let x = Polynomial::new(vec![BigInt::from(1), BigInt::from(0)]);
        let mut acc = Polynomial::new(vec![BigInt::from(1)]);
        for d in domain {
            acc = acc * (x.clone() - Polynomial::new(vec![d.clone()])); 
        }

        acc 
    }

    // zeroifier over domain
    pub fn zeroifier_domain(domain: Vec<FieldElement>) -> Polynomial {
        let x = Polynomial{coeffs: vec![FieldElement::one(), FieldElement::zero()]};
        let mut acc = Polynomial{coeffs: vec![FieldElement::one()]};
        for d in domain {

            // (x - d)
            acc = acc * (x.clone() - Polynomial{coeffs: vec![d]});
        }

        acc 
    }

    // test for colinearity
    pub fn test_colinearity(points: Vec<(BigInt, BigInt)>) -> bool {

        // collect intepolation domains and values
        let mut domain: Vec<FieldElement> = vec![];
        let mut values: Vec<FieldElement> = vec![];
        for (x, y) in points { 
            domain.push(FieldElement::new(x));
            values.push(FieldElement::new(y));
        } 

        // lagrange interpolation
        let poly = Polynomial::lagrange(domain, values);

        // ensure colinearity
        if poly.degree() == 1 { true } else { false }
    }

    // pow
    pub fn pow(&self, exponent: u128) -> Self {
        // special cases
        if self.is_zero() { return Polynomial::new(vec![BigInt::from(0)]); }
        if exponent == 0 { return Polynomial::new(vec![BigInt::from(1)]); }

        let mut acc = Polynomial::new(vec![BigInt::from(1)]); // Start with the identity polynomial
        let base = self.clone();

        // Process each bit of the exponent, from least significant to most significant
        for i in (0..exponent.leading_zeros() as u128).rev() {
            acc = acc.clone() * acc.clone();  // Square the result

            if (exponent & (1 << i)) != 0 {
                acc = acc.clone() * base.clone();  // Multiply by base if the i-th bit is 1
            }
        }

        acc
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
                res.push(max_poly.coeffs.get(i).unwrap().clone());
            }
        }
        
        // then add together the remaining shared degrees elementwise
        for i in 0..min_poly.coeffs.len() {
            let coeff_1: FieldElement = max_poly.coeffs.get(diff_len+i).unwrap().clone();
            let coeff_2: FieldElement = min_poly.coeffs.get(i).unwrap().clone();
            res.push(coeff_1 + coeff_2);
        }

        Polynomial::new(res.iter().map(|x| x.value.clone()).collect())
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
            for i in 0..diff_len {res.push(FieldElement::zero() - rhs.coeffs.get(i).unwrap().clone()); }

            // set offsets for degree union subtraction
            rhs_offset = diff_len; lhs_offset = 0;
        }
        // if minuend is longer, push lhs terms until poly degree union
        else if self.coeffs.len() > rhs.coeffs.len() {

            // set max and min polynomials, and difference in length
            max_poly = self.clone(); min_poly = rhs.clone();
            diff_len = max_poly.coeffs.len() - min_poly.coeffs.len();

            // push lhs terms until poly degree union
            for i in 0..diff_len {res.push(max_poly.coeffs.get(i).unwrap().clone()); }

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
            let rhs_coeff: FieldElement = rhs.coeffs.get(rhs_offset+i).unwrap().clone();
            let lhs_coeff: FieldElement = self.coeffs.get(lhs_offset+i).unwrap().clone();
            res.push(lhs_coeff - rhs_coeff);
        }

        Polynomial::new(res.iter().map(|x| x.value.clone()).collect())
    }
}

// polynomial multiplication
impl Mul for Polynomial {
    type Output = Polynomial;

    fn mul(self, rhs: Polynomial) -> Polynomial {
       
        // init result coeffs vec with (combined poly coef len - 2) + 1 zero() elements
        let mut res: Vec<FieldElement> = vec![FieldElement::zero(); self.coeffs.len() + rhs.coeffs.len() - 1];    

        // multiply polynomials by convolution
        for (i, c1) in self.coeffs.iter().enumerate() {
            for (j, c2) in rhs.coeffs.iter().enumerate() {
                res[i + j] = res[i + j].clone() + c1.clone() * c2.clone();
            }
        }

        Polynomial::new(res.iter().map(|x| x.value.clone()).collect())
    }
}

// negation
impl Neg for Polynomial {
    type Output = Polynomial;

    // negate by subtracting from zero
    fn neg(self) -> Polynomial {   
        let mut res: Vec<FieldElement> = vec![];
        for coef in self.coeffs.iter() {
            res.push(FieldElement::zero() - coef.clone());
        }
        Polynomial::new(res.iter().map(|x| x.value.clone()).collect())
    }
}

// division
impl Div for Polynomial {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        // Ensure denominator is not the zero polynomial
        if rhs.is_zero() { panic!("Attempted to divide by zero polynomial."); }

        // Initialize the numerator and denominator  
        let mut num = self.coeffs;
        let mut denom = rhs.coeffs;

        // Reverse for easier handling of indices (from highest to lowest degree)
        num.reverse(); denom.reverse();

        // Initialize the quotient with the correct degree
        let num_degree = num.len() - 1;
        let denom_degree = denom.len() - 1;
        let mut quotient = vec![FieldElement::zero(); num_degree - denom_degree + 1];

        // perform polynomial long division
        for i in 0..=num_degree - denom_degree {
            if num[i] != FieldElement::zero() {
                // Calculate the coefficient for the current term of the quotient
                let coeff = num[i].clone() / denom[0].clone();
                quotient[i] = coeff.clone();

                // Subtract the appropriate multiple of the denominator from the numerator
                for (j, denom_coeff) in denom.iter().enumerate() {
                    if i + j <= num_degree {
                        num[i + j] = num[i+j].clone() - coeff.clone() * denom_coeff.clone();
                    }
                }
            }
        }

        // Remove leading zeros from the quotient and reverse back to normal order
        while let Some(last) = quotient.last() {
            if last.clone() == FieldElement::zero() {
                quotient.pop();
            } else {
                break;
            }
        }

        quotient.reverse();

        // Return the resulting polynomial
        Polynomial { coeffs: quotient }
    }
}

// prints polynomial formatted polynomial
impl fmt::Display for Polynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let len = self.coeffs.len();

        // display coeefs from highest term to lowest
        for (i, coeff) in self.coeffs.iter().enumerate() {
            if i > 0 && coeff.value != BigInt::zero() { 
                write!(f, " + ")?; 
                write!(f, "{}x^{}", coeff, len - i - 1)?;
            }
        }
        Ok(())
    }
}

// tests
#[cfg(test)]
mod tests {
    use super::*;
    use rand::Rng;

    #[test]
    fn test_distributivity() {
        let zero = BigInt::from(0);
        let one = BigInt::from(1);
        let two = BigInt::from(2);
        let five = BigInt::from(5);

        let a = Polynomial::new(vec![zero.clone(), one.clone(), two.clone()]);
        let b = Polynomial::new(vec![two.clone(), two.clone(), one.clone()]);
        let c = Polynomial::new(vec![zero.clone(), five.clone(), two.clone(), five.clone(), five.clone(), one.clone()]);

        let lhs = a.clone() * (b.clone() + c.clone());
        let rhs = a.clone() * b.clone() + a.clone() * c.clone();

        assert_eq!(lhs, rhs);
    }

    #[test] 
    fn test_vanishing_polynomial_fuzz(){

        // randomly set vanishing points
        let mut rng = rand::thread_rng();
        let mut vanish_at: Vec<BigInt> = vec![];
        for i in 0..10 {
            let upper_bound = FieldElement::modulus();
            let random_value = rng.gen_bigint_range(&BigInt::zero(), &upper_bound);
            vanish_at.push(random_value);
        }

        // create vanishing polynomial
        let vanish_poly = Polynomial::zeroifier(vanish_at.clone());

        // check all zero
        for x in vanish_at{
            let eval = vanish_poly.eval(FieldElement::new(x.clone()));
            assert_eq!(eval.value, BigInt::from(0));
        } 
    }

    #[test] 
    fn test_test_colinearity(){

        let points = vec![
            (BigInt::from(0), BigInt::from(0)),
            (BigInt::from(1), BigInt::from(1)),
            (BigInt::from(2), BigInt::from(2)),
        ];


        let test = Polynomial::test_colinearity(points);
        assert_eq!(test, true);
    }

    #[test]
    fn test_add_same_degree() {

        // add polynomials
        let poly_1: Polynomial = Polynomial::new(vec![BigInt::from(10), BigInt::from(3), BigInt::from(1)]);
        let poly_2: Polynomial = Polynomial::new(vec![BigInt::from(90), BigInt::from(3), BigInt::from(1)]);
        let poly = poly_1.clone() + poly_2.clone();

        // ensure correct
        assert_eq!(poly.coeffs.get(0).unwrap().value, BigInt::from(100));
        assert_eq!(poly.coeffs.get(1).unwrap().value, BigInt::from(6));
        assert_eq!(poly.coeffs.get(2).unwrap().value, BigInt::from(2));
        }

    #[test]
    fn test_sub_same_degree() {

        // add polynomials
        let poly_1: Polynomial = Polynomial::new(vec![BigInt::from(80), BigInt::from(6), BigInt::from(1)]);
        let poly_2: Polynomial = Polynomial::new(vec![BigInt::from(40), BigInt::from(3), BigInt::from(6)]);
        let poly = poly_1.clone() - poly_2.clone();

        println!("coeffs len: {}", poly.coeffs.len());

        // ensure correct
        assert_eq!(poly.coeffs.get(0).unwrap().value, BigInt::from(40));
        assert_eq!(poly.coeffs.get(1).unwrap().value, BigInt::from(3));
        assert_eq!(poly.coeffs.get(2).unwrap().value, FieldElement::modulus()-BigInt::from(5));
    }

    #[test]
    fn test_sub_larger_lhs() {

        // add polynomials
        let lhs: Polynomial = Polynomial::new(vec![BigInt::from(80), BigInt::from(6), BigInt::from(1), BigInt::from(80), BigInt::from(6), BigInt::from(1)]);
        let rhs: Polynomial = Polynomial::new(vec![BigInt::from(40), BigInt::from(3), BigInt::from(6)]);
        let poly = lhs.clone() - rhs.clone();

        // ensure correct
        assert_eq!(poly.coeffs.len(), 6);
        assert_eq!(poly.coeffs.get(0).unwrap().value, BigInt::from(80));
        assert_eq!(poly.coeffs.get(1).unwrap().value, BigInt::from(6));
        assert_eq!(poly.coeffs.get(2).unwrap().value, BigInt::from(1));
        assert_eq!(poly.coeffs.get(3).unwrap().value, BigInt::from(40));
        assert_eq!(poly.coeffs.get(4).unwrap().value, BigInt::from(3));
        assert_eq!(poly.coeffs.get(5).unwrap().value, FieldElement::modulus()-BigInt::from(5));
    }

    #[test]
    fn test_sub_larger_rhs() {

        // add polynomials
        let lhs: Polynomial = Polynomial::new(vec![BigInt::from(40), BigInt::from(3), BigInt::from(6)]);
        let rhs: Polynomial = Polynomial::new(vec![BigInt::from(80), BigInt::from(6),BigInt::from(1), BigInt::from(80), BigInt::from(6), BigInt::from(1)]);
        let poly = lhs.clone() - rhs.clone();

        // ensure correct
        assert_eq!(poly.coeffs.len(), 6);
        assert_eq!(poly.coeffs.get(0).unwrap().value, FieldElement::modulus()-BigInt::from(80));
        assert_eq!(poly.coeffs.get(1).unwrap().value, FieldElement::modulus()-BigInt::from(6));
        assert_eq!(poly.coeffs.get(2).unwrap().value, FieldElement::modulus()-BigInt::from(1));
        assert_eq!(poly.coeffs.get(3).unwrap().value, FieldElement::modulus()-BigInt::from(40));
        assert_eq!(poly.coeffs.get(4).unwrap().value, FieldElement::modulus()-BigInt::from(3));
        assert_eq!(poly.coeffs.get(5).unwrap().value, BigInt::from(5));
    }

    #[test]
    fn test_add_diff_degree() {

        // polynomials
        let poly_1: Polynomial = Polynomial::new(vec![BigInt::from(10), BigInt::from(3), BigInt::from(1)]);
        let poly_2: Polynomial = Polynomial::new(vec![BigInt::from(1), BigInt::from(2), BigInt::from(3), BigInt::from(90), BigInt::from(3), BigInt::from(1)]);

        // add them
        let poly = poly_1.clone() + poly_2.clone();

        // ensure correct
        assert_eq!(poly.coeffs.get(0).unwrap().value, BigInt::from(1));
        assert_eq!(poly.coeffs.get(1).unwrap().value, BigInt::from(2));
        assert_eq!(poly.coeffs.get(2).unwrap().value, BigInt::from(3));
        assert_eq!(poly.coeffs.get(3).unwrap().value, BigInt::from(100));
        assert_eq!(poly.coeffs.get(4).unwrap().value, BigInt::from(6));
        assert_eq!(poly.coeffs.get(5).unwrap().value, BigInt::from(2));
    }

    #[test]
    fn test_mul() {
        // multiply, evaluate, and ensure correct
        let poly_1: Polynomial = Polynomial::new(vec![BigInt::from(10), BigInt::from(1), BigInt::from(0), BigInt::from(1)]);
        let poly_2: Polynomial = Polynomial::new(vec![BigInt::from(3), BigInt::from(1), BigInt::from(17)]);
        let poly = poly_1.clone() * poly_2.clone();
        let eval_poly = poly.eval(FieldElement::new(BigInt::from(2)));

        // ensure correct
        assert_eq!(eval_poly.value,  BigInt::from(2635));
    }

    #[test]
    fn test_negate() {
        // negate, evaluate, and ensure correct
        let poly: Polynomial = Polynomial::new(vec![BigInt::from(10), BigInt::from(3), BigInt::from(1)]);
        let neg_poly = -poly.clone();
        let eval_neg_poly = neg_poly.eval(FieldElement::new(BigInt::from(9)));
        assert_eq!(eval_neg_poly.value,  FieldElement::modulus()-838);
    }

    #[test]
    fn test_div_1() {
        //setup
        let a = Polynomial::new(vec![BigInt::from(1), BigInt::from(2)]);
        let b = Polynomial::new(vec![BigInt::from(1), BigInt::from(1)]);
        let c = a.clone() * b.clone();

        // division
        let d = c.clone() / a.clone();
        let e = c.clone() / b.clone();

        assert_eq!(b.coeffs, d.coeffs);
        assert_eq!(a.coeffs, e.coeffs);
    }

    #[test]
    fn test_div_2() {
        //setup
        let a = Polynomial::new(vec![BigInt::from(1), BigInt::from(20), BigInt::from(35), BigInt::from(-1460), BigInt::from(-9396), BigInt::from(-12960), BigInt::from(0)]);
        let b = Polynomial::new(vec![BigInt::from(1), BigInt::from(-7), BigInt::from(-18)]);
        let c = a.clone() / b.clone();

        println!("{}", c); // checksout in printout
    }
   
    #[test]
    fn test_eval(){
        let poly = Polynomial::new(vec![BigInt::from(10), BigInt::from(3), BigInt::from(1)]);
        let out = poly.eval(FieldElement::new(BigInt::from(2)));
        assert_eq!(out.value,  BigInt::from(47));
    }

    #[test]
    fn test_fuzz_lagrange() {

        // get random x and y values
        let mut rng = rand::thread_rng();
        let mut x: Vec<FieldElement> = vec![];
        let mut y: Vec<FieldElement> = vec![];
        for i in 0..3 {
            x.push(FieldElement::new(BigInt::from(i)));

            let upper_bound = FieldElement::modulus();
            let random_value = rng.gen_bigint_range(&BigInt::zero(), &upper_bound);
            y.push(FieldElement::new(random_value));
        }

        // lagrange interpolation
        let poly = Polynomial::lagrange(x.clone(), y.clone());

        // eval 
        let out_1 = poly.eval(FieldElement::new(x[0].value.clone()));        
        let out_2 = poly.eval(FieldElement::new(x[1].value.clone()));
        let out_3 = poly.eval(FieldElement::new(x[2].value.clone()));

        // ensure correct
        assert_eq!(out_1.value, y[0].value);
        assert_eq!(out_2.value, y[1].value);
        assert_eq!(out_3.value, y[2].value);   
    }

    #[test]
    fn test_polynomial_coefficient_order() {
        // Test 1: Create polynomial x (coefficients [1, 0] if highest to lowest)
        let poly_x = Polynomial::new(vec![BigInt::from(1), BigInt::from(0)]);
        
        // Evaluate at x=2
        let x = FieldElement::new(BigInt::from(2));
        let result = poly_x.eval(x.clone());
        
        println!("Test 1: Polynomial [1, 0] evaluated at x=2 gives: {}", result.value);
        println!("If coefficient order is highest to lowest, this should be 2 (since 1*x + 0 = x)");
        println!("If coefficient order is lowest to highest, this should be 1 (since 1 + 0*x = 1)");
        assert_eq!(result.value, BigInt::from(2), "Coefficients are stored highest to lowest degree");
        
        // Test 2: Create polynomial 3x + 5 (coefficients [3, 5] if highest to lowest)
        let poly_3x_plus_5 = Polynomial::new(vec![BigInt::from(3), BigInt::from(5)]);
        let result2 = poly_3x_plus_5.eval(x.clone());
        
        println!("\nTest 2: Polynomial [3, 5] evaluated at x=2 gives: {}", result2.value);
        println!("If coefficient order is highest to lowest, this should be 11 (since 3*x + 5 = 3*2 + 5 = 11)");
        println!("If coefficient order is lowest to highest, this should be 13 (since 3 + 5*x = 3 + 5*2 = 13)");
        assert_eq!(result2.value, BigInt::from(11), "Coefficients are stored highest to lowest degree");
        
        // Test 3: Create polynomial x^2 + 2x + 3 (coefficients [1, 2, 3] if highest to lowest)
        let poly_quadratic = Polynomial::new(vec![BigInt::from(1), BigInt::from(2), BigInt::from(3)]);
        let result3 = poly_quadratic.eval(x);
        
        println!("\nTest 3: Polynomial [1, 2, 3] evaluated at x=2 gives: {}", result3.value);
        println!("If coefficient order is highest to lowest, this should be 11 (since x^2 + 2*x + 3 = 4 + 4 + 3 = 11)");
        println!("If coefficient order is lowest to highest, this should be 17 (since 1 + 2*x + 3*x^2 = 1 + 4 + 12 = 17)");
        assert_eq!(result3.value, BigInt::from(11), "Coefficients are stored highest to lowest degree");
    }
}


