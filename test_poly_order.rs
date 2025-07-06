use zk_stark::modules::field::FieldElement;
use zk_stark::modules::univariate_poly::Polynomial;
use num_bigint::BigInt;

fn main() {
    // Test 1: Create polynomial x (coefficients [1, 0] if highest to lowest)
    let poly_x = Polynomial::new(vec![BigInt::from(1), BigInt::from(0)]);
    
    // Evaluate at x=2
    let x = FieldElement::new(BigInt::from(2));
    let result = poly_x.eval(x.clone());
    
    println!("Test 1: Polynomial [1, 0] evaluated at x=2 gives: {}", result.value);
    println!("If coefficient order is highest to lowest, this should be 2 (since 1*x + 0 = x)");
    println!("If coefficient order is lowest to highest, this should be 1 (since 1 + 0*x = 1)");
    println!();
    
    // Test 2: Create polynomial 3x + 5 (coefficients [3, 5] if highest to lowest)
    let poly_3x_plus_5 = Polynomial::new(vec![BigInt::from(3), BigInt::from(5)]);
    let result2 = poly_3x_plus_5.eval(x.clone());
    
    println!("Test 2: Polynomial [3, 5] evaluated at x=2 gives: {}", result2.value);
    println!("If coefficient order is highest to lowest, this should be 11 (since 3*x + 5 = 3*2 + 5 = 11)");
    println!("If coefficient order is lowest to highest, this should be 13 (since 3 + 5*x = 3 + 5*2 = 13)");
    println!();
    
    // Test 3: Create polynomial x^2 + 2x + 3 (coefficients [1, 2, 3] if highest to lowest)
    let poly_quadratic = Polynomial::new(vec![BigInt::from(1), BigInt::from(2), BigInt::from(3)]);
    let result3 = poly_quadratic.eval(x);
    
    println!("Test 3: Polynomial [1, 2, 3] evaluated at x=2 gives: {}", result3.value);
    println!("If coefficient order is highest to lowest, this should be 11 (since x^2 + 2*x + 3 = 4 + 4 + 3 = 11)");
    println!("If coefficient order is lowest to highest, this should be 11 (since 1 + 2*x + 3*x^2 = 1 + 4 + 12 = 17)");
}
