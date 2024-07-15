use crate::modules::field::{*};
use std::vec;


// agreed upon fibonacci square sequence program a_n+2 = (a_n+1)**2 + (a_n)**2
pub fn fib_square_trace (x: i128) -> Vec<FieldElement> {

    // create vec adn add starting point
    let mut a: Vec<FieldElement> = Vec::new();
    a.push(FieldElement::new(1));
    a.push(FieldElement::new(x));

    while a.len() < 1023 {
        let last: FieldElement = *a.get(a.len()-1).unwrap();
        let second_last: FieldElement = *a.get(a.len()-2).unwrap();
        a.push(second_last * second_last + last * last);
    }
    a
}

// generate domain of squared generators for interpolated polynomial
pub fn poly_domain() -> Vec<FieldElement> {

    let mut domain: Vec<FieldElement> = Vec::new();
    let gen: FieldElement = FieldElement::generator().pow( 3 * 2_i128.pow(20));

    for i in 0..1024 {
        domain.push(gen.pow(i));
    }
    domain
}

// tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trace_sanity_check() {
        let a = fib_square_trace(3141592);
        let first = *a.get(0).unwrap();
        let last = *a.get(a.len()-1).unwrap();
        assert_eq!(1, first.value);
        assert_eq!(2338775057, last.value);
    }

    #[test]
    fn test_domain() {
        
        let domain: Vec<FieldElement> = poly_domain();

        assert_eq!(domain.get(0).unwrap().value, 1);
        assert_eq!(domain.get(2).unwrap().value, 764652596);
        assert_eq!(domain.get(domain.len()-1).unwrap().value, 532203874)
    }
}