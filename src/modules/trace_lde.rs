use crate::modules::field::{*};
use crate::modules::polynomial::{*};


// agreed upon fibonacci square sequence program a_n+2 = (a_n+1)**2 + (a_n)**2
// this function accepts the second term and generates up until 1023 terms
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
    for i in 0..1024 { domain.push(gen.pow(i)); }
    domain
}

// generate coset for evaluation 8x larger than interpolation domain
// this is reed-solomon error correction, which introduces redundancy
// in order to allow errors in the interpolated polynomial when evaluated
pub fn extended_domain() -> Vec<FieldElement> {

    // get generator 
    let w = FieldElement::generator();
    println!("h: {}", w.value);

    println!("power: {}", 3 * 2_i128.pow(30) / 8192);

    let h = w.pow((3 * 2_i128.pow(30)) / 8192); // 8192 / 8 = 1024
    println!("h: {}", h.value);
    
    let mut H: Vec<FieldElement> = Vec::new();
    for i in 0..8192 { H.push(h.pow(i)); }

    let mut domain: Vec<FieldElement> = Vec::new();
    for i in 0..8192 { domain.push(w * *H.get(i).unwrap()); }


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

    #[test] 
    fn test_extend_domain() {
        let extended_domain: Vec<FieldElement> = extended_domain();
        assert_eq!(extended_domain.get(0).unwrap().value, 5);
        assert_eq!(extended_domain.get(extended_domain.len()-1).unwrap().value, 1375380442);
    }





}