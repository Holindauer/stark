use crate::modules::field::{*};
use std::vec;


// agreed upon fibonacci square sequence program a_n+2 = (a_n+1)**2 + (a_n)**2
pub fn fib_square (x: i128) -> Vec<FieldElement> {

    // create vec adn add starting point
    let mut a: Vec<FieldElement> = Vec::new();
    a.push(new_field_element(1));
    a.push(new_field_element(x));

    while a.len() < 1023 {
        let last: FieldElement = *a.get(a.len()-1).unwrap();
        let second_last: FieldElement = *a.get(a.len()-2).unwrap();
        a.push(second_last * second_last + last * last);
    }
    a
}



// tests
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sanity_check() {
        let a = fib_square(3141592);
        let last = *a.get(a.len()-1).unwrap();
        assert_eq!(2338775057, last.value);
    }
}