use std::f64;

pub fn modulus(a: f64, b: f64) -> f64 {
    ((a % b) + b) % b
}

#[cfg(test)]
pub fn assert_almost_eq<T>(a:T,b:T)
where T: std::ops::Sub+PartialOrd<f64>+std::fmt::Debug+PartialOrd<T>+Copy,
<T as std::ops::Sub>::Output: PartialOrd<f64>
{
    if a-b<1e-4{
        assert!(true)
    } else{
        assert_eq!(a,b)
    }
}