pub fn modulus(a: f64, b: f64) -> f64 {
    ((a % b) + b) % b
}