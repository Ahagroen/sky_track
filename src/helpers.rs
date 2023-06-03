pub fn modulus(a: f64, b: f64) -> f64 {
    ((a % b) + b) % b
}

pub fn get_tle_data() -> String {
    //Modify for variable satellite
    //https://celestrak.org/NORAD/elements/gp.php?NAME=Delfi&FORMAT=TLE - produces all 3 delfi sats
    let mut res =
        reqwest::blocking::get("https://celestrak.org/NORAD/elements/gp.php?CATNR=51074").unwrap();
    let mut body = String::new();
    std::io::Read::read_to_string(&mut res, &mut body).unwrap(); //Pattern matching to handle errors
    body
}
