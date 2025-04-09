pub const F: f64 = 1. / 298.257223563;
pub const A: f64 = 6378.137;//meters
#[derive(Clone, Copy)]
pub struct Eci {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}
#[derive(Clone, Copy)]
pub struct SatAngle {
    pub elevation: f64,
    pub azimuth: f64,
    pub range: f64,
}
#[derive(Clone, Copy)]
pub struct SubPoint {
    pub lat: f64,
    pub long: f64,
    pub alt: f64,
}

