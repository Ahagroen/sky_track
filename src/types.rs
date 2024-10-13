use serde::{Deserialize, Serialize};

pub const F: f64 = 1. / 298.257223563;
pub const A: f64 = 6378.135;
#[derive(Clone, Copy)]
pub struct Eci {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}
pub struct SatAngle {
    pub elevation: f64,
    pub azimuth: f64,
    pub range: f64,
}
pub struct SubPoint {
    pub lat: f64,
    pub long: f64,
    pub alt: f64,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct GroundStation {
    pub lat: f64,
    pub long: f64,
    pub alt: f64,
    pub name: String,
}
impl GroundStation {
    pub fn new(point: [f64; 3], name: &str) -> GroundStation {
        GroundStation {
            lat: point[0],
            long: point[1],
            alt: point[2],
            name: name.to_string(),
        }
    }
}