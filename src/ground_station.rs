use serde::{Deserialize, Serialize};
#[derive(Clone, Serialize, Deserialize)]
pub struct GroundStation {
    pub lat: f64,
    pub long: f64,
    pub alt: f64,
    pub name: String,
}
impl GroundStation {
    ///Point is Lat/long/alt
    pub fn new(point: [f64; 3], name: &str) -> GroundStation {
        GroundStation {
            lat: point[0],
            long: point[1],
            alt: point[2],
            name: name.to_string(),
        }
    }
}