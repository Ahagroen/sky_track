use serde::{Deserialize, Serialize};

mod types;
mod helpers;
pub mod satellite;
//DataSets
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


#[cfg(test)]
mod tests {
    use crate::satellite::Satellite;
    #[test]
    fn test_satellite_gen() {
        let sat = Satellite::new_from_tle(
            "DELFI-PQ                
            1 51074U 22002CU  23120.77859283  .00033391  00000+0  10673-2 0  9997
            2 51074  97.4622 192.5713 0010271  72.5102 287.7261 15.32323264 71737",
        );
        assert_eq!(sat.get_name(), "DELFI-PQ");
    }
}
