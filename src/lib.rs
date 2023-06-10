use chrono::NaiveDateTime;
use drivers::Eci;
use serde::{Deserialize, Serialize};
use sgp4::Constants;

mod drivers;
pub mod ground_track;
mod helpers;
pub mod pass_list;
//DataSets

#[derive(Clone, Serialize, Deserialize)]
pub struct GroundStation {
    lat: f64,
    long: f64,
    alt: f64,
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
pub struct Satellite {
    constants: Constants,
    pub name: String,
    pub epoch: NaiveDateTime,
    pub time_since_launch: u64,
}
impl Satellite {
    pub fn find_tle(_sat_name: &str) -> Satellite {
        let tle = helpers::get_tle_data(); //Will take sat_name once implemented
        Self::create_from_tle(&tle)
    }
    pub fn create_from_tle(tle: &str) -> Satellite {
        let mut line_count = 0;
        let mut name: &str = "";
        let mut first_line: &str = "";
        let mut second_line: &str = "";
        for line in tle.lines() {
            if line_count == 2 {
                //println!("second line: {}",line);
                second_line = line.trim();
                line_count = -1;
            } else if line_count == 1 {
                //println!("first line: {}",line);
                first_line = line.trim();
                line_count = 2;
            } else if line_count == 0 {
                //println!("name: {}",line);
                name = line.trim();
                line_count = 1;
            }
        }
        let elem = sgp4::Elements::from_tle(
            Some(name.to_string()),
            first_line.as_bytes(),
            second_line.as_bytes(),
        )
        .unwrap();
        let epoch = elem.datetime;
        let constants = sgp4::Constants::from_elements(&elem).unwrap();
        let orbit_days = (elem.revolution_number as f64 / elem.mean_motion).floor() as u64;
        Satellite {
            constants,
            name: name.to_string(),
            epoch,
            time_since_launch: orbit_days,
        }
    }
    pub fn get_speed(&self, time: f64) -> f64 {
        let set = Eci::new_from_point(self, time);
        (set.vx.powf(2.) + set.vy.powf(2.) + set.vz.powf(2.)).sqrt()
    }
    pub fn get_name(&self) -> String {
        self.name.to_string()
    }
}

#[cfg(test)]
mod tests {
    use crate::ground_track::Track;
    use crate::pass_list::PassTime;
    use crate::pass_list::Passes;
    use crate::GroundStation;
    use crate::Satellite;
    use chrono::NaiveDate;
    use chrono::NaiveDateTime;
    use chrono::NaiveTime;
    #[test]
    fn test_find_passes_number() {
        let sat = Satellite::create_from_tle(
            "DELFI-PQ                
            1 51074U 22002CU  23120.77859283  .00033391  00000+0  10673-2 0  9997
            2 51074  97.4622 192.5713 0010271  72.5102 287.7261 15.32323264 71737",
        );
        let delft = GroundStation::new([52.0022, 4.3736, 0.], "Test");
        let stations = vec![&delft];
        let passes = Passes::find_passes(
            stations,
            &sat,
            12 * 60,
            NaiveDateTime::new(
                NaiveDate::from_ymd_opt(2023, 05, 1).unwrap(),
                NaiveTime::from_num_seconds_from_midnight_opt(54900, 0).unwrap(),
            ),
        );
        for i in &passes.pass_list {
            println!("{} Pass", i.aos)
        }
        assert_eq!(3, passes.pass_list.len())
    }
    #[test]
    fn test_satellite_gen() {
        let sat = Satellite::create_from_tle(
            "DELFI-PQ                
            1 51074U 22002CU  23120.77859283  .00033391  00000+0  10673-2 0  9997
            2 51074  97.4622 192.5713 0010271  72.5102 287.7261 15.32323264 71737",
        );
        assert_eq!(sat.name, "DELFI-PQ");
    }
    #[test]
    fn test_track_creation() {
        let sat = Satellite::create_from_tle(
            "DELFI-PQ                
            1 51074U 22002CU  23120.77859283  .00033391  00000+0  10673-2 0  9997
            2 51074  97.4622 192.5713 0010271  72.5102 287.7261 15.32323264 71737",
        );
        let _epoch = NaiveDateTime::new(
            NaiveDate::from_ymd_opt(2023, 5, 1).unwrap(),
            NaiveTime::from_num_seconds_from_midnight_opt(45960, 0).unwrap(),
        );
        let track = Track::get_track(25., &sat, 0.);
        assert_eq!(track.points.len(), 25)
    }
    #[test]
    fn test_track_offset() {
        let sat = Satellite::create_from_tle(
            "DELFI-PQ                
            1 51074U 22002CU  23120.77859283  .00033391  00000+0  10673-2 0  9997
            2 51074  97.4622 192.5713 0010271  72.5102 287.7261 15.32323264 71737",
        );
        let _epoch = NaiveDateTime::new(
            NaiveDate::from_ymd_opt(2023, 5, 1).unwrap(),
            NaiveTime::from_num_seconds_from_midnight_opt(45960, 0).unwrap(),
        );
        let track = Track::get_track(25., &sat, -200.);
        assert_eq!(track.points.len(), 25)
    }
    #[test]
    fn test_track_last_point() {
        let sat = Satellite::create_from_tle(
            "DELFI-PQ                
            1 51074U 22002CU  23120.77859283  .00033391  00000+0  10673-2 0  9997
            2 51074  97.4622 192.5713 0010271  72.5102 287.7261 15.32323264 71737",
        );
        let _epoch = NaiveDateTime::new(
            NaiveDate::from_ymd_opt(2023, 5, 1).unwrap(),
            NaiveTime::from_num_seconds_from_midnight_opt(45960, 0).unwrap(),
        );
        let track = Track::get_track(500., &sat, -250.);
        assert_eq!(track.points.last().unwrap().time, 249)
    }
    #[test]
    fn test_pass_generation() {
        let sat = Satellite::create_from_tle(
            "DELFI-PQ                
            1 51074U 22002CU  23120.77859283  .00033391  00000+0  10673-2 0  9997
            2 51074  97.4622 192.5713 0010271  72.5102 287.7261 15.32323264 71737",
        );
        let station = GroundStation::new([52.0022, 4.3736, 0.], "Test");
        let pass = PassTime::new(&sat, 1608, &station);
        println!(
            "Pass at: {}\nTill: {}\nfor: {} seconds with max elevation: {}",
            pass.aos, pass.los, pass.length, pass.max_elevation
        );
        assert_eq!(pass.max_elevation, 63.955533284304856)
    }
}
