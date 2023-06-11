use crate::drivers::SatAngle;
use crate::GroundStation;
use crate::Satellite;
use chrono::NaiveDateTime;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize, Clone)]
pub struct PassTime {
    pub satellite: Satellite,
    pub ground_station: GroundStation,
    pub aos: NaiveDateTime,
    pub los: NaiveDateTime,
    pub length: f64,
    pub max_elevation: f64,
}
impl PassTime {
    pub fn search_for_passes(
        satellite: &Satellite,
        ground: &GroundStation,
        range: i64,
        start_time: NaiveDateTime,
    ) -> Vec<PassTime> {
        let start_delta = start_time.signed_duration_since(satellite.epoch);
        let base_offset = start_delta.num_seconds() / 60;
        let mut last_pass_time: u16 = 6000;
        let mut usable_pass = Vec::new();
        let mut pass = false;
        for i in base_offset..(range + base_offset) {
            if !(120..=480).contains(&last_pass_time) {
                let angle = SatAngle::get_angle_offset(satellite, ground, i as f64);
                if angle.elevation > 5.0 && !pass {
                    //can change depending on constraints - make variable?
                    pass = true;
                    //println!("{}pass",i);
                    let pass = PassTime::new(satellite, i, ground);
                    usable_pass.push(pass);
                } else if pass && angle.elevation < 5.0 {
                    //println!("{}End",i);
                    pass = !pass;
                    last_pass_time = 0;
                }
            } else {
                last_pass_time += 1
            }
        }
        usable_pass
    }

    pub fn new(sat: &Satellite, guess: i64, station: &GroundStation) -> PassTime {
        //println!("stop");
        let los: NaiveDateTime = Self::get_stop_time(sat, station, guess);
        //println!("start");
        let aos: NaiveDateTime = Self::get_start_time(sat, station, guess);
        let length = los.signed_duration_since(aos).num_milliseconds() as f64 / 1000.;
        //println!("{} Duration in seconds",length);
        let max_elevation = Self::get_highest_elevation(sat, station, length, aos);
        PassTime {
            satellite: sat.clone(),
            aos,
            los,
            length,
            max_elevation,
            ground_station: station.clone(),
        }
    }
    fn get_stop_time(sat: &Satellite, station: &GroundStation, guess: i64) -> NaiveDateTime {
        //works
        //Search Algorithm - Recursive Bisection
        let end_min = Self::angle_bisection(guess as f64, (guess + 10) as f64, sat, station);
        let end_time = sat.epoch.timestamp() + (end_min * 60.).ceil() as i64;
        NaiveDateTime::from_timestamp_opt(end_time, 0).unwrap()
    }
    fn get_start_time(sat: &Satellite, station: &GroundStation, guess: i64) -> NaiveDateTime {
        //Search Algorithm - Recursive Bisection
        let start_min = Self::angle_bisection((guess - 5) as f64, guess as f64, sat, station);
        let start_time = sat.epoch.timestamp() + (start_min * 60.).floor() as i64;
        NaiveDateTime::from_timestamp_opt(start_time, 0).unwrap()
    }
    fn get_highest_elevation(
        sat: &Satellite,
        station: &GroundStation,
        length: f64,
        aos: NaiveDateTime,
    ) -> f64 {
        //Taken at the midpoint time of the pass
        //println!("Mid time - Minutes: {}",mid_time);
        let mid_time =
            aos.signed_duration_since(sat.epoch).num_seconds() as f64 / 60. + length / 120.;
        let max_elevation = SatAngle::get_angle_offset(sat, station, mid_time); //Always at the midpoint??
        max_elevation.elevation
    }
    fn angle_bisection(left: f64, right: f64, sat: &Satellite, station: &GroundStation) -> f64 {
        //1. Compute midpoint elevation
        //2. set that to the new bound based on +- 5 degrees
        //3. re-compute
        let mut rs_carry = right;
        let mut ls_carry = left;
        let left_side_angle = SatAngle::get_angle_offset(sat, station, ls_carry);
        let low_to_high: bool = left_side_angle.elevation < 0.;
        loop {
            let midpoint = (rs_carry + ls_carry) / 2.;
            let angle = SatAngle::get_angle_offset(sat, station, midpoint);
            if (angle.elevation - 5.0).abs() < 0.1 || (ls_carry - rs_carry).abs() < 0.001 {
                return midpoint;
            } else if angle.elevation - 5.0 < 0. {
                //FIX: Refactor
                if low_to_high {
                    ls_carry = midpoint
                } else {
                    rs_carry = midpoint;
                }
            } else {
                if low_to_high {
                    rs_carry = midpoint
                } else {
                    ls_carry = midpoint;
                }
            }
        }
    }
}
