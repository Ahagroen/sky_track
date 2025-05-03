use std::fmt::Display;

use chrono::DateTime;
use chrono::Utc;
pub use ground_station::GroundStation;
pub use satellite::Satellite;
pub use types::Eci;
pub use types::SatAngle;
pub use types::SubPoint;
mod ground_station;
mod helpers;
mod satellite;
mod types;

pub struct Pass {
    aos: i64,
    tme: i64,
    max_elevation: f64,
    los: i64,
}

impl Display for Pass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "aos: {}, los: {}, tme: {}, max_el: {:.2}",
            self.get_aos_datetime(),
            self.get_los_datetime(),
            self.get_tme_datetime(),
            self.max_elevation
        )
    }
}

impl Pass {
    pub fn get_tme(&self) -> i64 {
        self.tme
    }

    pub fn get_max_elevation(&self) -> f64 {
        self.max_elevation
    }

    pub fn get_los(&self) -> i64 {
        self.los
    }

    pub fn get_aos(&self) -> i64 {
        self.aos
    }
    pub fn get_aos_datetime(&self) -> DateTime<Utc> {
        DateTime::from_timestamp(self.aos, 0).unwrap()
    }
    pub fn get_tme_datetime(&self) -> DateTime<Utc> {
        DateTime::from_timestamp(self.tme, 0).unwrap()
    }
    pub fn get_los_datetime(&self) -> DateTime<Utc> {
        DateTime::from_timestamp(self.los, 0).unwrap()
    }
}

pub fn find_passes_datetime(
    sat: &Satellite,
    gs: &GroundStation,
    start_datetime: &DateTime<Utc>,
    end_datetime: &DateTime<Utc>,
) -> Vec<Pass> {
    find_passes_offset(
        sat,
        gs,
        sat.seconds_since_epoch(start_datetime),
        sat.seconds_since_epoch(end_datetime),
    )
}

pub fn find_passes_offset(
    sat: &Satellite,
    gs: &GroundStation,
    start_offset: i64,
    end_offset: i64,
) -> Vec<Pass> {
    //https://celestrak.org/columns/v03n04/ - Reference
    println!("starting");
    let orbital_period = sat.get_orbital_period();
    let mut timestamp = start_offset - orbital_period.floor() as i64;
    let mut mins_maxes: Vec<(i64, bool)> = vec![]; //heapless vec!!
    let mut next_search: bool = false;
    while timestamp < end_offset {
        let next_data = find_min_max(
            sat,
            gs,
            timestamp as f64,
            timestamp as f64 + orbital_period,
            next_search,
        );
        next_search = !next_data.1;
        timestamp = next_data.0;
        mins_maxes.push(next_data);
    }
    let mut passes: Vec<Pass> = vec![];
    for i in (0..mins_maxes.len()).step_by(2) {
        if i + 1 == mins_maxes.len() {
            break;
        }
        if !mins_maxes[i].1 && mins_maxes[i + 1].1 {
            // println!("{},{}",sat.get_look_angle(gs,mins_maxes[i].0.floor() as i64).elevation,sat.get_look_angle(gs,mins_maxes[i+1].0.floor() as i64).elevation);
            if sat.get_look_angle(gs, mins_maxes[i].0).elevation < 0.
                && sat.get_look_angle(gs, mins_maxes[i + 1].0).elevation > 0.
            {
                //Pass with an above horizon maxima. i is the minima before and i+1 is the maxima after
                let aos = find_zero(sat, gs, mins_maxes[i].0, mins_maxes[i + 1].0);
                if aos < start_offset {
                    continue;
                }
                let los = find_zero(sat, gs, mins_maxes[i + 2].0, mins_maxes[i + 1].0);
                let max_el = sat.get_look_angle(gs, mins_maxes[i + 1].0 as i64);
                passes.push(Pass {
                    aos: compute_timestamp(sat, aos),
                    tme: compute_timestamp(sat, mins_maxes[i + 1].0),
                    max_elevation: max_el.elevation,
                    los: compute_timestamp(sat, los),
                })
            }
        } else {
            println!("something went wrong")
        }
    }
    passes
}

fn find_zero(sat: &Satellite, gs: &GroundStation, lower_bound: i64, upper_bound: i64) -> i64 {
    let mut a = lower_bound;
    let mut b = upper_bound;
    let mut c;
    while (b - a).abs() > 1 {
        c = (a + b) / 2;
        let computed_c = sat.get_look_angle_refraction(gs, c);
        if computed_c.elevation < 0. {
            a = c
        } else {
            b = c
        }
    }
    b as i64
}

fn compute_timestamp(sat: &Satellite, time: i64) -> i64 {
    time + sat.get_epoch().timestamp()
}

const GOLDEN_RATIO: f64 = 1.6180339887498948482045868343656; //Should be enough precision 
fn find_min_max(
    sat: &Satellite,
    gs: &GroundStation,
    start_offset: f64,
    end_offset: f64,
    next_search: bool,
) -> (i64, bool) {
    //returns time and if max/min (max=true)
    let mut ls = start_offset;
    let mut rs = end_offset;
    let mut guess = end_offset - (end_offset - start_offset) / GOLDEN_RATIO;
    let maxima = next_search;
    while (rs - ls).abs() > 1.0 {
        // ls_value = sat.get_look_angle(gs, ls.floor() as i64);
        // rs_value = sat.get_look_angle(gs, rs.floor() as i64);
        guess = rs - (rs - ls) / GOLDEN_RATIO;
        let second_guess = ls + (rs - ls) / GOLDEN_RATIO; //look between A and C
        if maxima {
            //Then guess is a minima
            if sat
                .get_look_angle(gs, second_guess.floor() as i64)
                .elevation
                < sat.get_look_angle(gs, guess.floor() as i64).elevation
            {
                rs = second_guess;
            } else {
                ls = guess;
            }
        } else {
            if sat
                .get_look_angle(gs, second_guess.floor() as i64)
                .elevation
                > sat.get_look_angle(gs, guess.floor() as i64).elevation
            {
                rs = second_guess;
            } else {
                ls = guess;
            }
        }
    }
    (ls.round() as i64, maxima)
}

#[cfg(test)]
mod tests {
    use crate::{GroundStation, Satellite, find_passes_datetime, helpers::quick_gen_datetime};
    #[test]
    fn test_find_one_pass() {
        let sat = Satellite::new_from_tle(
            "ISS (ZARYA)             
1 25544U 98067A   25122.54440123  .00015063  00000+0  27814-3 0  9994
2 25544  51.6345 173.1350 0002187  74.2134 285.9096 15.49297959508085",
        );
        let gs = GroundStation::new([51.9861, 4.3876, 0.], "Test");
        let start_date_time = quick_gen_datetime(2025, 05, 02, 21, 55, 14);
        let end_date_time = quick_gen_datetime(2025, 05, 02, 23, 00, 00);
        let passes = find_passes_datetime(&sat, &gs, &start_date_time, &end_date_time);
        for i in &passes {
            println!(
                "AOS: {}, azimuth: {:.2}, elevation: {:.2}",
                i.get_aos_datetime(),
                sat.get_look_angle_refraction(&gs, sat.offset_timestamp(i.get_aos()))
                    .azimuth,
                sat.get_look_angle_refraction(&gs, sat.offset_timestamp(i.get_aos()))
                    .elevation
            );
            println!(
                "TME: {}, azimuth: {:.2}, elevation: {:.2}",
                i.get_tme_datetime(),
                sat.get_look_angle_refraction(&gs, sat.offset_timestamp(i.get_tme()))
                    .azimuth,
                sat.get_look_angle_refraction(&gs, sat.offset_timestamp(i.get_tme()))
                    .elevation
            );
            println!(
                "LOS: {}, azimuth: {:.2}, elevation: {:.2}",
                i.get_los_datetime(),
                sat.get_look_angle_refraction(&gs, sat.offset_timestamp(i.get_los()))
                    .azimuth,
                sat.get_look_angle_refraction(&gs, sat.offset_timestamp(i.get_los()))
                    .elevation
            );
        }
        assert!(false)
    }
    #[test]
    fn test_find_passes() {
        let sat = Satellite::new_from_tle(
            "ISS (ZARYA)             
1 25544U 98067A   25122.54440123  .00015063  00000+0  27814-3 0  9994
2 25544  51.6345 173.1350 0002187  74.2134 285.9096 15.49297959508085",
        );
        let gs = GroundStation::new([51.9861, 4.3876, 0.], "Test");
        let start_date_time = quick_gen_datetime(2025, 05, 02, 21, 55, 14);
        let end_date_time = quick_gen_datetime(2025, 05, 05, 21, 55, 14);
        let passes = find_passes_datetime(&sat, &gs, &start_date_time, &end_date_time);
        for i in &passes {
            println!(
                "AOS: {}, azimuth: {:.2}, elevation: {:.2}",
                i.get_aos_datetime(),
                sat.get_look_angle_refraction(&gs, sat.offset_timestamp(i.get_aos()))
                    .azimuth,
                sat.get_look_angle_refraction(&gs, sat.offset_timestamp(i.get_aos()))
                    .elevation
            );
            println!(
                "TME: {}, azimuth: {:.2}, elevation: {:.2}",
                i.get_tme_datetime(),
                sat.get_look_angle_refraction(&gs, sat.offset_timestamp(i.get_tme()))
                    .azimuth,
                sat.get_look_angle_refraction(&gs, sat.offset_timestamp(i.get_tme()))
                    .elevation
            );
            println!(
                "LOS: {}, azimuth: {:.2}, elevation: {:.2}",
                i.get_los_datetime(),
                sat.get_look_angle_refraction(&gs, sat.offset_timestamp(i.get_los()))
                    .azimuth,
                sat.get_look_angle_refraction(&gs, sat.offset_timestamp(i.get_los()))
                    .elevation
            );
            println!()
        }
        println!(
            "{},{}",
            sat.get_look_angle_refraction(
                &gs,
                sat.offset_timestamp(quick_gen_datetime(2025, 5, 2, 22, 40, 43).timestamp())
            )
            .azimuth,
            sat.get_look_angle_refraction(
                &gs,
                sat.offset_timestamp(quick_gen_datetime(2025, 5, 2, 22, 40, 43).timestamp())
            )
            .elevation
        );
        assert_eq!(passes.len(), 19);
        println!(
            "{}",
            sat.get_look_angle(&gs, sat.offset_timestamp(passes[0].get_aos()))
                .elevation
        );
        assert_eq!(
            passes[0].get_aos_datetime(),
            quick_gen_datetime(2025, 5, 2, 22, 31, 54)
        );
        assert_eq!(
            passes[0].get_tme_datetime(),
            quick_gen_datetime(2025, 5, 2, 22, 36, 18)
        );
        assert_eq!(
            passes[0].get_aos_datetime(),
            quick_gen_datetime(2025, 5, 2, 22, 40, 43)
        );
    }
}
