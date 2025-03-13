use chrono::NaiveDateTime;
pub use satellite::Satellite;
pub use types::Eci;
pub use types::SatAngle;
pub use types::SubPoint;
pub use ground_station::GroundStation;
mod types;
mod helpers;
mod satellite;
mod ground_station;

struct Pass{
    aos:i64,
    tme:i64,
    max_elevation:f64,
    los:i64
}

pub fn find_passes_datetime(sat:&Satellite,gs:&GroundStation,datetime:&NaiveDateTime)->Vec<Pass>{
    find_passes_offset(sat, gs, sat.seconds_since_epoch(datetime))
}

pub fn find_passes_offset(sat:&Satellite,gs:&GroundStation,offset:i64)->Vec<Pass>{

}

#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn test_find_passes(){
        let sat = Satellite::new_from_tle(
            "ISS (ZARYA)             
            1 25544U 98067A   25072.43808874  .00018974  00000+0  33994-3 0  9997
            2 25544  51.6354  61.2721 0006420  16.6184 343.5014 15.49959635500318",
        );

    }
}