use std::fmt::Display;

use chrono::DateTime;
use chrono::Utc;
pub use satellite::Satellite;
pub use types::Eci;
pub use types::SatAngle;
pub use types::SubPoint;
pub use ground_station::GroundStation;
mod types;
mod helpers;
mod satellite;
mod ground_station;

pub struct Pass{
    aos:i64,
    tme:i64,
    max_elevation:f64,
    los:i64
}

impl Display for Pass{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f,"aos: {}, los: {}, tme: {}, max_el: {:.2}",self.get_aos_datetime(),self.get_los_datetime(),self.get_tme_datetime(),self.max_elevation)
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
    pub fn get_aos_datetime(&self)->DateTime<Utc>{
        DateTime::from_timestamp(self.aos, 0).unwrap()
    }
    pub fn get_tme_datetime(&self)->DateTime<Utc>{
        DateTime::from_timestamp(self.tme, 0).unwrap()
    }
    pub fn get_los_datetime(&self)->DateTime<Utc>{
        DateTime::from_timestamp(self.los, 0).unwrap()
    }
}

pub fn find_passes_datetime(sat:&Satellite,gs:&GroundStation,start_datetime:&DateTime<Utc>,end_datetime:&DateTime<Utc>)->Vec<Pass>{
    find_passes_offset(sat, gs, sat.seconds_since_epoch(start_datetime),sat.seconds_since_epoch(end_datetime))
}

pub fn find_passes_offset(sat:&Satellite,gs:&GroundStation,start_offset:i64,end_offset:i64)->Vec<Pass>{
    //https://celestrak.org/columns/v03n04/ - Reference
    println!("starting");
    let orbital_period = sat.get_orbital_period();
    let mut timestamp = start_offset-orbital_period.floor() as i64;
    let mut mins_maxes:Vec<(i64,bool)> = vec![];//heapless vec!!
    let mut next_search:bool = false;
    while timestamp < end_offset{
        let next_data = find_min_max(sat, gs, timestamp as f64, timestamp as f64 + orbital_period,next_search);
        next_search = !next_data.1;
        timestamp = next_data.0;
        mins_maxes.push(next_data);
    }
    let mut passes:Vec<Pass> = vec![];
    println!("minmaxes:{}",mins_maxes.len());
    for i in (0..mins_maxes.len()).step_by(2){
        if i+1 == mins_maxes.len(){
            break
        }
        if !mins_maxes[i].1 && mins_maxes[i+1].1{
            // println!("{},{}",sat.get_look_angle(gs,mins_maxes[i].0.floor() as i64).elevation,sat.get_look_angle(gs,mins_maxes[i+1].0.floor() as i64).elevation);
            if sat.get_look_angle(gs,mins_maxes[i].0).elevation < 0. && sat.get_look_angle(gs,mins_maxes[i+1].0).elevation > 0.{//Pass with an above horizon maxima. i is the minima before and i+1 is the maxima after
                let aos = find_zero(sat, gs, mins_maxes[i].0, mins_maxes[i+1].0);
                if aos<start_offset{
                    continue
                }
                let los = find_zero(sat, gs, mins_maxes[i+2].0, mins_maxes[i+1].0);
                let max_el = sat.get_look_angle(gs, mins_maxes[i+1].0 as i64);
                passes.push(Pass{ aos: compute_timestamp(sat, aos), tme: compute_timestamp(sat,mins_maxes[i+1].0), max_elevation: max_el.elevation, los: compute_timestamp(sat,los) })
            }
        }
        else{
            println!("something went wrong")
        }
    }
    passes
}

fn find_zero(sat: &Satellite, gs: &GroundStation, lower_bound:i64,upper_bound:i64) ->i64{
    let mut a = lower_bound;
    let mut b =  upper_bound;
    let mut c;
    while (b-a).abs() > 1{
        c = (a+b)/2;
        let computed_c = sat.get_look_angle(gs, c);
        if computed_c.elevation<0.{
            a = c
        } else {
            b = c
        }
    }
    let elev_a = sat.get_look_angle(gs, a).elevation;
    let elev_b = sat.get_look_angle(gs, b).elevation;
    let fraction = elev_a / (elev_a - elev_b);
    let refined_c = a as f64 + fraction * (b - a) as f64;
    return refined_c as i64
}

fn compute_timestamp(sat:&Satellite,time:i64)->i64{
    time+sat.get_epoch().timestamp()
}


const GOLDEN_RATIO:f64 = 1.6180339887498948482045868343656381177203091798057;//Should be enough precision 
fn find_min_max(sat:&Satellite,gs:&GroundStation,start_offset:f64,end_offset:f64,next_search:bool)->(i64,bool){//returns time and if max/min (max=true)
    let mut ls = start_offset;
    let mut rs = end_offset;
    let mut guess = end_offset-(end_offset-start_offset)/GOLDEN_RATIO;
    let maxima= next_search;
    while (rs-ls).abs() > 1.0{
        // ls_value = sat.get_look_angle(gs, ls.floor() as i64);
        // rs_value = sat.get_look_angle(gs, rs.floor() as i64);
        guess = rs-(rs-ls)/GOLDEN_RATIO;
        let second_guess = ls+(rs-ls)/GOLDEN_RATIO;//look between A and C
        if maxima{//Then guess is a minima
            if sat.get_look_angle(gs, second_guess.floor() as i64).elevation <sat.get_look_angle(gs, guess.floor() as i64).elevation{
                rs = second_guess;
            } else {
                ls = guess;
            }
        } else {
            if sat.get_look_angle(gs, second_guess.floor() as i64).elevation > sat.get_look_angle(gs, guess.floor() as i64).elevation{
                rs = second_guess;
            } else {
                ls = guess;
            }
        }
    }
    (guess.round() as i64,maxima)
    
}

#[cfg(test)]
mod tests{
    use crate::{find_passes_datetime, helpers::quick_gen_datetime, GroundStation, Satellite};
    #[test]
    fn test_find_passes(){
        let sat = Satellite::new_from_tle(
            "ISS (ZARYA)             
1 25544U 98067A   25078.36999458  .00023040  00000+0  41584-3 0  9998
2 25544  51.6365  31.8868 0003892  28.0409 332.0788 15.49628144501233",
        );
        println!("{}",quick_gen_datetime(2025, 3, 20, 11, 10, 0).timestamp());
        let gs = GroundStation::new([51.9861,4.3876,74.4], "Test");
        let start_date_time = quick_gen_datetime(2025, 03, 21, 00, 22, 0);
        let end_date_time = quick_gen_datetime(2025, 03, 23, 23, 22, 0);
        let passes = find_passes_datetime(&sat, &gs, &start_date_time, &end_date_time);
        for i in &passes{
            println!("{}",i);
        }
        assert_eq!(passes.len(),18);
        println!("{}",sat.get_look_angle(&gs,sat.offset_timestamp(quick_gen_datetime(2025, 3, 19, 21, 37, 48).timestamp())).elevation);
        assert_eq!(passes[0].get_aos_datetime(),quick_gen_datetime(2025, 3, 21, 15, 13, 34));
        assert_eq!(passes[0].get_tme_datetime(),quick_gen_datetime(2025, 3, 21, 15, 17, 53));
        assert_eq!(passes[0].get_aos_datetime(),quick_gen_datetime(2025, 3, 21, 15, 22, 13));
    }
}