use std::f64::consts::PI;

use serde::{Serialize,Deserialize};
use crate::{types::{Eci, SatAngle, SubPoint, A, F}, helpers::modulus,GroundStation};
use chrono::{Datelike, NaiveDateTime, Timelike};
use sgp4::Constants;


#[derive(Serialize, Deserialize, Clone)]
pub struct Satellite {
    constants: Constants,
    name: String,
    epoch: NaiveDateTime,//Epoch of the structs TLE
}
impl Satellite {
    pub fn new_from_tle(tle: &str) -> Satellite {
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
        let constants = sgp4::Constants::from_elements(&elem).unwrap();//need to get rid of this?
        Satellite {
            constants,
            name: name.to_string(),
            epoch,
        }
    }
    pub fn get_speed(&self, offset: i64) -> f64 {
        let set = self.get_point_eci(offset);
        (set.vx.powf(2.) + set.vy.powf(2.) + set.vz.powf(2.)).sqrt()
    }
    pub fn get_name(&self) -> String {
        self.name.to_string()
    }
    pub fn get_point_eci(&self,offset:i64)->Eci{
        let consts = &self.constants;
        let prop = consts.propagate(offset as f64/60.).unwrap();
        Eci {
            x: prop.position[0],
            y: prop.position[1],
            z: prop.position[2],
            vx: prop.velocity[0],
            vy: prop.velocity[1],
            vz: prop.velocity[2],
        }
    }
    pub fn get_look_angle(&self, station: &GroundStation, offset: i64) -> SatAngle {
        fn xyz_from_lla(loc: &GroundStation, sidereal_angle: f64) -> [f64; 3] {
            //Alt in meters, lat/long in decimal degrees
            let true_alt = 6378.137 + loc.alt / 1000.;
            let c = 1. / (1. + F * (F - 2.) * (loc.lat.to_radians().sin()).powf(2.)).sqrt();
            let s = (1. - F).powf(2.) * c;
            let z = true_alt * loc.lat.to_radians().sin() * s;
            let r = true_alt * loc.lat.to_radians().cos() * c;
            let theta_true = sidereal_angle.to_degrees() + loc.long;
            let x = r * theta_true.to_radians().cos();
            let y = r * theta_true.to_radians().sin();
            [x, y, z]
        }
        let time_stamp= self.epoch.timestamp()+offset;
        let sidereal = Self::sidereal_angle(time_stamp);
        let satellite = self.get_point_eci(offset);
        let observer = xyz_from_lla(station, sidereal);
        let rx = satellite.x - observer[0];
        let ry = satellite.y - observer[1];
        let rz = satellite.z - observer[2];
        let rad_lat = station.lat.to_radians();
        let side_angle = modulus(
            (sidereal.to_degrees() + station.long).to_radians(),
            2. * PI,
        );
        let rs = rad_lat.sin() * side_angle.cos() * rx + rad_lat.sin() * side_angle.sin() * ry
            - rad_lat.cos() * rz;
        let re = -side_angle.sin() * rx + side_angle.cos() * ry;
        let rz = rad_lat.cos() * side_angle.cos() * rx
            + rad_lat.cos() * side_angle.sin() * ry
            + rad_lat.sin() * rz;
        let range = (rs * rs + re * re + rz * rz).sqrt();
        let elevation = ((rz / range).asin()).to_degrees();
        let mut azimuth = (-re).atan2(rs);
        if rs > 0. {
            azimuth += PI
        } else if rs < 0. {
            azimuth += 2. * PI
        }
        azimuth = modulus(azimuth.to_degrees(), 360.);
        SatAngle{elevation, azimuth, range}
    }
    ///compute the sub_point at a given time since epoch in seconds. 
    pub fn get_sub_point(&self,offset:i64) -> SubPoint {
        fn get_lat_and_alt(satellite: &Eci) -> [f64; 2] {
            let mut guess = (satellite.z).atan2((satellite.x.powf(2.) + satellite.y.powf(2.)).sqrt());
            let e2 = 2. * F - F * F;
            let mut c = 1. / (1. - e2 * guess.sin().powf(2.)).sqrt();
            let mut last_guess: f64 = 0.;
            //let r = A*c*(guess.cos());
            let r = (satellite.x.powf(2.) + satellite.y.powf(2.)).sqrt();
            //println!("{} R",r);
            while (guess.to_degrees() - last_guess.to_degrees()).abs() > 0.00001 {
                last_guess = guess;
                c = 1. / ((1. - e2 * last_guess.sin().powf(2.)).sqrt());
                guess = (satellite.z + A * c * e2 * last_guess.sin()).atan2(r);
                //println!("Angle = {}",guess.to_degrees())
            }
            let alt = (r / guess.cos()) - A * c;
            [guess.to_degrees(), alt]
        }
        fn get_long(satellite: &Eci, sidereal_angle:f64) -> f64 {
            -1.*(modulus((((satellite.y/satellite.x).atan()) - sidereal_angle).to_degrees(),360.)-180.)
            //println!("Angle = {}",constrained_angle);
        }
        let time_stamp= self.epoch.timestamp()+offset;
        let sidereal_angle = Self::sidereal_angle(time_stamp);
        let eci_point = self.get_point_eci(offset);
        let longitude = get_long(&eci_point, sidereal_angle);
        let lat_alt = get_lat_and_alt(&eci_point);
        let latitude = lat_alt[0];
        let altitude = lat_alt[1];
        SubPoint {
            lat: latitude,
            long: longitude,
            alt: altitude,
        }
    }
    fn sidereal_angle(time_stamp: i64)->f64{
        //Meeus' approach from celestrak https://celestrak.org/columns/v02n02/
        let base_time = NaiveDateTime::from_timestamp_opt(time_stamp,0).unwrap();
        let years = base_time.year() as f64 - 1.;
        let a = (years / 100.).trunc();
        let b = 2. - a + (a / 4.).trunc();
        let julian_year =
            (365.25 * years).trunc() + ((30.6001 * 14.) as f64).trunc() + 1720994.5 + b;
        let day = base_time.ordinal();
        let julian_day = julian_year + day as f64;
        let j2000_day = julian_day - 2451545.0;
        let offset_time = base_time.num_seconds_from_midnight() as f64;
        //Adjust for UT1 off of nutation etc table
        //based on
        // https://celestrak.org/columns/v02n02/
        //From https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev2.pdf
        let offset = offset_time / 86400.;
        let t = j2000_day / 36525.0;
        let theta_0 = 24110.54841 + 8640184.812866 * t + 0.093104 * t * t
            - t * t * t * 6.2 * 10_f64.powf(-6.);
        let side_time = modulus(theta_0 + 1.00273790934 * offset * 86400., 86400.); //Why the 1.00?
                                                                                    //println!{"{}",side_time} //Make constant
        2. * PI * side_time / 86400.
    }
    pub fn seconds_since(&self,other:&NaiveDateTime)->i64{
        other.timestamp()-self.epoch.timestamp()
    }
}


#[cfg(test)]
mod tests{
    use chrono::{NaiveDate, NaiveTime};

    use crate::helpers::assert_almost_eq;

    use super::*;
    #[test]
    fn test_satellite_gen() {
        let sat = Satellite::new_from_tle(
            "DELFI-PQ                
            1 51074U 22002CU  23120.77859283  .00033391  00000+0  10673-2 0  9997
            2 51074  97.4622 192.5713 0010271  72.5102 287.7261 15.32323264 71737",
        );
        assert_eq!(sat.get_name(), "DELFI-PQ");
    }
    #[test]
    fn test_eci(){
        let sat = Satellite::new_from_tle(
            "RUST_SGP4(TESTING)               
            1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753
            2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667",
        );
        let eci_point = sat.get_point_eci(0);
        assert_almost_eq(eci_point.x,7022.46647);
        assert_almost_eq(eci_point.y,-1400.06656);
        assert_almost_eq(eci_point.z,0.05106);
    }
    #[test]
    fn test_sub_point(){
        let sat = Satellite::new_from_tle(
            "ISS (ZARYA)             
            1 25544U 98067A   25072.43808874  .00018974  00000+0  33994-3 0  9997
            2 25544  51.6354  61.2721 0006420  16.6184 343.5014 15.49959635500318",
        );
        let date = NaiveDate::from_ymd_opt(2025, 03, 13).unwrap();
        let time = NaiveTime::from_hms_opt(21, 54, 42).unwrap();
        let reference_sub_point = sat.get_sub_point(sat.seconds_since(&NaiveDateTime::new(date, time)));
        assert_almost_eq(reference_sub_point.long,-63.7546);
        assert_almost_eq(reference_sub_point.lat,35.8798);
        assert_almost_eq(reference_sub_point.alt,421.1125)
    }
    #[test]
    fn test_look_angle(){
        let sat = Satellite::new_from_tle(
            "ISS (ZARYA)             
            1 25544U 98067A   25072.43808874  .00018974  00000+0  33994-3 0  9997
            2 25544  51.6354  61.2721 0006420  16.6184 343.5014 15.49959635500318",
        );
        let date = NaiveDate::from_ymd_opt(2025, 03, 13).unwrap();
        let time = NaiveTime::from_hms_opt(23, 19, 42).unwrap();
        let refence_look_angle = sat.get_look_angle(&GroundStation::new([51.9861,4.3876,74.4], "Delft"), sat.seconds_since(&NaiveDateTime::new(date, time)));
        assert_almost_eq(refence_look_angle.azimuth, 168.7);
        assert_almost_eq(refence_look_angle.elevation, 65.2);
        assert_almost_eq(refence_look_angle.range, 463.006316);//off by 1km
    }
}