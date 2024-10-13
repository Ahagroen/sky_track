use std::f64::consts::PI;

use serde::{Serialize,Deserialize};
use crate::{types::{Eci, SatAngle, SubPoint, A, F}, helpers::modulus,GroundStation};
use chrono::{DateTime, Datelike, NaiveDateTime, Timelike};
use sgp4::Constants;


#[derive(Serialize, Deserialize, Clone)]
pub struct Satellite {
    constants: Constants,
    name: String,
    epoch: NaiveDateTime,
    time_since_launch: u64,//is this needed
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
        let constants = sgp4::Constants::from_elements(&elem).unwrap();
        let orbit_days = (elem.revolution_number as f64 / elem.mean_motion).floor() as u64;
        Satellite {
            constants,
            name: name.to_string(),
            epoch,
            time_since_launch: orbit_days,
        }
    }
    pub fn get_speed(&self, time: u64) -> f64 {
        let set = self.get_point_eci(time);
        (set.vx.powf(2.) + set.vy.powf(2.) + set.vz.powf(2.)).sqrt()
    }
    pub fn get_name(&self) -> String {
        self.name.to_string()
    }
    pub fn get_point_eci(&self,timestamp:u64)->Eci{
        let consts = &self.constants;
        let offset = (timestamp - self.epoch.and_utc().timestamp() as u64) as f64/60.;
        let prop = consts.propagate(offset).unwrap();
        Eci {
            x: prop.position[0],
            y: prop.position[1],
            z: prop.position[2],
            vx: prop.velocity[0],
            vy: prop.velocity[1],
            vz: prop.velocity[2],
        }
    }
    pub fn get_look_angle(&self, station: &GroundStation, time_stamp: u64) -> SatAngle {
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
        let sidereal = Self::sidereal_angle(time_stamp);
        let satellite = self.get_point_eci(time_stamp);
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
    pub fn get_sub_point(&self,time_stamp:u64) -> SubPoint {
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
            let angle = ((satellite.y.atan2(satellite.x)) - sidereal_angle).to_degrees() + 180.;
            modulus(angle, 360.) - 180.
            //println!("Angle = {}",constrained_angle);
        }
        let sidereal_angle = Self::sidereal_angle(time_stamp);
        let eci_point = self.get_point_eci(time_stamp);
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
    fn sidereal_angle(time_stamp: u64)->f64{
        //Meeus' approach from celestrak https://celestrak.org/columns/v02n02/
        let base_time = DateTime::from_timestamp(time_stamp.try_into().expect("Couldn't convert U64 into I64"),0).unwrap();
        let years = base_time.year() as f64 - 1.;
        let a = (years / 100.).trunc();
        let b = 2. - a + (a / 4.).trunc();
        let julian_year =
            (365.25 * years).trunc() + ((30.6001 * 14.) as f64).trunc() + 1720994.5 + b;
        let day = base_time.ordinal();
        let julian_day = julian_year + day as f64;
        let j2000_day = julian_day - 2451545.0;
        let offset_time = base_time.num_seconds_from_midnight() as f64;
        //Adjust for UT1 off of nutation etc tablee
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
}



