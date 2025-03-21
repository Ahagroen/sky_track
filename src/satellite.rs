use std::f64::consts::PI;

use serde::{Serialize,Deserialize};
use crate::{helpers::modulus, types::{Eci, SatAngle, SubPoint, A, F}, GroundStation};
use chrono::{DateTime, Datelike, NaiveDateTime, Timelike, Utc};
use sgp4::Constants;


#[derive(Serialize, Deserialize, Clone)]
pub struct Satellite {
    constants: Constants,
    name: String,
    orbital_period:f64,
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
        let period = 1./elem.mean_motion*86400.;
        let constants = sgp4::Constants::from_elements(&elem).unwrap();//need to get rid of this?
        Satellite {
            constants,
            orbital_period:period,
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

    fn xyz_from_lla(loc: &GroundStation, sidereal_angle: f64) -> [f64; 3] {
        let e2 = 2.*F - F.powf(2.0);
        let lat_true = (1.0 / (1.0 - e2)*loc.lat.to_radians().tan()).atan();
        let theta_true = sidereal_angle.to_degrees() + loc.long;
        let n = A/(1.-e2*(loc.lat.to_radians().sin()).powf(2.)).sqrt();
        let n_z = n*(1.-e2);
        let x = (n+loc.alt/1000.)*lat_true.cos()*theta_true.to_radians().cos();
        let y = (n+loc.alt/1000.)*lat_true.cos()*theta_true.to_radians().sin();
        let z = (n_z+loc.alt/1000.)*loc.lat.to_radians().sin();
        [x, y, z]
    }

    pub fn get_look_angle(&self, station: &GroundStation, offset: i64) -> SatAngle {
        let time_stamp= self.epoch.and_utc().timestamp()+offset;
        let satellite = self.get_point_eci(offset);
        Self::get_look_angle_direct(station, time_stamp, satellite)
    }

    fn get_look_angle_direct(station: &GroundStation, time_stamp: i64, satellite: Eci) -> SatAngle {
        let sidereal = Self::sidereal_angle(time_stamp);
        let observer = Self::xyz_from_lla(station, sidereal);
        let rx = satellite.x - observer[0];
        let ry = satellite.y - observer[1];
        let rz = satellite.z - observer[2];
        let rad_lat = station.lat.to_radians(); 
        let side_angle = modulus((sidereal.to_degrees() + station.long).to_radians(),2.*PI);
        let rs = rad_lat.sin() * side_angle.cos() * rx + rad_lat.sin() * side_angle.sin() * ry
            - rad_lat.cos() * rz;
        let re = -side_angle.sin() * rx + side_angle.cos() * ry;
        let rz = rad_lat.cos() * side_angle.cos() * rx
            + rad_lat.cos() * side_angle.sin() * ry
            + rad_lat.sin() * rz;
        let range = (rs * rs + re * re + rz * rz).sqrt();
        let elevation = ((rz / range).asin()).to_degrees();
        let mut azimuth = (-re/rs).atan();
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
        let time_stamp= self.epoch.and_utc().timestamp()+offset;
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
        let base_time = DateTime::from_timestamp(time_stamp,0).unwrap();
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
    pub fn seconds_since_epoch(&self,other:&DateTime<Utc>)->i64{
        other.timestamp()-self.epoch.and_utc().timestamp()
    }
    pub fn offset_timestamp(&self,other:i64)->i64{
        other-self.epoch.and_utc().timestamp()
    }
    pub fn get_orbital_period(&self)->f64{
        return self.orbital_period;
    }
    pub fn get_epoch(&self) -> DateTime<Utc> {
        self.epoch.and_utc()
    }
}


#[cfg(test)]
mod tests{
    use crate::helpers::{assert_almost_eq, quick_gen_datetime};

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
    fn test_sidereal_angle(){
        let time_stamp= quick_gen_datetime(1995, 10, 1, 9, 0, 0).timestamp();
        let sidereal_angle = Satellite::sidereal_angle(time_stamp);
        assert_almost_eq(sidereal_angle,2.524218);
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
        let date =quick_gen_datetime(2025, 03, 13,21, 54, 42);
        let reference_sub_point = sat.get_sub_point(sat.seconds_since_epoch(&date));
        assert_almost_eq(reference_sub_point.long,-63.7546);
        assert_almost_eq(reference_sub_point.lat,35.8798);
        assert_almost_eq(reference_sub_point.alt,421.1125)
    }
    #[test]
    fn test_gs_position(){
        let refence_gs = &GroundStation::new([40.,-75.,0.], "Test");
        let time_stamp= quick_gen_datetime(1995, 10, 1, 9, 0, 0).timestamp();
        let sidereal = Satellite::sidereal_angle(time_stamp);
        assert_almost_eq(Satellite::xyz_from_lla(refence_gs, sidereal)[2],4077.9855);
        assert_almost_eq(Satellite::xyz_from_lla(refence_gs, sidereal)[0],1703.295618);
        assert_almost_eq(Satellite::xyz_from_lla(refence_gs, sidereal)[1],4586.6514);
    }
    #[test]
    fn test_look_angle_direct(){
        let refence_gs = &GroundStation::new([40.,-93.,0.], "Test");
        let time_stamp= quick_gen_datetime(1995, 11, 18, 12, 46, 0).timestamp();
        let look_angles = Satellite::get_look_angle_direct(&refence_gs, time_stamp, Eci { x: -4400.594, y: 1932.870, z: 4760.712, vx: 0., vy: 0., vz: 0. });
        assert_almost_eq(look_angles.azimuth, 100.3600);
        assert_almost_eq(look_angles.elevation, 81.5200);
    }
    #[test]
    fn test_look_angle(){
        let sat = Satellite::new_from_tle(
            "ISS (ZARYA)             
1 25544U 98067A   25078.36999458  .00023040  00000+0  41584-3 0  9998
2 25544  51.6365  31.8868 0003892  28.0409 332.0788 15.49628144501233",
        );
        let date = quick_gen_datetime(2025, 03, 19,21, 00, 04);
        let refence_look_angle = sat.get_look_angle(&GroundStation::new([51.9861,4.3876,74.4], "Delft"), sat.seconds_since_epoch(&date));
        assert_almost_eq(refence_look_angle.azimuth, 289.23516);
        assert_almost_eq(refence_look_angle.elevation, -80.286);
        assert_almost_eq(refence_look_angle.range, 12991.4);//off by 1km
    }
}