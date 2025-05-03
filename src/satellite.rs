use crate::{
    GroundStation,
    helpers::modulus,
    types::{A, B, Eci, F, SatAngle, SubPoint},
};
use chrono::{DateTime, Datelike, NaiveDateTime, Timelike, Utc};
use serde::{Deserialize, Serialize};
use sgp4::Constants;
use std::{cell::LazyCell, collections::HashMap, f64::consts::PI, io::BufRead};

#[derive(Serialize, Deserialize, Clone)]
pub struct Satellite {
    constants: Constants,
    name: String,
    orbital_period: f64,
    epoch: NaiveDateTime, //Epoch of the structs TLE
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
        let period = 1. / elem.mean_motion * 86400.;
        let constants = sgp4::Constants::from_elements(&elem).unwrap(); //need to get rid of this?
        Satellite {
            constants,
            orbital_period: period,
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
    pub fn get_point_eci(&self, offset: i64) -> Eci {
        let consts = &self.constants;
        let prop = consts.propagate(offset as f64 / 60.).unwrap();
        Eci {
            x: prop.position[0],
            y: prop.position[1],
            z: prop.position[2],
            vx: prop.velocity[0],
            vy: prop.velocity[1],
            vz: prop.velocity[2],
        }
    }

    fn lla_to_ecef(loc: [f64; 3], sidereal_angle: f64, time_stamp: i64) -> [f64; 3] {
        let initial_condition = Self::ecef_xyz_from_lla([loc[0], loc[1], loc[2]], sidereal_angle);
        let adjustments = Self::compute_precession_nutation(time_stamp);
        Self::rotate_precession_nuation(initial_condition, adjustments)
    }
    fn ecef_xyz_from_lla(loc: [f64; 3], sidereal_angle: f64) -> [f64; 3] {
        let c = 1. / (1. + F * (F - 2.) * loc[0].to_radians().sin().powi(2)).sqrt();
        let s = (1. - F).powi(2) * c;
        let local_sidereal = loc[1].to_radians() + sidereal_angle;
        let ecef_x =
            (A + loc[2] / 1000.0) * c * (loc[0].to_radians().cos()) * (local_sidereal.cos());
        let ecef_y =
            (A + loc[2] / 1000.0) * c * (loc[0].to_radians().cos()) * (local_sidereal.sin());
        let ecef_z = (A + loc[2] / 1000.0) * s * (loc[0].to_radians().sin());
        [ecef_x, ecef_y, ecef_z]
    }
    fn compute_precession_nutation(time_stamp: i64) -> (f64, f64) {
        //https://www2.mps.mpg.de/homes/fraenz/systems/systems2art/node3.html
        let (j2000_day, _) = Self::get_julian_date(time_stamp);
        let t0 = j2000_day / 36525.0;
        let epsilon_j2k = 23.439291111;
        let epsilon_date =
            epsilon_j2k - 0.013004167 * t0 - 0.000000164 * t0.powi(2) - 0.000000504 * t0.powi(3);
        let del_epsilon = 0.0026 * (125. - 0.05295 * j2000_day).to_radians().cos()
            + 0.0002 * (200.9 + 1.97129 * j2000_day).to_radians().cos();
        let obliquity = epsilon_date + del_epsilon;
        let long_nutation = -0.0048 * (125.0 - 0.05295 * j2000_day).to_radians().sin()
            - 0.0004 * (200.9 + 1.97129 * j2000_day).to_radians().sin();
        (obliquity.to_radians(), long_nutation.to_radians())
    }
    fn rotate_precession_nuation(
        initial_state: [f64; 3],
        precession_nutation: (f64, f64),
    ) -> [f64; 3] {
        let first_x = initial_state[0];
        let first_y = initial_state[1] * (precession_nutation.0).cos()
            + initial_state[2] * (precession_nutation.0).sin();
        let first_z = -initial_state[1] * (precession_nutation.0).sin()
            + initial_state[2] * (precession_nutation.0).cos();
        let second_x = first_x * (-1. * precession_nutation.1).cos()
            + first_y * (-1. * precession_nutation.1).sin();
        let second_y = -1. * first_x * (-1. * precession_nutation.1).sin()
            + first_y * (-1. * precession_nutation.1).cos();
        let second_z = first_z;
        let third_x = second_x;
        let third_y =
            second_y * (-precession_nutation.0).cos() + second_z * (-precession_nutation.0).sin();
        let third_z =
            -second_y * (-precession_nutation.0).sin() + second_z * (-precession_nutation.0).cos();
        [third_x, third_y, third_z]
    }

    ///Computes the look angle from the station to the satellite, accounting for refractive effects. Returns None if the station cannot actually see the satellite.
    pub fn get_look_angle_refraction(&self, station: &GroundStation, offset: i64) -> SatAngle {
        let base_look_angle = self.get_look_angle(station, offset);
        let min_elv = -0.875 * (station.alt).sqrt();
        if base_look_angle.elevation > 10. {
            return base_look_angle;
        } else {
            let adjustment_factor = 1.
                / (1.728
                    + 0.5411 * base_look_angle.elevation
                    + 0.03723 * base_look_angle.elevation.powi(2)
                    + station.alt
                        * (0.1815
                            + 0.06272 * base_look_angle.elevation
                            + 0.011380 * base_look_angle.elevation.powi(2))
                    + station.alt.powi(2) * (0.01727 + 0.008288 * base_look_angle.elevation));
            SatAngle {
                azimuth: base_look_angle.azimuth,
                elevation: base_look_angle.elevation + adjustment_factor,
                range: base_look_angle.range,
            }
        }
    }

    /// Compute the current look angle between the satellite and the ground. NOTE: computes the free space look angle, and does not account for refraction.
    /// Use the functions concerning pass propagation to include atmospheric refraction
    ///
    pub fn get_look_angle(&self, station: &GroundStation, offset: i64) -> SatAngle {
        let time_stamp = self.epoch.and_utc().timestamp() + offset;
        let updated_ts = time_stamp; // + Self::get_delta_t(time_stamp) as i64;
        let satellite = self.get_point_eci(offset);
        Self::get_look_angle_direct(station, updated_ts, [satellite.x, satellite.y, satellite.z])
    }

    fn get_look_angle_direct(
        station: &GroundStation,
        time_stamp: i64,
        satellite: [f64; 3],
    ) -> SatAngle {
        let sidereal = Self::sidereal_angle(time_stamp);
        let observer = Self::lla_to_ecef(
            [station.lat, station.long, station.alt],
            sidereal,
            time_stamp,
        );
        let local_sidereal = modulus(station.long.to_radians() + sidereal, 2.0 * PI);
        let rx = satellite[0] - observer[0];
        let ry = satellite[1] - observer[1];
        let rz = satellite[2] - observer[2];
        let rs = station.lat.to_radians().sin() * local_sidereal.cos() * rx
            + station.lat.to_radians().sin() * local_sidereal.sin() * ry
            - station.lat.to_radians().cos() * rz;
        let re = -1. * local_sidereal.sin() * rx + local_sidereal.cos() * ry;
        let rd = station.lat.to_radians().cos() * local_sidereal.cos() * rx
            + station.lat.to_radians().cos() * local_sidereal.sin() * ry
            + station.lat.to_radians().sin() * rz;
        let range = (rx.powi(2) + ry.powi(2) + rz.powi(2)).sqrt();
        let elevation = ((rd / range).asin()).to_degrees();
        let mut azimuth = (-re / rs).atan();
        if rs > 0. {
            azimuth += PI
        }
        if azimuth < 0. {
            azimuth += 2. * PI
        }
        azimuth = azimuth.to_degrees();
        SatAngle {
            elevation,
            azimuth,
            range,
        }
    }
    ///compute the sub_point at a given time since epoch in seconds.
    pub fn get_sub_point(&self, offset: i64) -> SubPoint {
        fn get_lat_and_alt(satellite: &Eci) -> [f64; 2] {
            let mut guess =
                (satellite.z).atan2((satellite.x.powf(2.) + satellite.y.powf(2.)).sqrt());
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
        fn get_long(satellite: &Eci, sidereal_angle: f64) -> f64 {
            let data = (satellite.y.atan2(satellite.x) - sidereal_angle).to_degrees();
            println!("data: {}", data);
            if data < -180. {
                return data + 360.;
            } else if data > 180. {
                return data - 360.;
            } else {
                data
            }
            //println!("Angle = {}",constrained_angle);
        }
        let time_stamp = self.epoch.and_utc().timestamp() + offset;
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
    fn get_julian_date(time_stamp: i64) -> (f64, f64) {
        let base_time = DateTime::from_timestamp(time_stamp, 0).unwrap();
        let years = base_time.year() as f64 - 1.;
        let a = (years / 100.).trunc();
        let b = 2. - a + (a / 4.).trunc();
        let julian_year =
            (365.25 * years).trunc() + ((30.6001 * 14.) as f64).trunc() + 1720994.5 + b;
        let day = base_time.ordinal();
        let julian_day = julian_year + day as f64;
        let j2000_day = julian_day - 2451545.0;
        let offset_time = base_time.num_seconds_from_midnight() as f64;
        (j2000_day, offset_time)
    }
    fn sidereal_angle(time_stamp: i64) -> f64 {
        //Meeus' approach from celestrak https://celestrak.org/columns/v02n02/
        //Adjust for UT1 off of nutation etc table
        //based on
        // https://celestrak.org/columns/v02n02/
        //From https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev2.pdf
        let (j2000_day, offset_time) = Self::get_julian_date(time_stamp);
        let offset = offset_time / 86400.;
        let t = j2000_day / 36525.0;
        let theta_0 = 24110.54841 + 8640184.812866 * t + 0.093104 * t * t
            - t * t * t * 6.2 * 10_f64.powf(-6.);
        let side_time = (theta_0 + 1.00273790934 * offset * 86400.) % 86400.; //Why the 1.00?
        //println!{"{}",side_time} //Make constant
        2. * PI * side_time / 86400.
    }
    pub fn seconds_since_epoch(&self, other: &DateTime<Utc>) -> i64 {
        other.timestamp() - self.epoch.and_utc().timestamp()
    }
    pub fn offset_timestamp(&self, other: i64) -> i64 {
        other - self.epoch.and_utc().timestamp()
    }
    pub fn get_orbital_period(&self) -> f64 {
        self.orbital_period
    }
    pub fn get_epoch(&self) -> DateTime<Utc> {
        self.epoch.and_utc()
    }
}

fn load_deltat() -> HashMap<i64, f64> {
    let delta_t_data = include_bytes!("deltat.txt");
    let mut output: HashMap<i64, f64> = Default::default();
    for i in delta_t_data.lines() {
        let line_data = i.unwrap();
        let line: Vec<&str> = line_data.split_whitespace().collect();
        output.insert(
            line[0].parse::<i64>().unwrap() * 100 + line[1].parse::<i64>().unwrap(),
            line[3].parse::<f64>().unwrap(),
        );
    }
    output
}
#[cfg(test)]
mod tests {

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
    fn test_sidereal_angle() {
        let time_stamp = quick_gen_datetime(1995, 10, 1, 9, 0, 0).timestamp();
        let sidereal_angle = Satellite::sidereal_angle(time_stamp);
        assert_almost_eq(sidereal_angle, 2.524218);
    }
    #[test]
    fn test_eci() {
        let sat = Satellite::new_from_tle(
            "RUST_SGP4(TESTING)               
            1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753
            2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667",
        );
        let eci_point = sat.get_point_eci(0);
        assert_almost_eq(eci_point.x, 7022.46647);
        assert_almost_eq(eci_point.y, -1400.06656);
        assert_almost_eq(eci_point.z, 0.05106);
    }
    #[test]
    fn test_sub_point() {
        let sat = Satellite::new_from_tle(
            "ISS (ZARYA)             
            1 25544U 98067A   25072.43808874  .00018974  00000+0  33994-3 0  9997
            2 25544  51.6354  61.2721 0006420  16.6184 343.5014 15.49959635500318",
        );
        let date = quick_gen_datetime(2025, 03, 13, 21, 54, 42);
        let reference_sub_point = sat.get_sub_point(sat.seconds_since_epoch(&date));
        assert_almost_eq(reference_sub_point.long, -63.7546);
        assert_almost_eq(reference_sub_point.lat, 35.8798);
        assert_almost_eq(reference_sub_point.alt, 421.112511)
    }
    #[test]
    fn test_gs_position() {
        let time_stamp = quick_gen_datetime(1995, 10, 1, 9, 0, 0).timestamp();
        let sidereal = Satellite::sidereal_angle(time_stamp);
        assert_almost_eq(
            Satellite::ecef_xyz_from_lla([40., -75., 0.], sidereal)[2],
            4077.985572,
        );
        assert_almost_eq(
            Satellite::ecef_xyz_from_lla([40., -75., 0.], sidereal)[0],
            1703.295618,
        );
        assert_almost_eq(
            Satellite::ecef_xyz_from_lla([40., -75., 0.], sidereal)[1],
            4586.651468,
        );
    }
    #[test]
    fn test_eci_from_lla() {
        let refence_gs = &GroundStation::new([40., -93., 0.], "Test");
        let time_stamp = quick_gen_datetime(1995, 11, 18, 12, 46, 0).timestamp();
        let sidereal_angle = Satellite::sidereal_angle(time_stamp);
        let output = Satellite::ecef_xyz_from_lla(
            [refence_gs.lat, refence_gs.long, refence_gs.alt],
            sidereal_angle,
        );
        assert_almost_eq(output[0], 1703.295);
        assert_almost_eq(output[1], 4586.650);
        assert_almost_eq(output[2], 4077.98557);
    }
    #[test]
    fn test_look_angle_direct() {
        let refence_gs = &GroundStation::new([40., -93., 0.], "Test");
        let time_stamp = quick_gen_datetime(1995, 11, 18, 12, 46, 0).timestamp();
        let updated_ts = time_stamp;
        let look_angles = Satellite::get_look_angle_direct(
            &refence_gs,
            updated_ts,
            [-4400.594, 1932.870, 4760.712],
        );
        assert_almost_eq(look_angles.azimuth, 100.3600);
        assert_almost_eq(look_angles.elevation, 81.5200);
    }
    #[test]
    fn test_look_angle() {
        let sat = Satellite::new_from_tle(
            "ISS (ZARYA)             
1 25544U 98067A   25122.54440123  .00015063  00000+0  27814-3 0  9994
2 25544  51.6345 173.1350 0002187  74.2134 285.9096 15.49297959508085",
        );
        let date = quick_gen_datetime(2025, 05, 2, 21, 0, 0);
        let refence_look_angle = sat.get_look_angle(
            &GroundStation::new([51.9861, 4.3876, 0.], "Delft"),
            sat.seconds_since_epoch(&date),
        );
        assert_almost_eq(refence_look_angle.azimuth, 128.3);
        assert_almost_eq(refence_look_angle.elevation, -3.896);
        assert_almost_eq(refence_look_angle.range, 2824.966); //off by 3km
    }
    #[test]
    fn test_look_angle_2() {
        let sat = Satellite::new_from_tle(
            "ISS (ZARYA)             
1 25544U 98067A   25122.54440123  .00015063  00000+0  27814-3 0  9994
2 25544  51.6345 173.1350 0002187  74.2134 285.9096 15.49297959508085",
        );
        let date = quick_gen_datetime(2025, 05, 2, 22, 0, 0);
        let refence_look_angle = sat.get_look_angle(
            &GroundStation::new([51.9861, 4.3876, 0.], "Delft"),
            sat.seconds_since_epoch(&date),
        );
        println!("{:?}", sat.get_sub_point(sat.seconds_since_epoch(&date)));
        assert_almost_eq(refence_look_angle.azimuth, 238.8);
        assert_almost_eq(refence_look_angle.elevation, -66.);
        assert_almost_eq(refence_look_angle.range, 12104.405); //off by 3km
    }
}
