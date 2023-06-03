//Helpers
use crate::{helpers::modulus, GroundStation, Satellite};
use chrono::{Datelike, NaiveDateTime, Timelike, Utc};
use std::f64::consts::PI;

pub const F: f64 = 1. / 298.257223563;
#[derive(Clone, Copy)]
pub struct Eci {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
}
impl Eci {
    pub fn new_from_point(sat: &Satellite, offset: f64) -> Eci {
        let consts = &sat.constants;
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
}

pub struct SatAngle {
    pub elevation: f64,
    pub azimuth: f64,
    pub range: f64,
}
impl SatAngle {
    pub fn get_current(ground: &GroundStation, sat: &Satellite) -> SatAngle {
        //For GS driver
        let offset = Utc::now()
            .naive_utc()
            .signed_duration_since(sat.epoch)
            .num_seconds() as f64
            / 60.;
        let eci = Eci::new_from_point(sat, offset);
        Self::compute_angle(ground, eci, &sat.epoch, offset)
    }
    pub fn get_angle_offset(sat: &Satellite, station: &GroundStation, time_step: f64) -> SatAngle {
        let point = Eci::new_from_point(sat, time_step);
        SatAngle::compute_angle(station, point, &sat.epoch, time_step)
    }
    fn compute_angle(
        ground: &GroundStation,
        loc: Eci,
        epoch: &NaiveDateTime,
        offset: f64,
    ) -> SatAngle {
        let sidereal = SiderealAngle::find_sidereal_angle(epoch, offset);
        let observer = Self::xyz_from_lla(ground, &sidereal);
        let angle_set = Self::compute_look(observer, loc, &sidereal, ground.lat, ground.long);
        SatAngle {
            elevation: angle_set[0],
            azimuth: angle_set[1],
            range: angle_set[2],
        }
    }
    fn xyz_from_lla(loc: &GroundStation, sidereal_angle: &SiderealAngle) -> [f64; 3] {
        //Alt in meters, lat/long in decimal degrees
        let true_alt = 6378.137 + loc.alt / 1000.;
        let c = 1. / (1. + F * (F - 2.) * (loc.lat.to_radians().sin()).powf(2.)).sqrt();
        let s = (1. - F).powf(2.) * c;
        let z = true_alt * loc.lat.to_radians().sin() * s;
        let r = true_alt * loc.lat.to_radians().cos() * c;
        let theta_true = sidereal_angle.angle.to_degrees() + loc.long;
        let x = r * theta_true.to_radians().cos();
        let y = r * theta_true.to_radians().sin();
        [x, y, z]
    }

    fn compute_look(
        observer: [f64; 3],
        satellite: Eci,
        base_side_angle: &SiderealAngle,
        latitude: f64,
        longitude: f64,
    ) -> [f64; 3] {
        let rx = satellite.x - observer[0];
        let ry = satellite.y - observer[1];
        let rz = satellite.z - observer[2];
        let rad_lat = latitude.to_radians();
        let side_angle = modulus(
            (base_side_angle.angle.to_degrees() + longitude).to_radians(),
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
        [elevation, azimuth, range]
    }
}
pub struct SiderealAngle {
    pub angle: f64,
}
impl SiderealAngle {
    pub fn find_sidereal_angle(base_time: &NaiveDateTime, offset: f64) -> SiderealAngle {
        let time = Self::time_to_jd(*base_time, offset);
        SiderealAngle {
            angle: Self::sidereal_angle(time),
        }
    }
    fn time_to_jd(base_time: NaiveDateTime, offset: f64) -> [f64; 2] {
        //Meeus' approach from celestrak https://celestrak.org/columns/v02n02/
        let years = base_time.date().year() as f64 - 1.;
        let a = (years / 100.).trunc();
        let b = 2. - a + (a / 4.).trunc();
        let julian_year =
            (365.25 * years).trunc() + ((30.6001 * 14.) as f64).trunc() + 1720994.5 + b;
        let day = base_time.date().ordinal();
        let julian_day = julian_year + day as f64;
        let j2000_day = julian_day - 2451545.0;
        let offset_time = base_time.num_seconds_from_midnight() as f64 + offset * 60.;
        [j2000_day, offset_time]
        //Adjust for UT1 off of nutation etc table
    }

    fn sidereal_angle(julian: [f64; 2]) -> f64 {
        //based on
        // https://celestrak.org/columns/v02n02/
        //From https://celestrak.org/publications/AIAA/2006-6753/AIAA-2006-6753-Rev2.pdf
        let julian_day = julian[0];
        let offset = julian[1] / 86400.;
        let t = julian_day / 36525.0;
        let theta_0 = 24110.54841 + 8640184.812866 * t + 0.093104 * t * t
            - t * t * t * 6.2 * 10_f64.powf(-6.);
        let side_time = modulus(theta_0 + 1.00273790934 * offset * 86400., 86400.); //Why the 1.00?
                                                                                    //println!{"{}",side_time} //Make constant
        2. * PI * side_time / 86400. //In radian - returns the angle
                                     //println!("{}",theta);
    }
}
#[cfg(test)]
mod tests {
    use crate::{
        drivers::{Eci, SatAngle},
        GroundStation,
    };
    use chrono::{NaiveDate, NaiveDateTime, NaiveTime};
    #[test]
    fn test_look_angle() {
        let angle = SatAngle::compute_angle(
            &GroundStation {
                lat: 45.,
                long: -93.,
                alt: 0.,
                name: "test".to_string(),
            },
            Eci {
                x: -4400.594,
                y: 1932.870,
                z: 4760.712,
                vx: 0.,
                vy: 0.,
                vz: 0.,
            },
            &NaiveDateTime::new(
                NaiveDate::from_ymd_opt(1995, 11, 18).unwrap(),
                NaiveTime::from_num_seconds_from_midnight_opt(45960, 0).unwrap(),
            ),
            0.,
        );
        assert_eq!(81.5182305805702, angle.elevation);
        assert_eq!(100.35899093320114, angle.azimuth);
    }
}
