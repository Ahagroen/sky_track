use crate::helpers::modulus;
use chrono::Datelike;
use chrono::{NaiveDateTime, Timelike, Utc};
use sgp4::Constants;
use std::f64::consts::PI;
mod helpers;

const F: f64 = 1. / 298.257223563;
const A: f64 = 6378.135;

pub struct Passes {
    pub pass_list: Vec<PassTime>,
}
impl Passes {
    pub fn find_passes(
        gs: Vec<&GroundStation>,
        sat: &Satellite,
        max_offset: i64,
        start_time: NaiveDateTime,
    ) -> Passes {
        let mut passes = Vec::new();
        let start_delta = start_time.signed_duration_since(sat.epoch);
        let offset = start_delta.num_seconds() / 60;
        //println!("Offset {}",offset);
        for station in gs {
            passes.append(&mut Self::search_for_pass(sat, station, max_offset, offset))
        }
        Passes { pass_list: passes }
    }
    fn search_for_pass(
        satellite: &Satellite,
        ground: &GroundStation,
        range: i64,
        base_offset: i64,
    ) -> Vec<PassTime> {
        let mut last_pass_time: u16 = 6000;
        let mut usable_pass = Vec::new();
        let mut pass = false;
        for i in base_offset..(range + base_offset) {
            if !(120..=480).contains(&last_pass_time) {
                let point = Eci::new_from_point(satellite, i as f64);
                let angle = SatAngle::compute_angle(ground, point, &satellite.epoch, i as f64);
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

    pub fn get_next_pass(&self) -> &PassTime {
        self.pass_list.first().unwrap()
    }
    pub fn propagate(&mut self) {
        if Utc::now()
            .naive_utc()
            .signed_duration_since(self.pass_list[0].los)
            > chrono::Duration::zero()
        {
            self.pass_list.remove(0);
        }
    }
}
pub struct PassTime {
    pub aos: NaiveDateTime,
    pub los: NaiveDateTime,
    pub length: f64,
    pub max_elevation: f64,
    pub ground_station: GroundStation,
}
impl PassTime {
    fn new(sat: &Satellite, guess: i64, station: &GroundStation) -> PassTime {
        //println!("stop");
        let los: NaiveDateTime = Self::get_stop_time(sat, station, guess);
        //println!("start");
        let aos: NaiveDateTime = Self::get_start_time(sat, station, guess);
        let length = los.signed_duration_since(aos).num_milliseconds() as f64 / 1000.;
        //println!("{} Duration in seconds",length);
        let max_elevation = Self::get_highest_elevation(sat, station, length, aos);
        PassTime {
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
        let max_elevation = Self::get_angle(sat, station, mid_time); //Always at the midpoint??
        max_elevation.elevation
    }
    fn angle_bisection(left: f64, right: f64, sat: &Satellite, station: &GroundStation) -> f64 {
        //1. Compute midpoint elevation
        //2. set that to the new bound based on +- 5 degrees
        //3. re-compute
        let mut rs_carry = right;
        let mut ls_carry = left;
        let left_side_angle = Self::get_angle(sat, station, ls_carry);
        let low_to_high: bool = left_side_angle.elevation < 0.;
        loop {
            let midpoint = (rs_carry + ls_carry) / 2.;
            let angle = Self::get_angle(sat, station, midpoint);
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
    fn get_angle(sat: &Satellite, station: &GroundStation, time_step: f64) -> SatAngle {
        let point = Eci::new_from_point(sat, time_step);
        SatAngle::compute_angle(station, point, &sat.epoch, time_step)
    }
}

//DataSets
#[derive(Clone)]
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

pub struct Track {
    points: Vec<SubPoint>, //Second field is minutes after epoch, data is polled per minute
    base_offset: i64,
    pub last_gen: NaiveDateTime,
}
impl Track {
    pub fn get_track(duration: f64, satellite: &Satellite, start_offset: f64) -> Track {
        //Duration in seconds
        //Start offset is offset from current time
        let mut points = Vec::new();
        let mut delta = start_offset; //Starts at "base_time"
        let base_offset = Utc::now()
            .naive_utc()
            .signed_duration_since(satellite.epoch)
            .num_seconds();
        delta += base_offset as f64;
        let start_val = base_offset as f64;
        while delta < duration + base_offset as f64 + start_offset {
            //In seconds since base_time
            let loc = Eci::new_from_point(satellite, delta / 60.);
            let sub_point =
                SubPoint::compute_sub_point(loc, satellite.epoch, delta / 60., start_val / 60.);
            delta += 1.;
            points.push(sub_point)
        }
        let last_gen = Utc::now().naive_utc();
        Track {
            points,
            base_offset,
            last_gen,
        }
    }
    pub fn current_point(&self, offset: usize) -> &SubPoint {
        &self.points[offset]
    }
    pub fn get_point_array(&self) -> Vec<[f64; 4]> {
        let mut output = Vec::new();
        for i in &self.points {
            output.push([i.long, i.lat, i.alt, i.time as f64])
        }
        output
    }
    pub fn propagate(&mut self, duration: f64, satellite: &Satellite) {
        self.points.remove(0);
        self.base_offset = Utc::now()
            .naive_utc()
            .signed_duration_since(satellite.epoch)
            .num_seconds();
        let loc = Eci::new_from_point(satellite, (duration + self.base_offset as f64) / 60.);
        let sub_point = SubPoint::compute_sub_point(
            loc,
            satellite.epoch,
            (duration + self.base_offset as f64) / 60.,
            self.base_offset as f64,
        );
        self.last_gen = Utc::now().naive_utc();
        self.points.push(sub_point);
    }
}
pub struct SubPoint {
    pub lat: f64,
    pub long: f64,
    pub alt: f64,
    pub time: i64,
}
impl SubPoint {
    fn compute_sub_point(
        loc: Eci,
        base_time: NaiveDateTime,
        offset: f64,
        start_offset: f64,
    ) -> SubPoint {
        let longitude = Self::get_long(&loc, &base_time, offset);
        let lat_alt = Self::get_lat_and_alt(&loc);
        let latitude = lat_alt[0];
        let altitude = lat_alt[1];
        let time = offset * 60. - start_offset * 60.;
        SubPoint {
            lat: latitude,
            long: longitude,
            alt: altitude,
            time: time.round() as i64,
        }
    }
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
    fn get_long(satellite: &Eci, base_time: &NaiveDateTime, offset: f64) -> f64 {
        let sidereal_angle = SiderealAngle::find_sidereal_angle(base_time, offset);
        let angle = ((satellite.y.atan2(satellite.x)) - sidereal_angle.angle).to_degrees() + 180.;
        modulus(angle, 360.) - 180.
        //println!("Angle = {}",constrained_angle);
    }
}

//Helpers
#[derive(Clone, Copy)]
struct Eci {
    x: f64,
    y: f64,
    z: f64,
    vx: f64,
    vy: f64,
    vz: f64,
}
impl Eci {
    fn new_from_point(sat: &Satellite, offset: f64) -> Eci {
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
pub struct SatAngle {
    pub elevation: f64,
    pub azimuth: f64,
    pub range: f64,
}
impl SatAngle {
    pub fn get_current(ground: &GroundStation, sat: &Satellite) -> SatAngle {
        let offset = Utc::now()
            .naive_utc()
            .signed_duration_since(sat.epoch)
            .num_seconds() as f64
            / 60.;
        let eci = Eci::new_from_point(sat, offset);
        Self::compute_angle(ground, eci, &sat.epoch, offset)
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
struct SiderealAngle {
    angle: f64,
}
impl SiderealAngle {
    fn find_sidereal_angle(base_time: &NaiveDateTime, offset: f64) -> SiderealAngle {
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
    use super::Track;
    use crate::Eci;
    use crate::GroundStation;
    use crate::PassTime;
    use crate::Passes;
    use crate::SatAngle;
    use crate::Satellite;
    use crate::SubPoint;
    use chrono::NaiveDate;
    use chrono::NaiveDateTime;
    use chrono::NaiveTime;
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
    #[test]
    fn test_sub_point() {
        let sub_point = SubPoint::compute_sub_point(
            Eci {
                x: -4400.594,
                y: 1932.870,
                z: 4760.712,
                vx: 0.,
                vy: 0.,
                vz: 0.,
            },
            NaiveDateTime::new(
                NaiveDate::from_ymd_opt(1995, 11, 18).unwrap(),
                NaiveTime::from_num_seconds_from_midnight_opt(45960, 0).unwrap(),
            ),
            0.,
            0.,
        );
        assert_eq!(44.90766377931492, sub_point.lat); //Should be 44.91
        assert_eq!(-92.30530912969857, sub_point.long);
        assert_eq!(397.5072802741115, sub_point.alt);
    }
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
