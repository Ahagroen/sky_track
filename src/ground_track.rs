use crate::{
    drivers::{Eci, SiderealAngle, F},
    helpers::modulus,
    Satellite,
};
use chrono::{NaiveDateTime, Utc};
const A: f64 = 6378.135;
pub struct Track {
    pub points: Vec<SubPoint>, //Second field is minutes after epoch, data is polled per minute
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
#[cfg(test)]
mod tests {
    use chrono::{NaiveDate, NaiveDateTime, NaiveTime};

    use crate::{drivers::Eci, ground_track::SubPoint};

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
}
