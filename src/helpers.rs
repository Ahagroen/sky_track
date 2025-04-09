use std::f64;

pub fn modulus(a: f64, b: f64) -> f64 {
    ((a % b) + b) % b
}

#[cfg(test)]
pub fn assert_almost_eq<T>(a:T,b:T)
where T: std::ops::Sub+PartialOrd<f64>+std::fmt::Debug+PartialOrd<T>+Copy,
<T as std::ops::Sub>::Output: PartialOrd<f64>
{
    if a-b<1e-5{
        assert!(true)
    } else{
        assert_eq!(a,b)
    }
}

#[cfg(test)]
use chrono::{DateTime,Utc};
#[cfg(test)]
pub fn quick_gen_datetime(y:i32,m:u32,d:u32,h:u32,min:u32,s:u32)->DateTime<Utc>{
    use chrono::{NaiveDate, NaiveDateTime, NaiveTime};

    let date = NaiveDate::from_ymd_opt(y, m, d).unwrap();
    let time = NaiveTime::from_hms_opt(h, min, s).unwrap();
    DateTime::from_naive_utc_and_offset(NaiveDateTime::new(date, time),Utc)
}