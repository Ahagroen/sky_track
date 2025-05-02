mod astro_adjust{
    pub fn load_precession_nutation(time_stamp:i64)->(f64,f64,f64){
        let base_data = include_bytes!("EOP_20_C04_0h_dPsi_dEps_1962-now.txt");
        let reader = Cursor::new(base_data);
        let target = DateTime::from_timestamp(time_stamp, 0).unwrap();
        let initial_epsilon = -23.4392911; //this is the J2000 epsilon value (I think)
        let mut carry_epsilon:f64 = initial_epsilon;
        let mut dpsi:f64 = 0.;
        let mut deps:f64 = 0.;
        for line in reader.lines() {
            let line = line.unwrap();
            if line.chars().peekable().peek().unwrap().to_owned() == '#'{
                continue;
            }
            let columns: Vec<&str> = line.split_whitespace().collect();
        // Parse each column into the desired data type
            let year: i32 = columns[0].parse().unwrap();
            let month: i32 = columns[1].parse().unwrap();
            let day: i32 = columns[2].parse().unwrap();
            if year < target.year(){
                carry_epsilon += columns[9].parse::<f64>().unwrap();
            } else if year == target.year(){
                if month < target.month() as i32{
                    carry_epsilon += columns[9].parse::<f64>().unwrap();
                } else if month == target.month() as i32{
                    if day < target.day() as i32{
                        carry_epsilon += columns[9].parse::<f64>().unwrap();
                    } else if day == target.day() as i32{
                        dpsi = columns[8].parse().unwrap();
                        deps = columns[9].parse().unwrap();
                        break
                    }
                }
            }
        }
        (dpsi,deps,carry_epsilon)
    }

}
