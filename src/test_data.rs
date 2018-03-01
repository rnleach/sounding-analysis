//! Data used in tests.

use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use std::str::FromStr;

use metfor;
use sounding_base::{Profile, Sounding, StationInfo, Surface};
use sounding_validate::validate;

fn load_test_csv_sounding(
    location: &PathBuf,
) -> (Sounding, HashMap<String, Vec<f64>>, HashMap<String, i64>) {
    let mut f = File::open(location).expect(&format!("Error opening file: {:?}", location));

    let mut contents = String::new();
    f.read_to_string(&mut contents)
        .expect(&format!("Error reading file: {:?}", location));

    let lines: Vec<&str> = contents.split('\n').collect();
    let mut line_iter = lines.iter();

    //
    // Parse profile data
    //
    let mut height: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut temp: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut dp: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut press: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut wspd: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut wdir: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut wet_bulb: Vec<Option<f64>> = Vec::with_capacity(lines.len());

    for line in line_iter.by_ref() {
        if line.starts_with("### Surface Data ###") {
            break;
        }

        let tokens: Vec<&str> = line.split(',').collect();
        if tokens.len() < 6 {
            continue;
        }
        let t_c = f64::from_str(tokens[1]).ok();
        let dp_c = f64::from_str(tokens[2]).ok();
        let press_hpa = f64::from_str(tokens[3]).ok();

        height.push(f64::from_str(tokens[0]).ok());
        temp.push(t_c);
        dp.push(dp_c);
        press.push(press_hpa);
        wspd.push(f64::from_str(tokens[4]).ok());
        wdir.push(f64::from_str(tokens[5]).ok());

        if let (Some(t_c), Some(dp_c), Some(press_hpa)) = (t_c, dp_c, press_hpa) {
            wet_bulb.push(metfor::wet_bulb_c(t_c, dp_c, press_hpa).ok());
        } else {
            wet_bulb.push(None);
        }
    }

    let mut snd = Sounding::new()
        .set_profile(Profile::GeopotentialHeight, height)
        .set_profile(Profile::Temperature, temp)
        .set_profile(Profile::DewPoint, dp)
        .set_profile(Profile::Pressure, press)
        .set_profile(Profile::WetBulb, wet_bulb)
        .set_profile(Profile::WindSpeed, wspd)
        .set_profile(Profile::WindDirection, wdir);

    //
    // Surface data
    //
    for line in line_iter.by_ref() {
        if line.starts_with("### Analysis Int Section ###") {
            break;
        }

        let tokens: Vec<&str> = line.split(',').collect();
        if tokens.len() < 6 {
            continue;
        }

        let height = f64::from_str(tokens[0]).ok();
        let t_c = f64::from_str(tokens[1]).ok();
        let dp_c = f64::from_str(tokens[2]).ok();
        let press_hpa = f64::from_str(tokens[3]).ok();
        let _wspd = f64::from_str(tokens[4]).ok();
        let _wdir = f64::from_str(tokens[5]).ok();

        snd = snd
            .set_surface_value(Surface::Temperature, t_c)
            .set_station_info(StationInfo::new().with_elevation(height))
            .set_surface_value(Surface::DewPoint, dp_c)
            .set_surface_value(Surface::StationPressure, press_hpa)
            // FIXME: set these!
            .set_surface_value(Surface::UWind, None)
            .set_surface_value(Surface::VWind, None)
    }

    //
    // Integer values.
    //
    let mut target_int_vals = HashMap::new();
    for line in line_iter.by_ref() {
        if line.starts_with("### Analysis Float Section ###") {
            break;
        }

        let tokens: Vec<String> = line.split(',')
            .filter_map(|val| {
                let v = val.trim();
                if v != "" {
                    Some(v.to_owned())
                } else {
                    None
                }
            })
            .collect();

        if tokens.len() < 2 {
            continue;
        }

        let key = tokens[0].to_owned();
        let value = i64::from_str(&tokens[1]).unwrap();

        target_int_vals.insert(key, value);
    }

    // Float values.
    let mut target_float_vals = HashMap::new();
    for line in line_iter.by_ref() {
        let tokens: Vec<String> = line.split(',')
            .filter_map(|val| {
                let v = val.trim();
                if v != "" {
                    Some(v.to_owned())
                } else {
                    None
                }
            })
            .collect();
        if tokens.len() < 2 {
            continue;
        }
        let key = tokens[0].to_owned();
        let mut values = vec![];
        for token in tokens.iter().skip(1) {
            values.push(f64::from_str(token).unwrap());
        }
        target_float_vals.insert(key, values);
    }

    (snd, target_float_vals, target_int_vals)
}

pub fn load_test_file(fname: &str) -> (Sounding, HashMap<String, Vec<f64>>, HashMap<String, i64>) {
    let mut test_path = PathBuf::new();
    test_path.push("test_data");
    test_path.push(fname);
    load_test_csv_sounding(&test_path)
}

pub fn approx_equal(tgt: f64, guess: f64, tol: f64) -> bool {
    use std::f64;

    assert!(tol > 0.0);

    f64::abs(tgt - guess) <= tol
}

#[test]
fn test_load_test_csv_sounding() {
    let (snd, _, _) = load_test_file("standard.csv");
    assert!(validate(&snd).is_ok(), "Failed validation.");
}
