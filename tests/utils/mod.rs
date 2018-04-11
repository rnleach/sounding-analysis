use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use std::str::FromStr;

use metfor;
use sounding_base::{Profile, Sounding, StationInfo, Surface};

pub mod layer_tests;
pub mod level_tests;

#[allow(unused_macros)] // False alarm
macro_rules! check_file_complete {
    ($test_name:ident, $fname:expr) => {
        #[test]
        fn $test_name() {
            let (snd, ivals, fvals) = utils::load_test_file($fname);

            assert!(sounding_validate::validate(&snd).is_ok(), "Failed validation.");

            let ival_keys = [
                "num dendritic zones",
                "num warm dry bulb aloft",
                "num warm wet bulb aloft",
                "num inversions",
                "num freezing level",
                "num wet bulb zeros",
            ];

            let fval_keys = [
                "dendritic zone pressures",
                "warm dry bulb layer pressures",
                "cold surface layer pressures",
                "warm wet bulb layer pressures",
                "6km agl layer pressures",
                "700-500 hPa layer heights",
                "inversion layer pressures",
                "freezing level pressures",
                "wet bulb zero pressures",
                "max wet bulb aloft",
                "max wet bulb pressure",
            ];

            // Make sure all of these keys are in the hashmaps
            for key in ival_keys.iter() {
                assert!(ivals.contains_key(*key), *key);
            }

            for key in fval_keys.iter() {
                assert!(fvals.contains_key(*key), *key);
            }

            // Make sure there are no extra keys in there being ignored.
            for key in ivals.keys() {
                assert!(ival_keys.contains(&key.as_str()), "extra ival key found");
            }

            // Make sure there are no extra keys in there being ignored.
            for key in fvals.keys() {
                assert!(fval_keys.contains(&key.as_str()), "extra fval key found");
            }
        }
    };
}

#[allow(unused_macros)] // False alarm
macro_rules! test_file {
    ($test_mod_name:ident, $fname:expr) => {

        mod $test_mod_name {

            use std::collections::HashMap;

            use ::utils;
            use ::sounding_base::{Sounding};

            fn load_data() -> (Sounding, HashMap<String, i64>, HashMap<String, Vec<f64>>) {
                utils::load_test_file($fname)
            }

            mod levels {
                use ::$test_mod_name::load_data;
                use ::utils::level_tests;

                #[test]
                fn freezing_level() {
                    let (snd, ivals, fvals) = load_data();
                    level_tests::test_freezing_levels(&snd, &ivals, &fvals);
                }

                #[test]
                fn wet_bulb_zero(){
                    let (snd, ivals, fvals) = load_data();
                    level_tests::test_wet_bulb_zero_levels(&snd, &ivals, &fvals);
                }

                #[test]
                fn max_wet_bulb_in_profile(){
                    let (snd, _, fvals) = load_data();
                    level_tests::test_max_wet_bulb_in_profile(&snd, &fvals);
                }
            }

            mod layers {
                use ::$test_mod_name::load_data;
                use ::utils::layer_tests;

                #[test]
                fn dendritic_layers() {
                    let (snd, ivals, fvals) = load_data();
                    layer_tests::test_dendritic_layers(&snd, &ivals, &fvals);
                }

                #[test]
                fn warm_dry_bulb_aloft_and_cold_sfc_layers(){
                    let (snd, ivals, fvals) = load_data();
                    layer_tests::test_warm_dry_bulb_aloft_and_cold_sfc_layers(
                        &snd,
                        &ivals,
                        &fvals,
                    );
                }

                #[test]
                fn warm_wet_bulb_aloft(){
                    let (snd, ivals, fvals) = load_data();
                    layer_tests::test_warm_wet_bulb_aloft(&snd, &ivals, &fvals);
                }

                #[test]
                fn layer_agl(){
                    let (snd, _, fvals) = load_data();
                    layer_tests::test_layer_agl(&snd, &fvals);
                }

                #[test]
                fn pressure_layer() {
                    let (snd, _, fvals) = load_data();
                    layer_tests::test_pressure_layer(&snd, &fvals);
                }

                #[test]
                fn inversions(){
                    let (snd, ivals, fvals) = load_data();
                    layer_tests::test_inversion_layers(&snd, &ivals, &fvals);
                }
            }

        }
    };
}

pub fn load_test_file(fname: &str) -> (Sounding, HashMap<String, i64>, HashMap<String, Vec<f64>>) {
    let mut test_path = PathBuf::new();
    test_path.push("test_data");
    test_path.push(fname);
    load_test_csv_sounding(&test_path)
}

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
fn approx_equal(tgt: f64, guess: f64, tol: f64) -> bool {
    use std::f64;

    assert!(tol > 0.0);

    f64::abs(tgt - guess) <= tol
}

fn load_test_csv_sounding(
    location: &PathBuf,
) -> (Sounding, HashMap<String, i64>, HashMap<String, Vec<f64>>) {
    let mut f = File::open(location).expect(&format!("Error opening file: {:#?}", location));

    let mut contents = String::new();
    f.read_to_string(&mut contents)
        .expect(&format!("Error reading file: {:#?}", location));

    let lines: Vec<&str> = contents.split('\n').collect();
    let mut line_iter = lines.iter();

    //
    // Parse profile data
    //
    let mut height: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut temp: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut wb: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut dp: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut press: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut wspd: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut wdir: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut wet_bulb: Vec<Option<f64>> = Vec::with_capacity(lines.len());

    for line in line_iter.by_ref() {
        if line.starts_with("### Surface Data ###")
            || line.starts_with("### Analysis Int Section ###")
            || line.starts_with("### Analysis Float Section ###")
        {
            break;
        }

        let tokens: Vec<&str> = line.split(',').collect();
        if tokens.len() < 6 {
            continue;
        }
        let t_c = f64::from_str(tokens[1]).ok();
        let dp_c = f64::from_str(tokens[2]).ok();
        let press_hpa = f64::from_str(tokens[3]).ok();
        let wb_c = t_c.and_then(|t| {
            dp_c.and_then(|dp| press_hpa.and_then(|p| metfor::wet_bulb_c(t, dp, p).ok()))
        });

        height.push(f64::from_str(tokens[0]).ok());
        temp.push(t_c);
        wb.push(wb_c);
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
        .set_profile(Profile::WetBulb, wb)
        .set_profile(Profile::DewPoint, dp)
        .set_profile(Profile::Pressure, press)
        .set_profile(Profile::WetBulb, wet_bulb)
        .set_profile(Profile::WindSpeed, wspd)
        .set_profile(Profile::WindDirection, wdir);

    //
    // Surface data
    //
    for line in line_iter.by_ref() {
        if line.starts_with("### Analysis Int Section ###")
            || line.starts_with("### Analysis Float Section ###")
        {
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
        let wspd = f64::from_str(tokens[4]).ok();
        let wdir = f64::from_str(tokens[5]).ok();

        snd = snd.set_surface_value(Surface::Temperature, t_c)
            .set_station_info(StationInfo::new().with_elevation(height))
            .set_surface_value(Surface::DewPoint, dp_c)
            .set_surface_value(Surface::StationPressure, press_hpa)
            .set_surface_value(Surface::WindSpeed, wspd)
            .set_surface_value(Surface::WindDirection, wdir)
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
        if tokens.len() < 1 {
            continue;
        }
        let key = tokens[0].to_owned();
        let mut values = vec![];
        for token in tokens.iter().skip(1) {
            values.push(f64::from_str(token).unwrap());
        }
        target_float_vals.insert(key, values);
    }

    (snd, target_int_vals, target_float_vals)
}
