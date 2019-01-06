use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use std::str::FromStr;

use metfor::{self, Celsius, HectoPascal, WindSpdDir, Knots, Meters};
use optional::{Optioned};
use sounding_base::{Sounding, StationInfo};

pub mod index_tests;
pub mod layer_tests;
pub mod level_tests;

#[allow(unused_macros)] // False alarm
macro_rules! check_file_complete {
    ($test_name:ident, $fname:expr) => {
        #[test]
        fn $test_name() {
            let (snd, ivals, fvals) = utils::load_test_file($fname);

            assert!(
                sounding_validate::validate(&snd).is_ok(),
                "Failed validation."
            );

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
                "max temperature aloft",
                "max temperature pressure",
                "warm layer max t",
                "haines_low",
                "haines_mid",
                "haines_high",
                "haines",
                "hdw",
                "kindex",
                "swet",
                "total_totals",
                "precipitable_water",
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

            use crate::utils;
            use sounding_base::Sounding;

            fn load_data() -> (Sounding, HashMap<String, i64>, HashMap<String, Vec<f64>>) {
                utils::load_test_file($fname)
            }

            mod levels {
                use crate::utils::level_tests;
                use crate::$test_mod_name::load_data;

                #[test]
                fn freezing_level() {
                    let (snd, ivals, fvals) = load_data();
                    level_tests::test_freezing_levels(&snd, &ivals, &fvals);
                }

                #[test]
                fn wet_bulb_zero() {
                    let (snd, ivals, fvals) = load_data();
                    level_tests::test_wet_bulb_zero_levels(&snd, &ivals, &fvals);
                }

                #[test]
                fn max_wet_bulb_in_profile() {
                    let (snd, _, fvals) = load_data();
                    level_tests::test_max_wet_bulb_in_profile(&snd, &fvals);
                }

                #[test]
                fn max_temperature_in_profile() {
                    let (snd, _, fvals) = load_data();
                    level_tests::test_max_temperature(&snd, &fvals);
                }
            }

            mod layers {
                use crate::utils::layer_tests;
                use crate::$test_mod_name::load_data;

                #[test]
                fn dendritic_layers() {
                    let (snd, ivals, fvals) = load_data();
                    layer_tests::test_dendritic_layers(&snd, &ivals, &fvals);
                }

                #[test]
                fn warm_dry_bulb_aloft_and_cold_sfc_layers() {
                    let (snd, ivals, fvals) = load_data();
                    layer_tests::test_warm_dry_bulb_aloft_and_cold_sfc_layers(&snd, &ivals, &fvals);
                }

                #[test]
                fn warm_wet_bulb_aloft() {
                    let (snd, ivals, fvals) = load_data();
                    layer_tests::test_warm_wet_bulb_aloft(&snd, &ivals, &fvals);
                }

                #[test]
                fn layer_agl() {
                    let (snd, _, fvals) = load_data();
                    layer_tests::test_layer_agl(&snd, &fvals);
                }

                #[test]
                fn pressure_layer() {
                    let (snd, _, fvals) = load_data();
                    layer_tests::test_pressure_layer(&snd, &fvals);
                }

                #[test]
                fn inversions() {
                    let (snd, ivals, fvals) = load_data();
                    layer_tests::test_inversion_layers(&snd, &ivals, &fvals);
                }
            }

            mod indexes {
                use crate::utils::index_tests;
                use crate::$test_mod_name::load_data;
                use sounding_analysis;

                #[test]
                fn test_haines_low() {
                    let (snd, _, fvals) = load_data();
                    index_tests::test_index(
                        &snd,
                        &fvals,
                        sounding_analysis::haines_low,
                        "haines_low",
                        0,
                        0,
                    );
                }

                #[test]
                fn test_haines_mid() {
                    let (snd, _, fvals) = load_data();
                    index_tests::test_index(
                        &snd,
                        &fvals,
                        sounding_analysis::haines_mid,
                        "haines_mid",
                        0,
                        0,
                    );
                }

                #[test]
                fn test_haines_high() {
                    let (snd, _, fvals) = load_data();
                    index_tests::test_index(
                        &snd,
                        &fvals,
                        sounding_analysis::haines_high,
                        "haines_high",
                        0,
                        0,
                    );
                }

                #[test]
                fn test_haines() {
                    let (snd, _, fvals) = load_data();
                    index_tests::test_index(
                        &snd,
                        &fvals,
                        sounding_analysis::haines,
                        "haines",
                        0,
                        0,
                    );
                }

                #[test]
                fn test_hdw() {
                    let (snd, _, fvals) = load_data();
                    index_tests::test_index(
                        &snd,
                        &fvals,
                        sounding_analysis::hot_dry_windy,
                        "hdw",
                        5.0,
                        0.0,
                    );
                }

                #[test]
                fn test_kindex() {
                    let (snd, _, fvals) = load_data();
                    index_tests::test_index(
                        &snd,
                        &fvals,
                        sounding_analysis::kindex,
                        "kindex",
                        metfor::Celsius(0.1),
                        metfor::Celsius(-1000.0),
                    );
                }

                #[test]
                fn test_swet() {
                    let (snd, _, fvals) = load_data();
                    index_tests::test_index(
                        &snd,
                        &fvals,
                        sounding_analysis::swet,
                        "swet",
                        1.0,
                        -1000.0,
                    );
                }

                #[test]
                fn test_total_totals() {
                    let (snd, _, fvals) = load_data();
                    index_tests::test_index(
                        &snd,
                        &fvals,
                        sounding_analysis::total_totals,
                        "total_totals",
                        1.0,
                        0.0,
                    );
                }

                #[test]
                fn test_precipitable_water() {
                    let (snd, _, fvals) = load_data();
                    index_tests::test_index(
                        &snd,
                        &fvals,
                        sounding_analysis::precipitable_water,
                        "precipitable_water",
                        metfor::Mm(0.5),
                        metfor::Mm(0.0),
                    );
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

    let passed = f64::abs(tgt - guess) <= tol;

    if !passed {
        println!("left = {} and right = {}", tgt, guess);
    }

    passed
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
    let mut height: Vec<Optioned<Meters>> = Vec::with_capacity(lines.len());
    let mut temp: Vec<Optioned<Celsius>> = Vec::with_capacity(lines.len());
    let mut wb: Vec<Optioned<Celsius>> = Vec::with_capacity(lines.len());
    let mut dp: Vec<Optioned<Celsius>> = Vec::with_capacity(lines.len());
    let mut press: Vec<Optioned<HectoPascal>> = Vec::with_capacity(lines.len());
    let mut wind: Vec<Optioned<WindSpdDir<Knots>>> = Vec::with_capacity(lines.len());

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
        let t_c: Option<Celsius> = f64::from_str(tokens[1]).ok().map(Celsius);
        let dp_c: Option<Celsius> = f64::from_str(tokens[2]).ok().map(Celsius);
        let press_hpa: Option<HectoPascal> = f64::from_str(tokens[3]).ok().map(HectoPascal).into();
        let wb_c: Option<Celsius> = t_c.and_then(|t| {
            dp_c.and_then(|dp| press_hpa.and_then(|p| metfor::wet_bulb(t, dp, p)))
        });
        let wspd = f64::from_str(tokens[4]).ok();
        let wdir = f64::from_str(tokens[5]).ok();
        let wind_val: Option<WindSpdDir<Knots>> = wspd.and_then(|wspd| wdir.map(|wdir| WindSpdDir{speed: Knots(wspd), direction: wdir}));

        height.push(f64::from_str(tokens[0]).ok().map(Meters).into());
        temp.push(t_c.into());
        wb.push(wb_c.into());
        dp.push(dp_c.into());
        press.push(press_hpa.into());
        wind.push(wind_val.into());
        
    }

    let mut snd = Sounding::new()
        .with_height_profile(height)
        .with_temperature_profile(temp)
        .with_wet_bulb_profile(wb)
        .with_dew_point_profile(dp)
        .with_pressure_profile(press)
        .with_wind_profile(wind);

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

        let height = f64::from_str(tokens[0]).ok().map(Meters);
        let t_c = f64::from_str(tokens[1]).ok().map(Celsius);
        let dp_c = f64::from_str(tokens[2]).ok().map(Celsius);
        let press_hpa = f64::from_str(tokens[3]).ok().map(HectoPascal);
        let wspd = f64::from_str(tokens[4]).ok();
        let wdir = f64::from_str(tokens[5]).ok();
        let wind = wspd
            .and_then(|wspd| wdir
                .map(|wdir| WindSpdDir{speed: Knots(wspd), direction: wdir}));

        snd = snd
            .with_sfc_temperature(t_c)
            .with_station_info(StationInfo::new().with_elevation(height))
            .with_sfc_dew_point(dp_c)
            .with_station_pressure(press_hpa)
            .with_sfc_wind(wind)
    }

    //
    // Integer values.
    //
    let mut target_int_vals = HashMap::new();
    for line in line_iter.by_ref() {
        if line.starts_with("### Analysis Float Section ###") {
            break;
        }

        let tokens: Vec<String> = line
            .split(',')
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
        let tokens: Vec<String> = line
            .split(',')
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
