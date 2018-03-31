//! Data used in tests.

use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use std::str::FromStr;

use metfor;
use sounding_base::{Profile, Sounding, StationInfo, Surface};
use sounding_validate::validate;

pub fn approx_equal(tgt: f64, guess: f64, tol: f64) -> bool {
    use std::f64;

    assert!(tol > 0.0);

    f64::abs(tgt - guess) <= tol
}

// Load a sounding file, and assert the target values in the file match the values analyzed.
pub fn test_sounding(fname: &str) {
    let (snd, ivals, fvals) = load_test_file(fname);
    assert!(validate(&snd).is_ok(), "Failed validation.");

    test_layers(&snd, &ivals, &fvals);
}

fn load_test_file(fname: &str) -> (Sounding, HashMap<String, i64>, HashMap<String, Vec<f64>>) {
    let mut test_path = PathBuf::new();
    test_path.push("test_data");
    test_path.push(fname);
    load_test_csv_sounding(&test_path)
}

fn load_test_csv_sounding(
    location: &PathBuf,
) -> (Sounding, HashMap<String, i64>, HashMap<String, Vec<f64>>) {
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

    (snd, target_int_vals, target_float_vals)
}

fn test_layers(
    snd: &Sounding,
    tgt_int_vals: &HashMap<String, i64>,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
) {
    test_dendritic_layers(snd, tgt_int_vals, tgt_float_vals);
    test_warm_dry_bulb_aloft_and_cold_sfc_layers(snd, tgt_int_vals, tgt_float_vals);
    test_warm_wet_bulb_aloft(snd, tgt_int_vals, tgt_float_vals);
    test_layer_agl(snd, tgt_float_vals);
}

fn test_dendritic_layers(
    snd: &Sounding,
    tgt_int_vals: &HashMap<String, i64>,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
) {
    use layers::dendritic_snow_zone;

    if let Some(num_dendritic_layers) = tgt_int_vals.get("num dendritic zones") {
        let num_dendritic_layers = *num_dendritic_layers as usize;
        let analysis = dendritic_snow_zone(&snd).unwrap();
        let analyzed_num = analysis.len();
        assert_eq!(num_dendritic_layers, analyzed_num);

        // Check the pressure levels of those growth zones.
        if let Some(dendritic_zone_pressures) = tgt_float_vals.get("dendritic zone pressures") {
            let dendritic_zone_pressures = dendritic_zone_pressures.chunks(2);
            for (lyr, it) in analysis.iter().zip(dendritic_zone_pressures) {
                println!(
                    "\nbottom {:?}  ---  {:?}",
                    lyr.bottom.pressure.unwrap(),
                    it[0]
                );
                assert!(approx_equal(lyr.bottom.pressure.unwrap(), it[0], 0.1));

                println!("top {:?}  ---  {:?}", lyr.top.pressure.unwrap(), it[1]);
                assert!(approx_equal(lyr.top.pressure.unwrap(), it[1], 0.1));
            }
        } else if num_dendritic_layers > 0 {
            panic!("No pressure levels given for dendtritic growth zones.");
        }
    } else {
        panic!("No num dendritic zones analysis in test file.")
    }
}

fn test_warm_dry_bulb_aloft_and_cold_sfc_layers(
    snd: &Sounding,
    tgt_int_vals: &HashMap<String, i64>,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
) {
    use layers::{cold_surface_temperature_layer, warm_temperature_layer_aloft};

    if let Some(num_warm_layers) = tgt_int_vals.get("num warm dry bulb aloft") {
        let num_warm_layers = *num_warm_layers as usize;
        let analysis = warm_temperature_layer_aloft(&snd).unwrap();
        let analyzed_num = analysis.len();
        assert_eq!(num_warm_layers, analyzed_num);

        // Check the pressure levels of those warm layers.
        if let Some(warm_layer_pressures) = tgt_float_vals.get("warm dry bulb layer pressures") {
            let warm_layer_pressures = warm_layer_pressures.chunks(2);
            for (lyr, it) in analysis.iter().zip(warm_layer_pressures) {
                println!(
                    "\nbottom {:?}  ---  {:?}",
                    lyr.bottom.pressure.unwrap(),
                    it[0]
                );
                assert!(approx_equal(lyr.bottom.pressure.unwrap(), it[0], 0.1));

                println!("top {:?}  ---  {:?}", lyr.top.pressure.unwrap(), it[1]);
                assert!(approx_equal(lyr.top.pressure.unwrap(), it[1], 0.1));
            }
        } else if num_warm_layers > 0 {
            panic!("No pressure levels given for warm dry bulb temperature layers.");
        }

        // Check out the cold surface layer
        if let Some(cold_surface_layer_pressures) =
            tgt_float_vals.get("cold surface layer pressures")
        {
            let cold_sfc_analysis = cold_surface_temperature_layer(&snd, &analysis).unwrap();
            let cold_surface_layer_pressures = cold_surface_layer_pressures.chunks(2);
            for (lyr, it) in [cold_sfc_analysis].iter().zip(cold_surface_layer_pressures) {
                println!(
                    "\nbottom {:?}  ---  {:?}",
                    lyr.bottom.pressure.unwrap(),
                    it[0]
                );
                assert!(approx_equal(lyr.bottom.pressure.unwrap(), it[0], 0.1));

                println!("top {:?}  ---  {:?}", lyr.top.pressure.unwrap(), it[1]);
                assert!(approx_equal(lyr.top.pressure.unwrap(), it[1], 0.1));
            }
        } else if num_warm_layers > 0 {
            panic!("No pressure levels given for cold surface layer.");
        }
    } else {
        panic!("No num warm dry bulb layer info in test file.")
    }
}

fn test_warm_wet_bulb_aloft(
    snd: &Sounding,
    tgt_int_vals: &HashMap<String, i64>,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
) {
    use layers::warm_wet_bulb_layer_aloft;

    if let Some(num_warm_layers) = tgt_int_vals.get("num warm wet bulb aloft") {
        let num_warm_layers = *num_warm_layers as usize;
        let analysis = warm_wet_bulb_layer_aloft(&snd).unwrap();
        let analyzed_num = analysis.len();
        assert_eq!(num_warm_layers, analyzed_num);

        // Check the pressure levels of these layers.
        if let Some(warm_layer_pressures) = tgt_float_vals.get("warm wet bulb layer pressures") {
            let warm_layer_pressures = warm_layer_pressures.chunks(2);
            for (lyr, it) in analysis.iter().zip(warm_layer_pressures) {
                println!(
                    "\nbottom {:?}  ---  {:?}",
                    lyr.bottom.pressure.unwrap(),
                    it[0]
                );
                assert!(approx_equal(lyr.bottom.pressure.unwrap(), it[0], 0.1));

                println!("top {:?}  ---  {:?}", lyr.top.pressure.unwrap(), it[1]);
                assert!(approx_equal(lyr.top.pressure.unwrap(), it[1], 0.1));
            }
        } else if num_warm_layers > 0 {
            panic!("No pressure levels given for warm dry bulb temperature layers.");
        }
    } else {
        panic!("No num warm wet bulb layer info in test file.")
    }
}

fn test_layer_agl(
    snd: &Sounding,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
){
    use ::layers::layer_agl;

    let analysis = layer_agl(snd, 6000.0).unwrap();
    println!("6km AGL layer: {:?}", analysis);
    if let Some(agl_layer_pressures) = tgt_float_vals.get("6km agl layer pressures") {
            let agl_layer_pressures = agl_layer_pressures.chunks(2);
            for (lyr, it) in [analysis].iter().zip(agl_layer_pressures) {
                println!(
                    "\nbottom {:?}  ---  {:?}",
                    lyr.bottom.pressure.unwrap(),
                    it[0]
                );
                assert!(approx_equal(lyr.bottom.pressure.unwrap(), it[0], 0.1));

                println!("top {:?}  ---  {:?}", lyr.top.pressure.unwrap(), it[1]);
                assert!(approx_equal(lyr.top.pressure.unwrap(), it[1], 0.1));
            }
        } else {
            panic!("No pressure levels given for agl layers.");
        }
}

#[test]
fn test_load_test_csv_sounding() {
    let (snd, ivals, fvals) = load_test_file("standard.csv");
    assert!(validate(&snd).is_ok(), "Failed validation.");

    assert!(ivals.contains_key("num dendritic zones"));
    assert!(Some(&1) == ivals.get("num dendritic zones"));
    assert!(ivals.contains_key("num warm dry bulb aloft"));
    assert!(Some(&0) == ivals.get("num warm dry bulb aloft"));
    assert!(ivals.contains_key("num warm wet bulb aloft"));
    assert!(Some(&0) == ivals.get("num warm wet bulb aloft"));
    assert!(ivals.contains_key("num inversions"));
    assert!(Some(&0) == ivals.get("num inversions"));

    assert!(fvals.contains_key("dendritic zone pressures"));
    assert!(Some(&vec![604.2, 534.7]) == fvals.get("dendritic zone pressures"));
    assert!(fvals.contains_key("6km agl layer pressures"));
    assert!(Some(&vec![1080.0, 508.5]) == fvals.get("6km agl layer pressures"));

}

#[test]
fn test_standard() {
    test_sounding("standard.csv");
}
