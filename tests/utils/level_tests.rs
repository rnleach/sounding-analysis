use std::collections::HashMap;

use sounding_analysis::error::*;
use sounding_analysis::Levels;
use sounding_base::Sounding;
use super::*;

fn test_levels<F: FnOnce(&Sounding) -> Result<Levels>>(
    snd: &Sounding,
    tgt_int_vals: &HashMap<String, i64>,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
    anal_func: F,
    num_key: &str,
    level_key: &str,
) {
    if let Some(num_levels) = tgt_int_vals.get(num_key) {
        let num_levels = *num_levels as usize;
        let analysis = anal_func(&snd).unwrap();
        println!("\nanalysis = [");
        for lvl in &analysis {
            println!("{:#?}", lvl);
        }
        println!("]");

        let analyzed_num = analysis.len();
        assert_eq!(num_levels, analyzed_num);

        // Check the pressure level of these levels.
        if num_levels > 0 {
            if let Some(level_pressures) = tgt_float_vals.get(level_key) {
                assert!(level_pressures.len() >= 1);
                let mut count_levels_compared = 0;
                for (lvl, it) in analysis.iter().zip(level_pressures) {
                    println!("\nLevel {:#?}  ---  {:#?}", lvl.pressure.unwrap(), it,);
                    assert!(approx_equal(lvl.pressure.unwrap(), *it, 1.0));

                    count_levels_compared += 1;
                }
                assert_eq!(count_levels_compared, num_levels);
            } else {
                panic!("No pressure levels given for analysis target.");
            }
        }
    } else {
        panic!("No num value given..")
    }
}

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
pub fn test_freezing_levels(
    snd: &Sounding,
    tgt_int_vals: &HashMap<String, i64>,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
) {
    use sounding_analysis::levels::freezing_melting_levels;

    test_levels(
        snd,
        tgt_int_vals,
        tgt_float_vals,
        freezing_melting_levels,
        "num freezing level",
        "freezing level pressures",
    );
}

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
pub fn test_wet_bulb_zero_levels(
    snd: &Sounding,
    tgt_int_vals: &HashMap<String, i64>,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
) {
    use sounding_analysis::levels::wet_bulb_zero_levels;

    test_levels(
        snd,
        tgt_int_vals,
        tgt_float_vals,
        wet_bulb_zero_levels,
        "num wet bulb zeros",
        "wet bulb zero pressures",
    );
}

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
pub fn test_max_wet_bulb_in_profile(snd: &Sounding, tgt_float_vals: &HashMap<String, Vec<f64>>) {
    use sounding_analysis::levels::max_wet_bulb_in_profile;

    let analysis = max_wet_bulb_in_profile(snd).unwrap();

    println!("\nanalysis = {:#?}", analysis);

    if let (Some(mwb_pressures), Some(mwb)) = (
        tgt_float_vals.get("max wet bulb pressure"),
        tgt_float_vals.get("max wet bulb aloft"),
    ) {
        assert_eq!(mwb_pressures.len(), 1);
        assert_eq!(mwb.len(), 1);

        for (lvl, it) in [analysis].iter().zip(mwb_pressures) {
            println!("\nLevel {:#?}  ---  {:#?}", lvl.pressure.unwrap(), it,);
            assert!(approx_equal(lvl.pressure.unwrap(), *it, 1.0));
        }

        for (lvl, it) in [analysis].iter().zip(mwb) {
            println!(
                "\nMax Wet Bulb {:#?}  ---  {:#?}",
                lvl.wet_bulb.unwrap(),
                it,
            );
            assert!(approx_equal(lvl.wet_bulb.unwrap(), *it, 1.0));
        }
    } else {
        panic!("Missing max wet bulb aloft level or value.");
    }
}

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
pub fn test_max_temperature(snd: &Sounding,
    tgt_float_vals: &HashMap<String, Vec<f64>>) {
    use sounding_analysis::levels::{max_temperature_in_profile, max_temperature_in_layer};

    let analysis = max_temperature_in_profile(snd).unwrap();

    println!("\nanalysis = {:#?}", analysis);

    if let (Some(mt_pressures), Some(mt)) = (
        tgt_float_vals.get("max temperature pressure"),
        tgt_float_vals.get("max temperature aloft"),
    ) {
        assert_eq!(mt_pressures.len(), 1);
        assert_eq!(mt.len(), 1);

        for (lvl, it) in [analysis].iter().zip(mt_pressures) {
            println!("\nLevel {:#?}  ---  {:#?}", lvl.pressure.unwrap(), it,);
            assert!(approx_equal(lvl.pressure.unwrap(), *it, 1.0));
        }

        for (lvl, it) in [analysis].iter().zip(mt) {
            println!(
                "\nMax Temperature {:#?}  ---  {:#?}",
                lvl.temperature.unwrap(),
                it,
            );
            assert!(approx_equal(lvl.temperature.unwrap(), *it, 1.0));
        }
    } else {
        panic!("Missing max temperature level or value.");
    }

    let warm_layers = ::sounding_analysis::layers::warm_temperature_layer_aloft(snd).unwrap();
    let num_warm_layers = warm_layers.len();
    if num_warm_layers > 0 {
        if let Some(max_ts) = tgt_float_vals.get("warm layer max t") {
            assert_eq!(max_ts.len(), num_warm_layers);

            for (lyr, tgt_temp) in warm_layers.iter().zip(max_ts) {
                let layer_anal = max_temperature_in_layer(snd, lyr).unwrap();
                let max_t_in_this_layer = layer_anal.temperature.unwrap();

                println!("tgt_temp = {} and found value = {} in layer {:#?}", tgt_temp, max_t_in_this_layer, lyr);
                assert!(approx_equal(layer_anal.temperature.unwrap(), *tgt_temp, 0.5));
            }
        } else {   
            panic!("Missing max t value in warm layers aloft.");
        }
    }

}
