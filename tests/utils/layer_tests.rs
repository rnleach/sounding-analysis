use std::collections::HashMap;

use sounding_base::Sounding;
use sounding_analysis::{Layer, error::*, VEC_SIZE};
use super::*;

extern crate smallvec;
use self::smallvec::SmallVec;

fn test_layers<F: FnOnce(&Sounding) -> Result<SmallVec<[Layer; VEC_SIZE]>>>(
    snd: &Sounding,
    tgt_int_vals: &HashMap<String, i64>,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
    anal_func: F,
    num_key: &str,
    levels_key: &str,
) {
    if let Some(num_layers) = tgt_int_vals.get(num_key) {
        let num_layers = *num_layers as usize;
        let analysis = anal_func(&snd).unwrap();
        println!("\nanalysis = {:#?}", analysis);

        let analyzed_num = analysis.len();
        assert_eq!(num_layers, analyzed_num);

        // Check the pressure levels of these layers.
        if num_layers > 0 {
            if let Some(layer_pressures) = tgt_float_vals.get(levels_key) {
                assert!(layer_pressures.len() >= 2);
                let layer_pressures = layer_pressures.chunks(2);
                let mut count_layers_compared = 0;
                for (lyr, it) in analysis.iter().zip(layer_pressures) {
                    println!(
                        "\nbottom {:#?}  ---  {:#?}",
                        lyr.bottom.pressure.unwrap(),
                        it[0]
                    );
                    assert!(approx_equal(lyr.bottom.pressure.unwrap(), it[0], 0.1));

                    println!("top {:#?}  ---  {:#?}", lyr.top.pressure.unwrap(), it[1]);
                    assert!(approx_equal(lyr.top.pressure.unwrap(), it[1], 0.1));
                    count_layers_compared += 1;
                }
                assert_eq!(count_layers_compared, num_layers);
            } else {
                panic!("No pressure levels given for analysis target.");
            }
        }
    } else {
        panic!("No num value given..")
    }
}

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
pub fn test_dendritic_layers(
    snd: &Sounding,
    tgt_int_vals: &HashMap<String, i64>,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
) {
    use sounding_analysis::layers::dendritic_snow_zone;

    test_layers(
        snd,
        tgt_int_vals,
        tgt_float_vals,
        dendritic_snow_zone,
        "num dendritic zones",
        "dendritic zone pressures",
    );
}

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
pub fn test_warm_dry_bulb_aloft_and_cold_sfc_layers(
    snd: &Sounding,
    tgt_int_vals: &HashMap<String, i64>,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
) {
    use sounding_analysis::layers::{cold_surface_temperature_layer, warm_temperature_layer_aloft};

    // Test warm layers aloft.
    test_layers(
        snd,
        tgt_int_vals,
        tgt_float_vals,
        warm_temperature_layer_aloft,
        "num warm dry bulb aloft",
        "warm dry bulb layer pressures",
    );

    // Test the cold surface layer.
    if let Some(num_warm_layers) = tgt_int_vals.get("num warm dry bulb aloft") {
        let num_warm_layers = *num_warm_layers as usize;
        let analysis = warm_temperature_layer_aloft(&snd).unwrap();

        // Check the pressure levels of those warm layers.
        if num_warm_layers > 0 {
            // Check out the cold surface layer
            if let Some(cold_surface_layer_pressures) =
                tgt_float_vals.get("cold surface layer pressures")
            {
                let cold_sfc_analysis = cold_surface_temperature_layer(&snd, &analysis).unwrap();
                let cold_surface_layer_pressures = cold_surface_layer_pressures.chunks(2);
                for (lyr, it) in [cold_sfc_analysis].iter().zip(cold_surface_layer_pressures) {
                    println!(
                        "\nbottom {:#?}  ---  {:#?}",
                        lyr.bottom.pressure.unwrap(),
                        it[0]
                    );
                    assert!(approx_equal(lyr.bottom.pressure.unwrap(), it[0], 0.1));

                    println!("top {:#?}  ---  {:#?}", lyr.top.pressure.unwrap(), it[1]);
                    assert!(approx_equal(lyr.top.pressure.unwrap(), it[1], 0.1));
                }
            } else {
                panic!("No pressure levels given for cold surface layer.");
            }
        }
    }
}

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
pub fn test_warm_wet_bulb_aloft(
    snd: &Sounding,
    tgt_int_vals: &HashMap<String, i64>,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
) {
    use sounding_analysis::layers::warm_wet_bulb_layer_aloft;

    test_layers(
        snd,
        tgt_int_vals,
        tgt_float_vals,
        warm_wet_bulb_layer_aloft,
        "num warm wet bulb aloft",
        "warm wet bulb layer pressures",
    );
}

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
pub fn test_layer_agl(snd: &Sounding, tgt_float_vals: &HashMap<String, Vec<f64>>) {
    use sounding_analysis::layers::layer_agl;

    let analysis = layer_agl(snd, 6000.0).unwrap();
    println!("\n6km AGL layer: {:#?}", analysis);
    if let Some(agl_layer_pressures) = tgt_float_vals.get("6km agl layer pressures") {
        assert_eq!(agl_layer_pressures.len(), 2);
        let agl_layer_pressures = agl_layer_pressures.chunks(2);
        for (lyr, it) in [analysis].iter().zip(agl_layer_pressures) {
            println!(
                "bottom {:#?}  ---  {:#?}",
                lyr.bottom.pressure.unwrap(),
                it[0]
            );
            assert!(approx_equal(lyr.bottom.pressure.unwrap(), it[0], 1.0));

            println!("top {:#?}  ---  {:#?}", lyr.top.pressure.unwrap(), it[1]);
            assert!(approx_equal(lyr.top.pressure.unwrap(), it[1], 1.0));
        }
    } else {
        panic!("No pressure levels given for agl layers.");
    }
}

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
pub fn test_pressure_layer(snd: &Sounding, tgt_float_vals: &HashMap<String, Vec<f64>>) {
    use sounding_analysis::layers::pressure_layer;

    let analysis = pressure_layer(snd, 700.0, 500.0).unwrap();
    println!("\n700-500 hPa layer: {:#?}", analysis);
    if let Some(pressure_layer_heights) = tgt_float_vals.get("700-500 hPa layer heights") {
        assert_eq!(pressure_layer_heights.len(), 2);
        let pressure_layer_heights = pressure_layer_heights.chunks(2);
        for (lyr, it) in [analysis].iter().zip(pressure_layer_heights) {
            println!(
                "bottom {:#?}  ---  {:#?}",
                lyr.bottom.height.unwrap(),
                it[0]
            );
            assert!(approx_equal(lyr.bottom.height.unwrap(), it[0], 1.0));

            println!("top {:#?}  ---  {:#?}", lyr.top.height.unwrap(), it[1]);
            assert!(approx_equal(lyr.top.height.unwrap(), it[1], 1.0));
        }
    } else {
        panic!("No heights given for pressure layer.");
    }
}

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
pub fn test_inversion_layers(
    snd: &Sounding,
    tgt_int_vals: &HashMap<String, i64>,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
) {
    use sounding_analysis::layers::inversions;

    test_layers(
        snd,
        tgt_int_vals,
        tgt_float_vals,
        inversions,
        "num inversions",
        "inversion layer pressures",
    );
}