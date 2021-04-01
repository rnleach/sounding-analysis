//! Run these benches with `cargo bench --bench layers -- --verbose`

use criterion::{criterion_group, criterion_main, Criterion};

mod utils;

criterion_main!(
    temperature_layers_benches,
    height_pressure_layers_benches,
    inversion_layers_benches,
    convective_layers_benches
);

/**************************************************************************************************
 *                                     Temperature Layers
 **************************************************************************************************/
criterion_group!(
    temperature_layers_benches,
    hail_growth_zone_bench,
    warm_temperature_layer_aloft_bench,
    cold_surface_temperature_layer_bench,
    warm_surface_temperature_layer_bench,
    melting_freeing_enery_area_bench
);

// No bench for dendritic snow zone since it uses the same inner function as this one does.
fn hail_growth_zone_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("hail_growth_zone", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::hail_growth_zone(&snd).expect("oops");
            }
        });
    });
}

// No bench for warm_wet_bulb_layer_aloft since it uses the same inner function as this one does.
fn warm_temperature_layer_aloft_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("warm_temperature_layer_aloft", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::warm_temperature_layer_aloft(&snd).expect("oops");
            }
        });
    });
}

fn cold_surface_temperature_layer_bench(c: &mut Criterion) {
    use sounding_analysis::{Layer, Sounding};

    let snds = utils::load_all_test_files().to_vec();
    let pairs: Vec<(Sounding, Vec<Layer>)> = snds
        .into_iter()
        .map(|snd| {
            let warm_layers = sounding_analysis::warm_temperature_layer_aloft(&snd).unwrap();
            (snd, warm_layers)
        })
        .collect();

    c.bench_function("cold_surface_temperature_layer", |b| {
        b.iter(|| {
            for (snd, warm_layers) in &pairs {
                let _x = sounding_analysis::cold_surface_temperature_layer(&snd, &warm_layers)
                    .expect("oops");
            }
        });
    });
}

fn warm_surface_temperature_layer_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("warm_surface_temperature_layer", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::warm_surface_temperature_layer(&snd).expect("oops");
            }
        });
    });
}

fn melting_freeing_enery_area_bench(c: &mut Criterion) {
    use sounding_analysis::{Layer, Sounding};

    let snds: Vec<Sounding> = utils::load_all_test_files().to_vec();

    let pairs: Vec<(Sounding, Vec<Layer>)> = snds
        .into_iter()
        .filter_map(|snd| {
            let warm_layers = sounding_analysis::warm_temperature_layer_aloft(&snd).unwrap();
            let warm_surface_layer =
                sounding_analysis::warm_surface_temperature_layer(&snd).unwrap();
            let cold_surface_layer =
                sounding_analysis::cold_surface_temperature_layer(&snd, &warm_layers).unwrap();

            let layers: Vec<sounding_analysis::Layer> = cold_surface_layer
                .into_iter()
                .chain(warm_surface_layer.into_iter())
                .chain(warm_layers.into_iter())
                .filter(|lyr| lyr.bottom.height.is_some() && lyr.top.height.is_some())
                .collect();

            if layers.is_empty() {
                None
            } else {
                Some((snd, layers))
            }
        })
        .collect();

    c.bench_function("melting_freezing_energy_area", |b| {
        b.iter(|| {
            for (snd, layers) in &pairs {
                for lyr in layers {
                    let _x =
                        sounding_analysis::melting_freezing_energy_area(&snd, &lyr).expect("oops");
                }
            }
        });
    });
}

/**************************************************************************************************
 *                                     Height/Pressure Layers
 **************************************************************************************************/
criterion_group!(
    height_pressure_layers_benches,
    layer_hgt_bench,
    pressure_layer_bench
);

fn layer_hgt_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("layer_agl", |b| {
        b.iter(|| {
            for snd in &snds {
                for hgt in &[
                    metfor::Meters(500.0),
                    metfor::Meters(3_000.0),
                    metfor::Meters(6_000.0),
                ] {
                    let _x = sounding_analysis::layer_agl(&snd, *hgt).expect("oops");
                }
            }
        });
    });
}

fn pressure_layer_bench(c: &mut Criterion) {
    use metfor::HectoPascal;
    use optional::some;
    use sounding_analysis::{DataRow, Layer};

    let snds = utils::load_all_test_files();

    let layers: Vec<_> = [
        (1000.0, 500.0),
        (850.0, 500.0),
        (700.0, 500.0),
        (1000.0, 700.0),
        (850.0, 700.0),
        (1000.0, 850.0),
    ]
    .iter()
    .map(|&(bottom, top)| (HectoPascal(bottom), HectoPascal(top)))
    .map(|(bottom, top)| (some(bottom), some(top)))
    .map(|(bottom, top)| {
        (
            DataRow {
                pressure: bottom,
                ..DataRow::default()
            },
            DataRow {
                pressure: top,
                ..DataRow::default()
            },
        )
    })
    .map(|(bottom, top)| Layer { bottom, top })
    .collect();

    c.bench_function("pressure_layer", |b| {
        b.iter(|| {
            for snd in &snds {
                for layer in &layers {
                    let _x = sounding_analysis::pressure_layer(
                        &snd,
                        layer.bottom.pressure.unwrap(),
                        layer.top.pressure.unwrap(),
                    )
                    .expect("oops");
                }
            }
        });
    });
}

/**************************************************************************************************
 *                                     Inversion Layers
 **************************************************************************************************/
criterion_group!(
    inversion_layers_benches,
    inversions_bench,
    sfc_based_inversion_bench
);

fn inversions_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("inversions", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x =
                    sounding_analysis::inversions(&snd, metfor::HectoPascal(400.0)).expect("oops");
            }
        });
    });
}

fn sfc_based_inversion_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("sfc_based_inversion", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::sfc_based_inversion(&snd).expect("oops");
            }
        });
    });
}

/**************************************************************************************************
 *                                     Convective Layers
 **************************************************************************************************/
criterion_group!(convective_layers_benches, effective_inflow_layer_bench);

fn effective_inflow_layer_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("effective_inflow_layer", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::effective_inflow_layer(&snd);
            }
        });
    });
}
