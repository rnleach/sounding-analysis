//! Run these benches with `cargo bench --bench wind -- --verbose`
use criterion::{criterion_group, criterion_main, Criterion};

mod utils;

fn build_tester() -> Criterion {
    Criterion::default()
        .sample_size(200)
        .measurement_time(std::time::Duration::from_secs(10))
        .noise_threshold(0.03)
        .significance_level(0.01)
}

criterion_main!(wind_benches);

criterion_group!(
    name = wind_benches;
    config = build_tester();
    targets = mean_wind_bench, sr_helicity_bench, bunkers_storm_motion_bench
);

fn mean_wind_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files().to_vec();
    let pairs = snds
        .into_iter()
        .filter_map(|snd| {
            sounding_analysis::pressure_layer(
                &snd,
                metfor::HectoPascal(850.0),
                metfor::HectoPascal(700.0),
            )
            .ok()
            .map(|lyr| (snd, lyr))
        })
        .filter(|(_, lyr)| lyr.bottom.height.is_some())
        .collect::<Vec<_>>();

    assert!(!pairs.is_empty());

    c.bench_function("mean_wind", |b| {
        b.iter(|| {
            for (snd, layer) in &pairs {
                let _x = sounding_analysis::mean_wind(layer, &snd).expect("oops");
            }
        });
    });
}

fn sr_helicity_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files().to_vec();
    let pairs = snds
        .into_iter()
        .filter_map(|snd| {
            sounding_analysis::pressure_layer(
                &snd,
                metfor::HectoPascal(850.0),
                metfor::HectoPascal(700.0),
            )
            .ok()
            .map(|lyr| (snd, lyr))
        })
        .filter(|(_, lyr)| lyr.bottom.height.is_some())
        .collect::<Vec<_>>();

    let storm_motion: metfor::WindUV<metfor::MetersPSec> = metfor::WindUV {
        u: metfor::MetersPSec(6.0),
        v: metfor::MetersPSec(6.0),
    };

    assert!(!pairs.is_empty());

    c.bench_function("sr_helicity", |b| {
        b.iter(|| {
            for (snd, layer) in &pairs {
                let _x = sounding_analysis::sr_helicity(layer, storm_motion, &snd).expect("oops");
            }
        });
    });
}

fn bunkers_storm_motion_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("bunkers_storm_motion", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::bunkers_storm_motion(&snd);
            }
        });
    });
}
