//! Run these benches with `cargo bench --bench parcel -- --verbose`

use criterion::{criterion_group, criterion_main, Criterion};

mod utils;

fn build_tester() -> Criterion {
    Criterion::default()
        .sample_size(200)
        .measurement_time(std::time::Duration::from_secs(10))
        .noise_threshold(0.03)
        .significance_level(0.01)
}

criterion_main!(parcel_benches);

criterion_group!(
    name = parcel_benches;
    config = build_tester();
    targets = mixed_layer_parcel_bench, convective_parcel_bench, most_unstable_parcel_bench,
              average_parcel_bench
);

fn mixed_layer_parcel_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("mixed_layer_parcel", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::mixed_layer_parcel(&snd).expect("oops");
            }
        });
    });
}

// No bench for surface_parcel, lowest_level_parcel, pressure_parcel, or effective_layer_parcel
// because those are so simple.  If there's room for improvement in those functions it's probably
// buried in the code they call from some other module.

fn most_unstable_parcel_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("most_unstable_parcel", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::most_unstable_parcel(&snd).expect("oops");
            }
        });
    });
}

fn convective_parcel_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("convective_parcel", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::convective_parcel(&snd).expect("oops");
            }
        });
    });
}

fn average_parcel_bench(c: &mut Criterion) {
    use metfor::HectoPascal;
    use sounding_analysis::{Layer, Sounding};

    let snds = utils::load_all_test_files().to_vec();
    let pairs: Vec<(Sounding, Layer)> = snds
        .into_iter()
        .map(|snd| {
            let lyr =
                sounding_analysis::pressure_layer(&snd, HectoPascal(850.0), HectoPascal(700.0))
                    .unwrap();
            (snd, lyr)
        })
        .collect();

    c.bench_function("average_parcel", |b| {
        b.iter(|| {
            for (snd, lyr) in &pairs {
                let _x = sounding_analysis::average_parcel(&snd, lyr).expect("oops");
            }
        });
    });
}
