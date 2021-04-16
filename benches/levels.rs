//! Run these benches with `cargo bench --bench levels -- --verbose`

use criterion::{criterion_group, criterion_main, Criterion};

mod utils;

fn build_tester() -> Criterion {
    Criterion::default()
        .sample_size(200)
        .measurement_time(std::time::Duration::from_secs(10))
        .noise_threshold(0.03)
        .significance_level(0.01)
}

criterion_main!(levels_benches);

criterion_group!(
    name = levels_benches;
    config = build_tester();
    targets = freezing_levels_bench, max_temperature_in_profile_bench,
              max_temperature_in_layer_bench
);

fn freezing_levels_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("freezing_levels", |b| {
        b.iter(|| {
            for snd in &snds {
                sounding_analysis::freezing_levels(&snd).expect("oops");
            }
        });
    });
}

// Don't need to bench wet_bulb_zero_levels because it uses the same inner function as
// freezing_levels.

fn max_temperature_in_profile_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("max_temperature_in_profile", |b| {
        b.iter(|| {
            for snd in &snds {
                sounding_analysis::max_temperature_in_profile(&snd).expect("oops");
            }
        });
    });
}

// Don't need to bench max_wet_bulb_in_profile because it uses the same inner function as
// max_temperature_in_profile.

fn max_temperature_in_layer_bench(c: &mut Criterion) {
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

    c.bench_function("max_wet_bulb_in_profile", |b| {
        b.iter(|| {
            for snd in &snds {
                for layer in &layers {
                    sounding_analysis::max_temperature_in_layer(&snd, layer).expect("oops");
                }
            }
        });
    });
}

// Don't need to bench max_wet_bulb_in_layer because it uses the same inner function as
// max_temperature_in_layer.
