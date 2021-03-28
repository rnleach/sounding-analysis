//! Run these benches with `cargo bench --bench funtions -- --verbose`
//!
//! Run with `cargo bench --bench funtions -- --verbose vapor_pressure_over_ice` to select the
//! single benchmark.

use criterion::{criterion_group, criterion_main, Criterion};

mod utils;

criterion_main!(levels_benches);

criterion_group!(
    levels_benches,
    freezing_levels_bench,
    max_temperature_in_profile_bench
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
