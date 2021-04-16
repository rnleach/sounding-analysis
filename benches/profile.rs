//! Run these benches with `cargo bench --bench profile -- --verbose`
use criterion::{criterion_group, criterion_main, Criterion};

mod utils;

fn build_tester() -> Criterion {
    Criterion::default()
        .sample_size(200)
        .measurement_time(std::time::Duration::from_secs(10))
        .noise_threshold(0.03)
        .significance_level(0.01)
}

criterion_main!(profile_benches);

criterion_group!(
    name = profile_benches;
    config = build_tester();
    targets = wet_bulb_bench, relative_humidity_bench, relative_humidity_ice_bench,
              potential_temperature_bench, equivalent_potential_temperature_bench,
              temperature_lapse_rate_bench, sfc_to_level_temperature_lapse_rate_bench,
              hydrolapse_bench
);

fn wet_bulb_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("wet_bulb", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::wet_bulb(&snd);
            }
        });
    });
}

fn relative_humidity_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("relative_humidity", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::relative_humidity(&snd);
            }
        });
    });
}

fn relative_humidity_ice_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("relative_humidity_ice", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::relative_humidity_ice(&snd);
            }
        });
    });
}

fn potential_temperature_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("potential_temperature", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::potential_temperature(&snd);
            }
        });
    });
}

fn equivalent_potential_temperature_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("equivalent_potential_temperature", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::equivalent_potential_temperature(&snd);
            }
        });
    });
}

fn temperature_lapse_rate_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("temperature_lapse_rate", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::temperature_lapse_rate(&snd);
            }
        });
    });
}

fn sfc_to_level_temperature_lapse_rate_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("sfc_to_level_temperature_lapse_rate", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::sfc_to_level_temperature_lapse_rate(&snd);
            }
        });
    });
}

fn hydrolapse_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("hydrolapse", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::hydrolapse(&snd);
            }
        });
    });
}
