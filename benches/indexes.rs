//! Run these benches with `cargo bench --bench indexes -- --verbose`

use criterion::{criterion_group, criterion_main, Criterion};

mod utils;

criterion_main!(inexes_benches);

criterion_group!(
    inexes_benches,
    precipitable_water_bench,
    haines_mid_bench,
    hot_dry_windy_bench
);

fn precipitable_water_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("precipitable_water", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::precipitable_water(&snd).expect("oops");
            }
        });
    });
}

fn haines_mid_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("haines_mid", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::haines_mid(&snd).expect("oops");
            }
        });
    });
}

fn hot_dry_windy_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("hot_dry_windy", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::hot_dry_windy(&snd).expect("oops");
            }
        });
    });
}
