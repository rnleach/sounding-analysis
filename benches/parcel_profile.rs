//! Run these benches with `cargo bench --bench parcel_profile -- --verbose`
use criterion::{criterion_group, criterion_main, Criterion};

mod utils;

fn build_tester() -> Criterion {
    Criterion::default()
        .sample_size(200)
        .measurement_time(std::time::Duration::from_secs(10))
        .noise_threshold(0.03)
        .significance_level(0.01)
}

criterion_main!(parcel_profile_benches);

criterion_group!(
    name = parcel_profile_benches;
    config = build_tester();
    targets = lift_parcel_bench, robust_convective_parcel_ascent_bench, mix_down_bench, dcape_bench
);

fn lift_parcel_bench(c: &mut Criterion) {
    use sounding_analysis::{Parcel, Sounding};

    let snds = utils::load_all_test_files().to_vec();

    let pairs: Vec<(Sounding, Parcel)> = snds
        .into_iter()
        .map(|snd| {
            let parcel = sounding_analysis::mixed_layer_parcel(&snd).unwrap();
            (snd, parcel)
        })
        .collect();

    c.bench_function("lift_parcel", |b| {
        b.iter(|| {
            for (snd, parcel) in &pairs {
                let _x = sounding_analysis::lift_parcel(*parcel, snd).expect("oops");
            }
        });
    });
}

fn robust_convective_parcel_ascent_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("robust_convective_parcel_ascent", |b| {
        b.iter(|| {
            for snd in &snds {
                let _x = sounding_analysis::robust_convective_parcel_ascent(snd).expect("oops");
            }
        });
    });
}

fn mix_down_bench(c: &mut Criterion) {
    use metfor::HectoPascal;
    use sounding_analysis::{Parcel, Sounding};

    let snds = utils::load_all_test_files().to_vec();

    let pairs: Vec<(Sounding, Parcel)> = snds
        .into_iter()
        .map(|snd| {
            let parcel = sounding_analysis::pressure_parcel(&snd, HectoPascal(700.0)).unwrap();
            (snd, parcel)
        })
        .collect();

    c.bench_function("mix_down", |b| {
        b.iter(|| {
            for (snd, parcel) in &pairs {
                let _x = sounding_analysis::mix_down(*parcel, snd).expect("oops");
            }
        });
    });
}

fn dcape_bench(c: &mut Criterion) {
    let snds = utils::load_all_test_files();

    c.bench_function("dcape", |b| {
        b.iter(|| {
            for (i, snd) in snds.iter().enumerate() {
                let x = sounding_analysis::dcape(snd);
                match i {
                    1 | 2 => continue,
                    _ => x.expect("oops"),
                };
            }
        });
    });
}
