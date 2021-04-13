//! Run these benches with `cargo bench --bench fire -- --verbose`
use criterion::{criterion_group, criterion_main, Criterion};

mod utils;

criterion_main!(fire_benches);

criterion_group!(
    fire_benches,
    calc_plumes_bench,
    blow_up_bench,
    plume_heating_analysis_bench,
    analyze_plume_parcel_bench,
    lift_plume_parcel_bench
);

fn calc_plumes_bench(c: &mut Criterion) {
    use metfor::CelsiusDiff;
    use sounding_analysis::Sounding;

    const MIN_DT: CelsiusDiff = CelsiusDiff(-2.0);
    const MAX_DT: CelsiusDiff = CelsiusDiff(20.0);
    const INCREMENT: CelsiusDiff = CelsiusDiff(0.1);

    let snds = utils::load_all_test_files().to_vec();

    let args: Vec<(Sounding, Option<f64>)> = snds
        .into_iter()
        .flat_map(|snd| {
            [None, Some(8.0), Some(15.0)]
                .iter()
                .map(move |moisture_ratio| (snd.clone(), *moisture_ratio))
        })
        .collect();

    c.bench_function("calc_plumes", |b| {
        b.iter(|| {
            for (snd, moisture_ratio) in &args {
                let _x = sounding_analysis::experimental::fire::calc_plumes(
                    snd,
                    INCREMENT,
                    MIN_DT,
                    MAX_DT,
                    *moisture_ratio,
                )
                .expect("oops");
            }
        });
    });
}

fn blow_up_bench(c: &mut Criterion) {
    use sounding_analysis::Sounding;

    let snds = utils::load_all_test_files().to_vec();

    let args: Vec<(Sounding, Option<f64>)> = snds
        .into_iter()
        .flat_map(|snd| {
            [None, Some(8.0), Some(15.0)]
                .iter()
                .map(move |moisture_ratio| (snd.clone(), *moisture_ratio))
        })
        .collect();

    c.bench_function("blow_up", |b| {
        b.iter(|| {
            for (snd, moisture_ratio) in &args {
                let _x = sounding_analysis::experimental::fire::blow_up(snd, *moisture_ratio)
                    .expect("oops");
            }
        });
    });
}

fn plume_heating_analysis_bench(c: &mut Criterion) {
    use sounding_analysis::Sounding;

    let snds = utils::load_all_test_files().to_vec();

    let args: Vec<(Sounding, Option<f64>)> = snds
        .into_iter()
        .flat_map(|snd| {
            [None, Some(8.0), Some(15.0)]
                .iter()
                .map(move |moisture_ratio| (snd.clone(), *moisture_ratio))
        })
        .collect();

    c.bench_function("plume_heating_analysis", |b| {
        b.iter(|| {
            for (snd, moisture_ratio) in &args {
                let _x = sounding_analysis::experimental::fire::plume_heating_analysis(
                    snd,
                    *moisture_ratio,
                )
                .expect("oops");
            }
        });
    });
}

fn analyze_plume_parcel_bench(c: &mut Criterion) {
    use sounding_analysis::{Parcel, Sounding};
    const HEATING: metfor::CelsiusDiff = metfor::CelsiusDiff(5.5);

    let snds = utils::load_all_test_files().to_vec();

    let args: Vec<(Sounding, Parcel)> = snds
        .into_iter()
        .flat_map(|snd| {
            [None, Some(8.0), Some(15.0)]
                .iter()
                .map(move |moisture_ratio| (snd.clone(), *moisture_ratio))
        })
        .filter_map(|(snd, moisture_ratio)| {
            sounding_analysis::mixed_layer_parcel(&snd)
                .ok()
                .map(|pcl| (snd, moisture_ratio, pcl))
        })
        .map(|(snd, moisture_ratio, env_pcl)| {
            (
                snd,
                sounding_analysis::experimental::fire::create_plume_parcel_from(
                    env_pcl,
                    HEATING,
                    moisture_ratio,
                ),
            )
        })
        .collect();

    c.bench_function("analyze_plume_parcel", |b| {
        b.iter(|| {
            for (snd, pcl) in &args {
                let _x = sounding_analysis::experimental::fire::analyze_plume_parcel(*pcl, snd)
                    .expect("oops");
            }
        });
    });
}

fn lift_plume_parcel_bench(c: &mut Criterion) {
    use sounding_analysis::{Parcel, Sounding};
    const HEATING: metfor::CelsiusDiff = metfor::CelsiusDiff(5.5);

    let snds = utils::load_all_test_files().to_vec();

    let args: Vec<(Sounding, Parcel)> = snds
        .into_iter()
        .flat_map(|snd| {
            [None, Some(8.0), Some(15.0)]
                .iter()
                .map(move |moisture_ratio| (snd.clone(), *moisture_ratio))
        })
        .filter_map(|(snd, moisture_ratio)| {
            sounding_analysis::mixed_layer_parcel(&snd)
                .ok()
                .map(|pcl| (snd, moisture_ratio, pcl))
        })
        .map(|(snd, moisture_ratio, env_pcl)| {
            (
                snd,
                sounding_analysis::experimental::fire::create_plume_parcel_from(
                    env_pcl,
                    HEATING,
                    moisture_ratio,
                ),
            )
        })
        .collect();

    c.bench_function("lift_plume_parcel", |b| {
        b.iter(|| {
            for (snd, pcl) in &args {
                let _x = sounding_analysis::experimental::fire::lift_plume_parcel(*pcl, snd)
                    .expect("oops");
            }
        });
    });
}
