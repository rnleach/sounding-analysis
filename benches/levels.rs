#![feature(test)]

extern crate test;
use test::Bencher;

mod utils;

#[bench]
fn test_freezing_levels(b: &mut Bencher) {
    let snds = utils::load_all_test_files();

    let test = || {
        for snd in &snds {
            sounding_analysis::freezing_levels(&snd).expect("oops");
        }
    };

    b.iter(test)
}

#[bench]
fn test_max_temperature_in_profile(b: &mut Bencher) {
    let snds = utils::load_all_test_files();

    let test = || {
        for snd in &snds {
            sounding_analysis::max_temperature_in_profile(&snd).expect("oops");
        }
    };

    b.iter(test)
}
