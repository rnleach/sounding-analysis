//! Data used in tests.

use std::fs::File;
use std::io::Read;
use std::path::PathBuf;
use std::str::FromStr;

use sounding_base::Profile;
use sounding_base::{Sounding, StationInfo};
use sounding_validate::validate;

pub fn create_test_sounding() -> Sounding {
    use sounding_base::Index::*;
    use sounding_base::Surface;

    let snd = Sounding::new()
        .set_station_info(StationInfo::new_with_values(1, (45.0, -115.0), 600.0))
        .set_valid_time(None)
        .set_lead_time(0)
        .set_index(Showalter, -2.0)
        .set_index(LI, -2.0)
        .set_index(SWeT, 35.0)
        .set_index(K, 45.0)
        .set_index(LCL, 850.0)
        .set_index(PWAT, 2.0)
        .set_index(TotalTotals, 55.0)
        .set_index(CAPE, 852.0)
        .set_index(LCLTemperature, 290.0)
        .set_index(CIN, -200.0)
        .set_index(EquilibrimLevel, 222.0)
        .set_index(LFC, 800.0)
        .set_index(BulkRichardsonNumber, 1.2)
        .set_profile(
            Profile::Pressure,
            vec![
                Option::from(1000.0),
                Option::from(975.0),
                Option::from(925.0),
                Option::from(900.0),
                Option::from(875.0),
                Option::from(850.0),
                Option::from(800.0),
                Option::from(700.0),
                Option::from(500.0),
                Option::from(300.0),
                Option::from(250.0),
                Option::from(200.0),
                Option::from(100.0),
            ],
        )
        .set_profile(
            Profile::Temperature,
            vec![
                Option::from(30.0),
                Option::from(29.0),
                Option::from(27.0),
                Option::from(25.0),
                Option::from(22.0),
                Option::from(20.0),
                Option::from(15.0),
                Option::from(2.0),
                Option::from(-10.0),
                Option::from(-20.0),
                Option::from(-30.0),
                Option::from(-50.0),
                Option::from(-45.0),
            ],
        )
        .set_profile(
            Profile::WetBulb,
            // FIXME: Actually calculate these to be correct
            vec![
                Option::from(20.0),
                Option::from(20.0),
                Option::from(20.0),
                Option::from(20.0),
                Option::from(20.0),
                Option::from(20.0),
                Option::from(14.0),
                Option::from(1.0),
                Option::from(-11.0),
                Option::from(-25.0),
                Option::from(-39.0),
                Option::from(-58.0),
                Option::from(-60.0),
            ],
        )
        .set_profile(
            Profile::DewPoint,
            vec![
                Option::from(18.0),
                Option::from(17.0),
                Option::from(16.0),
                Option::from(17.0),
                Option::from(19.0),
                Option::from(20.0),
                Option::from(13.0),
                Option::from(0.0),
                Option::from(-12.0),
                Option::from(-27.0),
                Option::from(-45.0),
                Option::from(-62.0),
                Option::from(-80.0),
            ],
        )
        .set_profile(
            Profile::WindDirection,
            vec![
                Option::from(0.0),
                Option::from(40.0),
                Option::from(80.0),
                Option::from(120.0),
                Option::from(160.0),
                Option::from(0.0),
                Option::from(40.0),
                Option::from(80.0),
                Option::from(120.0),
                Option::from(160.0),
                Option::from(200.0),
                Option::from(240.0),
                Option::from(280.0),
            ],
        )
        .set_profile(
            Profile::WindSpeed,
            vec![
                Option::from(5.0),
                Option::from(10.0),
                Option::from(15.0),
                Option::from(12.0),
                Option::from(27.0),
                Option::from(5.0),
                Option::from(10.0),
                Option::from(15.0),
                Option::from(12.0),
                Option::from(27.0),
                Option::from(45.0),
                Option::from(62.0),
                Option::from(80.0),
            ],
        )
        .set_profile(
            Profile::GeopotentialHeight,
            vec![
                Option::from(650.0),
                Option::from(700.0),
                Option::from(800.0),
                Option::from(900.0),
                Option::from(1000.0),
                Option::from(1050.0),
                Option::from(2000.0),
                Option::from(3000.0),
                Option::from(4000.0),
                Option::from(5000.0),
                Option::from(6500.0),
                Option::from(7000.0),
                Option::from(8000.0),
            ],
        )
        .set_profile(
            Profile::CloudFraction,
            vec![
                Option::from(0.0),
                Option::from(20.0),
                Option::from(40.0),
                Option::from(60.0),
                Option::from(80.0),
                Option::from(100.0),
                Option::from(85.0),
                Option::from(70.0),
                Option::from(50.0),
                Option::from(30.0),
                Option::from(25.0),
                Option::from(20.0),
                Option::from(10.0),
            ],
        )
        .set_surface_value(Surface::MSLP, 1014.0)
        .set_surface_value(Surface::StationPressure, 1010.0)
        .set_surface_value(Surface::UWind, 0.0)
        .set_surface_value(Surface::VWind, 0.0);

    assert!(validate(&snd).is_ok(), "Test data failed validation.");

    snd
}

/// Single dendritic layer, basic crossing
pub fn create_simple_dendtritic_test_sounding() -> Sounding {
    let snd = create_test_sounding()
        .set_profile(
            Profile::Temperature,
            vec![
                Option::from(-8.1),
                Option::from(-10.0),
                Option::from(-12.0),
                Option::from(-14.0),
                Option::from(-16.0),
                Option::from(-18.0),
                Option::from(-20.0),
                Option::from(-22.0),
                Option::from(-24.0),
                Option::from(-26.0),
                Option::from(-30.0),
                Option::from(-50.0),
                Option::from(-45.0),
            ],
        )
        .set_profile(
            Profile::WetBulb,
            // FIXME: Actually calculate these to be correct
            vec![
                Option::from(-8.5),
                Option::from(-10.5),
                Option::from(-12.5),
                Option::from(-14.5),
                Option::from(-16.5),
                Option::from(-18.5),
                Option::from(-20.5),
                Option::from(-22.5),
                Option::from(-24.5),
                Option::from(-27.0),
                Option::from(-39.0),
                Option::from(-58.0),
                Option::from(-60.0),
            ],
        )
        .set_profile(
            Profile::DewPoint,
            vec![
                Option::from(-9.0),
                Option::from(-11.00),
                Option::from(-13.0),
                Option::from(-15.0),
                Option::from(-17.0),
                Option::from(-19.0),
                Option::from(-21.0),
                Option::from(-23.0),
                Option::from(-25.0),
                Option::from(-28.0),
                Option::from(-45.0),
                Option::from(-62.0),
                Option::from(-80.0),
            ],
        );

    assert!(
        validate(&snd).is_ok(),
        "Simple dendritic test data failed validation."
    );

    snd
}

/// Multiple dendritic layer, jumps over, 3-layers
pub fn create_complex_dendtritic_test_sounding() -> Sounding {
    let snd = create_test_sounding()
        .set_profile(
            Profile::Temperature,
            vec![
                Option::from(-12.1),
                Option::from(-10.0),
                Option::from(-20.0),
                Option::from(-14.0),
                Option::from(-16.0),
                Option::from(-18.0),
                Option::from(-20.0),
                Option::from(-22.0),
                Option::from(-24.0),
                Option::from(-26.0),
                Option::from(-30.0),
                Option::from(-50.0),
                Option::from(-45.0),
            ],
        )
        .set_profile(
            Profile::WetBulb,
            // FIXME: Actually calculate these to be correct
            vec![
                Option::from(-12.5),
                Option::from(-10.5),
                Option::from(-20.5),
                Option::from(-14.5),
                Option::from(-16.5),
                Option::from(-18.5),
                Option::from(-20.5),
                Option::from(-22.5),
                Option::from(-24.5),
                Option::from(-27.0),
                Option::from(-39.0),
                Option::from(-58.0),
                Option::from(-60.0),
            ],
        )
        .set_profile(
            Profile::DewPoint,
            vec![
                Option::from(-13.0),
                Option::from(-11.00),
                Option::from(-21.0),
                Option::from(-15.0),
                Option::from(-17.0),
                Option::from(-19.0),
                Option::from(-21.0),
                Option::from(-23.0),
                Option::from(-25.0),
                Option::from(-28.0),
                Option::from(-45.0),
                Option::from(-62.0),
                Option::from(-80.0),
            ],
        );

    assert!(
        validate(&snd).is_ok(),
        "Complex dendritic test data failed validation."
    );

    snd
}

fn load_test_csv_sounding(location: &PathBuf) -> Sounding {
    let mut f = File::open(location).expect(&format!("Error opening file: {:?}", location));

    let mut contents = String::new();
    f.read_to_string(&mut contents)
        .expect(&format!("Error reading file: {:?}", location));

    let lines: Vec<&str> = contents.split('\n').skip(1).collect();
    let mut height: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut temp: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut dp: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut press: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut wspd: Vec<Option<f64>> = Vec::with_capacity(lines.len());
    let mut wdir: Vec<Option<f64>> = Vec::with_capacity(lines.len());

    for line in lines {
        let tokens: Vec<&str> = line.split(',').collect();
        if tokens.len() < 6 {
            continue;
        }
        height.push(f64::from_str(tokens[0]).ok());
        temp.push(f64::from_str(tokens[1]).ok());
        dp.push(f64::from_str(tokens[2]).ok());
        press.push(f64::from_str(tokens[3]).ok());
        wspd.push(f64::from_str(tokens[4]).ok());
        wdir.push(f64::from_str(tokens[5]).ok());
    }

    Sounding::new()
        .set_profile(Profile::GeopotentialHeight, height)
        .set_profile(Profile::Temperature, temp)
        .set_profile(Profile::DewPoint, dp)
        .set_profile(Profile::Pressure, press)
        .set_profile(Profile::WindSpeed, wspd)
        .set_profile(Profile::WindDirection, wdir)
}

#[test]
fn test_load_test_csv_sounding() {
    let mut test_path = PathBuf::new();
    test_path.push("test_data");
    test_path.push("standard.csv");
    let snd = load_test_csv_sounding(&test_path);
    assert!(validate(&snd).is_ok(), "Failed validation.");
}
