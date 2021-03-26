use metfor::{Celsius, HectoPascal, Knots, Meters, WindSpdDir};
use optional::Optioned;
use sounding_analysis::{Sounding, StationInfo};
use std::{fs::File, io::Read, path::PathBuf, str::FromStr};

pub fn load_all_test_files() -> [Sounding; 4] {
    let snd1 = load_test_file("standard.csv");
    let snd2 = load_test_file("complex_dendritic.csv");
    let snd3 = load_test_file("multiple_warm_layers_aloft.csv");
    let snd4 = load_test_file("multiple_inversions_aloft.csv");

    [snd1, snd2, snd3, snd4]
}

fn load_test_file(fname: &str) -> Sounding {
    let mut test_path = PathBuf::new();
    test_path.push("test_data");
    test_path.push(fname);
    load_test_csv_sounding(&test_path)
}

fn load_test_csv_sounding(location: &PathBuf) -> Sounding {
    let mut f = File::open(location).expect(&format!("Error opening file: {:#?}", location));

    let mut contents = String::new();
    f.read_to_string(&mut contents)
        .expect(&format!("Error reading file: {:#?}", location));

    let lines: Vec<&str> = contents.split('\n').collect();
    let mut line_iter = lines.iter();

    //
    // Parse profile data
    //
    let mut height: Vec<Optioned<Meters>> = Vec::with_capacity(lines.len());
    let mut temp: Vec<Optioned<Celsius>> = Vec::with_capacity(lines.len());
    let mut wb: Vec<Optioned<Celsius>> = Vec::with_capacity(lines.len());
    let mut dp: Vec<Optioned<Celsius>> = Vec::with_capacity(lines.len());
    let mut press: Vec<Optioned<HectoPascal>> = Vec::with_capacity(lines.len());
    let mut wind: Vec<Optioned<WindSpdDir<Knots>>> = Vec::with_capacity(lines.len());

    for line in line_iter.by_ref() {
        if line.starts_with("### Surface Data ###")
            || line.starts_with("### Analysis Int Section ###")
            || line.starts_with("### Analysis Float Section ###")
        {
            break;
        }

        let tokens: Vec<&str> = line.split(',').collect();
        if tokens.len() < 6 {
            continue;
        }
        let t_c: Option<Celsius> = f64::from_str(tokens[1]).ok().map(Celsius);
        let dp_c: Option<Celsius> = f64::from_str(tokens[2]).ok().map(Celsius);
        let press_hpa: Option<HectoPascal> = f64::from_str(tokens[3]).ok().map(HectoPascal).into();
        let wb_c: Option<Celsius> = t_c
            .and_then(|t| dp_c.and_then(|dp| press_hpa.and_then(|p| metfor::wet_bulb(t, dp, p))));
        let wspd = f64::from_str(tokens[4]).ok();
        let wdir = f64::from_str(tokens[5]).ok();
        let wind_val: Option<WindSpdDir<Knots>> = wspd.and_then(|wspd| {
            wdir.map(|wdir| WindSpdDir {
                speed: Knots(wspd),
                direction: wdir,
            })
        });

        height.push(f64::from_str(tokens[0]).ok().map(Meters).into());
        temp.push(t_c.into());
        wb.push(wb_c.into());
        dp.push(dp_c.into());
        press.push(press_hpa.into());
        wind.push(wind_val.into());
    }

    let mut snd = Sounding::new()
        .with_height_profile(height)
        .with_temperature_profile(temp)
        .with_wet_bulb_profile(wb)
        .with_dew_point_profile(dp)
        .with_pressure_profile(press)
        .with_wind_profile(wind);

    //
    // Surface data
    //
    for line in line_iter.by_ref() {
        if line.starts_with("### Analysis Int Section ###")
            || line.starts_with("### Analysis Float Section ###")
        {
            break;
        }

        let tokens: Vec<&str> = line.split(',').collect();
        if tokens.len() < 6 {
            continue;
        }

        let height = f64::from_str(tokens[0]).ok().map(Meters);
        let t_c = f64::from_str(tokens[1]).ok().map(Celsius);
        let dp_c = f64::from_str(tokens[2]).ok().map(Celsius);
        let press_hpa = f64::from_str(tokens[3]).ok().map(HectoPascal);
        let wspd = f64::from_str(tokens[4]).ok();
        let wdir = f64::from_str(tokens[5]).ok();
        let wind = wspd.and_then(|wspd| {
            wdir.map(|wdir| WindSpdDir {
                speed: Knots(wspd),
                direction: wdir,
            })
        });

        snd = snd
            .with_sfc_temperature(t_c)
            .with_station_info(StationInfo::new().with_elevation(height))
            .with_sfc_dew_point(dp_c)
            .with_station_pressure(press_hpa)
            .with_sfc_wind(wind)
    }

    snd
}
