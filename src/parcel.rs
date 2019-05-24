//! Functions for doing parcel analysis on a sounding, specifically related to convection.
//!
use crate::{
    error::{AnalysisError, Result},
    interpolation::{linear_interp, linear_interpolate_sounding},
    profile::equivalent_potential_temperature,
};
use itertools::izip;
use metfor::{self, Celsius, HectoPascal, Kelvin, Quantity};
use optional::Optioned;
use sounding_base::{DataRow, Sounding};

/// Variables defining a parcel as used in parcel analysis.
#[derive(Debug, Clone, Copy)]
pub struct Parcel {
    /// Temperature in C
    pub temperature: Celsius,
    /// Pressure in hPa
    pub pressure: HectoPascal,
    /// Dew point in C
    pub dew_point: Celsius,
}

impl Parcel {
    /// Get the potential temperatures of the parcel
    pub fn theta(&self) -> Kelvin {
        metfor::theta(self.pressure, self.temperature)
    }

    /// Get the equivalent potential temperature of the parcel
    pub fn theta_e(&self) -> Result<Kelvin> {
        metfor::theta_e(self.temperature, self.dew_point, self.pressure)
            .ok_or(AnalysisError::MetForError)
    }

    /// Get the mixing ratio of the parcel
    pub fn mixing_ratio(&self) -> Result<f64> {
        metfor::mixing_ratio(self.dew_point, self.pressure).ok_or(AnalysisError::MetForError)
    }

    /// Get the virtual temperature of the parcel
    pub fn virtual_temperature(&self) -> Result<Kelvin> {
        metfor::virtual_temperature(self.temperature, self.dew_point, self.pressure)
            .ok_or(AnalysisError::MetForError)
    }

    /// Try to convert a `DataRow` to a `Parcel`.
    pub fn from_datarow(mut dr: DataRow) -> Option<Self> {
        let temperature = dr.temperature.take()?;
        let pressure = dr.pressure.take()?;
        let dew_point = dr.dew_point.take()?;

        Some(Parcel {
            temperature,
            pressure,
            dew_point,
        })
    }
}

/// Create a mixed layer parcel.
///
/// Calculated by averaging potential temperature and mixing ratio in the lowest 100 hPa of the
/// sounding and calculating a temperature and dew point at the lowest pressure level using those
/// averaged values.
pub fn mixed_layer_parcel(snd: &Sounding) -> Result<Parcel> {
    let press = snd.pressure_profile();
    let t = snd.temperature_profile();
    let dp = snd.dew_point_profile();

    if press.is_empty() || t.is_empty() || dp.is_empty() {
        return Err(AnalysisError::MissingProfile);
    }

    let bottom_p = press
        .iter()
        // Remove None values
        .filter_map(|&p| p.into_option())
        // Get the first of the non-None values
        .nth(0)
        // Check to make sure we got one, return error if not.
        .ok_or(AnalysisError::NoDataProfile)?;

    let (sum_p, sum_t, sum_mw) = izip!(press, t, dp)
        // remove levels with missing values
        .filter_map(|(p, t, dp)| {
            if p.is_some() && t.is_some() && dp.is_some() {
                Some((p.unpack(), t.unpack(), dp.unpack()))
            } else {
                None
            }
        })
        // only go up 100 hPa above the lowest level
        .take_while(|&(p, _, _)| p >= bottom_p - HectoPascal(100.0))
        // Get theta
        .map(|(p, t, dp)| (p, dp, metfor::theta(p, t)))
        // convert to mw
        .filter_map(|(p, dp, th)| metfor::mixing_ratio(dp, p).and_then(|mw| Some((p, th, mw))))
        // calculate the sums and count needed for the average
        .fold((HectoPascal(0.0), 0.0f64, 0.0f64), |acc, (p, theta, mw)| {
            let (sum_p, sum_t, sum_mw) = acc;
            (
                sum_p + p,
                sum_t + theta.unpack() * p.unpack(),
                sum_mw + mw * p.unpack(),
            )
        });

    if sum_p == HectoPascal(0.0) {
        return Err(AnalysisError::NotEnoughData);
    }

    // convert back to temperature and dew point at the lowest pressure level.
    let pressure = bottom_p;
    let temperature = Celsius::from(metfor::temperature_from_theta(
        Kelvin(sum_t / sum_p.unpack()),
        bottom_p,
    ));
    let dew_point = metfor::dew_point_from_p_and_mw(bottom_p, sum_mw / sum_p.unpack())
        .ok_or(AnalysisError::InvalidInput)?;

    Ok(Parcel {
        temperature,
        pressure,
        dew_point,
    })
}

/// Get a surface parcel.
pub fn surface_parcel(snd: &Sounding) -> Result<Parcel> {
    let pressure = snd.station_pressure().ok_or(AnalysisError::MissingValue)?;
    let temperature = snd.sfc_temperature().ok_or(AnalysisError::MissingValue)?;
    let dew_point = snd.sfc_dew_point().ok_or(AnalysisError::MissingValue)?;

    Ok(Parcel {
        temperature,
        pressure,
        dew_point,
    })
}

/// Get the lowest level parcel. This should be the surface parcel, but some files do not have
/// complete information at the surface, so the first level above the ground is best you can do.
pub fn lowest_level_parcel(snd: &Sounding) -> Result<Parcel> {
    snd.bottom_up()
        .filter_map(Parcel::from_datarow)
        .nth(0)
        .ok_or(AnalysisError::NotEnoughData)
}

/// Get the parcel at a specific pressure.
pub fn pressure_parcel(snd: &Sounding, pressure: HectoPascal) -> Result<Parcel> {
    let row = linear_interpolate_sounding(snd, pressure)?;

    Parcel::from_datarow(row).ok_or(AnalysisError::MissingValue)
}

/// Get the most unstable parcel.
///
/// This is defined as the parcel in the lowest 300 hPa of the sounding with the highest equivalent
/// potential temperature.
pub fn most_unstable_parcel(snd: &Sounding) -> Result<Parcel> {
    let temp_vec: Vec<Optioned<Kelvin>>;

    let theta_e = if !snd.theta_e_profile().is_empty() {
        snd.theta_e_profile()
    } else {
        temp_vec = equivalent_potential_temperature(snd);
        &temp_vec
    };
    let press = snd.pressure_profile();

    let top_pressure = press
        .iter()
        .filter_map(|p| p.into_option())
        .nth(0)
        .ok_or(AnalysisError::NotEnoughData)?
        - HectoPascal(300.0);

    let (idx, _) = izip!(0.., press, theta_e)
        .take_while(|&(_, press_opt, _)| {
            if press_opt.is_some() {
                press_opt.unpack() >= top_pressure
            } else {
                true // Don't stop just because one is missing or the none value
            }
        })
        .filter_map(|(idx, _, theta_e_opt)| {
            if theta_e_opt.is_some() {
                Some((idx, theta_e_opt.unpack()))
            } else {
                None
            }
        })
        .fold(
            (::std::usize::MAX, metfor::ABSOLUTE_ZERO),
            |(max_idx, max_val), (i, theta_e)| {
                if theta_e > max_val {
                    (i, theta_e)
                } else {
                    (max_idx, max_val)
                }
            },
        );

    let row = snd.data_row(idx).ok_or(AnalysisError::NoDataProfile)?;

    Parcel::from_datarow(row).ok_or(AnalysisError::MissingValue)
}

/// Get the convective parcel - this is the surface parcel that will rise without any lifting.
///
/// This is based off the mixing ratio of the mixed layer parcel.
pub fn convective_parcel(snd: &Sounding) -> Result<Parcel> {
    let Parcel {
        dew_point: dp,
        pressure: initial_p,
        temperature: initial_t,
    } = mixed_layer_parcel(snd)?;

    let tgt_mw = metfor::mixing_ratio(dp, initial_p).ok_or(AnalysisError::MetForError)?;

    let pressure = snd.pressure_profile();
    let temperature = snd.temperature_profile();

    let (tgt_p, tgt_t) = izip!(pressure, temperature)
        .filter_map(|(p, t)| {
            let p = p.into_option()?;
            let t = t.into_option()?;
            let mw = metfor::mixing_ratio(t, p)?;

            Some((p, t, mw))
        })
        // Get past a cool surface layer where mw < tgt_mw possible due to using a mixed parcel for
        // the target mw in combination with a surface based inversion.
        .skip_while(|(_, _, mw)| *mw <= tgt_mw)
        // Scan to find the level where the environmental mixing ratio becomes less than our parcel
        // (surface or mixed layer parcel) mixing ratio. This is the CCL, the level where a cloud
        // first forms.
        .scan((0.0, 0.0, 0.0), |(old_p, old_t, old_mw), (p, t, mw)| {
            let result = if mw <= tgt_mw && *old_mw >= tgt_mw {
                // found the crossing
                let tgt_p = HectoPascal(linear_interp(tgt_mw, *old_mw, mw, *old_p, p.unpack()));
                let tgt_t = Celsius(linear_interp(tgt_mw, *old_mw, mw, *old_t, t.unpack()));

                Some(Some((tgt_p, tgt_t)))
            } else {
                Some(None)
            };

            // save these as the new ones
            *old_p = p.unpack();
            *old_t = t.unpack();
            *old_mw = mw;

            result
        })
        // Each layer that is not the CCL will be a Some(None), so skip past it
        .filter_map(|opt_opt| opt_opt)
        // Grab the first (and only) layer where the mixing ratios meet.
        .nth(0)
        // Probably a bad choice to use an error to signal this, but it is impossible to tell if
        // we never found it because there wasn't enough data, or it just didn't exist.
        .ok_or(AnalysisError::NotEnoughData)?;

    // Extrapolate dry adiabatically back to the parcel level.
    let tgt_theta = metfor::theta(tgt_p, tgt_t);
    let tgt_t = Celsius::from(metfor::temperature_from_theta(tgt_theta, initial_p)).max(initial_t);

    Ok(Parcel {
        pressure: initial_p,
        temperature: tgt_t,
        dew_point: dp,
    })
}
