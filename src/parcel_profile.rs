//! Create and analyze a profile from lifting or descending a parcel.
use crate::{
    error::{AnalysisError, Result},
    interpolation::linear_interpolate_sounding,
    parcel::Parcel,
    sounding::{DataRow, Sounding},
};
use itertools::izip;
use metfor::{self, Celsius, CelsiusDiff, HectoPascal, JpKg, Kelvin, Meters, MetersPSec, Quantity};
use optional::Optioned;

/// Hold profiles for a parcel and it's environment.
#[derive(Debug, Clone)]
pub struct ParcelProfile {
    /// Pressure profile
    pub pressure: Vec<HectoPascal>,
    /// Height profile
    pub height: Vec<Meters>,
    /// Parcel virtual temperature profile
    pub parcel_t: Vec<Celsius>,
    /// Environment virtual temperature profile
    pub environment_t: Vec<Celsius>,
}

pub(crate) mod lift;

/// Parcel analysis, this is a way to package the analysis of a parcel.
///
/// These are done by converting the profiles to virtual temperature. It is assumed the reason for
/// lifting the parcel and doing the analysis is related to bouyancy and some kind of convection
/// or stability analysis.
#[derive(Debug, Clone)]
pub struct ParcelAscentAnalysis {
    // The orginal parcel and profile
    parcel: Parcel,
    profile: ParcelProfile,

    // Indicies from analysis
    cape: Optioned<JpKg>,
    hail_cape: Optioned<JpKg>,
    ncape: Optioned<f64>,
    lcl_height_agl: Optioned<Meters>,    // cloud base for aviation
    lcl_pressure: Optioned<HectoPascal>, // plotting on skew-t
    lcl_temperature: Optioned<Celsius>,  // ice or ice/water cloud?
    cin: Optioned<JpKg>,
    el_pressure: Optioned<HectoPascal>,      // plotting on skew-t
    el_height_asl: Optioned<Meters>,         // Calculating convective cloud tops for aviation
    el_temperature: Optioned<Celsius>,       // useful for comparing to satellite
    lfc_pressure: Optioned<HectoPascal>,     // plotting on skew-t
    lfc_virt_temperature: Optioned<Celsius>, // plotting on skew-t
}

impl ParcelAscentAnalysis {
    /// Get the CAPE.
    pub fn cape(&self) -> Optioned<JpKg> {
        self.cape
    }

    /// Get the CAPE in the hail growth zone.
    pub fn hail_cape(&self) -> Optioned<JpKg> {
        self.hail_cape
    }

    /// Get the normalized cape.
    pub fn ncape(&self) -> Optioned<f64> {
        self.ncape
    }

    /// Get the LCL height AGL.
    pub fn lcl_height_agl(&self) -> Optioned<Meters> {
        self.lcl_height_agl
    }

    /// Get the LCL pressrue level.
    pub fn lcl_pressure(&self) -> Optioned<HectoPascal> {
        self.lcl_pressure
    }

    /// Get the temperature at the LCL.
    pub fn lcl_temperature(&self) -> Optioned<Celsius> {
        self.lcl_temperature
    }

    /// Get the CIN.
    pub fn cin(&self) -> Optioned<JpKg> {
        self.cin
    }

    /// Get the pressure at the equilibrium level.
    pub fn el_pressure(&self) -> Optioned<HectoPascal> {
        self.el_pressure
    }

    /// Get the height ASL of the equilibrium level.
    pub fn el_height_asl(&self) -> Optioned<Meters> {
        self.el_height_asl
    }

    /// Get the temperature at the equilibrium level.
    pub fn el_temperature(&self) -> Optioned<Celsius> {
        self.el_temperature
    }

    /// Get the pressure at the LFC.
    pub fn lfc_pressure(&self) -> Optioned<HectoPascal> {
        self.lfc_pressure
    }

    /// Get the virtual temperature at the LFC.
    pub fn lfc_virt_temperature(&self) -> Optioned<Celsius> {
        self.lfc_virt_temperature
    }

    /// Retrieve the parcel's profile
    #[inline]
    pub fn profile(&self) -> &ParcelProfile {
        &self.profile
    }

    /// Retrieve the original parcel.
    #[inline]
    pub fn parcel(&self) -> &Parcel {
        &self.parcel
    }

    /// Calculate the parcel vertical speed at the equilibrium level. Note that this is an over
    /// estimate of updraft speed due to the effects of entrainment and water/ice loading.
    #[inline]
    pub fn calculate_cape_speed(&self) -> Option<MetersPSec> {
        self.cape
            .map(|cape| MetersPSec::pack(f64::sqrt(2.0 * cape.unpack())))
    }
}

/// Lift a parcel for a convective parcel analysis.
///
/// The resulting `ParcelProfile` and analysis are based off of virtual temperatures and the idea
/// that if there is no *moist* convection, or convective cloud, then there is no CAPE or CIN.
pub fn lift_parcel(parcel: Parcel, snd: &Sounding) -> Result<ParcelAscentAnalysis> {
    lift::lift_parcel(parcel, snd)
}

/// In order for parcel lifting to work and create a parallel environmental profile, we need to
/// start at a level in the sounding with pressure, height, temperature, and dew point. Otherwise
/// we end up with too much missing data in the sounding.
pub(crate) fn find_parcel_start_data(snd: &Sounding, parcel: &Parcel) -> Result<(DataRow, Parcel)> {
    let good_row = |row: &DataRow| -> bool {
        row.temperature.is_some()
            && row.dew_point.is_some()
            && row.pressure.is_some()
            && row.height.is_some()
    };

    let first_guess = linear_interpolate_sounding(snd, parcel.pressure)?;
    if good_row(&first_guess) {
        return Ok((first_guess, *parcel));
    }

    let second_guess = snd
        .bottom_up()
        .find(good_row)
        .ok_or(AnalysisError::NotEnoughData)?;

    let pressure = second_guess.pressure.ok_or(AnalysisError::InvalidInput)?;
    let theta = parcel.theta();
    let temperature = Celsius::from(metfor::temperature_from_theta(theta, pressure));
    let mw = parcel.mixing_ratio()?;
    let dew_point =
        metfor::dew_point_from_p_and_mw(pressure, mw).ok_or(AnalysisError::MetForError)?;
    let new_parcel = Parcel {
        pressure,
        temperature,
        dew_point,
    };

    Ok((second_guess, new_parcel))
}

/// A more robust convective parcel analysis.
///
/// Some approximations are used in many algorithms which are usually good enough. However,
/// sometimes they are close but miss the mark. Not to mention we are using linear interpolation
/// in so many places. Convective parcel analysis is one of those areas where sometimes this comes
/// up and we have a "convective parcel" with the equilibrium level below the lifting condensation
/// level. It's almost always very close though.
///
/// This algorithm finds the convective parcel the fast way, and if it is good, then it just uses
/// that parcel. Otherwise it tweaks the parcel to find a better convective parcel. Better meaning
/// that the EL is above or equal to the LCL. This algorithm  is MUCH slower in cases where the
/// 'fast way' doesn't work.
pub fn robust_convective_parcel_ascent(snd: &Sounding) -> Result<ParcelAscentAnalysis> {
    let mut start_parcel = crate::parcel::convective_parcel(snd)?;
    let mut analysis = lift_parcel(start_parcel, snd)?;

    if analysis.lcl_pressure <= analysis.el_pressure {
        // We've got a problem, so refine our parcel
        let mut warmer_t = start_parcel.temperature + CelsiusDiff(1.0);
        let mut warmer_parcel = Parcel {
            temperature: warmer_t,
            ..start_parcel
        };
        let mut anal = lift_parcel(warmer_parcel, snd)?;

        // Bracket the convective t in a 1 C range
        while anal.lcl_pressure <= anal.el_pressure {
            start_parcel = warmer_parcel;
            warmer_t += CelsiusDiff(1.0);
            warmer_parcel = Parcel {
                temperature: warmer_t,
                ..start_parcel
            };
            anal = lift_parcel(warmer_parcel, snd)?;
        }
        analysis = anal;

        let mut diff = 1.0;
        while diff > 0.1 {
            // Within 0.1 is probably overkill
            diff = (warmer_parcel.temperature - start_parcel.temperature).unpack() / 2.0;
            let mid_t = start_parcel.temperature + CelsiusDiff(diff);
            let mid_parcel = Parcel {
                temperature: mid_t,
                ..start_parcel
            };
            anal = lift_parcel(mid_parcel, snd)?;
            if anal.lcl_pressure <= anal.el_pressure {
                start_parcel = mid_parcel;
            } else {
                warmer_parcel = mid_parcel;
                analysis = anal;
            }
        }
    }

    Ok(analysis)
}

/// Descend a parcel dry adiabatically.
///
/// The resulting `ParcelProfile` has actual temperatures and not virtual temperatures. This is for
/// analyzing inversions and visualizing what a sounding would look like if deep, dry mixing were
/// to occur from surface heating alone.
pub fn mix_down(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let theta = parcel.theta();
    let theta_func = |theta_val, press| {
        Some(Celsius::from(metfor::temperature_from_theta(
            theta_val, press,
        )))
    };

    descend_parcel(parcel, snd, theta, theta_func, false, false)
}

/// Descend a parcel moist adiabatically.
///
/// The resulting `ParcelProfile` has virtual temperatures and is intended for calculating
/// DCAPE.
fn descend_moist(parcel: Parcel, snd: &Sounding) -> Result<ParcelProfile> {
    let theta = parcel.theta_e()?;

    let theta_func =
        |theta_e, press| metfor::temperature_from_theta_e_saturated_and_pressure(press, theta_e);

    descend_parcel(parcel, snd, theta, theta_func, true, true)
}

#[inline]
fn descend_parcel<F>(
    parcel: Parcel,
    snd: &Sounding,
    theta: Kelvin,
    theta_func: F,
    saturated: bool,
    virtual_t: bool,
) -> Result<ParcelProfile>
where
    F: Fn(Kelvin, HectoPascal) -> Option<Celsius>,
{
    let mut pressure = Vec::new();
    let mut height = Vec::new();
    let mut parcel_t = Vec::new();
    let mut environment_t = Vec::new();

    // Actually start at the bottom and work up.
    let press = snd.pressure_profile();
    let env_t = snd.temperature_profile();
    let env_dp = snd.dew_point_profile();
    let hght = snd.height_profile();

    let pcl_mw = parcel.mixing_ratio()?;

    // Nested scope to limit borrows
    {
        // Helper function to add row to parcel profile
        let mut add_row = |pp, hh, pcl_tt, env_tt| {
            pressure.push(pp);
            height.push(hh);
            parcel_t.push(pcl_tt);
            environment_t.push(env_tt);
        };

        izip!(press, hght, env_t, env_dp)
            .take_while(|(p_opt, _, _, _)| {
                if p_opt.is_some() {
                    p_opt.unpack() >= parcel.pressure
                } else {
                    true // Just skip over levels with missing data
                }
            })
            // Remove levels with missing data
            .filter_map(|(p_opt, h_opt, e_t_opt, e_dp_opt)| {
                if p_opt.is_some() && h_opt.is_some() && e_t_opt.is_some() && e_dp_opt.is_some() {
                    Some((
                        p_opt.unpack(),
                        h_opt.unpack(),
                        e_t_opt.unpack(),
                        e_dp_opt.unpack(),
                    ))
                } else {
                    None
                }
            })
            // Get the parcel temperature
            .filter_map(|(p, h, e_t, e_dp)| {
                theta_func(theta, p).map(|pcl_t| (p, h, pcl_t, e_t, e_dp))
            })
            // Get the parcel dew point
            .filter_map(|(p, h, pcl_t, e_t, e_dp)| {
                let p_dp = if saturated {
                    pcl_t
                } else {
                    metfor::dew_point_from_p_and_mw(p, pcl_mw)?
                };

                Some((p, h, pcl_t, p_dp, e_t, e_dp))
            })
            // Convert to virtual temperature if needed.
            .filter_map(|(p, h, pcl_t, p_dp, e_t, e_dp)| {
                let pcl_t = if virtual_t {
                    Celsius::from(metfor::virtual_temperature(pcl_t, p_dp, p)?)
                } else {
                    pcl_t
                };

                let e_t = if virtual_t {
                    Celsius::from(metfor::virtual_temperature(e_t, e_dp, p)?)
                } else {
                    e_t
                };

                Some((p, h, pcl_t, e_t))
            })
            .for_each(|(p, h, pt, et)| {
                add_row(p, h, pt, et);
            });

        // Add the parcel layer also
        let parcel_level = linear_interpolate_sounding(snd, parcel.pressure)?;
        let parcel_height = parcel_level.height.ok_or(AnalysisError::MissingValue)?;
        let env_t = parcel_level
            .temperature
            .ok_or(AnalysisError::MissingValue)?;
        add_row(parcel.pressure, parcel_height, parcel.temperature, env_t);
    }

    Ok(ParcelProfile {
        pressure,
        height,
        parcel_t,
        environment_t,
    })
}

/// Downdraft CAPE.
///
/// Defined as the net area between a parcel descended moist adiabatically from the level of the
/// lowest theta-e in the lowest 400 hPa of the sounding.
///
/// Returns the profile, downdraft cape, and downrush temperature in a tuple.
pub fn dcape(snd: &Sounding) -> Result<(ParcelProfile, JpKg, Celsius)> {
    let t = snd.temperature_profile();
    let dp = snd.dew_point_profile();
    let p = snd.pressure_profile();

    // Find the lowest pressure, 400 mb above the surface (or starting level)
    let top_p = p
        .iter()
        .filter_map(|p| p.into_option())
        .next()
        .ok_or(AnalysisError::NotEnoughData)?
        - HectoPascal(400.0);

    // Check that pressure is positive. This is theoretically possible if the surface pressure is
    // less than 400 hPa. Maybe in the Himalayas? Moisture is insignficant at 100 hPa anyway. This
    // is a really gross error check.
    debug_assert!(
        top_p.into_option().is_some() && top_p > HectoPascal(100.0),
        "surface pressure is below 500 mb"
    );

    // Find the starting parcel.
    let pcl = izip!(p, t, dp)
        .filter_map(|(p, t, dp)| {
            if p.is_some() && t.is_some() && dp.is_some() {
                Some((p.unpack(), t.unpack(), dp.unpack()))
            } else {
                None
            }
        })
        .take_while(|(p, _, _)| *p >= top_p)
        .fold(
            Err(AnalysisError::NotEnoughData),
            |acc: Result<Parcel>, (p, t, dp)| match metfor::theta_e(t, dp, p) {
                Some(th_e) => match acc {
                    Ok(parcel) => {
                        if let Ok(old_theta) = parcel.theta_e() {
                            if old_theta < th_e {
                                Ok(parcel)
                            } else {
                                Ok(Parcel {
                                    temperature: t,
                                    dew_point: dp,
                                    pressure: p,
                                })
                            }
                        } else {
                            Ok(Parcel {
                                temperature: t,
                                dew_point: dp,
                                pressure: p,
                            })
                        }
                    }
                    Err(_) => Ok(Parcel {
                        temperature: t,
                        dew_point: dp,
                        pressure: p,
                    }),
                },
                None => acc,
            },
        )?;

    let profile = descend_moist(pcl, snd)?;

    let mut dcape = 0.0;
    let mut h0 = Meters(std::f64::MAX); // Big number
    let mut pt0 = Kelvin(0.0);
    let mut et0 = Kelvin(0.0);
    for (&p, &h, &pt, &et) in izip!(
        &profile.pressure,
        &profile.height,
        &profile.parcel_t,
        &profile.environment_t
    ) {
        let pt = metfor::theta(p, pt);
        let et = metfor::theta(p, et);

        let dz = h - h0;
        // we must be starting out, becuase h0 starts as a large positive number
        if dz <= Meters(0.0) {
            h0 = h;
            pt0 = pt;
            et0 = et;
            continue;
        }

        dcape +=
            ((pt - et).unpack() / et.unpack() + (pt0 - et0).unpack() / et0.unpack()) * dz.unpack();

        h0 = h;
        pt0 = pt;
        et0 = et;
    }

    // - for integration direction should be top down, 9.81 for gravity, and 2.0 for trapezoid rule.
    dcape *= metfor::g / 2.0;

    let downrush_t = *profile.parcel_t.get(0).ok_or(AnalysisError::MissingValue)?;

    Ok((profile, JpKg(dcape), downrush_t))
}
