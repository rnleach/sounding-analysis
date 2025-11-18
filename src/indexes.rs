//! Indexes that are specific to a sounding, but not a particular parcel analysis of that sounding.

use crate::{
    error::Result,
    sounding::Sounding,
};
use itertools::{izip, Itertools};
use metfor::{mixing_ratio, Mm, Quantity};

/// Precipitable water (mm)
#[inline]
pub fn precipitable_water(snd: &Sounding) -> Result<Mm> {
    let p_profile = snd.pressure_profile();
    let dp_profile = snd.dew_point_profile();

    let integrated_mw = izip!(p_profile, dp_profile)
        // Remove levels with missing data
        .filter(|(p, dp)| p.is_some() && dp.is_some())
        // Unpack from the Optioned type
        .map(|(p, dp)| (p.unpack(), dp.unpack()))
        // Converte dew point to mixing ratio, removing failed levels.
        .filter_map(|(p, dp)| mixing_ratio(dp, p).map(|mw| (p, mw)))
        // View them as pairs for integration with the trapezoid method
        .tuple_windows::<(_, _)>()
        // Do the sum for integrating
        .fold(0.0, |mut acc_mw, ((p0, mw0), (p1, mw1))| {
            let dp = p0 - p1;
            acc_mw += (mw0 + mw1) * dp.unpack();

            acc_mw
        });

    Ok(Mm(integrated_mw / 9.81 / 997.0 * 100_000.0 / 2.0))
}

