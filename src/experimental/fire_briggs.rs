//! Experimental sounding analysis for fire plumes utilizing the Briggs plume model for the
//! subcloud plume.

use crate::{
    error::{AnalysisError, Result},
    experimental::fire::{analyze_plume_parcel, PlumeAscentAnalysis},
    fire::{
        entrained_layer_mean_wind_speed, entrained_mixed_layer,
        potential_t_and_specific_humidity_to_pressure_and_temperature,
    },
    parcel::Parcel,
    sounding::Sounding,
};
use itertools::izip;
use metfor::{self, Celsius, GigaWatts, JpKg, Kelvin, Meters, Quantity};
use optional::Optioned;

/// Various analysis results of lifting plumes parcels vs the heating supplied. The subcloud
/// portion of the plume is modeled with the Briggs plume model.
#[derive(Debug, Clone)]
pub struct BriggsPlumeHeatingAnalysis {
    /// The Briggs Plume mixed subcloud layers average parcel.
    pub starting_parcel: Parcel,
    /// How much heating does the fire need to produce for these results.
    pub fire_power: Vec<GigaWatts>,
    /// The parcels lifted for each step.
    pub fp_parcels: Vec<Parcel>,
    /// Basically CAPE, but it takes the heating of the fire into account as well.
    pub max_int_buoyancies: Vec<Optioned<JpKg>>,
    /// The ratio of the CAPE that was from latent heat release.
    pub wet_ratio: Vec<Optioned<f64>>,
    /// LCL - cloud base.
    pub lcl_heights: Vec<Optioned<Meters>>,
    /// Equiplibrium level (leve of maximum integrated buoyancy.
    pub el_heights: Vec<Optioned<Meters>>,
    /// If we keep integrating above the EL, how far up before the buoyancy gets to zero. The max
    /// height minus the EL would be the depth of an overshooting top.
    pub max_heights: Vec<Optioned<Meters>>,
}

/// Do a PlumeHeatingAnalysis.
pub fn briggs_plume_heating_analysis(
    snd: &Sounding,
    moisture_ratio: Option<f64>,
) -> Result<BriggsPlumeHeatingAnalysis> {
    let (starting_parcel, anal_iter) = plumes_heating_iter(snd, moisture_ratio)?;

    let mut fire_power: Vec<GigaWatts> = vec![];
    let mut fp_parcels: Vec<Parcel> = vec![];
    let mut max_int_buoyancies: Vec<Optioned<JpKg>> = vec![];
    let mut wet_ratio: Vec<Optioned<f64>> = vec![];
    let mut lcl_heights: Vec<Optioned<Meters>> = vec![];
    let mut el_heights: Vec<Optioned<Meters>> = vec![];
    let mut max_heights: Vec<Optioned<Meters>> = vec![];

    anal_iter.for_each(|(fp, anal)| {
        fire_power.push(fp);
        fp_parcels.push(anal.parcel);
        max_int_buoyancies.push(anal.max_int_buoyancy);

        let a_wet_ratio = anal.max_dry_int_buoyancy.and_then(|dry| {
            anal.max_int_buoyancy.map_t(|total| {
                if total > JpKg(0.0) {
                    (total - dry) / total
                } else {
                    0.0
                }
            })
        });
        wet_ratio.push(a_wet_ratio);
        lcl_heights.push(anal.lcl_height);
        el_heights.push(anal.el_height);
        max_heights.push(anal.max_height);
    });

    Ok(BriggsPlumeHeatingAnalysis {
        starting_parcel,
        fire_power,
        fp_parcels,
        max_int_buoyancies,
        wet_ratio,
        lcl_heights,
        el_heights,
        max_heights,
    })
}

fn plumes_heating_iter(
    snd: &Sounding,
    moisture_ratio: Option<f64>,
) -> Result<(
    Parcel,
    impl Iterator<Item = (GigaWatts, PlumeAscentAnalysis)> + '_,
)> {
    const INCREMENT: f64 = 0.000_1;
    const MAX_BETA: f64 = 0.20;

    let (starting_parcel, parcel_iter) = plume_parcels(snd, MAX_BETA, INCREMENT, moisture_ratio)?;

    let anal_iter = parcel_iter
        // Do the analysis, ignore errors.
        .filter_map(move |(fp, pcl)| analyze_plume_parcel(pcl, snd).ok().map(|anal| (fp, anal)))
        // Skip levels with no useful information.
        .skip_while(|(_fp, anal)| {
            anal.el_height.is_none() && anal.max_height.is_none() && anal.lcl_height.is_none()
        })
        // Take while there is some useful information.
        .take_while(|(_fp, anal)| {
            anal.el_height.is_some() || anal.max_height.is_some() || anal.lcl_height.is_some()
        });

    Ok((starting_parcel, anal_iter))
}

/// Given a sounding, return an iterator that creates parcels starting with the mixed layer parcel
/// and then incrementing along the SP curve.
///
/// Arguments:
/// * snd - the environmental sounding.
/// * max_beta - the maximum beta along the SP curve
/// * increment - the difference in beta along hte SP curve for each parcel.
/// * moisture_ratio - a value of 10 means for every 10C of heating (dt), add 1 g/kg of moisture
///   to the parcel. If it is `None`, don't add any moisture.
///
fn plume_parcels(
    snd: &Sounding,
    max_beta: f64,
    increment: f64,
    moisture_ratio: Option<f64>,
) -> Result<(Parcel, impl Iterator<Item = (GigaWatts, Parcel)> + use<'_>)> {
    let (starting_theta, starting_sh, _z_ml, bottom_p, _p_ml) = entrained_mixed_layer(snd)?;
    let starting_t = Celsius::from(metfor::temperature_from_pot_temp(starting_theta, bottom_p));
    let starting_dp = Celsius::from(
        metfor::dew_point_from_p_and_specific_humidity(bottom_p, starting_sh)
            .ok_or(AnalysisError::InvalidInput)?,
    );

    let parcel = Parcel {
        temperature: starting_t,
        pressure: bottom_p,
        dew_point: starting_dp,
    };

    Ok((
        parcel,
        PlumeParcelIterator {
            starting_pcl: parcel,
            next_beta: 0.0f64 - increment,
            max_beta,
            increment,
            moisture_ratio,
            snd,
        },
    ))
}

/// Iterator for `plume_parcels` function that generates increasingly warmer parcels along an SP
/// curve.
struct PlumeParcelIterator<'a> {
    starting_pcl: Parcel,
    next_beta: f64,
    max_beta: f64,
    increment: f64,
    moisture_ratio: Option<f64>,
    snd: &'a Sounding,
}

impl Iterator for PlumeParcelIterator<'_> {
    type Item = (GigaWatts, Parcel);

    fn next(&mut self) -> Option<Self::Item> {
        let (sfc_height, p_sfc) =
            match izip!(self.snd.height_profile(), self.snd.pressure_profile())
                .filter(|(h, p)| h.is_some() && p.is_some())
                .map(|(h, p)| (h.unpack(), p.unpack()))
                .next()
                .ok_or(AnalysisError::NotEnoughData)
            {
                Ok((h, p)) => (h, p),
                Err(_) => return None,
            };

        loop {
            self.next_beta += self.increment;
            if self.next_beta > self.max_beta {
                break;
            }

            // Calculate the position along the SP curve.
            let theta = metfor::potential_temperature(
                self.starting_pcl.pressure,
                self.starting_pcl.temperature,
            );
            let sh = match metfor::specific_humidity(
                self.starting_pcl.dew_point,
                self.starting_pcl.pressure,
            ) {
                Some(sh) => sh,
                None => continue,
            };

            let theta_fc = Kelvin((1.0 + self.next_beta) * theta.unpack());
            let sh_fc =
                sh + self.next_beta / self.moisture_ratio.unwrap_or(1.0) / 1000.0 * theta.unpack();

            let t_parcel = metfor::temperature_from_pot_temp(theta_fc, self.starting_pcl.pressure);
            let dp_parcel = match metfor::dew_point_from_p_and_specific_humidity(
                self.starting_pcl.pressure,
                sh_fc,
            ) {
                Some(dp) => dp,
                None => continue,
            };

            let (p_lfc, _t_lfc) =
                match potential_t_and_specific_humidity_to_pressure_and_temperature(theta_fc, sh_fc)
                {
                    Ok((p, t)) => (p, t),
                    Err(_) => continue,
                };

            let dtheta = theta_fc - theta;

            let heights = self.snd.height_profile();
            let pressures = self.snd.pressure_profile();
            let height_asl_lfc =
                match crate::interpolation::linear_interpolate(pressures, heights, p_lfc)
                    .into_option()
                {
                    Some(h) => h,
                    None => continue,
                };

            let z_fc = height_asl_lfc - sfc_height;

            let u = match entrained_layer_mean_wind_speed(z_fc, self.snd) {
                Ok(u) => u,
                Err(_) => continue,
            };

            let fp = metfor::pft(z_fc, p_lfc, u, dtheta, theta_fc, p_sfc);

            return Some((
                fp,
                Parcel {
                    temperature: t_parcel.into(),
                    dew_point: dp_parcel.into(),
                    ..self.starting_pcl
                },
            ));
        }

        None
    }
}
