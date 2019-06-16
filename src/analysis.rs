//! Data type and methods for building and describing an analysis.
//!
//! Not every possible analysis is in this data.
use crate::{
    indexes::{
        haines, haines_high, haines_low, haines_mid, hot_dry_windy, kindex, precipitable_water,
        swet, total_totals,
    },
    layers::{effective_inflow_layer, Layer},
    parcel::{mixed_layer_parcel, most_unstable_parcel, surface_parcel},
    parcel_profile::{
        dcape, lift_parcel, partition_cape, robust_convective_parcel, ParcelAnalysis, ParcelProfile,
    },
    wind::{self, bunkers_storm_motion, mean_wind},
};
use metfor::{
    Celsius, CelsiusDiff, IntHelicityM2pS2, JpKg, Length, Meters, MetersPSec, Mm, Quantity, WindUV,
};
use optional::{none, some, Noned, Optioned};
use sounding_base::Sounding;
use std::collections::HashMap;

/// Convenient package for commonly requested analysis values.
///
/// All parcel related values are assumed to be for the 100hPa mixed layer at the surface.
#[derive(Debug, Clone)]
pub struct Analysis {
    // Sounding used to make the analysis
    sounding: Sounding,

    // Profile specific indicies
    swet: Optioned<f64>,
    k_index: Optioned<Celsius>,
    total_totals: Optioned<f64>,
    precipitable_water: Optioned<Mm>,
    convective_t: Optioned<Celsius>,
    right_mover: Optioned<WindUV<MetersPSec>>,
    left_mover: Optioned<WindUV<MetersPSec>>,
    mean_wind: Optioned<WindUV<MetersPSec>>,
    sr_helicity_3k_rm: Optioned<IntHelicityM2pS2>,
    sr_helicity_3k_lm: Optioned<IntHelicityM2pS2>,
    effective_inflow_layer: Option<Layer>,
    sr_helicity_eff_rm: Optioned<IntHelicityM2pS2>,
    sr_helicity_eff_lm: Optioned<IntHelicityM2pS2>,

    // Fire weather indicies
    haines: Optioned<u8>,
    haines_low: Optioned<u8>,
    haines_mid: Optioned<u8>,
    haines_high: Optioned<u8>,
    hdw: Optioned<f64>,
    convective_deficit: Optioned<CelsiusDiff>,
    cape_ratio: Optioned<f64>,

    // Downburst
    dcape: Optioned<JpKg>,
    downrush_t: Optioned<Celsius>,
    downburst_profile: Option<ParcelProfile>,

    // Parcel analysis
    mixed_layer: Option<ParcelAnalysis>,
    surface: Option<ParcelAnalysis>,
    most_unstable: Option<ParcelAnalysis>,
    convective: Option<ParcelAnalysis>,

    // Provider analysis
    provider_analysis: HashMap<&'static str, f64>,
}

impl Analysis {
    /// Create a new `Analysis`.
    pub fn new(snd: Sounding) -> Self {
        Analysis {
            sounding: snd,
            swet: none(),
            k_index: none(),
            total_totals: none(),
            precipitable_water: none(),
            convective_t: none(),
            right_mover: none(),
            left_mover: none(),
            mean_wind: none(),
            sr_helicity_3k_rm: none(),
            sr_helicity_3k_lm: none(),
            effective_inflow_layer: None,
            sr_helicity_eff_rm: none(),
            sr_helicity_eff_lm: none(),

            haines: none(),
            haines_low: none(),
            haines_mid: none(),
            haines_high: none(),
            hdw: none(),
            convective_deficit: none(),
            cape_ratio: none(),

            dcape: none(),
            downrush_t: none(),
            downburst_profile: None,

            mixed_layer: None,
            surface: None,
            most_unstable: None,
            convective: None,

            provider_analysis: HashMap::new(),
        }
    }

    /// Builder method to set the SWeT
    pub fn with_swet<T>(self, value: T) -> Self
    where
        Optioned<f64>: From<T>,
    {
        Self {
            swet: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the K index
    pub fn with_k_index<T>(self, value: T) -> Self
    where
        Optioned<Celsius>: From<T>,
    {
        Self {
            k_index: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the Total Totals
    pub fn with_total_totals<T>(self, value: T) -> Self
    where
        Optioned<f64>: From<T>,
    {
        Self {
            total_totals: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the precipitable water
    pub fn with_pwat<T, U>(self, value: T) -> Self
    where
        Optioned<U>: From<T>,
        U: Length + Noned,
        Mm: From<U>,
    {
        let u: Optioned<U> = Optioned::from(value);
        let pw: Optioned<Mm> = u.map_t(Mm::from);
        Self {
            precipitable_water: pw,
            ..self
        }
    }

    /// Builder method to set the convective temperature
    pub fn with_convective_t<T>(self, value: T) -> Self
    where
        Optioned<Celsius>: From<T>,
    {
        Self {
            convective_t: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to add right mover storm motion
    pub fn with_right_mover<T>(self, value: T) -> Self
    where
        Optioned<WindUV<MetersPSec>>: From<T>,
    {
        Self {
            right_mover: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to add left mover storm motion
    pub fn with_left_mover<T>(self, value: T) -> Self
    where
        Optioned<WindUV<MetersPSec>>: From<T>,
    {
        Self {
            left_mover: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to add the mean wind
    pub fn with_mean_wind<T>(self, value: T) -> Self
    where
        Optioned<WindUV<MetersPSec>>: From<T>,
    {
        Self {
            mean_wind: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to add the storm relative helicity
    pub fn with_sr_helicity_3k_rm<T>(self, value: T) -> Self
    where
        Optioned<IntHelicityM2pS2>: From<T>,
    {
        Self {
            sr_helicity_3k_rm: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to add the storm relative helicity
    pub fn with_sr_helicity_3k_lm<T>(self, value: T) -> Self
    where
        Optioned<IntHelicityM2pS2>: From<T>,
    {
        Self {
            sr_helicity_3k_lm: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to add an effective inflow layer
    pub fn with_effiective_inflow_layer<T>(self, value: T) -> Self
    where
        Option<Layer>: From<T>,
    {
        Self {
            effective_inflow_layer: Option::from(value),
            ..self
        }
    }

    /// Builder method to add the effective storm relative helicity
    pub fn with_sr_helicity_eff_rm<T>(self, value: T) -> Self
    where
        Optioned<IntHelicityM2pS2>: From<T>,
    {
        Self {
            sr_helicity_eff_rm: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to add the effective storm relative helicity
    pub fn with_sr_helicity_eff_lm<T>(self, value: T) -> Self
    where
        Optioned<IntHelicityM2pS2>: From<T>,
    {
        Self {
            sr_helicity_eff_lm: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the haines index. High, mid, or low Haines should be based on the
    /// station elevation of the sounding.
    pub fn with_haines<T>(self, value: T) -> Self
    where
        Optioned<u8>: From<T>,
    {
        Self {
            haines: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the low level haines index.
    pub fn with_haines_low<T>(self, value: T) -> Self
    where
        Optioned<u8>: From<T>,
    {
        Self {
            haines_low: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the mid level haines index.
    pub fn with_haines_mid<T>(self, value: T) -> Self
    where
        Optioned<u8>: From<T>,
    {
        Self {
            haines_mid: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the high level haines index.
    pub fn with_haines_high<T>(self, value: T) -> Self
    where
        Optioned<u8>: From<T>,
    {
        Self {
            haines_high: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the hot-dry-windy index.
    pub fn with_hdw<T>(self, value: T) -> Self
    where
        Optioned<f64>: From<T>,
    {
        Self {
            hdw: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the wet/dry cape ratio. EXPERIMENTAL
    pub fn with_cape_ratio<T>(self, value: T) -> Self
    where
        Optioned<f64>: From<T>,
    {
        Self {
            cape_ratio: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the convective temperature deficit. EXPERIMENTAL
    pub fn with_convective_deficit<T>(self, value: T) -> Self
    where
        Optioned<CelsiusDiff>: From<T>,
    {
        Self {
            convective_deficit: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the DCAPE.
    pub fn with_dcape<T>(self, value: T) -> Self
    where
        Optioned<JpKg>: From<T>,
    {
        Self {
            dcape: Optioned::from(value),
            ..self
        }
    }

    /// Builder method to set the downrush temperature of a wet micro-burst.
    pub fn with_downrush_t<T>(self, value: T) -> Self
    where
        Optioned<Celsius>: From<T>,
    {
        Self {
            downrush_t: Optioned::from(value),
            ..self
        }
    }

    /// Get the Swet
    pub fn swet(&self) -> Optioned<f64> {
        self.swet
    }

    /// Get the K index
    pub fn k_index(&self) -> Optioned<Celsius> {
        self.k_index
    }

    /// Get the Total Totals index
    pub fn total_totals(&self) -> Optioned<f64> {
        self.total_totals
    }

    /// Get the precipitable water.
    pub fn pwat(&self) -> Optioned<Mm> {
        self.precipitable_water
    }

    /// Get the convective temperature.
    pub fn convective_t(&self) -> Optioned<Celsius> {
        self.convective_t
    }

    /// Get the right mover.
    pub fn right_mover(&self) -> Optioned<WindUV<MetersPSec>> {
        self.right_mover
    }

    /// Get the left mover.
    pub fn left_mover(&self) -> Optioned<WindUV<MetersPSec>> {
        self.left_mover
    }

    /// Get the mean wind.
    pub fn mean_wind(&self) -> Optioned<WindUV<MetersPSec>> {
        self.mean_wind
    }

    /// Get the storm relative helicity for a right mover storm
    pub fn sr_helicity_3k_rm(&self) -> Optioned<IntHelicityM2pS2> {
        self.sr_helicity_3k_rm
    }

    /// Get the storm relative helicity for a left mover storm
    pub fn sr_helicity_3k_lm(&self) -> Optioned<IntHelicityM2pS2> {
        self.sr_helicity_3k_lm
    }

    /// Get the effective inflow layer
    pub fn effective_inflow_layer(&self) -> Option<Layer> {
        self.effective_inflow_layer
    }

    /// Get the effective storm relative helicity for a right mover storm
    pub fn sr_helicity_eff_rm(&self) -> Optioned<IntHelicityM2pS2> {
        self.sr_helicity_eff_rm
    }

    /// Get the effective storm relative helicity for a left mover storm
    pub fn sr_helicity_eff_lm(&self) -> Optioned<IntHelicityM2pS2> {
        self.sr_helicity_eff_lm
    }

    /// Get the downrush temperature from a microburst.
    pub fn downrush_t(&self) -> Optioned<Celsius> {
        self.downrush_t
    }

    /// Get the Haines Index.
    pub fn haines(&self) -> Optioned<u8> {
        self.haines
    }

    /// Get the low level Haines Index.
    pub fn haines_low(&self) -> Optioned<u8> {
        self.haines_low
    }

    /// Get the mid level Haines Index.
    pub fn haines_mid(&self) -> Optioned<u8> {
        self.haines_mid
    }

    /// Get the high level Haines Index.
    pub fn haines_high(&self) -> Optioned<u8> {
        self.haines_high
    }

    /// Get the hot-dry-windy index.
    pub fn hdw(&self) -> Optioned<f64> {
        self.hdw
    }

    /// Get the wet/dry CAPE ratio. EXPERIMENTAL.
    pub fn cape_ratio(&self) -> Optioned<f64> {
        self.cape_ratio
    }

    /// Get the convective temperature deficit.
    pub fn convective_deficit(&self) -> Optioned<CelsiusDiff> {
        self.convective_deficit
    }

    /// Get the DCAPE.
    pub fn dcape(&self) -> Optioned<JpKg> {
        self.dcape
    }

    /// Set the mixed layer parcel analysis.
    pub fn with_mixed_layer_parcel_analysis<T>(self, anal: T) -> Self
    where
        Option<ParcelAnalysis>: From<T>,
    {
        let mixed_layer = Option::from(anal);
        Analysis {
            mixed_layer,
            ..self
        }
    }

    /// Get the mixed layer parcel analysis
    pub fn mixed_layer_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
        self.mixed_layer.as_ref()
    }

    /// Set the surface parcel analysis.
    pub fn with_surface_parcel_analysis<T>(self, anal: T) -> Self
    where
        Option<ParcelAnalysis>: From<T>,
    {
        let surface = Option::from(anal);
        Analysis { surface, ..self }
    }

    /// Get the surface parcel analysis
    pub fn surface_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
        self.surface.as_ref()
    }

    /// Set the most unstable parcel analysis.
    pub fn with_most_unstable_parcel_analysis<T>(self, anal: T) -> Self
    where
        Option<ParcelAnalysis>: From<T>,
    {
        let most_unstable = Option::from(anal);
        Analysis {
            most_unstable,
            ..self
        }
    }

    /// Get the most unstable parcel analysis
    pub fn most_unstable_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
        self.most_unstable.as_ref()
    }

    /// Set the convective parcel analysis.
    pub fn with_convective_parcel_analysis<T>(self, anal: T) -> Self
    where
        Option<ParcelAnalysis>: From<T>,
    {
        let convective = Option::from(anal);
        Analysis { convective, ..self }
    }

    /// Get the convective parcel analysis
    pub fn convective_parcel_analysis(&self) -> Option<&ParcelAnalysis> {
        self.convective.as_ref()
    }

    /// Set the downburst profile
    pub fn with_downburst_profile<T>(self, parcel_profile: T) -> Self
    where
        Option<ParcelProfile>: From<T>,
    {
        let downburst_profile = Option::from(parcel_profile);
        Analysis {
            downburst_profile,
            ..self
        }
    }

    /// Get the downburst profile
    pub fn downburst_profile(&self) -> Option<&ParcelProfile> {
        self.downburst_profile.as_ref()
    }

    /// Set the provider analysis.
    ///
    /// This is just a table of what ever values you want to store, it may be empty.
    pub fn with_provider_analysis(self, provider_analysis: HashMap<&'static str, f64>) -> Self {
        Analysis {
            provider_analysis,
            ..self
        }
    }

    /// Get a reference to the provider analysis so you can query it.
    pub fn provider_analysis(&self) -> &HashMap<&'static str, f64> {
        &self.provider_analysis
    }

    /// Get a mutable reference to the provider analysis so you can modify it.
    pub fn provider_analysis_mut(&mut self) -> &mut HashMap<&'static str, f64> {
        &mut self.provider_analysis
    }

    /// Get a reference to the sounding.
    pub fn sounding(&self) -> &Sounding {
        &self.sounding
    }

    /// Analyze the sounding to get as much information as you can.
    pub fn fill_in_missing_analysis_mut(&mut self) {
        self.swet = self
            .swet
            .or_else(|| Optioned::from(swet(&self.sounding).ok()));
        self.total_totals = self
            .total_totals
            .or_else(|| Optioned::from(total_totals(&self.sounding).ok()));
        self.k_index = self
            .k_index
            .or_else(|| Optioned::from(kindex(&self.sounding).ok()));
        self.precipitable_water = self
            .precipitable_water
            .or_else(|| Optioned::from(precipitable_water(&self.sounding).ok()));

        self.haines = self
            .haines
            .or_else(|| Optioned::from(haines(&self.sounding).ok()));
        self.haines_low = self
            .haines_low
            .or_else(|| Optioned::from(haines_low(&self.sounding).ok()));
        self.haines_mid = self
            .haines_mid
            .or_else(|| Optioned::from(haines_mid(&self.sounding).ok()));
        self.haines_high = self
            .haines_high
            .or_else(|| Optioned::from(haines_high(&self.sounding).ok()));
        self.hdw = self
            .hdw
            .or_else(|| Optioned::from(hot_dry_windy(&self.sounding).ok()));

        if self.dcape.is_none() || self.downrush_t.is_none() || self.downburst_profile.is_none() {
            let result = dcape(&self.sounding);

            if let Ok((pp, dcape, down_t)) = result {
                self.dcape = some(dcape);
                self.downrush_t = some(down_t);
                self.downburst_profile = Some(pp);
            }
        }

        if self.mixed_layer.is_none() {
            self.mixed_layer = match mixed_layer_parcel(&self.sounding) {
                Ok(parcel) => lift_parcel(parcel, &self.sounding).ok(),
                Err(_) => None,
            };
        }
        if self.most_unstable.is_none() {
            self.most_unstable = match most_unstable_parcel(&self.sounding) {
                Ok(parcel) => lift_parcel(parcel, &self.sounding).ok(),
                Err(_) => None,
            };
        }
        if self.surface.is_none() {
            self.surface = match surface_parcel(&self.sounding) {
                Ok(parcel) => lift_parcel(parcel, &self.sounding).ok(),
                Err(_) => None,
            };
        }
        if self.convective.is_none() {
            self.convective = robust_convective_parcel(&self.sounding).ok();
        }

        // Convective T
        if self.convective_t.is_none() {
            self.convective_t = self
                .convective
                .as_ref()
                .map(|parcel_anal| parcel_anal.parcel().temperature)
                .into();
        }

        // Left and right mover storm motion
        if self.right_mover.is_none() || self.left_mover.is_none() {
            let (rm, lm) = match bunkers_storm_motion(&self.sounding) {
                Ok((rm, lm)) => (some(rm), some(lm)),
                Err(_) => (none(), none()),
            };

            self.right_mover = rm;
            self.left_mover = lm;
        }

        // Fill in the mean wind
        if self.mean_wind.is_none() {
            if let Some(layer) = &crate::layers::layer_agl(&self.sounding, Meters(6000.0)).ok() {
                self.mean_wind = Optioned::from(mean_wind(layer, &self.sounding).ok());
            }
        }

        // Fill in the storm relative helicity
        if self.sr_helicity_3k_rm.is_none() || self.sr_helicity_3k_lm.is_none() {
            if let (Some(layer), Some(sm), Some(lm)) = (
                &crate::layers::layer_agl(&self.sounding, Meters(3000.0)).ok(),
                self.right_mover.into_option(),
                self.left_mover.into_option(),
            ) {
                self.sr_helicity_3k_rm =
                    Optioned::from(wind::sr_helicity(layer, sm, &self.sounding()).ok());

                self.sr_helicity_3k_lm =
                    Optioned::from(wind::sr_helicity(layer, lm, &self.sounding()).ok());
            }
        }

        // Fill in the effective inflow layer
        if self.effective_inflow_layer.is_none() {
            self.effective_inflow_layer = effective_inflow_layer(&self.sounding());
        }

        // Fill in the effective storm relative helicity
        if self.sr_helicity_eff_rm.is_none() || self.sr_helicity_eff_lm.is_none() {
            if let (Some(layer), Some(sm), Some(lm)) = (
                &self.effective_inflow_layer,
                self.right_mover.into_option(),
                self.left_mover.into_option(),
            ) {
                self.sr_helicity_eff_rm =
                    Optioned::from(wind::sr_helicity(layer, sm, &self.sounding()).ok());

                self.sr_helicity_eff_lm =
                    Optioned::from(wind::sr_helicity(layer, lm, &self.sounding()).ok());
            }
        }

        // Convective deficit
        if self.convective_deficit.is_none() {
            self.convective_deficit = self.convective_t.and_then(|ct| {
                self.mixed_layer
                    .as_ref()
                    .map(|parcel_anal| ct - parcel_anal.parcel().temperature)
                    .into()
            });
        }

        // Cape ratio
        if self.cape_ratio.is_none() {
            self.cape_ratio = self
                .convective
                .as_ref()
                .and_then(|parcel_anal| partition_cape(parcel_anal).ok())
                .and_then(|(dry, wet)| {
                    if dry.unpack().abs() > 0.1 {
                        Some(wet / dry)
                    } else {
                        None
                    }
                })
                .into();
        }
    }

    /// Analyze the sounding to get as much information as you can.
    pub fn fill_in_missing_analysis(mut self) -> Self {
        self.fill_in_missing_analysis_mut();
        self
    }
}
