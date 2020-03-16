use super::Layer;
use crate::{parcel, parcel_profile, sounding::Sounding};
use metfor::JpKg;

/// Get the effective inflow layer.
pub fn effective_inflow_layer(snd: &Sounding) -> Option<Layer> {
    let mut vals_iter = snd
        .bottom_up()
        // Convert rows to parcels
        .filter_map(|row| parcel::Parcel::from_datarow(row).map(|pcl| (row, pcl)))
        // Lift the parcel, skip if there is an error
        .filter_map(|(row, pcl)| {
            parcel_profile::lift_parcel(pcl, snd)
                .ok()
                .map(|pcl_anal| (row, pcl_anal))
        })
        // Get the CAPE and CIN
        .filter_map(|(row, pcl_anal)| {
            pcl_anal
                .cape()
                .into_option()
                .and_then(|cape| pcl_anal.cin().map(|cin| (row, cape, cin)))
        })
        // Skip levels until we get one that meets criteria
        .skip_while(|(_, cape, cin)| *cape < JpKg(100.0) || *cin < JpKg(-250.0))
        // Take levels as long as they meet criteria
        .take_while(|(_, cape, cin)| *cape >= JpKg(100.0) && *cin >= JpKg(-250.0))
        // Discard the cape and cin values, we only need the rows
        .map(|(row, _, _)| row);

    let bottom = vals_iter.next();
    let top = vals_iter.last();

    if let (Some(bottom), Some(top)) = (bottom, top) {
        Some(Layer { bottom, top })
    } else {
        None
    }
}
