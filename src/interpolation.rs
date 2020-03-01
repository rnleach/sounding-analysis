use crate::{
    error::{AnalysisError, Result},
    sounding::{DataRow, Sounding},
};
use itertools::{izip, Itertools};
use metfor::{HectoPascal, Knots, Quantity, WindSpdDir, WindUV};
use optional::Optioned;
use std::ops::Sub;

/// Interpolate values from the vertical sounding using pressure as the primary coordinate.
///
/// Returns a `DataRow` struct with interpolated values.
pub fn linear_interpolate_sounding(snd: &Sounding, tgt_p: HectoPascal) -> Result<DataRow> {
    let pressure: &[Optioned<HectoPascal>] = snd.pressure_profile();

    // What kind of bracket is this?
    enum BracketType {
        Bracket(usize, usize),
        EndEquals(usize),
    }

    // Map this pair of slice index and pressure points to a BracketType
    let make_bracket = |pnt_0, pnt_1| -> Option<BracketType> {
        let (i0, p0): (_, HectoPascal) = pnt_0;
        let (i1, p1): (_, HectoPascal) = pnt_1;

        // Always assume pressure is sorted in descending order
        debug_assert!(p0 > p1);
        if p0 > tgt_p && p1 < tgt_p {
            Some(BracketType::Bracket(i0, i1))
        } else if (p0 - tgt_p).unpack().abs() < std::f64::EPSILON {
            Some(BracketType::EndEquals(i0))
        } else if (p1 - tgt_p).unpack().abs() < std::f64::EPSILON {
            Some(BracketType::EndEquals(i1))
        } else {
            None
        }
    };

    // Find the levels to interpolate between.
    pressure
        .iter()
        .enumerate()
        // Remove levels with missing pressure (SHOULD be none...but...) and then unwrap from the
        // Optioned type
        .filter_map(|(i, p_val_opt)| p_val_opt.map(|p_val| (i, p_val)))
        // Look at the levels two at a time...
        .tuple_windows::<(_, _)>()
        // Map these pairs to brackets and remove anything that isn't a bracket. Should leave
        // at most one bracket in the iterator!
        .filter_map(|(pnt_0, pnt_1)| make_bracket(pnt_0, pnt_1))
        // Get the first (and only) bracket
        .nth(0) // Option<BracketType>
        // Perform the interpolation!
        .and_then(|bracket| match bracket {
            BracketType::Bracket(i0, i1) => {
                let row0 = snd.data_row(i0)?;
                let row1 = snd.data_row(i1)?;
                linear_interp_data_rows(row0, row1, tgt_p)
            }
            BracketType::EndEquals(i) => snd.data_row(i),
        })
        // Map to error
        .ok_or(AnalysisError::InterpolationError)
}

/// Interpolate values given two parallel vectors of data and a target value.
///
/// Assumes that xs is monotonic.
#[inline]
pub fn linear_interpolate<X, Y>(xs: &[Optioned<X>], ys: &[Optioned<Y>], target_x: X) -> Optioned<Y>
where
    X: Quantity + optional::Noned + PartialOrd + Sub<X>,
    <X as Sub<X>>::Output: Quantity + optional::Noned,
    Y: Quantity + optional::Noned + Sub<Y>,
    <Y as Sub<Y>>::Output: Quantity,
{
    debug_assert_eq!(xs.len(), ys.len());

    enum BracketType<X, Y> {
        Bracket((X, Y), (X, Y)),
        EndEqual((X, Y)),
    }

    let make_bracket = |pnt_0, pnt_1| -> Option<BracketType<X, Y>> {
        let (x0, _) = pnt_0;
        let (x1, _) = pnt_1;

        if (x0 < target_x && x1 > target_x) || (x0 > target_x && x1 < target_x) {
            Some(BracketType::Bracket(pnt_0, pnt_1))
        } else if (x0 - target_x).unpack().abs() < std::f64::EPSILON {
            Some(BracketType::EndEqual(pnt_0))
        } else if (x1 - target_x).unpack().abs() < std::f64::EPSILON {
            Some(BracketType::EndEqual(pnt_1))
        } else {
            None
        }
    };

    let value_opt = izip!(xs, ys)
        // Filter out elements where one of the values is missing, this allows us to skip over
        // over a point with a missing value and use the points on either side of it for the
        // interpolation.
        .filter(|(x, y)| x.is_some() && y.is_some())
        // Unpack the values from the `Optioned` type
        .map(|(x, y)| (x.unpack(), y.unpack()))
        // Look at them in pairs.
        .tuple_windows::<(_, _)>()
        // Make a bracket and filter out all levels the don't create a bracket.
        .filter_map(|(pnt_0, pnt_1)| make_bracket(pnt_0, pnt_1))
        // Get the first (and only) one that brackets the target value
        .nth(0) // This is an Option<BracketType>
        // Map from the bracket type to the interpolated value
        .map(|val| match val {
            BracketType::Bracket(pnt_0, pnt_1) => {
                let (x0, y0) = pnt_0;
                let (x1, y1) = pnt_1;
                linear_interp(target_x, x0, x1, y0, y1)
            }
            BracketType::EndEqual(pnt) => pnt.1,
        });

    Optioned::from(value_opt)
}

#[inline]
pub(crate) fn linear_interp<X, Y>(x_val: X, x1: X, x2: X, y1: Y, y2: Y) -> Y
where
    X: Sub<X> + Copy + std::fmt::Debug + std::cmp::PartialEq,
    <X as Sub<X>>::Output: Quantity,
    Y: Quantity + Sub<Y>,
    <Y as Sub<Y>>::Output: Quantity,
{
    debug_assert_ne!(x1, x2);

    let run = (x2 - x1).unpack();
    let rise = (y2 - y1).unpack();
    let dx = (x_val - x1).unpack();

    Y::pack(y1.unpack() + dx * (rise / run))
}

#[inline]
fn linear_interp_data_rows(row0: DataRow, row1: DataRow, tgt_p: HectoPascal) -> Option<DataRow> {
    let p0 = row0.pressure.into_option()?;
    let p1 = row1.pressure.into_option()?;

    let run = p1 - p0;
    let dp = tgt_p - p0;

    let mut result = DataRow::default();
    result.pressure = Optioned::from(tgt_p);

    result.temperature = eval_linear_interp(row0.temperature, row1.temperature, run, dp);
    result.wet_bulb = eval_linear_interp(row0.wet_bulb, row1.wet_bulb, run, dp);
    result.dew_point = eval_linear_interp(row0.dew_point, row1.dew_point, run, dp);
    result.theta_e = eval_linear_interp(row0.theta_e, row1.theta_e, run, dp);

    // Special interpolation for vectors
    if let (Some(w_below), Some(w_above)) = (row0.wind.into_option(), row1.wind.into_option()) {
        let WindUV::<Knots> {
            u: x_below,
            v: y_below,
        } = WindUV::from(w_below);
        let WindUV::<Knots> {
            u: x_above,
            v: y_above,
        } = WindUV::from(w_above);
        let dp = dp.unpack();
        let run = run.unpack();

        let rise_x = x_above - x_below;
        let rise_y = y_above - y_below;

        let x = x_below + rise_x * (dp / run);
        let y = y_below + rise_y * (dp / run);

        let interped_wind = WindSpdDir::from(WindUV { u: x, v: y });

        result.wind = interped_wind.into();
    }

    result.pvv = eval_linear_interp(row0.pvv, row1.pvv, run, dp);
    result.height = eval_linear_interp(row0.height, row1.height, run, dp);
    result.cloud_fraction = eval_linear_interp(row0.cloud_fraction, row1.cloud_fraction, run, dp);

    Some(result)
}

#[inline]
fn eval_linear_interp<QX, Y>(
    low_val: Optioned<Y>,
    high_val: Optioned<Y>,
    run: QX,
    dp: QX,
) -> Optioned<Y>
where
    QX: Quantity + optional::Noned,
    Y: Quantity + optional::Noned,
{
    if low_val.is_some() && high_val.is_some() {
        let (val_below, val_above) = (low_val.unpack().unpack(), high_val.unpack().unpack());
        let rise: f64 = (val_above - val_below).unpack();
        let run: f64 = run.unpack();
        let dp: f64 = dp.unpack();
        Optioned::from(Y::pack(val_below + dp * rise / run))
    } else {
        Optioned::default()
    }
}
