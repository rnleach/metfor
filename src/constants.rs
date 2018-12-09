//! Meteorlogical constants.
//!
#![allow(non_upper_case_globals)]

use crate::types::*;
/// Acceleration due to gravity at the Earth's surface. (m s<sup>-2</sup>)
pub const g: f64 = 9.81;

/// The gas constant for dry air. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const Rd: f64 = 287.058;

/// The gas constant for water vapor. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const Rv: f64 = 461.5;

/// Specific heat of dry air at constant pressure. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const cp: f64 = 1005.0;

/// Specific heat of dry air at constant volume. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const cv: f64 = 718.0;

/// Ratio of R / Rv. (no units)
pub const epsilon: f64 = Rd / Rv;

/// Ratio of cp and cv. (unitless)
pub const gamma: f64 = cp / cv;

/// Absolute zero temperature.
pub const ABSOLUTE_ZERO: Kelvin = Kelvin(ABSOLUTE_ZERO_K);

pub(crate) const ABSOLUTE_ZERO_K: f64 = 0.0;
pub(crate) const ABSOLUTE_ZERO_C: f64 = -273.15;
// pub(crate) const ABSOLUTE_ZERO_F: f64 = -459.67;
