//! Meteorlogical constants.
//!
#![allow(non_upper_case_globals)]

use crate::types::*;

/// Acceleration due to gravity at the Earth's surface. (m s<sup>-2</sup>)
#[doc(hidden)]
pub const g: f64 = 9.806_65;

/// The gas constant for dry air. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const Rd: JpKgpK = JpKgpK(Rd_);
const Rd_: f64 = 287.058;

/// The gas constant for water vapor. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const Rv: JpKgpK = JpKgpK(Rv_);
const Rv_: f64 = 461.5;

/// Specific heat of dry air at constant pressure. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const cp: JpKgpK = JpKgpK(cp_);
const cp_: f64 = 1005.0;

/// Specific heat of dry air at constant volume. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const cv: JpKgpK = JpKgpK(cv_);
const cv_: f64 = 718.0;

/// Ratio of R / Rv. (no units)
pub const epsilon: f64 = Rd_ / Rv_;

/// Ratio of cp and cv. (unitless)
pub const gamma: f64 = cp_ / cv_;

/// Absolute zero temperature.
pub const ABSOLUTE_ZERO: Kelvin = Kelvin(ABSOLUTE_ZERO_K);

/// Freezing.
pub const FREEZING: Celsius = Celsius(0.0);

/// Dry adiabatic lapse rate
pub const DRY_ADIABATIC_LAPSE_RATE: CelsiusPKm = CelsiusPKm(-g);

pub(crate) const ABSOLUTE_ZERO_K: f64 = 0.0;
pub(crate) const ABSOLUTE_ZERO_C: f64 = -273.15;
// pub(crate) const ABSOLUTE_ZERO_F: f64 = -459.67;
