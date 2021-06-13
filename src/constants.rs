//! Meteorlogical constants.
//!
#![allow(non_upper_case_globals)]

use crate::types::*;

/// Acceleration due to gravity at the Earth's surface. (m s<sup>-2</sup>)
#[doc(hidden)]
pub const g: f64 = -9.806_65;

/// The gas constant for dry air. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const Rd: JpKgpK = JpKgpK(Rd_);
const Rd_: f64 = 287.058;

/// The gas constant for water vapor. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const Rv: JpKgpK = JpKgpK(Rv_);
const Rv_: f64 = 461.5;

/// Ratio of R / Rv. (no units)
pub const epsilon: f64 = Rd_ / Rv_;

/// Specific heat capacity of dry air at constant pressure. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const cpd: JpKgpK = JpKgpK(cpd_);
const cpd_: f64 = 1005.7;

/// Specific heat capacity of water vapor at constant pressure. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const cpv: JpKgpK = JpKgpK(cpv_);
const cpv_: f64 = 1870.0;

/// Specific heat capacity of liquid water. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const cl: JpKgpK = JpKgpK(cl_);
const cl_: f64 = 4190.0;

/// Specific heat capacity of dry air at constant volume. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const cvd: JpKgpK = JpKgpK(cvd_);
const cvd_: f64 = 718.0;

/// Ratio of cp and cv. (unitless)
pub const gamma: f64 = cpd_ / cvd_;

/// Absolute zero temperature.
pub const ABSOLUTE_ZERO: Kelvin = Kelvin(0.0);

/// Freezing.
pub const FREEZING: Celsius = Celsius(0.0);

/// Dry adiabatic lapse rate
pub const DRY_ADIABATIC_LAPSE_RATE: CelsiusPKm = CelsiusPKm(g);
