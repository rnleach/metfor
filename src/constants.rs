//! Meteorlogical constants.
//!
#![allow(non_upper_case_globals)]

/// Acceleration due to gravity at the Earth's surface. (m s<sup>-2</sup>)
pub const g: f64 = 9.81;

/// The gas constant for dry air. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const R: f64 = 287.04;

/// The gas constant for water vapor. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const Rv: f64 = 461.55;

/// Specific heat of dry air at constant pressure. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const cp: f64 = 1004.0;

/// Specific heat of dry air at constant volume. (J K<sup>-1</sup> kg<sup>-1</sup>)
pub const cv: f64 = 717.0;

/// Ratio of cp and cv. (unitless)
pub const gamma: f64 = cp / cv;

pub(crate) const ABSOLUTE_ZERO_K: f64 = 0.0;
pub(crate) const ABSOLUTE_ZERO_C: f64 = -273.15;
pub(crate) const ABSOLUTE_ZERO_F: f64 = -459.67;
