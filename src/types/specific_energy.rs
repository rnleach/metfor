//! Specific energy units (energy per unit mass)
use crate::types::{temperatures::Kelvin, Quantity};
use std::borrow::Borrow;
use std::cmp::Ordering;
use std::fmt::Display;
use std::ops::{Add, Div, Mul, Sub};

/// Marker trait for specific energy types.
pub trait SpecificEnergy: Quantity {}

/// Marker trait for specific energy per Kelvin types
pub trait SpecificEnergyPerKelvin: Quantity {}

/// Specific energy in J kg<sup>-1</sup> units. Used for CAPE, CIN, latent heat, etc.
#[derive(Clone, Copy, Debug)]
pub struct JpKg(pub f64);

/// Specific energy per Kelvin in J K<sup>-1</sup> kg<sup>-1</sup>. Used for gas constants and
/// specific heat valus.
#[derive(Clone, Copy, Debug)]
pub struct JpKgpK(pub f64);

impl SpecificEnergy for JpKg {}
macro_rules! implQuantity {
    ($t:tt) => {
        impl Quantity for $t {
            #[inline]
            fn pack(val: f64) -> Self {
                $t(val)
            }

            #[inline]
            fn unpack(self) -> f64 {
                self.0
            }

            #[inline]
            fn unwrap(self) -> f64 {
                self.0
            }

            #[inline]
            fn into_option(self) -> Option<f64> {
                Some(self.0)
            }
        }

        impl Borrow<f64> for $t {
            #[inline]
            fn borrow(&self) -> &f64 {
                &self.0
            }
        }

        implOpsForQuantity!($t);
    };
}

implQuantity!(JpKg);
implQuantity!(JpKgpK);

impl Display for JpKg {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{} J/kg", self.0)
    }
}

impl Display for JpKgpK {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{} J/kg/K", self.0)
    }
}

impl Div<JpKgpK> for JpKgpK {
    type Output = f64;

    fn div(self, rhs: JpKgpK) -> f64 {
        self.0 / rhs.0
    }
}

impl Mul<Kelvin> for JpKgpK {
    type Output = JpKg;

    fn mul(self, rhs: Kelvin) -> JpKg {
        JpKg(self.0 * rhs.unpack())
    }
}

impl Mul<JpKgpK> for Kelvin {
    type Output = JpKg;

    fn mul(self, rhs: JpKgpK) -> JpKg {
        rhs * self
    }
}
