//! Specific energy units (energy per unit mass)
use crate::error::*;
use crate::types::Quantity;
use std::cmp::Ordering;
use std::fmt::Display;
use std::ops::{Add, Sub};

/// Marker trait for specific energy types.
pub trait SpecificEnergy: Quantity {}

/// Specific energy in J/Kg units. Used for CAPE, CIN, latent heat, etc.
#[derive(Clone, Copy, Debug)]
pub struct JpKg(pub f64);

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
            fn into_result(self) -> Result<f64> {
                Ok(self.0)
            }

            #[inline]
            fn borrow_inner(&self) -> &f64 {
                &self.0
            }
        }

        implOpsForQuantity!($t);
    };
}

implQuantity!(JpKg);

impl Display for JpKg {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{} J/Kg", self.0)
    }
}
