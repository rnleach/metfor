//! Specific energy units (energy per unit mass)
use crate::types::{
    temperatures::{CelsiusDiff, Kelvin, KelvinDiff},
    Quantity,
};
use std::cmp::Ordering;
use std::ops::Mul;

/// Marker trait for specific energy types.
pub trait SpecificEnergy: Quantity {}

/// Marker trait for specific energy per Kelvin types
pub trait SpecificEnergyPerKelvin: Quantity {}

/// Specific energy in J kg<sup>-1</sup> units. Used for CAPE, CIN, latent heat, etc.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct JpKg(pub f64);

/// Specific energy per Kelvin in J K<sup>-1</sup> kg<sup>-1</sup>. Used for gas constants and
/// specific heat valus.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
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

        implOpsForQuantity!($t);
    };
}

implQuantity!(JpKg);
implQuantity!(JpKgpK);

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

impl Mul<CelsiusDiff> for JpKgpK {
    type Output = JpKg;

    fn mul(self, rhs: CelsiusDiff) -> JpKg {
        let diff: CelsiusDiff = rhs.into();
        JpKg(self.0 * diff.unpack())
    }
}

impl Mul<KelvinDiff> for JpKgpK {
    type Output = JpKg;

    fn mul(self, rhs: KelvinDiff) -> JpKg {
        let diff: KelvinDiff = rhs.into();
        JpKg(self.0 * diff.unpack())
    }
}

impl Mul<JpKgpK> for CelsiusDiff {
    type Output = JpKg;

    fn mul(self, rhs: JpKgpK) -> JpKg {
        rhs * self
    }
}
