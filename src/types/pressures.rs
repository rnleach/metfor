//! Pressure units
use crate::types::Quantity;
use std::cmp::Ordering;

/// Marker trait for Pressure types.
pub trait Pressure: Quantity {}

/// Pressure in hPa units.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct HectoPascal(pub f64);

/// Pressure in Pa units.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Pascal(pub f64);

/// Pressure in mb units.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Millibar(pub f64);

impl Pressure for HectoPascal {}
impl Pressure for Pascal {}
impl Pressure for Millibar {}

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
                if self.0 < 0.0 {
                    panic!("Negative Pressure");
                }

                self.0
            }

            #[inline]
            fn into_option(self) -> Option<f64> {
                if self.0 < 0.0 {
                    None
                } else {
                    Some(self.0)
                }
            }
        }

        impl optional::Noned for $t {
            fn is_none(&self) -> bool {
                self.0.is_none()
            }

            fn get_none() -> Self {
                $t(std::f64::NAN)
            }
        }

        implOpsForQuantity!($t);
    };
}

implQuantity!(HectoPascal);
implQuantity!(Pascal);
implQuantity!(Millibar);

double_conversion!(Millibar, HectoPascal, 1.0, 0.0, 1.0);
double_conversion!(HectoPascal, Pascal, 100.0, 0.0, 1.0);
double_conversion!(Millibar, Pascal, 100.0, 0.0, 1.0);
