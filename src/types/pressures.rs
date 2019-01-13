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

/// Pressure in mb units.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Millibar(pub f64);

impl Pressure for HectoPascal {}
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

        implOpsForQuantity!($t);
    };
}

implQuantity!(HectoPascal);
implQuantity!(Millibar);

impl From<Millibar> for HectoPascal {
    #[inline]
    fn from(mb: Millibar) -> Self {
        HectoPascal(mb.0)
    }
}

impl From<HectoPascal> for Millibar {
    #[inline]
    fn from(hpa: HectoPascal) -> Self {
        Millibar(hpa.0)
    }
}
