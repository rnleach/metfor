//! Helicity types
use crate::types::Quantity;
use std::cmp::Ordering;

/// Marker trait for helicity.
pub trait Helicity: Quantity + PartialEq + PartialOrd {}

/// Marker trait for vertically integrated helicity.
pub trait IntHelicity: Quantity + PartialOrd + PartialEq {}

/// Helicity in m s<sup>-2</sup> units
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct HelicityMpS2(pub f64);

/// Vertically integrated helicity in m<sup>2</sup> s<sup>-2</sup> units
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct IntHelicityM2pS2(pub f64);

impl Helicity for HelicityMpS2 {}
impl IntHelicity for IntHelicityM2pS2 {}

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

implQuantity!(HelicityMpS2);
implQuantity!(IntHelicityM2pS2);
