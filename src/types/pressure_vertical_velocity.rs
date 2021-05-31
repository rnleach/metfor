//! Pressure vertical velocity.
use crate::types::Quantity;
use std::cmp::Ordering;

/// Marker trait for pressure vertical veclocity types.
pub trait PVV: Quantity + PartialEq + PartialOrd {}

/// Pressure vertical velocity in Pa/s
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct PaPS(pub f64);

/// Pressure vertical velocity in microbar/s
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct MicroBarPS(pub f64);

impl PVV for PaPS {}
impl PVV for MicroBarPS {}

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

implQuantity!(PaPS);
implQuantity!(MicroBarPS);

impl From<MicroBarPS> for PaPS {
    #[inline]
    fn from(p: MicroBarPS) -> Self {
        PaPS(p.0 / 10.0)
    }
}

impl From<PaPS> for MicroBarPS {
    #[inline]
    fn from(p: PaPS) -> Self {
        MicroBarPS(p.0 * 10.0)
    }
}
