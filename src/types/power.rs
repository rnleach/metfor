//! Power units. (Used for some specialized fire plumes equations)
use crate::types::Quantity;
use std::cmp::Ordering;

/// Marker trait for power types.
pub trait Power: Quantity {}

/// Power.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct GigaWatts(pub f64);

impl Power for GigaWatts {}

impl Quantity for GigaWatts {
    #[inline]
    fn pack(val: f64) -> Self {
        Self(val)
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

implOpsForQuantity!(GigaWatts);
