//! Length units for elevation and geopotential height.
use crate::types::Quantity;
use std::cmp::Ordering;

/// Marker trait for speed types.
pub trait Speed: Quantity + PartialEq + PartialOrd {}

/// Speed in meters per second.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct MetersPSec(pub f64);

/// Speed in knots.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Knots(pub f64);

/// Speed in miles per hour
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct MilesPHour(pub f64);

/// Speed in Kilometers per hour
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct KmPHour(pub f64);

impl Speed for MetersPSec {}
impl Speed for Knots {}
impl Speed for MilesPHour {}
impl Speed for KmPHour {}

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
                // allow negative speeds for use as an element in vectors
                self.0
            }

            #[inline]
            fn into_option(self) -> Option<f64> {
                // allow negative speeds for use as an element in vectors
                Some(self.0)
            }
        }

        implOpsForQuantity!($t);
    };
}

implQuantity!(MetersPSec);
implQuantity!(Knots);
implQuantity!(MilesPHour);
implQuantity!(KmPHour);

double_conversion!(MetersPSec, KmPHour, 3.6, 0.0, 1.0);
double_conversion!(Knots, KmPHour, 1.852, 0.0, 1.0);
double_conversion!(MilesPHour, KmPHour, 1.609_34, 0.0, 1.0);

double_conversion!(MetersPSec, Knots, 1.943_844_494, 0.0, 1.0);
double_conversion!(MilesPHour, Knots, 0.868_976, 0.0, 1.0);

double_conversion!(MetersPSec, MilesPHour, 2.236_94, 0.0, 1.0);
