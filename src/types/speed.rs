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

//--------------------------------------------------------------------------------------------------
//                                    Kilometer per hour Conversions
//--------------------------------------------------------------------------------------------------
impl From<MetersPSec> for KmPHour {
    #[inline]
    fn from(s: MetersPSec) -> Self {
        KmPHour(s.0 / 1_000.0 * 3_600.0)
    }
}

impl From<Knots> for KmPHour {
    #[inline]
    fn from(s: Knots) -> Self {
        KmPHour(s.0 * 1.852)
    }
}

impl From<MilesPHour> for KmPHour {
    #[inline]
    fn from(s: MilesPHour) -> Self {
        KmPHour(s.0 * 1.609_34)
    }
}

//--------------------------------------------------------------------------------------------------
//                                      Knots Conversions
//--------------------------------------------------------------------------------------------------
impl From<MetersPSec> for Knots {
    #[inline]
    fn from(s: MetersPSec) -> Self {
        Knots(s.0 / 0.514_444_44)
    }
}

impl From<KmPHour> for Knots {
    #[inline]
    fn from(s: KmPHour) -> Self {
        Knots(s.0 * 0.5399_57)
    }
}

impl From<MilesPHour> for Knots {
    #[inline]
    fn from(s: MilesPHour) -> Self {
        Knots(s.0 * 0.868_976)
    }
}

//--------------------------------------------------------------------------------------------------
//                                      MilesPHour Conversions
//--------------------------------------------------------------------------------------------------
impl From<MetersPSec> for MilesPHour {
    #[inline]
    fn from(s: MetersPSec) -> Self {
        MilesPHour(s.0 * 2.236_94)
    }
}

impl From<KmPHour> for MilesPHour {
    #[inline]
    fn from(s: KmPHour) -> Self {
        MilesPHour(s.0 * 0.621_371)
    }
}

impl From<Knots> for MilesPHour {
    #[inline]
    fn from(s: Knots) -> Self {
        MilesPHour(s.0 * 1.150_78)
    }
}

//--------------------------------------------------------------------------------------------------
//                                    MeterPSec Conversions
//--------------------------------------------------------------------------------------------------
impl From<MilesPHour> for MetersPSec {
    #[inline]
    fn from(s: MilesPHour) -> Self {
        MetersPSec(s.0 * 0.447_04)
    }
}

impl From<KmPHour> for MetersPSec {
    #[inline]
    fn from(s: KmPHour) -> Self {
        MetersPSec(s.0 * 0.277_778)
    }
}

impl From<Knots> for MetersPSec {
    #[inline]
    fn from(s: Knots) -> Self {
        MetersPSec(s.0 * 0.514_444_44)
    }
}
