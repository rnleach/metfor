//! Length units for elevation and geopotential height.
use crate::types::Quantity;
use std::borrow::Borrow;
use std::cmp::Ordering;
use std::fmt::Display;

/// Marker trait for speed types.
pub trait Speed: Quantity + PartialEq + PartialOrd {}

/// Speed in meters per second.
#[derive(Clone, Copy, Debug)]
pub struct MetersPSec(pub f64);

/// Speed in knots.
#[derive(Clone, Copy, Debug)]
pub struct Knots(pub f64);

/// Speed in miles per hour
#[derive(Clone, Copy, Debug)]
pub struct MilesPHour(pub f64);

/// Speed in Kilometers per hour
#[derive(Clone, Copy, Debug)]
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
                if self.0 < 0.0 {
                    panic!("negative speed");
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

        impl Borrow<f64> for $t {
            #[inline]
            fn borrow(&self) -> &f64 {
                &self.0
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
        KmPHour(s.0 / 1000.0 * 3600.0)
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
        KmPHour(s.0 * 1.60934)
    }
}

//--------------------------------------------------------------------------------------------------
//                                      Knots Conversions
//--------------------------------------------------------------------------------------------------
impl From<MetersPSec> for Knots {
    #[inline]
    fn from(s: MetersPSec) -> Self {
        Knots(s.0 * 1.94384)
    }
}

impl From<KmPHour> for Knots {
    #[inline]
    fn from(s: KmPHour) -> Self {
        Knots(s.0 * 0.539957)
    }
}

impl From<MilesPHour> for Knots {
    #[inline]
    fn from(s: MilesPHour) -> Self {
        Knots(s.0 * 0.868976)
    }
}

//--------------------------------------------------------------------------------------------------
//                                      MilesPHour Conversions
//--------------------------------------------------------------------------------------------------
impl From<MetersPSec> for MilesPHour {
    #[inline]
    fn from(s: MetersPSec) -> Self {
        MilesPHour(s.0 * 2.23694)
    }
}

impl From<KmPHour> for MilesPHour {
    #[inline]
    fn from(s: KmPHour) -> Self {
        MilesPHour(s.0 * 0.621371)
    }
}

impl From<Knots> for MilesPHour {
    #[inline]
    fn from(s: Knots) -> Self {
        MilesPHour(s.0 * 1.15078)
    }
}

//--------------------------------------------------------------------------------------------------
//                                    MeterPSec Conversions
//--------------------------------------------------------------------------------------------------
impl From<MilesPHour> for MetersPSec {
    #[inline]
    fn from(s: MilesPHour) -> Self {
        MetersPSec(s.0 * 0.44704)
    }
}

impl From<KmPHour> for MetersPSec {
    #[inline]
    fn from(s: KmPHour) -> Self {
        MetersPSec(s.0 * 0.277778)
    }
}

impl From<Knots> for MetersPSec {
    #[inline]
    fn from(s: Knots) -> Self {
        MetersPSec(s.0 * 0.514444)
    }
}

//--------------------------------------------------------------------------------------------------
//                                     Display implementations
//--------------------------------------------------------------------------------------------------

impl Display for MetersPSec {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.1}m/s", self.0)
    }
}

impl Display for MilesPHour {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.0}mph", self.0)
    }
}

impl Display for KmPHour {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.0}kph", self.0)
    }
}

impl Display for Knots {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.0}kt", self.0)
    }
}