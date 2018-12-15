//! Length units for elevation and geopotential height.
use crate::types::Quantity;
use std::borrow::Borrow;
use std::cmp::Ordering;
use std::fmt::Display;
use std::ops::{Add, Sub};

/// Marker trait for elevation/height types.
pub trait GeoPotHeight: Quantity + PartialEq + PartialOrd {}

/// Geopotential Height in meters.
#[derive(Clone, Copy, Debug)]
pub struct Meters(pub f64);

/// Geopotential Height in feet.
#[derive(Clone, Copy, Debug)]
pub struct Feet(pub f64);

/// Geopotential Height in kilometers.
#[derive(Clone, Copy, Debug)]
pub struct KM(pub f64);

impl GeoPotHeight for Meters {}
impl GeoPotHeight for Feet {}
impl GeoPotHeight for KM {}

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
                // Allow negative values for negative elevation.
                self.0
            }

            #[inline]
            fn into_option(self) -> Option<f64> {
                Some(self.0)
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

implQuantity!(Meters);
implQuantity!(Feet);
implQuantity!(KM);

impl From<Feet> for Meters {
    #[inline]
    fn from(h: Feet) -> Self {
        Meters(h.0 / 3.2808)
    }
}

impl From<Feet> for KM {
    #[inline]
    fn from(h: Feet) -> Self {
        KM(h.0 / 3_280.8)
    }
}

impl From<Meters> for Feet {
    #[inline]
    fn from(h: Meters) -> Self {
        Feet(h.0 * 3.2808)
    }
}

impl From<Meters> for KM {
    #[inline]
    fn from(h: Meters) -> Self {
        KM(h.0 / 1000.0)
    }
}

impl From<KM> for Feet {
    #[inline]
    fn from(h: KM) -> Self {
        Feet(h.0 * 3_280.8)
    }
}

impl From<KM> for Meters {
    #[inline]
    fn from(h: KM) -> Self {
        Meters(h.0 * 1000.0)
    }
}

impl Display for Meters {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.0}m", self.0)
    }
}

impl Display for Feet {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.0}ft", self.0)
    }
}

impl Display for KM {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.2}km", self.0)
    }
}
