//! Pressure units
use error::*;
use std::cmp::Ordering;
use std::fmt::Display;
use std::ops::{Add, Sub};
use types::Quantity;

/// Marker trait for Pressure types.
pub trait Pressure: Quantity {}

/// Pressure in hPa units.
#[derive(Clone, Copy, Debug)]
pub struct HectoPascal(pub f64);

/// Pressure in mb units.
#[derive(Clone, Copy, Debug)]
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
            fn into_result(self) -> Result<f64> {
                if self.0 < 0.0 {
                    Err(MetForErr::NegativePressure)
                } else {
                    Ok(self.0)
                }
            }

            #[inline]
            fn borrow_inner(&self) -> &f64 {
                &self.0
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

impl Display for HectoPascal {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{} hPa", self.0)
    }
}

impl Display for Millibar {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{} mb", self.0)
    }
}
