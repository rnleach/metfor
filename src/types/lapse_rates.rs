//! Lapse rate units
use crate::types::{Feet, Km, Quantity};
use std::borrow::Borrow;
use std::cmp::Ordering;
use std::fmt::Display;

/// Marker trait for temperature lapse rate types.
pub trait TempLR: Quantity + PartialEq + PartialOrd {}

/// Temperature lapse rate in Fahrenheit per thousand feet (kft).
#[derive(Clone, Copy, Debug)]
pub struct FahrenheitPKft(pub f64);

/// Temperature lapse rate in Celsius per kilometer.
#[derive(Clone, Copy, Debug)]
pub struct CelsiusPKm(pub f64);

/// Temperature lapse rate in Kelvin per kilometer.
pub type KelvinPKm = CelsiusPKm;

impl TempLR for FahrenheitPKft {}
impl TempLR for CelsiusPKm {}

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

        impl Borrow<f64> for $t {
            #[inline]
            fn borrow(&self) -> &f64 {
                &self.0
            }
        }

        implOpsForQuantity!($t);
    };
}

implQuantity!(FahrenheitPKft);
implQuantity!(CelsiusPKm);

impl From<CelsiusPKm> for FahrenheitPKft {
    #[inline]
    fn from(c: CelsiusPKm) -> Self {
        let dt = c.unpack();
        let dz = Feet::from(Km(1.0)).unpack() / 1000.0;
        FahrenheitPKft(1.8 * dt / dz)
    }
}

impl From<FahrenheitPKft> for CelsiusPKm {
    #[inline]
    fn from(f: FahrenheitPKft) -> Self {
        let dt = f.unpack();
        let dz = Km::from(Feet(1000.0)).unpack();
        CelsiusPKm(dt / dz / 1.8)
    }
}

impl Display for CelsiusPKm {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.1}\u{00B0}C/km", self.0)
    }
}

impl Display for FahrenheitPKft {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.1}\u{00B0}F/kft", self.0)
    }
}
