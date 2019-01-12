//! Lapse rate units
use crate::types::{Feet, Km, Quantity};
use std::cmp::Ordering;
use std::fmt::Display;

/// Marker trait for temperature lapse rate types.
pub trait TempLR: Quantity + PartialEq + PartialOrd {}

/// Marker trait for hydrolapse.
pub trait Hydrolapse: Quantity + PartialOrd + PartialEq {}

/// Temperature lapse rate in Fahrenheit per thousand feet (kft).
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct FahrenheitPKft(pub f64);

/// Temperature lapse rate in Celsius per kilometer.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct CelsiusPKm(pub f64);

/// Temperature lapse rate in Kelvin per kilometer.
pub type KelvinPKm = CelsiusPKm;

/// Hydrolapse in for mixing ratio in km<sup>-1</sup>
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct HydrolapsePKm(pub f64);

/// Hydrolapse in for mixing ratio in g / kg/ km
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct HydrolapseGPKgPKm(pub f64);

impl TempLR for FahrenheitPKft {}
impl TempLR for CelsiusPKm {}

impl Hydrolapse for HydrolapsePKm {}

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

implQuantity!(FahrenheitPKft);
implQuantity!(CelsiusPKm);

implQuantity!(HydrolapsePKm);
implQuantity!(HydrolapseGPKgPKm);

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

impl From<HydrolapsePKm> for HydrolapseGPKgPKm {
    #[inline]
    fn from(hl: HydrolapsePKm) -> Self {
        HydrolapseGPKgPKm(1000.0 * hl.unpack())
    }
}

impl From<HydrolapseGPKgPKm> for HydrolapsePKm {
    #[inline]
    fn from(hl: HydrolapseGPKgPKm) -> Self {
        HydrolapsePKm(hl.unpack() / 1000.0)
    }
}

impl Display for CelsiusPKm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.1}\u{00B0}C/km", self.0)
    }
}

impl Display for FahrenheitPKft {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.1}\u{00B0}F/kft", self.0)
    }
}

impl Display for HydrolapsePKm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.4}km\u{207B}\u{2081}", self.0)
    }
}

impl Display for HydrolapseGPKgPKm {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.1}g kg\u{207B}\u{2081} km\u{207B}\u{2081}", self.0)
    }
}
