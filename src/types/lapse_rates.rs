//! Lapse rate units
use crate::types::{CelsiusDiff, FahrenheitDiff, Feet, Km, Quantity};
use std::cmp::Ordering;

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
        let dt = FahrenheitDiff::from(CelsiusDiff(c.unpack())).unpack();
        let dz = Feet::from(Km(1.0)).unpack() / 1000.0;
        FahrenheitPKft(dt / dz)
    }
}

impl From<FahrenheitPKft> for CelsiusPKm {
    #[inline]
    fn from(f: FahrenheitPKft) -> Self {
        let dt = CelsiusDiff::from(FahrenheitDiff(f.unpack())).unpack();
        let dz = Km::from(Feet(1000.0)).unpack();
        CelsiusPKm(dt / dz)
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

#[cfg(test)]
mod test {
    use super::*;
    use crate::test_utils::approx_equal;

    const TOL: f64 = 1.0e-5;

    #[test]
    fn test_there_and_back() {
        let mut lr: f64 = -20.0;
        while lr < 15.0 {
            let fpkft = FahrenheitPKft(lr);
            let cpkm = CelsiusPKm::from(fpkft);
            let back = FahrenheitPKft::from(cpkm);
            assert!(approx_equal(fpkft.unpack(), back.unpack(), TOL));

            lr += 0.01;
        }
    }

    #[test]
    fn test_standard() {
        const DRY_ADIABATIC_ENGLISH: FahrenheitPKft = FahrenheitPKft(-5.3803208);
        const DRY_ADIABATIC_METRIC: CelsiusPKm = CelsiusPKm(crate::constants::g);
        let dry_adiabatic_english = FahrenheitPKft::from(DRY_ADIABATIC_METRIC);
        let dry_adiabatic_metric = CelsiusPKm::from(DRY_ADIABATIC_ENGLISH);

        assert!(approx_equal(
            dry_adiabatic_english,
            DRY_ADIABATIC_ENGLISH,
            FahrenheitPKft(TOL)
        ));

        assert!(approx_equal(
            dry_adiabatic_metric,
            DRY_ADIABATIC_METRIC,
            CelsiusPKm(TOL)
        ));
    }
}
