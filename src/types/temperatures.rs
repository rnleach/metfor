//! Temperature units
use crate::constants::*;
use crate::types::Quantity;
use std::cmp::Ordering;
use std::fmt::Display;
use std::ops::{Add, Sub};

/// Marker trait for temperature types.
pub trait Temperature: Quantity + PartialEq + PartialOrd {}

/// Marker trait for temperature differences.
pub trait TempDiff: Quantity + PartialOrd + PartialEq {}

/// Temperature in Fahrenheit units.
#[derive(Clone, Copy, Debug)]
pub struct Fahrenheit(pub f64);

/// Temperature in Celsius units.
#[derive(Clone, Copy, Debug)]
pub struct Celsius(pub f64);

/// Temperature in Kelvin units.
#[derive(Clone, Copy, Debug)]
pub struct Kelvin(pub f64);

impl Temperature for Fahrenheit {}
impl Temperature for Celsius {}
impl Temperature for Kelvin {}

/// Temperature difference in Fahrenheit units.
#[derive(Clone, Copy, Debug)]
pub struct FahrenheitDiff(pub f64);

/// Temperature difference in Celsius units.
#[derive(Clone, Copy, Debug)]
pub struct CelsiusDiff(pub f64);

/// Temperature difference in Kelvin units.
pub type KelvinDiff = CelsiusDiff;

impl TempDiff for FahrenheitDiff {}
impl TempDiff for CelsiusDiff {}

macro_rules! implQuantityForT {
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
                if self < ABSOLUTE_ZERO {
                    panic!("Below absolute zero.");
                }

                self.0
            }

            #[inline]
            fn into_option(self) -> Option<f64> {
                if self < ABSOLUTE_ZERO {
                    None
                } else {
                    Some(self.0)
                }
            }
        }

        implNonedForQuantity!($t);
        implMulDivOpsForQuantity!($t);
        implOrdEqOpsForQuantity!($t);
    };
}

implQuantityForT!(Fahrenheit);
implQuantityForT!(Celsius);
implQuantityForT!(Kelvin);

macro_rules! implQuantityForDiffT {
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
implQuantityForDiffT!(FahrenheitDiff);
implQuantityForDiffT!(CelsiusDiff);

impl<T> Sub<T> for Celsius
where
    Celsius: From<T>,
{
    type Output = CelsiusDiff;

    #[inline]
    fn sub(self, rhs: T) -> Self::Output {
        let rhs = Celsius::from(rhs);
        Self::Output::pack(self.unpack() - rhs.unpack())
    }
}

impl<T> Add<T> for Celsius
where
    CelsiusDiff: From<T>,
{
    type Output = Celsius;

    #[inline]
    fn add(self, rhs: T) -> Self::Output {
        let rhs = CelsiusDiff::from(rhs);
        Self::Output::pack(self.unpack() + rhs.unpack())
    }
}

impl<T> Sub<T> for Kelvin
where
    Kelvin: From<T>,
{
    type Output = KelvinDiff;

    #[inline]
    fn sub(self, rhs: T) -> Self::Output {
        let rhs = Kelvin::from(rhs);
        Self::Output::pack(self.unpack() - rhs.unpack())
    }
}

impl<T> Add<T> for Kelvin
where
    CelsiusDiff: From<T>,
{
    type Output = Kelvin;

    #[inline]
    fn add(self, rhs: T) -> Self::Output {
        let rhs = CelsiusDiff::from(rhs);
        Self::Output::pack(self.unpack() + rhs.unpack())
    }
}

impl<T> Sub<T> for Fahrenheit
where
    Fahrenheit: From<T>,
{
    type Output = FahrenheitDiff;

    #[inline]
    fn sub(self, rhs: T) -> Self::Output {
        let rhs = Fahrenheit::from(rhs);
        Self::Output::pack(self.unpack() - rhs.unpack())
    }
}

impl<T> Add<T> for Fahrenheit
where
    FahrenheitDiff: From<T>,
{
    type Output = Fahrenheit;

    #[inline]
    fn add(self, rhs: T) -> Self::Output {
        let rhs = FahrenheitDiff::from(rhs);
        Self::Output::pack(self.unpack() + rhs.unpack())
    }
}

impl From<Celsius> for Fahrenheit {
    #[inline]
    fn from(c: Celsius) -> Self {
        Fahrenheit(1.8 * c.0 + 32.0)
    }
}

impl From<Kelvin> for Fahrenheit {
    #[inline]
    fn from(k: Kelvin) -> Self {
        Fahrenheit::from(Celsius::from(k))
    }
}

impl From<Fahrenheit> for Celsius {
    #[inline]
    fn from(f: Fahrenheit) -> Self {
        Celsius((f.0 - 32.0) / 1.8)
    }
}

impl From<Kelvin> for Celsius {
    #[inline]
    fn from(k: Kelvin) -> Self {
        Celsius(k.0 + ABSOLUTE_ZERO_C)
    }
}

impl From<Celsius> for Kelvin {
    #[inline]
    fn from(c: Celsius) -> Self {
        Kelvin(c.0 - ABSOLUTE_ZERO_C)
    }
}

impl From<Fahrenheit> for Kelvin {
    #[inline]
    fn from(f: Fahrenheit) -> Self {
        Kelvin::from(Celsius::from(f))
    }
}

impl From<FahrenheitDiff> for CelsiusDiff {
    #[inline]
    fn from(f: FahrenheitDiff) -> Self {
        CelsiusDiff(f.0 / 1.8)
    }
}

impl From<CelsiusDiff> for FahrenheitDiff {
    #[inline]
    fn from(c: CelsiusDiff) -> Self {
        FahrenheitDiff(c.0 * 1.8)
    }
}

impl Display for Celsius {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.1}\u{00B0}C", self.0)
    }
}

impl Display for Fahrenheit {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.1}\u{00B0}F", self.0)
    }
}

impl Display for Kelvin {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.1}\u{00B0}K", self.0)
    }
}

impl Display for CelsiusDiff {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.1}\u{0394}\u{00B0}(C or K)", self.0)
    }
}

impl Display for FahrenheitDiff {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.1}\u{0394}\u{00B0}F", self.0)
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::test_utils::*;

    pub const TOL: f64 = 1.0e-9;

    #[test]
    fn test_celsius_to_kelvin() {
        assert!(approx_equal(
            Kelvin::from(Celsius(-10.0)),
            Kelvin(263.15),
            CelsiusDiff(TOL)
        ));
        assert!(approx_equal(
            Kelvin::from(Celsius(0.0)),
            Kelvin(273.15),
            CelsiusDiff(TOL)
        ));
        assert!(approx_equal(
            Kelvin::from(Celsius(10.0)),
            Kelvin(283.15),
            CelsiusDiff(TOL)
        ));
        assert!(approx_equal(
            Kelvin::from(Celsius(-273.15)),
            Kelvin(0.0),
            CelsiusDiff(TOL)
        ));
        assert!(Kelvin::from(Celsius(-300.0)).into_option().is_none());
    }

    #[test]
    fn test_kelvin_to_celsius() {
        assert!(approx_equal(
            Celsius::from(Kelvin(263.15)),
            Celsius(-10.0),
            CelsiusDiff(TOL)
        ));
        assert!(approx_equal(
            Celsius::from(Kelvin(273.15)),
            Celsius(0.0),
            CelsiusDiff(TOL)
        ));
        assert!(approx_equal(
            Celsius::from(Kelvin(283.15)),
            Celsius(10.0),
            CelsiusDiff(TOL)
        ));
        assert!(approx_equal(
            Celsius::from(Kelvin(0.0)),
            Celsius(-273.15),
            CelsiusDiff(TOL)
        ));
        assert!(Celsius::from(Kelvin(-10.0)).into_option().is_none());
    }

    #[test]
    fn test_kelvin_celsius_conversions_are_inverses_of_each_other() {
        for celsius in temperatures() {
            let kelvin = Kelvin::from(celsius);
            let back_to_celsius = Celsius::from(kelvin);
            assert!(approx_equal(celsius, back_to_celsius, CelsiusDiff(TOL)));
        }
    }

    #[test]
    fn test_celsius_to_f() {
        assert!(approx_equal(
            Fahrenheit::from(Celsius(0.0)),
            Fahrenheit(32.0),
            FahrenheitDiff(TOL)
        ));
        assert!(approx_equal(
            Fahrenheit::from(Celsius(100.0)),
            Fahrenheit(212.0),
            FahrenheitDiff(TOL)
        ));
        assert!(approx_equal(
            Fahrenheit::from(Celsius(-40.0)),
            Fahrenheit(-40.0),
            FahrenheitDiff(TOL),
        ));
        assert!((ABSOLUTE_ZERO + -CelsiusDiff(1.0)).into_option().is_none());
    }

    #[test]
    fn test_f_to_celsius() {
        assert!(approx_equal(
            Celsius::from(Fahrenheit(32.0)),
            Celsius(0.0),
            CelsiusDiff(TOL)
        ));
        assert!(approx_equal(
            Celsius::from(Fahrenheit(212.0)),
            Celsius(100.0),
            CelsiusDiff(TOL)
        ));
        assert!(approx_equal(
            Celsius::from(Fahrenheit(-40.0)),
            Celsius(-40.0),
            CelsiusDiff(TOL)
        ));
        assert!((ABSOLUTE_ZERO + -FahrenheitDiff(1.0))
            .into_option()
            .is_none());
    }

    #[test]
    fn test_f_celsius_conversions_are_inverses_of_each_other() {
        for celsius in temperatures() {
            let fahrenheit = Fahrenheit::from(celsius);
            let back_to_celsius = Celsius::from(fahrenheit);
            assert!(approx_equal(celsius, back_to_celsius, CelsiusDiff(TOL)));
        }
    }
}
