//! Temperature units
use crate::constants::*;
use crate::types::Quantity;
use std::cmp::Ordering;
use std::ops::{Add, AddAssign, Sub};

/// Marker trait for temperature types.
pub trait Temperature: Quantity + PartialEq + PartialOrd {}

/// Marker trait for temperature differences.
pub trait TempDiff: Quantity + PartialOrd + PartialEq {}

/// Temperature in Fahrenheit units.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Fahrenheit(pub f64);

/// Temperature in Celsius units.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Celsius(pub f64);

/// Temperature in Kelvin units.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Kelvin(pub f64);

impl Temperature for Fahrenheit {}
impl Temperature for Celsius {}
impl Temperature for Kelvin {}

/// Temperature difference in Fahrenheit units.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct FahrenheitDiff(pub f64);

/// Temperature difference in Celsius units.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
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

impl<T> AddAssign<T> for Celsius
where
    CelsiusDiff: From<T>,
{
    #[inline]
    fn add_assign(&mut self, rhs: T) {
        let rhs = CelsiusDiff::from(rhs);
        *self = Self::pack(self.unpack() + rhs.unpack())
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

impl<T> AddAssign<T> for Kelvin
where
    KelvinDiff: From<T>,
{
    #[inline]
    fn add_assign(&mut self, rhs: T) {
        let rhs = KelvinDiff::from(rhs);
        *self = Self::pack(self.unpack() + rhs.unpack())
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

impl<T> AddAssign<T> for Fahrenheit
where
    FahrenheitDiff: From<T>,
{
    #[inline]
    fn add_assign(&mut self, rhs: T) {
        let rhs = FahrenheitDiff::from(rhs);
        *self = Self::pack(self.unpack() + rhs.unpack())
    }
}

double_conversion!(Celsius, Fahrenheit, 1.8, 32.0, 1.0);
double_conversion!(Kelvin, Fahrenheit, 1.8, 273.15 + 32.0, 1.0);
double_conversion!(Celsius, Kelvin, 1.0, 273.15, 1.0);

double_conversion!(CelsiusDiff, FahrenheitDiff, 1.8, 0.0, 1.0);

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

    #[test]
    fn test_temperature_diff_conversions() {
        let cd = CelsiusDiff(100.0);
        let fd = FahrenheitDiff::from(cd);

        assert!(approx_equal(fd, FahrenheitDiff(180.0), FahrenheitDiff(TOL)));
    }
}
