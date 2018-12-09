//! Temperature units
use crate::constants::*;
use crate::error::*;
use crate::types::Quantity;
use std::cmp::Ordering;
use std::fmt::Display;
use std::ops::{Add, Sub};

/// Marker trait for temperature types.
pub trait Temperature: Quantity + PartialEq + PartialOrd {}

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
                if self < ABSOLUTE_ZERO {
                    panic!("Below absolute zero.");
                }

                self.0
            }

            #[inline]
            fn into_result(self) -> Result<f64> {
                if self < ABSOLUTE_ZERO {
                    Err(MetForErr::BelowAbsoluteZero)
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

implQuantity!(Fahrenheit);
implQuantity!(Celsius);
implQuantity!(Kelvin);

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

#[cfg(test)]
mod test {
    use super::*;
    use crate::test_utils::*;

    pub const TOL: f64 = 1.0e-9;

    use crate::error::MetForErr::*;

    #[test]
    fn test_celsius_to_kelvin() {
        assert!(approx_equal(
            Kelvin::from(Celsius(-10.0)),
            Kelvin(263.15),
            Kelvin(TOL)
        ));
        assert!(approx_equal(
            Kelvin::from(Celsius(0.0)),
            Kelvin(273.15),
            Kelvin(TOL)
        ));
        assert!(approx_equal(
            Kelvin::from(Celsius(10.0)),
            Kelvin(283.15),
            Kelvin(TOL)
        ));
        assert!(approx_equal(
            Kelvin::from(Celsius(-273.15)),
            Kelvin(0.0),
            Kelvin(TOL)
        ));
        assert!(Kelvin::from(Celsius(-300.0)).into_result().unwrap_err() == BelowAbsoluteZero);
    }

    #[test]
    fn test_kelvin_to_celsius() {
        assert!(approx_equal(
            Celsius::from(Kelvin(263.15)),
            Celsius(-10.0),
            Celsius(TOL)
        ));
        assert!(approx_equal(
            Celsius::from(Kelvin(273.15)),
            Celsius(0.0),
            Celsius(TOL)
        ));
        assert!(approx_equal(
            Celsius::from(Kelvin(283.15)),
            Celsius(10.0),
            Celsius(TOL)
        ));
        assert!(approx_equal(
            Celsius::from(Kelvin(0.0)),
            Celsius(-273.15),
            Celsius(TOL)
        ));
        assert!(Celsius::from(Kelvin(-10.0)).into_result().unwrap_err() == BelowAbsoluteZero);
    }

    #[test]
    fn test_kelvin_celsius_conversions_are_inverses_of_each_other() {
        for celsius in temperatures() {
            let kelvin = Kelvin::from(celsius);
            let back_to_celsius = Celsius::from(kelvin);
            assert!(approx_equal(celsius, back_to_celsius, Celsius(TOL)));
        }
    }

    #[test]
    fn test_celsius_to_f() {
        assert!(approx_equal(
            Fahrenheit::from(Celsius(0.0)),
            Fahrenheit(32.0),
            Fahrenheit(TOL)
        ));
        assert!(approx_equal(
            Fahrenheit::from(Celsius(100.0)),
            Fahrenheit(212.0),
            Fahrenheit(TOL)
        ));
        assert!(approx_equal(
            Fahrenheit::from(Celsius(-40.0)),
            Fahrenheit(-40.0),
            Fahrenheit(TOL)
        ));
        assert!((ABSOLUTE_ZERO - Celsius(1.0)).into_result().unwrap_err() == BelowAbsoluteZero);
    }

    #[test]
    fn test_f_to_celsius() {
        assert!(approx_equal(
            Celsius::from(Fahrenheit(32.0)),
            Celsius(0.0),
            Celsius(TOL)
        ));
        assert!(approx_equal(
            Celsius::from(Fahrenheit(212.0)),
            Celsius(100.0),
            Celsius(TOL)
        ));
        assert!(approx_equal(
            Celsius::from(Fahrenheit(-40.0)),
            Celsius(-40.0),
            Celsius(TOL)
        ));
        assert!((ABSOLUTE_ZERO - Fahrenheit(1.0)).into_result().unwrap_err() == BelowAbsoluteZero);
    }

    #[test]
    fn test_f_celsius_conversions_are_inverses_of_each_other() {
        for celsius in temperatures() {
            let fahrenheit = Fahrenheit::from(celsius);
            let back_to_celsius = Celsius::from(fahrenheit);
            assert!(approx_equal(celsius, back_to_celsius, Celsius(TOL)));
        }
    }
}
