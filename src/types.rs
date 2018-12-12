#![macro_use]
//! New type wrappers for meteorlogical units.

pub use self::pressures::*;
pub use self::specific_energy::*;
pub use self::temperatures::*;
pub use self::unitless::*;
pub use self::winds::*;

// pub use self::pressure_vertical_velocity::*
// pub use self::geopotential_height::*

/// A quantity is a common super trait for types that represent units of measurement.
pub trait Quantity: Copy + Debug + Display + Sized + Borrow<f64> {
    /// Create a new instance of self by wrapping a value
    fn pack(_: f64) -> Self;

    /// Unpack a wrapped value without any error checking.
    fn unpack(self) -> f64;

    /// Unwrap the value from the new type and check for validity, panic if contents are invalid.
    fn unwrap(self) -> f64;

    /// Convert into an option that is `None` if the content is invalid.
    fn into_option(self) -> Option<f64>;
}

/// A version of `Quantity` for vectors.
pub trait VectorQuantity: Copy + Debug + Display + Sized {
    /// Create a new instance of self by wrapping some values.
    fn pack(_: (f64, f64)) -> Self;

    /// Unpack a wrapped value without any error checking.
    fn unpack(self) -> (f64, f64);

    /// Unwrap the values from the new type and check validity, panic if contents are invalid.
    fn unwrap(self) -> (f64, f64);

    /// Convert into an option that is `None` if the content is invalid.
    fn into_option(self) -> Option<(f64, f64)>;
}

// Not exported
use std::borrow::Borrow;
use std::fmt::{Debug, Display};

macro_rules! implOpsForQuantity {
    ($t:tt) => {
        impl<T> Add<T> for $t
        where
            $t: From<T> + Quantity,
            T: Quantity,
        {
            type Output = $t;

            #[inline]
            fn add(self, rhs: T) -> $t {
                let rhs = $t::from(rhs);
                Self::pack(self.unpack() + rhs.unpack())
            }
        }

        impl<T> Sub<T> for $t
        where
            $t: From<T> + Quantity,
            T: Quantity,
        {
            type Output = $t;

            #[inline]
            fn sub(self, rhs: T) -> $t {
                let rhs = $t::from(rhs);
                Self::pack(self.unpack() - rhs.unpack())
            }
        }

        impl<T> PartialEq<T> for $t
        where
            $t: From<T> + Quantity,
            T: Quantity,
        {
            #[inline]
            fn eq(&self, other: &T) -> bool {
                let other = $t::from(*other);
                self.unpack() == other.unpack()
            }
        }

        impl<T> PartialOrd<T> for $t
        where
            $t: From<T> + Quantity,
            T: Quantity,
        {
            #[inline]
            fn partial_cmp(&self, other: &T) -> Option<Ordering> {
                let other = $t::from(*other);
                let other: &f64 = other.borrow();
                let me: &f64 = self.borrow();
                std::cmp::PartialOrd::partial_cmp(me, other)
            }
        }

        #[cfg(feature = "use_optional")]
        impl optional::Noned for $t
        where
            $t: Quantity,
        {
            #[inline]
            fn is_none(&self) -> bool {
                optional::Noned::is_none(Borrow::<f64>::borrow(self))
            }

            #[inline]
            fn get_none() -> $t {
                $t::pack(optional::Noned::get_none())
            }
        }
    };
}

mod pressures;
mod specific_energy;
mod temperatures;
mod unitless;
mod winds;

#[cfg(all(test, feature = "use_optional"))]
mod test {
    use crate::types::temperatures::Celsius;
    use optional::*;

    #[test]
    fn test_optional() {
        let val: Optioned<Celsius> = none();

        assert!(val.is_none());
    }
}
