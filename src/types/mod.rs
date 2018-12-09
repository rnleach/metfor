#![macro_use]
//! New type wrappers for meteorlogical units.

pub use self::pressures::*;
pub use self::specific_energy::*;
pub use self::temperatures::*;
pub use self::winds::*;

/// A quantity is a common super trait for types that represent units of measurement.
pub trait Quantity: Copy + Debug + Display {
    /// Borrow the inner value
    fn borrow_inner(&self) -> &f64;

    /// Create a new instance of self by wrapping a value
    fn pack(_: f64) -> Self;

    /// Unpack a wrapped value without any error checking.
    fn unpack(self) -> f64;

    /// Unwrap the value from the new type and check for validity, panic if contents are invalid.
    fn unwrap(self) -> f64;

    /// Convert into a result that is an error if below absolute zero.
    fn into_result(self) -> Result<f64>;
}

// Not exported

use crate::error::*;
use std::fmt::{Debug, Display};

macro_rules! implOpsForQuantity {
    ($t:tt) => {
        impl<T> Add<T> for $t
        where
            $t: From<T> + Quantity,
            T: Quantity,
        {
            type Output = $t;

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
            fn partial_cmp(&self, other: &T) -> Option<Ordering> {
                let other = $t::from(*other);
                self.borrow_inner().partial_cmp(other.borrow_inner())
            }
        }
    };
}

mod pressures;
mod specific_energy;
mod temperatures;
mod winds;
