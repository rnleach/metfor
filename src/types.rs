#![macro_use]
//! New type wrappers for meteorlogical units.

pub use self::lapse_rates::*;
pub use self::length::*;
pub use self::pressure_vertical_velocity::*;
pub use self::pressures::*;
pub use self::specific_energy::*;
pub use self::temperatures::*;
pub use self::unitless::*;
pub use self::winds::*;

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
    /// Create a new instance of self by wrapping some values. This must be x-y coordinates from the
    /// standard cartesian coordinate system.
    fn pack_xy(_: (f64, f64)) -> Self;

    /// Unpack a wrapped value without any error checking. The returned values represent the vector
    /// in a standard x-y cartesian coordinate system.
    fn unpack_xy(self) -> (f64, f64);

    /// Unwrap the values from the new type and check validity, panic if contents are invalid. The
    /// returned values represent the vector in a standard x-y cartesian coordinate system.
    fn unwrap_xy(self) -> (f64, f64);

    /// Convert into an option that is `None` if the content is invalid. The returned values
    /// represent the vector in a standard x-y cartesian coordinate system.
    fn into_option(self) -> Option<(f64, f64)>;
}

//
// Not exported
//
use std::borrow::Borrow;
use std::fmt::{Debug, Display};

macro_rules! implAddSubOpsForQuantity {
    ($t:tt) => {
        impl<T> Add<T> for $t
        where
            $t: From<T> + Quantity,
            T: Quantity,
        {
            type Output = $t;

            #[inline]
            fn add(self, rhs: T) -> Self::Output {
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
            fn sub(self, rhs: T) -> Self::Output {
                let rhs = $t::from(rhs);
                Self::pack(self.unpack() - rhs.unpack())
            }
        }
    };
}

macro_rules! implMulDivOpsForQuantity {
    ($t:tt) => {
        impl<T> Div<T> for $t
        where
            $t: From<T> + Quantity,
            T: Quantity,
        {
            type Output = f64;

            #[inline]
            fn div(self, rhs: T) -> Self::Output {
                let rhs = $t::from(rhs);
                self.unpack() / rhs.unpack()
            }
        }

        impl Mul<f64> for $t {
            type Output = $t;

            #[inline]
            fn mul(self, rhs: f64) -> Self::Output {
                Self::pack(self.unpack() * rhs)
            }
        }
    };
}

macro_rules! implOrdEqOpsForQuantity {
    ($t:tt) => {
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
        impl optional::OptEq for $t
        where
            $t: Quantity,
        {
            #[inline]
            fn opt_eq(&self, other: &Self) -> bool {
                self == other
            }
        }

        #[cfg(feature = "use_optional")]
        impl optional::OptOrd for $t
        where
            $t: Quantity + optional::Noned,
        {
            #[inline]
            fn opt_cmp(&self, other: &Self) -> Ordering {
                use optional::Noned;

                if self.is_none() {
                    if other.is_none() {
                        Ordering::Equal
                    } else {
                        Ordering::Less
                    }
                } else if other.is_none() {
                    Ordering::Greater
                } else {
                    self.partial_cmp(other).unwrap()
                }
            }
        }
    };
}

macro_rules! implNonedForQuantity {
    ($t:tt) => {
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

macro_rules! implOpsForQuantity {
    ($t:tt) => {
        implNonedForQuantity!($t);
        implAddSubOpsForQuantity!($t);
        implMulDivOpsForQuantity!($t);
        implOrdEqOpsForQuantity!($t);
    };
}

macro_rules! implOpsForVectorQuantity {
    ($t:tt) => {
        impl<T> Add<T> for $t
        where
            $t: From<T> + VectorQuantity,
            T: VectorQuantity,
        {
            type Output = $t;

            #[inline]
            fn add(self, rhs: T) -> $t {
                let rhs = $t::from(rhs);
                let (x, y) = self.unpack_xy();
                let (rhs_x, rhs_y) = rhs.unpack_xy();

                let x_res = x + rhs_x;
                let y_res = y + rhs_y;

                Self::pack_xy((x_res, y_res))
            }
        }

        impl<T> Sub<T> for $t
        where
            $t: From<T> + VectorQuantity,
            T: VectorQuantity,
        {
            type Output = $t;

            #[inline]
            fn sub(self, rhs: T) -> $t {
                let rhs = $t::from(rhs);
                let (x, y) = self.unpack_xy();
                let (rhs_x, rhs_y) = rhs.unpack_xy();

                let x_res = x - rhs_x;
                let y_res = y - rhs_y;

                Self::pack_xy((x_res, y_res))
            }
        }

        impl<T> PartialEq<T> for $t
        where
            $t: From<T> + VectorQuantity,
            T: VectorQuantity,
        {
            #[inline]
            fn eq(&self, rhs: &T) -> bool {
                let rhs = $t::from(*rhs);
                let (x, y) = self.unpack_xy();
                let (rhs_x, rhs_y) = rhs.unpack_xy();

                x == rhs_x && y == rhs_y
            }
        }

        #[cfg(feature = "use_optional")]
        impl optional::OptEq for $t
        where
            $t: VectorQuantity,
        {
            #[inline]
            fn opt_eq(&self, other: &Self) -> bool {
                self == other
            }
        }
    };
}

mod lapse_rates;
mod length;
mod pressure_vertical_velocity;
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
