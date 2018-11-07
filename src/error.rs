//! Error reporting types for the `metfor` crate.
use std::error::Error;
use std::fmt::Display;

/// Error encountered during calculation.
#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum MetForErr {
    /// Colder than absolute zero.
    BelowAbsoluteZero,

    /// A negative pressure value was encountered.
    NegativePressure,

    /// A partial pressure was calculated as higher than the total atmospheric pressure.
    VaporPressureTooHigh,

    /// Invalid inputs, because they represent physically impossible values.
    UnphysicalInput,

    /// Input out of allowable range. This is usually for emperical methods which should not be
    /// allowed to extrapolate.
    InputOutOfRange,
}

/// Shorthand for a result type.
pub type Result<T> = ::std::result::Result<T, MetForErr>;

impl Display for MetForErr {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        use MetForErr::*;

        match self {
            BelowAbsoluteZero => write!(f, "temperature below absolute zero"),
            NegativePressure => write!(f, "negative pressure"),
            VaporPressureTooHigh => write!(f, "vapor pressure higher than total pressure"),
            UnphysicalInput => write!(f, "unspecified unphysical input"),
            InputOutOfRange => write!(f, "input out of allowable range"),
        }
    }
}

impl Error for MetForErr {}
