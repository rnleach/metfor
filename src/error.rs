//! Error reporting types for the `metfor` crate.

/// Error encountered during calculation.
#[derive(Fail, Debug, PartialEq, Eq, Clone, Copy)]
pub enum MetForErr {
    /// Colder than absolute zero.
    #[fail(display = "temperature below absolute zero encountered")]
    BelowAbsoluteZero,

    /// Invalid inputs, because they represent physically impossible values.
    #[fail(display = "unphysical input values")]
    UnphysicalInput,

    /// Input out of allowable range. This is usually for emperical methods which should not be
    /// allowed to extrapolate.
    #[fail(display = "input value out of allowable range")]
    InputOutOfRange,

    /// A negative pressure value was encountered.
    #[fail(display = "negative pressure value encountered")]
    NegativePressure,

    /// A partial pressure was calculated as higher than the total atmospheric pressure.
    #[fail(display = "vapor pressure greater than total pressure encountered")]
    VaporPressureTooHigh,
}

/// Shorthand for a result type.
pub type Result<T> = ::std::result::Result<T, MetForErr>;
