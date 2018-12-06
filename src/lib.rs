#![warn(missing_docs)]
//! Meteorological constants and formulas.
//!
//! I investigated using some sort of dimensional analysis via types with a crate like [uom][uom]
//! or [dimensioned][dimensioned]. However after experimentation, neither of these work well for a
//! library. Choosing to use one would force that library on the users of this library. In the
//! future I may make a feature to use one or the other of these crates.
//!
//! [uom]: https://crates.io/crates/uom
//! [dimensioned]: https://crates.io/crates/dimensioned

//
// API
//

pub mod constants;
pub use constants::*;
pub use error::*;
pub use functions::*;
pub use types::*;

//
// Internal use only
//
mod error;
mod functions;
mod types;

#[cfg(test)]
mod test_utils;
