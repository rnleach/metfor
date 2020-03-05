//! Length units for elevation, geopotential height, and visibility.
use crate::types::Quantity;
use std::cmp::Ordering;

/// Marker trait for elevation/height types.
pub trait Length: Quantity + PartialEq + PartialOrd {}

/// Length in kilometers, usually used for geopotential height above ground level.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Km(pub f64);

/// Length in decameters, usually used for geopotential height.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Decameters(pub f64);

/// Length in meters, usually used for geopotential height and elevation.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Meters(pub f64);

/// Length in centimeters, usually used for precipitation depth or precipitable water.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Cm(pub f64);

/// Length in millimeters, usually used for precipitation depth or precipitable water.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Mm(pub f64);

/// Length in statute miles, usually used for visibility.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct StatuteMiles(pub f64);

/// Length in feet, usually used for geopotential height and elevation.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Feet(pub f64);

/// Length in inches, usually used for precipitation depth or precipitable water.
#[derive(Clone, Copy, Debug)]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Serialize))]
#[cfg_attr(feature = "use_serde", derive(serde_derive::Deserialize))]
pub struct Inches(pub f64);

impl Length for Km {}
impl Length for Decameters {}
impl Length for Meters {}
impl Length for Cm {}
impl Length for StatuteMiles {}
impl Length for Mm {}
impl Length for Feet {}
impl Length for Inches {}

// Make this macro in each module, because temperatures may have different rules than
// length or energy.
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
                // Allow negative values for negative elevation.
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

implQuantity!(Km);
implQuantity!(Decameters);
implQuantity!(Meters);
implQuantity!(Cm);
implQuantity!(Mm);
implQuantity!(StatuteMiles);
implQuantity!(Feet);
implQuantity!(Inches);

//--------------------------------------------------------------------------------------------------
//                                       Test macros
//--------------------------------------------------------------------------------------------------

#[cfg(test)]
macro_rules! make_there_and_back_tests {
    ($test_name:tt, $from:tt, $to:tt, $from_start:expr, $from_end:expr, $tol:expr) => {
        #[test]
        fn $test_name() {
            let mut val = $from_start;
            while (val < $from_end) {
                let converted = $to::from(val);
                let back_again = $from::from(converted);
                assert!(approx_equal(val, back_again, $tol));
                val += $from(($from_end - $from_start).unpack() / 100.0);
            }

            let mut val = $to::from($from_start);
            while (val < $to::from($from_end)) {
                let converted = $from::from(val);
                let back_again = $to::from(converted);
                assert!(approx_equal(val, back_again, $to::from($tol)));
                val += $to(($from_end - $from_start).unpack() / 100.0);
            }
        }
    };
}

//--------------------------------------------------------------------------------------------------
//                                       Conversions
//--------------------------------------------------------------------------------------------------

double_conversion!(Km, Decameters, 100.0, 0.0, 1.0);
double_conversion!(Km, Meters, 1_000.0, 0.0, 1.0);
double_conversion!(Km, Cm, 100_000.0, 0.0, 1.0);
double_conversion!(Km, Mm, 1_000_000.0, 0.0, 1.0);
double_conversion!(Km, StatuteMiles, 0.621_371, 0.0, 1.0);
double_conversion!(Km, Feet, 3_280.84, 0.0, 1.0);
double_conversion!(Km, Inches, 39_370.1, 0.0, 1.0);

double_conversion!(Decameters, Meters, 10.0, 0.0, 1.0);
double_conversion!(Decameters, Cm, 1_000.0, 0.0, 1.0);
double_conversion!(Decameters, Mm, 10_000.0, 0.0, 1.0);
double_conversion!(Decameters, StatuteMiles, 0.006_213_71, 0.0, 1.0);
double_conversion!(Decameters, Feet, 32.808_4, 0.0, 1.0);
double_conversion!(Decameters, Inches, 393.701, 0.0, 1.0);

double_conversion!(Meters, Cm, 100.0, 0.0, 1.0);
double_conversion!(Meters, Mm, 1_000.0, 0.0, 1.0);
double_conversion!(Meters, StatuteMiles, 6.213_71e-4, 0.0, 1.0);
double_conversion!(Meters, Feet, 3.280_84, 0.0, 1.0);
double_conversion!(Meters, Inches, 39.370_1, 0.0, 1.0);

double_conversion!(Cm, Mm, 10.0, 0.0, 1.0);
double_conversion!(Cm, StatuteMiles, 6.213_71e-6, 0.0, 1.0);
double_conversion!(Cm, Feet, 0.032_808_4, 0.0, 1.0);
double_conversion!(Cm, Inches, 0.393_701, 0.0, 1.0);

double_conversion!(Mm, StatuteMiles, 6.213_71e-7, 0.0, 1.0);
double_conversion!(Mm, Feet, 0.003_280_84, 0.0, 1.0);
double_conversion!(Mm, Inches, 0.039_370_1, 0.0, 1.0);

double_conversion!(StatuteMiles, Feet, 5_280.0, 0.0, 1.0);
double_conversion!(StatuteMiles, Inches, 63_360.0, 0.0, 1.0);

double_conversion!(Feet, Inches, 12.0, 0.0, 1.0);

//--------------------------------------------------------------------------------------------------
//                                         Tests
//--------------------------------------------------------------------------------------------------
#[cfg(test)]
mod tests {
    use super::*;
    use crate::test_utils::approx_equal;

    const TOL: f64 = 1.0e-7;

    #[test]
    fn basics() {
        assert!(approx_equal(Km(1.0), Decameters(100.0), Km(TOL)));
        assert!(approx_equal(Km(1.0), Meters(1_000.0), Km(TOL)));
        assert!(approx_equal(Km(1.0), Cm(100_000.0), Km(TOL)));
        assert!(approx_equal(Km(1.0), Mm(1_000_000.0), Km(TOL)));
        assert!(approx_equal(Km(1.0), StatuteMiles(0.621_371), Km(TOL)));
        assert!(approx_equal(Km(1.0), Feet(3_280.84), Km(TOL)));
        assert!(approx_equal(Km(1.0), Inches(39_370.1), Km(TOL)));
    }

    make_there_and_back_tests!(km_decameters, Km, Decameters, Km(0.0), Km(100.0), Km(TOL));
    make_there_and_back_tests!(km_meters, Km, Meters, Km(0.0), Km(100.0), Km(TOL));
    make_there_and_back_tests!(km_cm, Km, Cm, Km(0.0), Km(100.0), Km(TOL));
    make_there_and_back_tests!(km_mm, Km, Mm, Km(0.0), Km(100.0), Km(TOL));
    make_there_and_back_tests!(km_miles, Km, StatuteMiles, Km(0.0), Km(100.0), Km(TOL));
    make_there_and_back_tests!(km_feet, Km, Feet, Km(0.0), Km(100.0), Km(TOL));
    make_there_and_back_tests!(km_inches, Km, Inches, Km(0.0), Km(100.0), Km(TOL));

    make_there_and_back_tests!(
        dm_meters,
        Decameters,
        Meters,
        Decameters(0.0),
        Decameters(100.0),
        Decameters(TOL)
    );
    make_there_and_back_tests!(
        dm_cm,
        Decameters,
        Cm,
        Decameters(0.0),
        Decameters(100.0),
        Decameters(TOL)
    );
    make_there_and_back_tests!(
        dm_mm,
        Decameters,
        Mm,
        Decameters(0.0),
        Decameters(100.0),
        Decameters(TOL)
    );
    make_there_and_back_tests!(
        dm_miles,
        Decameters,
        StatuteMiles,
        Decameters(0.0),
        Decameters(100.0),
        Decameters(TOL)
    );
    make_there_and_back_tests!(
        dm_feet,
        Decameters,
        Feet,
        Decameters(0.0),
        Decameters(100.0),
        Decameters(TOL)
    );
    make_there_and_back_tests!(
        dm_inches,
        Decameters,
        Inches,
        Decameters(0.0),
        Decameters(100.0),
        Decameters(TOL)
    );

    make_there_and_back_tests!(m_cm, Meters, Cm, Meters(0.0), Meters(100.0), Meters(TOL));
    make_there_and_back_tests!(m_mm, Meters, Mm, Meters(0.0), Meters(100.0), Meters(TOL));
    make_there_and_back_tests!(
        m_miles,
        Meters,
        StatuteMiles,
        Meters(0.0),
        Meters(100.0),
        Meters(TOL)
    );
    make_there_and_back_tests!(
        m_feet,
        Meters,
        Feet,
        Meters(0.0),
        Meters(100.0),
        Meters(TOL)
    );
    make_there_and_back_tests!(
        m_inches,
        Meters,
        Inches,
        Meters(0.0),
        Meters(100.0),
        Meters(TOL)
    );

    make_there_and_back_tests!(cm_mm, Cm, Mm, Cm(0.0), Cm(100.0), Cm(TOL));
    make_there_and_back_tests!(cm_miles, Cm, StatuteMiles, Cm(0.0), Cm(100.0), Cm(TOL));
    make_there_and_back_tests!(cm_feet, Cm, Feet, Cm(0.0), Cm(100.0), Cm(TOL));
    make_there_and_back_tests!(cm_inches, Cm, Inches, Cm(0.0), Cm(100.0), Cm(TOL));

    make_there_and_back_tests!(mm_miles, Mm, StatuteMiles, Mm(0.0), Mm(100.0), Mm(TOL));
    make_there_and_back_tests!(mm_feet, Mm, Feet, Mm(0.0), Mm(100.0), Mm(TOL));
    make_there_and_back_tests!(mm_inches, Mm, Inches, Mm(0.0), Mm(100.0), Mm(TOL));

    make_there_and_back_tests!(
        sm_feet,
        StatuteMiles,
        Feet,
        StatuteMiles(0.0),
        StatuteMiles(100.0),
        StatuteMiles(TOL)
    );
    make_there_and_back_tests!(
        sm_inches,
        StatuteMiles,
        Inches,
        StatuteMiles(0.0),
        StatuteMiles(100.0),
        StatuteMiles(TOL)
    );

    make_there_and_back_tests!(feet_inches, Feet, Inches, Feet(0.0), Feet(100.0), Feet(TOL));
}
