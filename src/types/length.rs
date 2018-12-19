//! Length units for elevation and geopotential height.
use crate::types::Quantity;
use std::borrow::Borrow;
use std::cmp::Ordering;
use std::fmt::Display;
use std::ops::{Add, Sub};

/// Marker trait for elevation/height types.
pub trait Length: Quantity + PartialEq + PartialOrd {}

/// Length in meters, usually used for geopotential height and elevation.
#[derive(Clone, Copy, Debug)]
pub struct Meters(pub f64);

/// Length in decameters, usually used for geopotential height.
#[derive(Clone, Copy, Debug)]
pub struct Decameters(pub f64);

/// Length in feet, usually used for geopotential height and elevation.
#[derive(Clone, Copy, Debug)]
pub struct Feet(pub f64);

/// Length in kilometers, usually used for geopotential height above ground level.
#[derive(Clone, Copy, Debug)]
pub struct Km(pub f64);

/// Length in inches, usually used for precipitation depth or precipitable water.
#[derive(Clone, Copy, Debug)]
pub struct Inches(pub f64);

/// Length in millimeters, usually used for precipitation depth or precipitable water.
#[derive(Clone, Copy, Debug)]
pub struct Mm(pub f64);

/// Length in centimeters, usually used for precipitation depth or precipitable water.
#[derive(Clone, Copy, Debug)]
pub struct Cm(pub f64);

impl Length for Meters {}
impl Length for Decameters {}
impl Length for Feet {}
impl Length for Km {}
impl Length for Inches {}
impl Length for Mm {}
impl Length for Cm {}

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

        impl Borrow<f64> for $t {
            #[inline]
            fn borrow(&self) -> &f64 {
                &self.0
            }
        }

        implOpsForQuantity!($t);
    };
}

implQuantity!(Meters);
implQuantity!(Decameters);
implQuantity!(Feet);
implQuantity!(Km);
implQuantity!(Inches);
implQuantity!(Mm);
implQuantity!(Cm);

//--------------------------------------------------------------------------------------------------
//                                    Kilometer Conversions
//--------------------------------------------------------------------------------------------------
impl From<Feet> for Km {
    #[inline]
    fn from(h: Feet) -> Self {
        Km(h.0 / 3_280.8)
    }
}

impl From<Inches> for Km {
    #[inline]
    fn from(h: Inches) -> Self {
        Km::from(Feet::from(h))
    }
}

impl From<Decameters> for Km {
    #[inline]
    fn from(h: Decameters) -> Self {
        Km(h.0 * 10.0 / 1000.0)
    }
}

impl From<Meters> for Km {
    #[inline]
    fn from(h: Meters) -> Self {
        Km(h.0 / 1000.0)
    }
}

impl From<Cm> for Km {
    #[inline]
    fn from(h: Cm) -> Self {
        Km(h.0 / 100.0 / 1000.0)
    }
}

impl From<Mm> for Km {
    #[inline]
    fn from(h: Mm) -> Self {
        Km(h.0 / 1000.0 / 1000.0)
    }
}

//--------------------------------------------------------------------------------------------------
//                                    Decameter Conversions
//--------------------------------------------------------------------------------------------------
impl From<Feet> for Decameters {
    #[inline]
    fn from(h: Feet) -> Self {
        Decameters(h.0 / 32.808)
    }
}

impl From<Inches> for Decameters {
    #[inline]
    fn from(h: Inches) -> Self {
        Decameters::from(Feet::from(h))
    }
}

impl From<Km> for Decameters {
    #[inline]
    fn from(h: Km) -> Self {
        Decameters(h.0 * 1000.0 / 10.0)
    }
}

impl From<Meters> for Decameters {
    #[inline]
    fn from(h: Meters) -> Self {
        Decameters(h.0 / 10.0)
    }
}

impl From<Cm> for Decameters {
    #[inline]
    fn from(h: Cm) -> Self {
        Decameters(h.0 / 100.0 / 10.0)
    }
}

impl From<Mm> for Decameters {
    #[inline]
    fn from(h: Mm) -> Self {
        Decameters(h.0 / 1000.0 / 10.0)
    }
}

//--------------------------------------------------------------------------------------------------
//                                      Meter Conversions
//--------------------------------------------------------------------------------------------------
impl From<Feet> for Meters {
    #[inline]
    fn from(h: Feet) -> Self {
        Meters(h.0 / 3.2808)
    }
}

impl From<Inches> for Meters {
    #[inline]
    fn from(h: Inches) -> Self {
        Meters::from(Feet::from(h))
    }
}

impl From<Km> for Meters {
    #[inline]
    fn from(h: Km) -> Self {
        Meters(h.0 * 1000.0)
    }
}

impl From<Decameters> for Meters {
    #[inline]
    fn from(h: Decameters) -> Self {
        Meters(h.0 * 10.0)
    }
}

impl From<Cm> for Meters {
    #[inline]
    fn from(h: Cm) -> Self {
        Meters(h.0 / 100.0)
    }
}

impl From<Mm> for Meters {
    #[inline]
    fn from(h: Mm) -> Self {
        Meters(h.0 / 1000.0)
    }
}

//--------------------------------------------------------------------------------------------------
//                                    Centimeter Conversions
//--------------------------------------------------------------------------------------------------
impl From<Feet> for Cm {
    #[inline]
    fn from(h: Feet) -> Self {
        Cm::from(Meters::from(h))
    }
}

impl From<Inches> for Cm {
    #[inline]
    fn from(h: Inches) -> Self {
        Cm::from(Meters::from(h))
    }
}

impl From<Km> for Cm {
    #[inline]
    fn from(h: Km) -> Self {
        Cm::from(Meters::from(h))
    }
}

impl From<Decameters> for Cm {
    #[inline]
    fn from(h: Decameters) -> Self {
        Cm::from(Meters::from(h))
    }
}

impl From<Meters> for Cm {
    #[inline]
    fn from(h: Meters) -> Self {
        Cm(h.0 * 100.0)
    }
}

impl From<Mm> for Cm {
    #[inline]
    fn from(h: Mm) -> Self {
        Cm::from(Meters::from(h))
    }
}

//--------------------------------------------------------------------------------------------------
//                                    Millimeter Conversions
//--------------------------------------------------------------------------------------------------
impl From<Feet> for Mm {
    #[inline]
    fn from(h: Feet) -> Self {
        Mm::from(Meters::from(h))
    }
}

impl From<Inches> for Mm {
    #[inline]
    fn from(h: Inches) -> Self {
        Mm::from(Meters::from(h))
    }
}

impl From<Km> for Mm {
    #[inline]
    fn from(h: Km) -> Self {
        Mm::from(Meters::from(h))
    }
}

impl From<Decameters> for Mm {
    #[inline]
    fn from(h: Decameters) -> Self {
        Mm::from(Meters::from(h))
    }
}

impl From<Meters> for Mm {
    #[inline]
    fn from(h: Meters) -> Self {
        Mm(h.0 * 1000.0)
    }
}

impl From<Cm> for Mm {
    #[inline]
    fn from(h: Cm) -> Self {
        Mm::from(Meters::from(h))
    }
}

//--------------------------------------------------------------------------------------------------
//                                       Feet Conversions
//--------------------------------------------------------------------------------------------------
impl From<Inches> for Feet {
    #[inline]
    fn from(h: Inches) -> Self {
        Feet(h.0 / 12.0)
    }
}

impl From<Km> for Feet {
    #[inline]
    fn from(h: Km) -> Self {
        Feet(h.0 * 3_280.8)
    }
}

impl From<Decameters> for Feet {
    #[inline]
    fn from(h: Decameters) -> Self {
        Feet(h.0 / 10.0 * 3.2808)
    }
}

impl From<Meters> for Feet {
    #[inline]
    fn from(h: Meters) -> Self {
        Feet(h.0 * 3.2808)
    }
}

impl From<Cm> for Feet {
    #[inline]
    fn from(h: Cm) -> Self {
        Feet(h.0 / 100.0 * 3.2808)
    }
}

impl From<Mm> for Feet {
    #[inline]
    fn from(h: Mm) -> Self {
        Feet(h.0 / 1000.0 * 3.2808)
    }
}

//--------------------------------------------------------------------------------------------------
//                                       Inch Conversions
//--------------------------------------------------------------------------------------------------
impl From<Feet> for Inches {
    #[inline]
    fn from(h: Feet) -> Self {
        Inches(h.0 * 12.0)
    }
}

impl From<Km> for Inches {
    #[inline]
    fn from(h: Km) -> Self {
        Inches::from(Feet::from(h))
    }
}

impl From<Decameters> for Inches {
    #[inline]
    fn from(h: Decameters) -> Self {
        Inches::from(Feet::from(h))
    }
}

impl From<Meters> for Inches {
    #[inline]
    fn from(h: Meters) -> Self {
        Inches::from(Feet::from(h))
    }
}

impl From<Cm> for Inches {
    #[inline]
    fn from(h: Cm) -> Self {
        Inches::from(Feet::from(h))
    }
}

impl From<Mm> for Inches {
    #[inline]
    fn from(h: Mm) -> Self {
        Inches::from(Feet::from(h))
    }
}

//--------------------------------------------------------------------------------------------------
//                                     Display implementations
//--------------------------------------------------------------------------------------------------

impl Display for Meters {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.0}m", self.0)
    }
}

impl Display for Decameters {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.0}dm", self.0)
    }
}

impl Display for Feet {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.0}ft", self.0)
    }
}

impl Display for Km {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.2}km", self.0)
    }
}

impl Display for Inches {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.2}in", self.0)
    }
}

impl Display for Cm {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.1}cm", self.0)
    }
}

impl Display for Mm {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:.0}mm", self.0)
    }
}
