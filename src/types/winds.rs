//! Wind units and vectors.
use crate::types::{speed::Speed, Quantity, VectorQuantity};
use std::fmt::Display;

/// Marker trait for Wind types.
pub trait Wind<S>: VectorQuantity<S>
where
    S: Speed,
{
}

/*--------------------------------------------------------------------------------------------------
                        Wind as speed in knots and the direction it is coming from.
--------------------------------------------------------------------------------------------------*/
/// Wind direction and speed in knots.
#[allow(missing_docs)]
#[derive(Clone, Copy, Debug)]
pub struct WindSpdDir<S>
where
    S: Quantity,
{
    pub speed: S,
    pub direction: f64,
}

impl<S> VectorQuantity<S> for WindSpdDir<S>
where
    S: Quantity + PartialOrd,
{
    #[inline]
    fn pack_xy(vals: (S, S)) -> Self {
        let (direction, speed) = from_cart_to_wind(vals.0, vals.1);

        WindSpdDir { speed, direction }
    }

    #[inline]
    fn unpack_xy(self) -> (S, S) {
        from_wind_to_cart(self.direction, self.speed)
    }

    #[inline]
    fn unwrap_xy(self) -> (S, S) {
        if self.speed < S::pack(0.0) {
            panic!("Speed cannot be less than 0.0!");
        }

        if self.direction < 0.0 || self.direction > 360.0 {
            panic!("Wind direction not in 0 - 360 range.")
        }

        self.unpack_xy()
    }

    #[inline]
    fn abs(&self) -> S {
        self.speed
    }

    #[inline]
    fn into_option(self) -> Option<(S, S)> {
        if self.speed < S::pack(0.0) || self.direction < 0.0 || self.direction > 360.0 {
            None
        } else {
            Some(self.unpack_xy())
        }
    }
}

#[cfg(feature = "use_optional")]
impl<S> optional::Noned for WindSpdDir<S>
where
    S: Quantity + optional::Noned,
{
    #[inline]
    fn is_none(&self) -> bool {
        optional::Noned::is_none(&self.speed) || optional::Noned::is_none(&self.direction)
    }

    #[inline]
    fn get_none() -> Self {
        Self {
            speed: optional::Noned::get_none(),
            direction: optional::Noned::get_none(),
        }
    }
}

impl<S> Display for WindSpdDir<S>
where
    S: Quantity,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:03.0} at {}", self.direction, self.speed)
    }
}

impl<S> Wind<S> for WindSpdDir<S> where S: Speed {}
implOpsForVectorQuantity!(WindSpdDir);

/*--------------------------------------------------------------------------------------------------
                                   Wind as U and V components.
--------------------------------------------------------------------------------------------------*/
/// Wind in U and V components in m/s.
#[allow(missing_docs)]
#[derive(Clone, Copy, Debug)]
pub struct WindUV<S> {
    pub u: S,
    pub v: S,
}

impl<S> VectorQuantity<S> for WindUV<S>
where
    S: Quantity,
{
    #[inline]
    fn pack_xy(vals: (S, S)) -> Self {
        WindUV {
            u: vals.0,
            v: vals.1,
        }
    }

    #[inline]
    fn unpack_xy(self) -> (S, S) {
        (self.u, self.v)
    }

    #[inline]
    fn unwrap_xy(self) -> (S, S) {
        (self.u, self.v)
    }

    #[inline]
    fn into_option(self) -> Option<(S, S)> {
        Some((self.u, self.v))
    }

    #[inline]
    fn abs(&self) -> S {
        let x = self.u.unpack();
        let y = self.v.unpack();
        S::pack(x.hypot(y))
    }
}

#[cfg(feature = "use_optional")]
impl<S> optional::Noned for WindUV<S>
where
    S: Quantity + optional::Noned,
{
    #[inline]
    fn is_none(&self) -> bool {
        optional::Noned::is_none(&self.u) || optional::Noned::is_none(&self.v)
    }

    #[inline]
    fn get_none() -> Self {
        Self::pack_xy((optional::Noned::get_none(), optional::Noned::get_none()))
    }
}

impl<S: Display> Display for WindUV<S> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "({}, {})", self.u, self.v)
    }
}

impl<S: Speed> Wind<S> for WindUV<S> {}
implOpsForVectorQuantity!(WindUV);

/*--------------------------------------------------------------------------------------------------
                                         Wind conversions.
--------------------------------------------------------------------------------------------------*/
impl<S1, S2> From<WindUV<S1>> for WindSpdDir<S2>
where
    S1: Speed,
    S2: From<S1>,
    S2: Speed,
{
    #[inline]
    fn from(wind: WindUV<S1>) -> Self {
        let (direction, speed) = from_cart_to_wind(wind.u, wind.v);
        let speed = S2::from(speed);

        WindSpdDir { direction, speed }
    }
}

impl<S1, S2> From<WindSpdDir<S1>> for WindUV<S2>
where
    S1: Speed,
    S2: From<S1>,
    S2: Speed,
{
    #[inline]
    fn from(wind: WindSpdDir<S1>) -> Self {
        let speed = S2::from(wind.speed);
        let (u, v) = from_wind_to_cart(wind.direction, speed);

        WindUV { u, v }
    }
}

/// Convert from standard cartesian coordinate vectors to meteorological wind from direction at
/// speed coordinates.
///
/// Returns a tuple (direction, speed)
#[inline]
fn from_cart_to_wind<S: Quantity>(x: S, y: S) -> (f64, S) {
    let spd = (x.unpack().powi(2) + y.unpack().powi(2)).sqrt();

    let mut direction = 180.0 + 90.0 - y.unpack().atan2(x.unpack()).to_degrees();
    while direction < 0.0 {
        direction += 360.0;
    }
    while direction >= 360.0 {
        direction -= 360.0;
    }

    (direction, S::pack(spd))
}

/// Convert from meteorological wind from direction at speed coordinates to standard x-y cartesian
/// coordinates.
///
/// Returns a tuple (x_velocity, y_velocity)
#[inline]
fn from_wind_to_cart<S: Quantity>(dir: f64, spd: S) -> (S, S) {
    let rads = dir.to_radians();

    let u = S::pack(-spd.unpack() * rads.sin());
    let v = S::pack(-spd.unpack() * rads.cos());

    (u, v)
}

/*--------------------------------------------------------------------------------------------------
                                              Unit Tests.
--------------------------------------------------------------------------------------------------*/
#[cfg(test)]
mod test {
    use crate::*;

    const TOL: f64 = 1.0e-6;

    macro_rules! assert_approx_eq {
        ($a:expr, $b:expr, $tol:expr, $msg:expr) => {{
            println!(
                "{:?} == {:?} with tolerance <= {:?}, required tolerance {:?}",
                $a,
                $b,
                ($a - $b).abs(),
                $tol
            );
            assert!($a.approx_eq($b, $tol), $msg);
        }};
    }
    fn get_test_winds() -> [((f64, Knots), (MetersPSec, MetersPSec)); 8] {
        [
            (
                (000.0, Knots(10.0)),
                (MetersPSec(0.0), MetersPSec(-5.14444444)),
            ),
            (
                (045.0, Knots(10.0)),
                (MetersPSec(-3.637671549), MetersPSec(-3.637671549)),
            ),
            (
                (090.0, Knots(10.0)),
                (MetersPSec(-5.14444444), MetersPSec(0.0)),
            ),
            (
                (135.0, Knots(10.0)),
                (MetersPSec(-3.637671549), MetersPSec(3.637671549)),
            ),
            (
                (180.0, Knots(10.0)),
                (MetersPSec(0.0), MetersPSec(5.14444444)),
            ),
            (
                (225.0, Knots(10.0)),
                (MetersPSec(3.637671549), MetersPSec(3.637671549)),
            ),
            (
                (270.0, Knots(10.0)),
                (MetersPSec(5.14444444), MetersPSec(0.0)),
            ),
            (
                (315.0, Knots(10.0)),
                (MetersPSec(3.637671549), MetersPSec(-3.637671549)),
            ),
        ]
    }

    #[test]
    fn test_spd_dir_to_uv() {
        let spd_dir_to_uv_data = get_test_winds();

        for &((dir, spd), (u, v)) in &spd_dir_to_uv_data {
            let dir_spd = WindSpdDir {
                speed: spd,
                direction: dir,
            };
            let WindUV::<MetersPSec> {
                u: calc_u,
                v: calc_v,
            } = WindUV::from(dir_spd);

            assert_approx_eq!(calc_u, u, MetersPSec(TOL), "U speed mismatch.");
            assert_approx_eq!(calc_v, v, MetersPSec(TOL), "V speed mismatch.");

            let wind_uv = WindUV::<MetersPSec>::from(dir_spd);
            assert_approx_eq!(dir_spd, wind_uv, MetersPSec(TOL), "Vector mismatch.")
        }
    }

    #[test]
    fn test_uv_to_spd_dir() {
        let spd_dir_to_uv_data = get_test_winds();

        for &((dir, spd), (u, v)) in &spd_dir_to_uv_data {
            let uv_wind = WindUV::<MetersPSec> { u, v };
            let WindSpdDir::<Knots> {
                speed: calc_spd,
                direction: calc_dir,
            } = WindSpdDir::from(uv_wind);

            assert_approx_eq!(calc_dir, dir, TOL, "Direction mismatch.");
            assert_approx_eq!(calc_spd, spd, Knots(TOL), "Speed mismatch.");

            let dir_spd = WindSpdDir::<Knots>::from(uv_wind);
            assert_approx_eq!(dir_spd, uv_wind, Knots(TOL), "Vector mismatch.")
        }
    }

}
