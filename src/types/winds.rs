//! Wind units and vectors.
use crate::types::VectorQuantity;
use std::fmt::Display;

/// Marker trait for Wind types.
pub trait Wind: VectorQuantity {}

/// Wind direction and speed in knots.
#[allow(missing_docs)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindSpdDir {
    pub speed: f64,
    pub direction: f64,
}

impl VectorQuantity for WindSpdDir {
    #[inline]
    fn pack(vals: (f64, f64)) -> Self {
        WindSpdDir {
            speed: vals.0,
            direction: vals.1,
        }
    }

    #[inline]
    fn unpack(self) -> (f64, f64) {
        (self.speed, self.direction)
    }

    #[inline]
    fn unwrap(self) -> (f64, f64) {
        if self.speed < 0.0 {
            panic!("Speed cannot be less than 0.0!");
        }

        if self.direction < 0.0 || self.direction > 360.0 {
            panic!("Wind direction not in 0 - 360 range.")
        }

        (self.speed, self.direction)
    }

    #[inline]
    fn into_option(self) -> Option<(f64, f64)> {
        if self.speed < 0.0 || self.direction < 0.0 || self.direction > 360.0 {
            None
        } else {
            Some((self.speed, self.direction))
        }
    }
}

#[cfg(features = "use_optional")]
impl optional::Noned for WindSpdDir {
    #[inline]
    fn is_none(&self) -> bool {
        optional::Noned::is_none(self.speed) || optional::Noned::is_none(self.direction)
    }

    #[inline]
    fn get_none() -> Self {
        Self::pack((optional::Noned::get_none(), optional::Noned::get_none()))
    }
}

impl Wind for WindSpdDir {}

/// Wind in U and V components in m/s.
#[allow(missing_docs)]
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WindUV {
    pub u: f64,
    pub v: f64,
}

impl VectorQuantity for WindUV {
    #[inline]
    fn pack(vals: (f64, f64)) -> Self {
        WindUV {
            u: vals.0,
            v: vals.1,
        }
    }

    #[inline]
    fn unpack(self) -> (f64, f64) {
        (self.u, self.v)
    }

    #[inline]
    fn unwrap(self) -> (f64, f64) {
        (self.u, self.v)
    }

    #[inline]
    fn into_option(self) -> Option<(f64, f64)> {
        Some((self.u, self.v))
    }
}

#[cfg(features = "use_optional")]
impl optional::Noned for WindUV {
    #[inline]
    fn is_none(&self) -> bool {
        optional::Noned::is_none(self.u) || optional::Noned::is_none(self.v)
    }

    #[inline]
    fn get_none() -> Self {
        Self::pack((optional::Noned::get_none(), optional::Noned::get_none()))
    }
}

impl Wind for WindUV {}

impl From<WindUV> for WindSpdDir {
    #[inline]
    fn from(wind: WindUV) -> Self {
        let spd_ms = (wind.u.powi(2) + wind.v.powi(2)).sqrt();
        let speed = mps_to_knots(spd_ms);

        let mut direction = 180.0 + 90.0 - wind.v.atan2(wind.u).to_degrees();
        while direction < 0.0 {
            direction += 360.0;
        }
        while direction >= 360.0 {
            direction -= 360.0;
        }

        WindSpdDir { direction, speed }
    }
}

impl From<WindSpdDir> for WindUV {
    #[inline]
    fn from(wind: WindSpdDir) -> Self {
        let rads = wind.direction.to_radians();
        let spd_ms = knots_to_mps(wind.speed);

        let u = -spd_ms * rads.sin();
        let v = -spd_ms * rads.cos();

        WindUV { u, v }
    }
}

impl Display for WindSpdDir {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "{:03.0} at {:02.0}kt", self.direction, self.speed)
    }
}

impl Display for WindUV {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::result::Result<(), std::fmt::Error> {
        write!(f, "({:.1}m/s, {:.1}m/s)", self.u, self.v)
    }
}

/// Convert knots to m/s.
#[inline]
pub fn knots_to_mps(spd: f64) -> f64 {
    spd * 0.514_444_444
}

/// Convert m/s to knots.
#[inline]
pub fn mps_to_knots(spd: f64) -> f64 {
    spd * 1.943_844_5
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::test_utils::*;

    const TOL: f64 = 1.0e-6;

    macro_rules! assert_approx_eq {
        ($a:expr, $b:expr) => {{
            println!("{} == {} with tolerance <= {}", $a, $b, f64::abs($a - $b));
            assert!(approx_equal($a, $b, TOL));
        }};
        ($a:expr, $b:expr, $msg:expr) => {{
            println!("{} == {} with tolerance <= {}", $a, $b, f64::abs($a - $b));
            assert!(approx_equal($a, $b, TOL), $msg);
        }};
        ($a:expr, $b:expr, $tol:expr, $msg:expr) => {{
            println!("{} == {} with tolerance <= {}", $a, $b, f64::abs($a - $b));
            assert!(approx_equal($a, $b, $tol), $msg);
        }};
    }
    fn get_test_winds() -> [((f64, f64), (f64, f64)); 8] {
        [
            ((000.0, 10.0), (0.0, -5.14444444)),
            ((045.0, 10.0), (-3.637671549, -3.637671549)),
            ((090.0, 10.0), (-5.14444444, 0.0)),
            ((135.0, 10.0), (-3.637671549, 3.637671549)),
            ((180.0, 10.0), (0.0, 5.14444444)),
            ((225.0, 10.0), (3.637671549, 3.637671549)),
            ((270.0, 10.0), (5.14444444, 0.0)),
            ((315.0, 10.0), (3.637671549, -3.637671549)),
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
            let WindUV {
                u: calc_u,
                v: calc_v,
            } = WindUV::from(dir_spd);

            assert_approx_eq!(calc_u, u, TOL, "U speed mismatch.");
            assert_approx_eq!(calc_v, v, TOL, "V speed mismatch.");
        }
    }

    #[test]
    fn test_uv_to_spd_dir() {
        let spd_dir_to_uv_data = get_test_winds();

        for &((dir, spd), (u, v)) in &spd_dir_to_uv_data {
            let uv_wind = WindUV { u, v };
            let WindSpdDir {
                speed: calc_spd,
                direction: calc_dir,
            } = WindSpdDir::from(uv_wind);

            assert_approx_eq!(calc_dir, dir, TOL, "Direction mismatch.");
            assert_approx_eq!(calc_spd, spd, TOL, "Speed mismatch.");
        }
    }

}
