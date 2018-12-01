use constants::*;
use error::MetForErr::*;
use error::*;

/// Convert temperature from Celsius to Kelvin.
#[inline(always)]
pub fn celsius_to_kelvin(celsius: f64) -> Result<f64> {
    if celsius < ABSOLUTE_ZERO_C {
        Err(BelowAbsoluteZero)
    } else {
        Ok(celsius - ABSOLUTE_ZERO_C)
    }
}

/// Convert Kelvin temperatures to Celsius.
#[inline(always)]
pub fn kelvin_to_celsius(kelvin: f64) -> Result<f64> {
    if kelvin < ABSOLUTE_ZERO_K {
        Err(BelowAbsoluteZero)
    } else {
        Ok(kelvin + ABSOLUTE_ZERO_C)
    }
}

/// Convert Celsius to Fahrenheit.
#[inline(always)]
pub fn celsius_to_f(temperature: f64) -> Result<f64> {
    if temperature < ABSOLUTE_ZERO_C {
        Err(BelowAbsoluteZero)
    } else {
        Ok(1.8 * temperature + 32.0)
    }
}

/// Convert Fahrenheit to Celsius
#[inline(always)]
pub fn f_to_celsius(temperature: f64) -> Result<f64> {
    if temperature < ABSOLUTE_ZERO_F {
        Err(BelowAbsoluteZero)
    } else {
        Ok((temperature - 32.0) / 1.8)
    }
}

/// Calculate potential temperature assuming a 1000 hPa reference level.
///
/// * `pressure_hpa` - Pressure in hPa
/// * `temperature_c` - Temperature in Celsius.
///
/// Returns: Potential temperature in Kelvin
#[inline]
pub fn theta_kelvin(pressure_hpa: f64, temperature_c: f64) -> Result<f64> {
    use std::f64;

    if pressure_hpa < 0.0 {
        Err(NegativePressure)
    } else if temperature_c < ABSOLUTE_ZERO_C {
        Err(BelowAbsoluteZero)
    } else {
        Ok(celsius_to_kelvin(temperature_c)? * f64::powf(1000.0 / pressure_hpa, Rd / cp))
    }
}

/// Given a potential temperature and pressure, calculate the temperature.
///
/// * `theta_kelvin` - Potential temperature in Kelvin.
/// * `pressure_hpa` - Pressure in hPa
///
/// Returns: The temperature in Celsius.
#[inline]
pub fn temperature_c_from_theta(theta_kelvin: f64, pressure_hpa: f64) -> Result<f64> {
    use std::f64;
    if pressure_hpa < 0.0 {
        Err(NegativePressure)
    } else if theta_kelvin < ABSOLUTE_ZERO_K {
        Err(BelowAbsoluteZero)
    } else {
        Ok(kelvin_to_celsius(
            theta_kelvin * f64::powf(pressure_hpa / 1000.0, Rd / cp),
        )?)
    }
}

// Constants used for the limits of applicability to the empirical relationships used for
// vapor pressure.
const MIN_T_VP_ICE: f64 = -80.0;
const MIN_T_VP_LIQUID: f64 = -80.0;
const MAX_T_VP_ICE: f64 = 0.0;
const MAX_T_VP_LIQUID: f64 = 50.0;

/// Get the vapor pressure over liquid water.
///
/// Alduchov, O.A., and Eskridge, R.E. Improved Magnus` form approximation of saturation vapor
/// pressure. United States: N. p., 1997. Web. doi:10.2172/548871.
///
///  * `dew_point_c` - the dew point in Celsius. If saturation is assumed this is also the
///                    temperature.
///
/// Returns: The vapor pressure of water vapor in hPa.
#[inline]
pub fn vapor_pressure_liquid_water(dew_point_c: f64) -> Result<f64> {
    if dew_point_c < MIN_T_VP_LIQUID || dew_point_c > MAX_T_VP_LIQUID {
        Err(InputOutOfRange)
    } else {
        Ok(6.1037 * f64::exp(17.641 * dew_point_c / (dew_point_c + 243.27)))
    }
}

/// Get the dew point given the vapor pressure of water over liquid water. This function is the
/// inverse of `vapor_pressure_liquid_water`.
///
/// `vp_hpa` - The vapor pressure of water in hPa.
///
/// Returns: The dew point in Celsius.
#[inline]
fn dew_point_from_vapor_pressure_over_liquid(vp_hpa: f64) -> Result<f64> {
    let a = f64::ln(vp_hpa / 6.1037) / 17.641;
    let dp = a * 243.27 / (1.0 - a);
    if dp < MIN_T_VP_LIQUID || dp > MAX_T_VP_LIQUID {
        Err(InputOutOfRange)
    } else {
        Ok(dp)
    }
}

/// Get the vapor pressure over ice.
///
/// Alduchov, O.A., and Eskridge, R.E. Improved Magnus` form approximation of saturation vapor
/// pressure. United States: N. p., 1997. Web. doi:10.2172/548871.
///
///  * `frost_point_c` - the frost point in Celsius. If saturation is assumed this is also the
///                      temperature.
///
/// Returns: The vapor pressure of water vapor in hPa.
#[inline]
pub fn vapor_pressure_ice(frost_point_c: f64) -> Result<f64> {
    if frost_point_c < MIN_T_VP_ICE || frost_point_c > MAX_T_VP_ICE {
        Err(InputOutOfRange)
    } else {
        Ok(6.1121 * f64::exp(22.587 * frost_point_c / (frost_point_c + 273.86)))
    }
}

/// Get the frost point given the vapor pressure of water over ice. This function is the inverse of
/// `vapor_pressure_ice`.
///
///
/// `vp_hpa` - The vapor pressure of water in hPa.
///
/// Returns: The frost point in Celsius.
#[inline]
pub fn frost_point_from_vapor_pressure_over_ice(vp_hpa: f64) -> Result<f64> {
    let a = f64::ln(vp_hpa / 6.1121) / 22.587;
    let fp = a * 273.86 / (1.0 - a);
    if fp < MIN_T_VP_ICE || fp > MAX_T_VP_ICE {
        Err(InputOutOfRange)
    } else {
        Ok(fp)
    }
}

/// Calculate the relative humidity with respect to liquid water.
///
/// * `temperature_c` - The temperature in Celsius.
/// * `dew_point_c` - The dew point in Celsius.
///
/// Returns: The relative humidity as a decimal, i.e. 0.95 instead of 95%.
#[inline]
pub fn rh(temperature_c: f64, dew_point_c: f64) -> Result<f64> {
    let (es, e) = (
        vapor_pressure_liquid_water(temperature_c)?,
        vapor_pressure_liquid_water(dew_point_c)?,
    );

    // Allow e > es for supersaturation.
    Ok(e / es)
}

/// Calculate the relative humidity with respect to ice.
///
/// * `temperature_c` - The temperature in Celsius.
/// * `frost_point_c` - The frost point in Celsius.
///
/// Returns: The relative humidity as a decimal, i.e. 0.95 instead of 95%.
#[inline]
pub fn rh_ice(temperature_c: f64, frost_point_c: f64) -> Result<f64> {
    let (es, e) = (
        vapor_pressure_ice(temperature_c)?,
        vapor_pressure_ice(frost_point_c)?,
    );

    // Allow e > es for supersaturation.
    Ok(e / es)
}

/// Calculate the mixing ratio of water.
///
/// Eq 5.9 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// * `dew_point_c` - the dew point in Celsius. If you want the saturation mixing ratio, just make
///                   this equal to the temperature.
/// * `pressure_hpa` - the pressure in hPa.
///
/// Returns: The mixing ratio as a unitless value. Note this is often reported as g/kg.
#[inline]
pub fn mixing_ratio(dew_point_c: f64, pressure_hpa: f64) -> Result<f64> {
    let vp = vapor_pressure_liquid_water(dew_point_c)?;
    if vp > pressure_hpa {
        Err(VaporPressureTooHigh)
    } else {
        Ok(epsilon * vp / (pressure_hpa - vp))
    }
}

/// Given a mixing ratio and pressure, calculate the dew point temperature. If saturation is
/// assumed this is also the temperature.
///
/// * `pressure_hpa` - the pressure in hPa
/// * `mw` - the mixing ratio.
///
/// Returns: The dew point in Celsius
#[inline]
pub fn dew_point_from_p_and_mw(pressure_hpa: f64, mw: f64) -> Result<f64> {
    let vp = mw * pressure_hpa / (mw + epsilon);
    dew_point_from_vapor_pressure_over_liquid(vp)
}

/// Calculate the temperature and pressure at the lifting condensation level (LCL).
///
/// Eqs 5.17 and 5.18 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// * `temperature_c` - the initial temperature in Celsius.
/// * `dew_point_c` - the initial dew point in Celsius.
/// * `pressure_hpa` - the initial pressure of the parcel in hPa.
///
/// Returns: The pressure at the LCL in hPa and the temperature at the LCL in Kelvin.
#[inline]
pub fn pressure_and_temperature_at_lcl(
    temperature_c: f64,
    dew_point_c: f64,
    pressure_hpa: f64,
) -> Result<(f64, f64)> {
    let plcl = pressure_hpa_at_lcl(temperature_c, dew_point_c, pressure_hpa)?;
    let tlcl = temperature_c_from_theta(theta_kelvin(pressure_hpa, temperature_c)?, plcl)?;
    let tlcl = celsius_to_kelvin(tlcl)?;

    Ok((plcl, tlcl))
}

/// Approximate pressure of the Lifting Condensation Level (LCL).
///
/// * `temperature_c` - the initial temperature in Celsius.
/// * `dew_point_c` - the initial dew point in Celsius.
/// * `pressure_hpa` - the initial pressure of the parcel in hPa.
///
/// Returns: The pressure at the LCL in hPa.
#[inline]
pub fn pressure_hpa_at_lcl(temperature_c: f64, dew_point_c: f64, pressure_hpa: f64) -> Result<f64> {
    if dew_point_c < ABSOLUTE_ZERO_C || temperature_c < ABSOLUTE_ZERO_C {
        Err(BelowAbsoluteZero)
    } else if dew_point_c >= temperature_c {
        Ok(pressure_hpa)
    } else {
        let theta = theta_kelvin(pressure_hpa, temperature_c)?;
        let mw = mixing_ratio(dew_point_c, pressure_hpa)?;
        let theta_e = theta_e_kelvin(temperature_c, dew_point_c, pressure_hpa)?;

        let eq = |p| -> Result<f64> {
            let t = temperature_c_from_theta(theta, p)?;
            let dp = dew_point_from_p_and_mw(p, mw)?;

            Ok(t - dp)
        };

        // Search between 1060 and 300 hPa. If it falls outside that range, give up!
        let first_guess = find_root(&eq, 300.0, 1060.0)?;

        let eq = |p| -> Result<f64> {
            let t1 = temperature_c_from_theta(theta, p)?;
            let t2 = temperature_c_from_theta_e_saturated_and_pressure(p, theta_e)?;

            Ok(t1 - t2)
        };

        const DELTA: f64 = 10.0;
        let lclp = if let Ok(lcl) = find_root(&eq, first_guess + DELTA, first_guess - DELTA) {
            lcl
        } else {
            first_guess
        };

        Ok(lclp)
    }
}

/// Calculate the specific humidity.
///
/// Eqs 5.11 and 5.12 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// * `dew_point_c` - the dew point in Celsius, if this is the same as the temperature then this
///                   calculates the saturation specific humidity.
/// * `pressure_hpa` - the pressure in hPa.
///
/// Returns the specific humidity. (no units)
#[inline]
pub fn specific_humidity(dew_point_c: f64, pressure_hpa: f64) -> Result<f64> {
    let vp = vapor_pressure_liquid_water(dew_point_c)?;

    if pressure_hpa < 0.0 {
        Err(NegativePressure)
    } else if vp > pressure_hpa {
        Err(VaporPressureTooHigh)
    } else {
        Ok(vp / pressure_hpa * epsilon)
    }
}

/// Calculate equivalent potential temperature.
///
/// Equation from ["The Glossary of Meteorology"]
/// (http://glossary.ametsoc.org/wiki/Equivalent_potential_temperature) online where the
/// approximation of ignoring the "total water mixing ratio" is used since most of the time we do
/// not have the necessary information to calculate that.
///
/// * `temperature_c` - the initial temperature in Celsius.
/// * `dew_point_c` - the initial dew point in Celsius.
/// * `pressure_hpa` - the initial pressure in hPa.
///
/// Returns: The equivalent potential temperature in Kelvin.
#[inline]
pub fn theta_e_kelvin(temperature_c: f64, dew_point_c: f64, pressure_hpa: f64) -> Result<f64> {
    let tk = celsius_to_kelvin(temperature_c)?;
    const P0: f64 = 1000.0; // Reference pressure for potential temperatures.
    let rv = mixing_ratio(dew_point_c, pressure_hpa)?;
    let pd = pressure_hpa - vapor_pressure_liquid_water(dew_point_c)?;

    if pd < 0.0 {
        return Err(VaporPressureTooHigh);
    }

    let h = rh(temperature_c, dew_point_c)?;

    let lv = latent_heat_of_condensation(temperature_c)?;

    Ok(
        tk * f64::powf(P0 / pd, Rd / cp)
            * f64::powf(h, -rv * Rv / cp)
            * f64::exp(lv * rv / cp / tk),
    )
}

/// Given the pressure and equivalent potential temperature, assume saturation and calculate the
/// temperature.
///
/// This function is useful if you were trying to calculate the temperature for plotting a
/// moist adiabat on a skew-t log-p diagram.
///
/// * `pressure_hpa` - the initial pressure in hPa.
/// * `theta_e_k` - the equivalent potential temperature in Kelvin.
///
/// Returns: The temperature in Celsius.
#[inline]
pub fn temperature_c_from_theta_e_saturated_and_pressure(
    pressure_hpa: f64,
    theta_e_k: f64,
) -> Result<f64> {
    find_root(
        &|t_c| Ok(theta_e_kelvin(t_c, t_c, pressure_hpa)? - theta_e_k),
        -80.0,
        50.0,
    )
}

/// Calculate the web bulb temperature.
///
/// * `temperature_c` - the temperature in Celsius.
/// * `dew_point_c` - the dew point in Celsius.
/// * `pressure_hpa` - the pressure in hPa
///
/// Returns: The wet bulb temperature in Celsius.
#[inline]
pub fn wet_bulb_c(temperature_c: f64, dew_point_c: f64, pressure_hpa: f64) -> Result<f64> {
    let (p_lcl, t_lcl) = pressure_and_temperature_at_lcl(temperature_c, dew_point_c, pressure_hpa)?;
    let t_lcl = kelvin_to_celsius(t_lcl)?;
    let theta_e = theta_e_kelvin(t_lcl, t_lcl, p_lcl)?;

    find_root(
        &|celsius| Ok(theta_e_kelvin(celsius, celsius, pressure_hpa)? - theta_e),
        dew_point_c,
        temperature_c,
    )
}

/// Latent heat of condensation for water.
///
/// Polynomial curve fit to Table 2.1. R. R. Rogers; M. K. Yau (1989). A Short Course in Cloud
/// Physics (3rd ed.). Pergamon Press. p. 16. ISBN 0-7506-3215-1.
///
/// `temperature_c` - temperature in Celsius.
///
/// Returns: the latent heat of condensation for water in J kg<sup>-1</sup>.
#[inline]
pub fn latent_heat_of_condensation(temperature_c: f64) -> Result<f64> {
    // The table has values from -40.0 to 40.0. So from -100.0 to -40.0 is actually an exrapolation.
    // I graphed the values from the extrapolation, and the curve looks good, and is approaching the
    // latent heat of sublimation, but does not exceed it. This seems very reasonable to me,
    // especially considering that a common approximation is to just us a constant value.
    if temperature_c < -100.0 || temperature_c > 60.0 {
        Err(InputOutOfRange)
    } else {
        let t = temperature_c;
        Ok((2500.8 - 2.36 * t + 0.0016 * t * t - 0.00006 * t * t * t) * 1000.0)
    }
}

/// Virtual temperature in Celsius
///
/// From the [Glossary of Meteorology].(http://glossary.ametsoc.org/wiki/Virtual_temperature)
///
/// * `temperature_c` - the temperature in Celsius.
/// * `dew_point_c` - the dew point in Celsius.
/// * `pressure_hpa` - the pressure in hPa
///
/// Returns the virtual temperature in Celsius.
#[inline]
pub fn virtual_temperature_c(
    temperature_c: f64,
    dew_point_c: f64,
    pressure_hpa: f64,
) -> Result<f64> {
    let rv = mixing_ratio(dew_point_c, pressure_hpa)?;
    let t_k = theta_kelvin(pressure_hpa, temperature_c)?;
    let vt_k = t_k * (1.0 + rv / epsilon) / (1.0 + rv);
    temperature_c_from_theta(vt_k, pressure_hpa)
}

/// Convert wind speed (wind direction is the direction the wind is blowing from) into
/// U-V components.
///
/// Returns a tuple `(u_mps, v_mps)`, which is the U and V wind components in meters per second.
#[inline]
pub fn spd_dir_to_uv(from_dir_in_degrees: f64, speed_in_knots: f64) -> (f64, f64) {
    let rads = from_dir_in_degrees.to_radians();
    let spd_ms = knots_to_mps(speed_in_knots);

    (-spd_ms * rads.sin(), -spd_ms * rads.cos())
}

/// Convert U-V wind speeds to speed and direction.
///
/// Returns a tuple `(direction_degrees, speed_in_knots)`.
#[inline]
pub fn uv_to_spd_dir(u_mps: f64, v_mps: f64) -> (f64, f64) {
    let spd_ms = (u_mps.powi(2) + v_mps.powi(2)).sqrt();
    let spd_knots = mps_to_knots(spd_ms);

    let mut degrees = 180.0 + 90.0 - v_mps.atan2(u_mps).to_degrees();
    while degrees < 0.0 {
        degrees += 360.0;
    }
    while degrees >= 360.0 {
        degrees -= 360.0;
    }

    (degrees, spd_knots)
}

/// Convert knots to m/s.
#[inline]
pub fn knots_to_mps(spd: f64) -> f64  {
    spd * 0.514_444_444
}

/// Convert m/s to knots.
#[inline]
pub fn mps_to_knots(spd: f64) -> f64 {
    spd * 1.943_844_5
}

/// Bisection algorithm for finding the root of an equation given values bracketing a root. Used
/// when finding wet bulb temperature.
// FIXME: Update to use Brent's method.
fn find_root(f: &Fn(f64) -> Result<f64>, mut low_val: f64, mut high_val: f64) -> Result<f64> {
    use std::f64;
    const MAX_IT: usize = 50;
    const EPS: f64 = 1.0e-10;

    if low_val > high_val {
        ::std::mem::swap(&mut low_val, &mut high_val);
    }

    let mut f_low = f(low_val)?;
    let f_high = f(high_val)?;

    // Check to make sure we have bracketed a root.
    if f_high * f_low > 0.0 {
        return Err(InputOutOfRange);
    }

    let mut mid_val = (high_val - low_val) / 2.0 + low_val;
    let mut f_mid = f(mid_val)?;
    for _ in 0..MAX_IT {
        if f_mid * f_low > 0.0 {
            low_val = mid_val;
            f_low = f_mid;
        } else {
            high_val = mid_val;
        }

        if (high_val - low_val).abs() < EPS {
            break;
        }

        mid_val = (high_val - low_val) / 2.0 + low_val;
        f_mid = f(mid_val)?;
    }

    Ok(mid_val)
}

#[cfg(test)]
mod test {
    use super::*;

    const TOL: f64 = 1.0e-9;

    fn approx_equal(left: f64, right: f64, tol: f64) -> bool {
        use std::f64;
        assert!(tol > 0.0);
        f64::abs(left - right) <= tol
    }

    struct DRange {
        start: f64,
        step: f64,
        stop: f64,
    }

    impl Iterator for DRange {
        type Item = f64;

        fn next(&mut self) -> Option<Self::Item> {
            if self.step > 0.0 && self.start > self.stop {
                None
            } else if self.step < 0.0 && self.start < self.stop {
                None
            } else {
                let next = self.start;
                self.start += self.step;
                Some(next)
            }
        }
    }

    fn pressure_levels() -> DRange {
        DRange {
            start: 1000.0,
            step: -10.0,
            stop: 100.0,
        }
    }

    fn temperatures() -> DRange {
        DRange {
            start: -100.0,
            step: 10.0,
            stop: 100.0,
        }
    }

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

    #[test]
    fn test_celsius_to_kelvin() {
        assert_approx_eq!(celsius_to_kelvin(-10.0).unwrap(), 263.15);
        assert_approx_eq!(celsius_to_kelvin(0.0).unwrap(), 273.15);
        assert_approx_eq!(celsius_to_kelvin(10.0).unwrap(), 283.15);
        assert_approx_eq!(celsius_to_kelvin(-273.15).unwrap(), 0.0);
        assert!(celsius_to_kelvin(-300.0).unwrap_err() == BelowAbsoluteZero);
    }

    #[test]
    fn test_kelvin_to_celsius() {
        assert_approx_eq!(kelvin_to_celsius(263.15).unwrap(), -10.0);
        assert_approx_eq!(kelvin_to_celsius(273.15).unwrap(), 0.0);
        assert_approx_eq!(kelvin_to_celsius(283.15).unwrap(), 10.0);
        assert_approx_eq!(kelvin_to_celsius(0.0).unwrap(), -273.15);
        assert!(kelvin_to_celsius(-10.0).unwrap_err() == BelowAbsoluteZero);
    }

    #[test]
    fn test_kelvin_celsius_conversions_are_inverses_of_each_other() {
        for celsius in temperatures() {
            let kelvin = celsius_to_kelvin(celsius).unwrap();
            let back_to_celsius = kelvin_to_celsius(kelvin).unwrap();
            assert_approx_eq!(celsius, back_to_celsius);
        }
    }

    #[test]
    fn test_celsius_to_f() {
        assert_approx_eq!(celsius_to_f(0.0).unwrap(), 32.0);
        assert_approx_eq!(celsius_to_f(100.0).unwrap(), 212.0);
        assert_approx_eq!(celsius_to_f(-40.0).unwrap(), -40.0);
        assert!(celsius_to_f(ABSOLUTE_ZERO_C - 1.0).unwrap_err() == BelowAbsoluteZero);
    }

    #[test]
    fn test_f_to_celsius() {
        assert_approx_eq!(f_to_celsius(32.0).unwrap(), 0.0);
        assert_approx_eq!(f_to_celsius(212.0).unwrap(), 100.0);
        assert_approx_eq!(f_to_celsius(-40.0).unwrap(), -40.0);
        assert!(f_to_celsius(ABSOLUTE_ZERO_F - 1.0).unwrap_err() == BelowAbsoluteZero);
    }

    #[test]
    fn test_f_celsius_conversions_are_inverses_of_each_other() {
        for celsius in temperatures() {
            let fahrenheit = celsius_to_f(celsius).unwrap();
            let back_to_celsius = f_to_celsius(fahrenheit).unwrap();
            assert_approx_eq!(celsius, back_to_celsius);
        }
    }

    #[test]
    fn test_theta_kelvin() {
        for celsius in temperatures() {
            let kelvin = celsius_to_kelvin(celsius).unwrap();
            let theta = theta_kelvin(1000.0, celsius).unwrap();

            assert_approx_eq!(kelvin, theta);
        }

        assert!(theta_kelvin(-1.0, 20.0).unwrap_err() == NegativePressure);
        assert!(theta_kelvin(100.0, ABSOLUTE_ZERO_C - 1.0).unwrap_err() == BelowAbsoluteZero);
    }

    #[test]
    fn test_temperature_c_from_theta() {
        for celsius in temperatures() {
            let kelvin = celsius_to_kelvin(celsius).unwrap();
            let back_to_celsius = temperature_c_from_theta(kelvin, 1000.0).unwrap();

            assert_approx_eq!(celsius, back_to_celsius);
        }

        assert!(temperature_c_from_theta(325.0, -100.0).unwrap_err() == NegativePressure);
        assert!(
            temperature_c_from_theta(ABSOLUTE_ZERO_K - 1.0, 100.0).unwrap_err()
                == BelowAbsoluteZero
        );
    }

    #[test]
    fn test_temperature_c_from_theta_and_back_by_theta_kelvin() {
        for p in pressure_levels() {
            for celsius in temperatures() {
                let theta = theta_kelvin(p, celsius).unwrap();
                let back_to_celsius = temperature_c_from_theta(theta, p).unwrap();
                assert_approx_eq!(celsius, back_to_celsius);
            }
        }
    }

    #[test]
    fn test_vapor_pressure_liquid_water_there_and_back() {
        for &dp in [
            -80.0, -60.0, -40.0, -20.0, -10.0, 0.0, 10.0, 20.0, 40.0, 49.0,
        ]
            .iter()
        {
            let forward = vapor_pressure_liquid_water(dp).unwrap();
            let back = dew_point_from_vapor_pressure_over_liquid(forward).unwrap();
            assert_approx_eq!(dp, back);
        }
    }

    #[test]
    fn test_vapor_pressure_ice_there_and_back() {
        for &dp in [-80.0, -60.0, -40.0, -20.0, -10.0, -5.0, 0.0].iter() {
            let forward = vapor_pressure_ice(dp).unwrap();
            let back = frost_point_from_vapor_pressure_over_ice(forward).unwrap();
            assert_approx_eq!(dp, back);
        }
    }

    #[test]
    fn test_rh() {
        for t in temperatures() {
            for dp in temperatures().filter(|dp| *dp <= t) {
                let rh_result = rh(t, dp);
                match rh_result {
                    Err(InputOutOfRange) => { /* Bubbled up from vapor pressure. */ }
                    Err(_) => panic!("Unexpected Error."),
                    Ok(rh) => {
                        assert!(rh > 0.0, "Negative RH!");
                        assert!(rh <= 1.0);
                        if t == dp {
                            assert!(rh == 1.0);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_rh_ice() {
        for t in temperatures() {
            for dp in temperatures().filter(|dp| *dp <= t) {
                let rh_result = rh_ice(t, dp);
                match rh_result {
                    Err(InputOutOfRange) => { /* Bubbled up from vapor pressure. */ }
                    Err(_) => panic!("Unexpected Error."),
                    Ok(rh) => {
                        assert!(rh > 0.0, "Negative RH!");
                        assert!(rh <= 1.0);
                        if t == dp {
                            assert!(rh == 1.0);
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn test_mixing_ratio() {
        for p in pressure_levels() {
            for dp in temperatures() {
                match vapor_pressure_liquid_water(dp) {
                    Ok(vp) => {
                        if vp >= p {
                            continue;
                        }

                        let mw = mixing_ratio(dp, p).unwrap();
                        assert!(mw > 0.0);
                    }
                    Err(VaporPressureTooHigh) => panic!("Vapor Pressure Too High."),
                    Err(_) => { /* Ignore error that bubbled up from other func*/ }
                }
            }
        }
    }

    #[test]
    fn test_mixing_ratio_and_back_by_dew_point_from_p_and_mw() {
        for press in pressure_levels() {
            for dp in temperatures() {
                let mw = mixing_ratio(dp, press);

                match mw {
                    Ok(mw) => {
                        let back = dew_point_from_p_and_mw(press, mw);
                        match back {
                            Ok(back) => {
                                println!(
                                    "{} == {} with tolerance <= {}",
                                    dp,
                                    back,
                                    f64::abs(dp - back)
                                );
                                assert!(approx_equal(dp, back, 1.0e-2));
                            }
                            Err(_) => { /* Ignore error that bubbled up. */ }
                        }
                    }
                    Err(VaporPressureTooHigh) => { /* Ignore this for now.*/ }
                    Err(_) => { /* Ignore error that bubbled up.*/ }
                }
            }
        }
    }

    #[test]
    fn test_pressure_and_temperature_at_lcl() {
        for p in pressure_levels() {
            for t in temperatures() {
                for dp in temperatures().filter(|dp| *dp <= t) {
                    match pressure_and_temperature_at_lcl(t, dp, p) {
                        Ok((p_lcl, t_lcl)) => {
                            let k = celsius_to_kelvin(t).unwrap();
                            assert!(p_lcl <= p);
                            assert!((t_lcl - k) < 0.001);
                            assert!(p_lcl > 0.0);
                            assert!(t_lcl > 0.0);
                        }
                        Err(_) => { /* Ignore these for now. */ }
                    }
                }
            }
        }
    }

    #[test]
    fn test_pressure_hpa_at_lcl() {
        for p in pressure_levels() {
            for t in temperatures() {
                for dp in temperatures().filter(|dp| *dp <= t) {
                    match pressure_hpa_at_lcl(t, dp, p) {
                        Ok(p_lcl) => {
                            assert!(p_lcl <= p);
                            assert!(p_lcl > 0.0);
                        }
                        Err(_) => { /* Ignore */ }
                    }
                }
            }
        }
    }

    #[test]
    fn test_specific_humidity() {
        for p in pressure_levels() {
            for dp in temperatures() {
                if let Ok(sh) = specific_humidity(dp, p) {
                    assert!(sh > 0.0 && sh <= 1.0);
                }
            }
        }
    }

    #[test]
    fn test_theta_e_kelvin() {
        for p in pressure_levels() {
            for t in temperatures() {
                for dp in temperatures() {
                    match theta_e_kelvin(t, dp, p) {
                        Ok(theta_e) => {
                            let theta = theta_kelvin(p, t).unwrap();
                            assert!(theta <= theta_e);
                        }
                        Err(_) => { /* Ignore forwarded errors. */ }
                    }
                }
            }
        }
    }

    #[test]
    fn test_theta_e_saturated_theta_e_theta_relationships() {
        for p in pressure_levels() {
            for t in temperatures() {
                for dp in temperatures().filter(|dp| *dp <= t) {
                    if let Ok((theta, theta_e, theta_es)) = theta_kelvin(p, t)
                        .and_then(|theta| Ok((theta, theta_e_kelvin(t, dp, p)?)))
                        .and_then(|pair| Ok((pair.0, pair.1, theta_e_kelvin(t, t, p)?)))
                    {
                        println!("{} <= {} <= {}", theta, theta_e, theta_es);
                        assert!(
                            approx_lte(theta, theta_e, TOL) && approx_lte(theta_e, theta_es, TOL)
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_theta_e_saturated_there_and_back() {
        for p in pressure_levels() {
            for t in temperatures() {
                if let Ok(t_back) = theta_e_kelvin(t, t, p).and_then(|theta_es| {
                    temperature_c_from_theta_e_saturated_and_pressure(p, theta_es)
                }) {
                    println!("{} {} {}", t, t_back, (t - t_back).abs());
                    assert!(approx_equal(t, t_back, 1.0e-7));
                }
            }
        }
    }

    #[test]
    fn test_wet_bulb_dp_temp_consistency() {
        for &press in [1000.0, 906.4, 850.0, 700.0, 314.6].iter() {
            for &t in [-45.46, -20.0, -10.0, -5.0, -1.46, 0.0, 5.0, 10.0, 20.0].iter() {
                for &dp in [-59.49, -21.0, -11.0, -5.1, -3.06, 0.0, 4.9, 9.9, 20.0].iter() {
                    if dp > t {
                        continue;
                    }
                    if let Ok(wb) = wet_bulb_c(t, dp, press) {
                        assert!(dp <= wb);
                        assert!(wb <= t);
                        assert!(dp <= t);
                    }
                }
            }
        }
    }

    #[test]
    fn test_virtual_temperature_c() {
        let t_c = [0.0, 10.0, 20.0, 25.0, 30.0];
        let dp_c = [-15.0, -15.0, 18.0, 19.0, 19.0];
        let p_hpa = [1000.0, 1000.0, 900.0, 875.0, 875.0];
        let vt = [0.2, 10.21, 22.57, 27.86, 32.91];

        for (&t, (&dp, (&p, &target))) in
            t_c.iter().zip(dp_c.iter().zip(p_hpa.iter().zip(vt.iter())))
        {
            println!(
                "target = {}, value = {}",
                target,
                virtual_temperature_c(t, dp, p).unwrap()
            );
            assert!(approx_equal(
                target,
                virtual_temperature_c(t, dp, p).unwrap(),
                0.05
            ));
        }
    }

    #[test]
    fn test_find_root() {
        assert!(approx_equal(
            1.0,
            find_root(&|x| Ok(x * x - 1.0), 2.0, 0.0).unwrap(),
            1.0e-10
        ));
        assert!(approx_equal(
            -1.0,
            find_root(&|x| Ok(x * x - 1.0), -2.0, 0.0).unwrap(),
            1.0e-10
        ));
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

        const LOCAL_TOL: f64 = 1.0e-6;

        for &((dir, spd), (u, v)) in &spd_dir_to_uv_data {
            let (calc_u, calc_v) = spd_dir_to_uv(dir, spd);
            assert_approx_eq!(calc_u, u, LOCAL_TOL, "U speed mismatch.");
            assert_approx_eq!(calc_v, v, LOCAL_TOL, "V speed mismatch.");
        }
    }

    #[test]
    fn test_uv_to_spd_dir() {
        let spd_dir_to_uv_data = get_test_winds();

        const LOCAL_TOL: f64 = 1.0e-6;

        for &((dir, spd), (u, v)) in &spd_dir_to_uv_data {
            let (calc_dir, calc_spd) = uv_to_spd_dir(u, v);
            assert_approx_eq!(calc_dir, dir, LOCAL_TOL, "Direction mismatch.");
            assert_approx_eq!(calc_spd, spd, LOCAL_TOL, "Speed mismatch.");
        }
    }

    fn approx_lte(a: f64, b: f64, tol: f64) -> bool {
        (a - b) <= tol
    }
}
