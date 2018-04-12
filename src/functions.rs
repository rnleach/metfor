use constants::*;
use error::*;
use error::MetForErr::*;

/// Convert temperatures Celsius to Kelvin.
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

/// Calculate the potential temperature of a parcel.
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
        Ok(celsius_to_kelvin(temperature_c)? * f64::powf(1000.0 / pressure_hpa, R / cp))
    }
}

/// Given a potential temperature and pressure, calculate the temperature of a parcel.
///
/// * `theta_kelvin` - The potential temperature in Kelvin.
/// * `pressure_hpa` - The pressure of the parcel in hPa
///
/// Returns: The temperature in
#[inline]
pub fn temperature_c_from_theta(theta_kelvin: f64, pressure_hpa: f64) -> Result<f64> {
    use std::f64;
    if pressure_hpa < 0.0 {
        Err(NegativePressure)
    } else if theta_kelvin < ABSOLUTE_ZERO_K {
        Err(BelowAbsoluteZero)
    } else {
        Ok(kelvin_to_celsius(
            theta_kelvin * f64::powf(pressure_hpa / 1000.0, R / cp),
        )?)
    }
}

// Constants used for the limits of applicability to the empirical relationships used for
// vapor pressure.
const MIN_T_VP_ICE: f64 = -80.0;
const MIN_T_VP_LIQUID: f64 = -40.0;
const MAX_T_VP_ICE: f64 = 0.0;
const MAX_T_VP_LIQUID: f64 = 50.0;

/// Get the vapor pressure over liquid water.
///
/// Alduchov, O.A., and Eskridge, R.E. Improved Magnus` form approximation of saturation vapor
/// pressure. United States: N. p., 1997. Web. doi:10.2172/548871.
///
///  * `dew_point_c` - the dew point in Celsius of the parcel. If saturation is assumed this is also
///                   the temperature of the parcel.
///
/// Returns: The vapor pressure of water vapor in hPa.
#[inline]
pub fn vapor_pressure_liquid_water(dew_point_c: f64) -> Result<f64> {
    use std::f64;
    if dew_point_c < MIN_T_VP_LIQUID || dew_point_c > MAX_T_VP_LIQUID {
        Err(InputOutOfRange)
    } else {
        Ok(6.1094 * f64::exp(17.625 * dew_point_c / (dew_point_c + 243.04)))
    }
}

/// Get the dew point given the vapor pressure of water over liquid water. This function is the
/// inverse of `vapor_pressure_liquid_water`.
///
/// The corresponding temperature limits on that function correspond to vapor pressures of
/// 0.018968 hPa and 123.60577 hPa. Values outside this range will be returned as an error for
/// inputs out of range.
///
/// `vp_hpa` - The vapor pressure of water in hPa.
///
/// Returns: The dew point in Celsius.
#[inline]
fn dew_point_from_vapor_pressure_over_liquid(vp_hpa: f64) -> Result<f64> {
    use std::f64;

    if vp_hpa < 0.018968 || vp_hpa > 123.60577 {
        Err(InputOutOfRange)
    } else {
        let a = f64::ln(vp_hpa / 6.1094) / 17.625;
        Ok(a * 243.04 / (1.0 - a))
    }
}

/// Get the vapor pressure over ice.
///
/// Alduchov, O.A., and Eskridge, R.E. Improved Magnus` form approximation of saturation vapor
/// pressure. United States: N. p., 1997. Web. doi:10.2172/548871.
///
///  * `dew_point_c` - the dew point in Celsius of the parcel. If saturation is assumed this is also
///                   the temperature of the parcel.
///
/// Returns: The vapor pressure of water vapor in hPa.
#[inline]
pub fn vapor_pressure_ice(dew_point_c: f64) -> Result<f64> {
    use std::f64;
    if dew_point_c < MIN_T_VP_ICE || dew_point_c > MAX_T_VP_ICE {
        Err(InputOutOfRange)
    } else {
        Ok(6.1121 * f64::exp(22.587 * dew_point_c / (dew_point_c + 273.86)))
    }
}

/// Get the dew point given the vapor pressure of water over ice. This function is the inverse of
/// `vapor_pressure_ice`.
///
/// The corresponding temperature limits on that function correspond to vapor pressures of
/// 0.0005472 hPa and 6.1121 hPa. Values outside this range will be returned as an error for inputs
/// out of range.
///
/// `vp_hpa` - The vapor pressure of water in hPa.
///
/// Returns: The dew point in Celsius.
#[inline]
fn dew_point_from_vapor_pressure_over_ice(vp_hpa: f64) -> Result<f64> {
    use std::f64;

    if vp_hpa < 0.0005472 || vp_hpa > 6.1121 {
        Err(InputOutOfRange)
    } else {
        let a = f64::ln(vp_hpa / 6.1121) / 22.587;
        Ok(a * 273.86 / (1.0 - a))
    }
}

/// Get the vapor pressure of water using the liquid equation above freezing and the ice equation
/// below freezing.
///
/// Uses a combination of the `vapor_pressure_ice` and `vapor_pressure_liquid_water` functions above
/// and chooses the most appropriate value.
///
/// * `dew_point_c` - the dew point in Celsius of the parcel. If saturation is assumed this is also
///                   the temperature of the parcel.
///
/// Returns: The vapor pressure of water vapor in hPa.
#[inline]
pub fn vapor_pressure_water(dew_point_c: f64) -> Result<f64> {
    if dew_point_c < MAX_T_VP_ICE {
        vapor_pressure_ice(dew_point_c)
    } else {
        vapor_pressure_liquid_water(dew_point_c)
    }
}

/// Get the dew point given the vapor pressure of water. This function is the inverse of
/// `vapor_pressure_water`.
///
/// `vp_hpa` - The vapor pressure of water in hPa.
///
/// Returns: The dew point in Celsius.
#[inline]
fn dew_point_from_vapor_pressure(vp_hpa: f64) -> Result<f64> {
    if vp_hpa < 6.1094 {
        dew_point_from_vapor_pressure_over_ice(vp_hpa)
    } else {
        dew_point_from_vapor_pressure_over_liquid(vp_hpa)
    }
}

/// Calculate the relative humidity.
///
/// * `temperature_c` - The temperature of the parcel in Celsius.
/// * `dew_point_c` - The dew point of the parcel in Celsius.
///
/// Returns: The relative humidity as a decimal, i.e. 0.95 instead of 95%.
#[inline]
pub fn rh(temperature_c: f64, dew_point_c: f64) -> Result<f64> {
    let (es, e) = if temperature_c < MAX_T_VP_ICE {
        (
            vapor_pressure_ice(temperature_c)?,
            vapor_pressure_ice(dew_point_c)?,
        )
    } else {
        (
            vapor_pressure_liquid_water(temperature_c)?,
            vapor_pressure_liquid_water(dew_point_c)?,
        )
    };

    // Allow e > es for supersaturation.
    Ok(e / es)
}

/// Calculate the relative humidity with respect to liquid water.
///
/// * `temperature_c` - The temperature of the parcel in Celsius.
/// * `dew_point_c` - The dew point of the parcel in Celsius.
///
/// Returns: The relative humidity as a decimal, i.e. 0.95 instead of 95%.
#[inline]
pub fn rh_water(temperature_c: f64, dew_point_c: f64) -> Result<f64> {
    let (es, e) = (
        vapor_pressure_liquid_water(temperature_c)?,
        vapor_pressure_liquid_water(dew_point_c)?,
    );

    // Allow e > es for supersaturation.
    Ok(e / es)
}

/// Calculate the relative humidity with respect to ice.
///
/// * `temperature_c` - The temperature of the parcel in Celsius.
/// * `dew_point_c` - The dew point of the parcel in Celsius.
///
/// Returns: The relative humidity as a decimal, i.e. 0.95 instead of 95%.
#[inline]
pub fn rh_ice(temperature_c: f64, dew_point_c: f64) -> Result<f64> {
    let (es, e) = (
        vapor_pressure_ice(temperature_c)?,
        vapor_pressure_ice(dew_point_c)?,
    );

    // Allow e > es for supersaturation.
    Ok(e / es)
}

/// Calculate the mixing ratio of water.
///
/// Eq 5.9 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// * `dew_point_c` - the dew point of the parcel in Celsius. If you want the saturation mixing
///                   ratio, just make this equal to the temperature.
/// * `pressure_hpa` - the pressure of the parcel in hPa.
///
/// Returns: The mixing ratio in kg/kg.
#[inline]
pub fn mixing_ratio(dew_point_c: f64, pressure_hpa: f64) -> Result<f64> {
    let vp = vapor_pressure_water(dew_point_c)?;
    if vp > pressure_hpa {
        Err(VaporPressureTooHigh)
    } else {
        Ok(R / Rv * vp / (pressure_hpa - vp))
    }
}

/// Given a mixing ratio and pressure, calculate the dew point temperature. If saturation is
/// assumed this is also the temperature.
///
/// * `pressure_hpa` - the pressure of the parcel in hPa
/// * `mw` - the mixing ratio in kg/kg.
///
/// Returns: The dew point in celsius
#[inline]
pub fn dew_point_from_p_and_mw(pressure_hpa: f64, mw: f64) -> Result<f64> {
    let vp = mw * pressure_hpa / (mw + R / Rv);
    dew_point_from_vapor_pressure(vp)
}

/// Approximate temperature at the Lifting Condensation Level (LCL).
///
/// Eq 5.17 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// * `temperature_c` - the initial temperature in Celsius.
/// * `dew_point_c` - the initial dew point in Celsius.
///
/// Returns: The temperature at the LCL in Kelvin.
#[inline]
pub fn temperature_kelvin_at_lcl(temperature_c: f64, dew_point_c: f64) -> Result<f64> {
    if dew_point_c < ABSOLUTE_ZERO_C || temperature_c < ABSOLUTE_ZERO_C {
        Err(BelowAbsoluteZero)
    } else if dew_point_c >= temperature_c {
        Ok(celsius_to_kelvin(temperature_c)?)
    } else {
        let celsius_lcl =
            dew_point_c - (0.001296 * dew_point_c + 0.1963) * (temperature_c - dew_point_c);
        celsius_to_kelvin(celsius_lcl)
    }
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
    use std::f64;

    if dew_point_c < ABSOLUTE_ZERO_C || temperature_c < ABSOLUTE_ZERO_C {
        Err(BelowAbsoluteZero)
    } else if dew_point_c >= temperature_c {
        Ok((pressure_hpa, celsius_to_kelvin(temperature_c)?))
    } else {
        let t_lcl = temperature_kelvin_at_lcl(temperature_c, dew_point_c)?;
        let t_kelvin = celsius_to_kelvin(temperature_c)?;

        let p_lcl = pressure_hpa * f64::powf(t_lcl / t_kelvin, cp / R);

        Ok((p_lcl, t_lcl))
    }
}

/// Approximate pressure of the Lifting Condensation Level (LCL).
///
/// Eq 5.18 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// * `temperature_c` - the initial temperature in Celsius.
/// * `dew_point_c` - the initial dew point in Celsius.
/// * `pressure_hpa` - the initial pressure of the parcel in hPa.
///
/// Returns: The pressure at the LCL in hPa.
#[inline]
pub fn pressure_hpa_at_lcl(temperature_c: f64, dew_point_c: f64, pressure_hpa: f64) -> Result<f64> {
    pressure_and_temperature_at_lcl(temperature_c, dew_point_c, pressure_hpa).map(|r| r.0)
}

/// Calculate the specific humidity.
///
/// Eqs 5.11 and 5.12 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// * `dew_point_c` - the dew point in Celsius, if this is the same as the temperature then this
///                   calculates the saturation specific humidity.
/// * `pressure_hpa` - thepressure of the parcel in hPa.
///
/// Returns the specific humidity. (no units)
#[inline]
pub fn specific_humidity(dew_point_c: f64, pressure_hpa: f64) -> Result<f64> {
    let vp = vapor_pressure_water(dew_point_c)?;

    if pressure_hpa < 0.0 {
        Err(NegativePressure)
    } else if vp > pressure_hpa {
        Err(UnphysicalInput)
    } else {
        Ok(vp / pressure_hpa * (R / Rv))
    }
}

/// Calculate equivalent potential temperature.
///
/// Eq 5.23 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// * `temperature_c` - the initial temperature in Celsius.
/// * `dew_point_c` - the initial dew point in Celsius.
/// * `pressure_hpa` - the initial pressure of the parcel in hPa.
///
/// Returns: The equivalent potential temperature in Kelvin.
#[inline]
pub fn theta_e_kelvin(temperature_c: f64, dew_point_c: f64, pressure_hpa: f64) -> Result<f64> {
    let theta = theta_kelvin(pressure_hpa, temperature_c)?;
    let t_lcl = temperature_kelvin_at_lcl(temperature_c, dew_point_c)?;
    let qs = specific_humidity(dew_point_c, pressure_hpa)?;
    let lc = latent_heat_of_condensation(temperature_c)?;

    Ok(theta * (1.0 + lc * qs / (cp * t_lcl)))
}

/// Calculate the equivalent potential temperature assuming a saturated parcel.
///
/// Eq 5.23 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// * `temperature_c` - the initial temperature in Celsius.
/// * `pressure_hpa` - the initial pressure of the parcel in hPa.
///
/// Returns: The saturated equivalent potential temperature in Kelvin.
#[inline]
pub fn theta_e_saturated_kelvin(pressure_hpa: f64, temperature_c: f64) -> Result<f64> {
    let theta = theta_kelvin(pressure_hpa, temperature_c)?;
    let qs = specific_humidity(temperature_c, pressure_hpa)?;
    let lc = latent_heat_of_condensation(temperature_c)?;

    Ok(theta * (1.0 + lc * qs / (cp * celsius_to_kelvin(temperature_c)?)))
}

/// Calculate the web bulb temperature.
///
/// * `temperature_c` - the temperature of the parcel in Celsius.
/// * `dew_point_c` - the dew point of the parcel in Celsius.
/// * `pressure_hpa` - the pressure of the parcel in hPa
///
/// Returns: The wet bulb temperature in Celsius.
#[inline]
pub fn wet_bulb_c(temperature_c: f64, dew_point_c: f64, pressure_hpa: f64) -> Result<f64> {
    let (p_lcl, t_lcl) = pressure_and_temperature_at_lcl(temperature_c, dew_point_c, pressure_hpa)?;
    let theta_e = theta_e_saturated_kelvin(p_lcl, kelvin_to_celsius(t_lcl)?)?;

    find_root(
        &|celsius| Ok(theta_e_saturated_kelvin(pressure_hpa, celsius)? - theta_e),
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
    // latent heat of sublimation, but does not exceed it. This seems very reasonable be to me,
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
/// From the Glossary of Meteorology. http://glossary.ametsoc.org/wiki/Virtual_temperature
///
/// * `temperature_c` - the temperature of the parcel in Celsius.
/// * `dew_point_c` - the dew point of the parcel in Celsius.
/// * `pressure_hpa` - the pressure of the parcel in hPa
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

/// Bisection algorithm for finding the root of an equation given values bracketing a root. Used
/// when finding wet bulb temperature.
fn find_root(f: &Fn(f64) -> Result<f64>, mut low_val: f64, mut high_val: f64) -> Result<f64> {
    use std::f64;
    const MAX_IT: usize = 50;
    const EPS: f64 = 1.0e-10;

    if low_val > high_val {
        ::std::mem::swap(&mut low_val, &mut high_val);
    }

    let mut f_low = f(low_val)?;

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
        ($a:expr, $b:expr) => {
            {
                println!("{} == {} with tolerance <= {}", $a, $b, f64::abs($a -$b));
                assert!(approx_equal($a,$b,TOL));
            }
        };
        ($a:expr, $b:expr, $msg:expr) => {
            {
                println!("{} == {} with tolerance <= {}", $a, $b, f64::abs($a -$b));
                assert!(approx_equal($a,$b,TOL), $msg);
            }
        };
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
        for &dp in [-40.0, -20.0, -10.0, 0.0, 10.0, 20.0, 40.0, 50.0].iter() {
            let forward = vapor_pressure_liquid_water(dp).unwrap();
            let back = dew_point_from_vapor_pressure_over_liquid(forward).unwrap();
            assert_approx_eq!(dp, back);
        }
    }

    #[test]
    fn test_vapor_pressure_ice_there_and_back() {
        for &dp in [-80.0, -60.0, -40.0, -20.0, -10.0, -5.0, 0.0].iter() {
            let forward = vapor_pressure_ice(dp).unwrap();
            let back = dew_point_from_vapor_pressure_over_ice(forward).unwrap();
            assert_approx_eq!(dp, back);
        }
    }

    #[test]
    fn test_vapor_pressure_water_and_back_by_dew_point_from_vapor_pressure() {
        for &dp in [-80.0, -20.0, -10.0, 0.0, 10.0, 20.0, 40.0, 50.0].iter() {
            let forward = vapor_pressure_water(dp).unwrap();
            let back = dew_point_from_vapor_pressure(forward).unwrap();
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
    fn test_rh_water() {
        for t in temperatures() {
            for dp in temperatures().filter(|dp| *dp <= t) {
                let rh_result = rh_water(t, dp);
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
                match vapor_pressure_water(dp) {
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
    fn test_temperature_kelvin_at_lcl() {
        for t in temperatures() {
            for dp in temperatures().filter(|dp| *dp <= t) {
                match temperature_kelvin_at_lcl(t, dp) {
                    Ok(lcl_t) => {
                        assert!(lcl_t >= 0.0);
                        assert!(lcl_t <= celsius_to_kelvin(t).unwrap());
                    }
                    Err(BelowAbsoluteZero) => { /* Ignore for now. */ }
                    Err(_) => { /* These errors bubbled up, ignore. */ }
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
                            assert!(t_lcl <= k);
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
    fn test_theta_e_saturated_kelvin() {
        for p in pressure_levels() {
            for t in temperatures() {
                match theta_e_saturated_kelvin(p, t) {
                    Ok(theta_es) => {
                        let theta = theta_kelvin(p, t).unwrap();
                        assert!(theta <= theta_es);
                    }
                    Err(_) => { /* Ignore errors, they were all passed up. */ }
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
                        .and_then(|pair| Ok((pair.0, pair.1, theta_e_saturated_kelvin(p, t)?)))
                    {
                        assert!(theta <= theta_e && theta_e <= theta_es);
                    }
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
}
