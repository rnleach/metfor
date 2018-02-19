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

pub mod constants;
pub use constants::*;

/// Convert temperatures Celsius to Kelvin.
#[inline(always)]
pub fn celsius_to_kelvin(celsius: f64) -> f64 {
    celsius + 273.15
}

/// Convert Kelvin temperatures to Celsius.
#[inline(always)]
pub fn kelvin_to_celsius(kelvin: f64) -> f64 {
    kelvin - 273.15
}

/// Convert Celsius to Fahrenheit.
#[inline(always)]
pub fn celsius_to_f(temperature: f64) -> f64 {
    1.8 * temperature + 32.0
}

/// Convert Fahrenheit to Celsius
#[inline(always)]
pub fn f_to_celsius(temperature: f64) -> f64 {
    (temperature - 32.0) / 1.8
}

/// Calculate the potential temperature of a parcel.
///
/// * `pressure_hpa` - Pressure in hPa
/// * `temperature_c` - Temperature in Celsius.
///
/// Returns: Potential temperature in Kelvin
#[inline]
pub fn theta_kelvin(pressure_hpa: f64, temperature_c: f64) -> f64 {
    use std::f64;

    celsius_to_kelvin(temperature_c) * f64::powf(1000.0 / pressure_hpa, R / cp)
}

/// Given a potential temperature and pressure, calculate the temperature of a parcel.
///
/// * `theta_kelvin` - The potential temperature in Kelvin.
/// * `pressure_hpa` - The pressure of the parcel in hPa
///
/// Returns: The temperature in
#[inline]
pub fn temperature_c_from_theta(theta_kelvin: f64, pressure_hpa: f64) -> f64 {
    use std::f64;

    kelvin_to_celsius(theta_kelvin * f64::powf(pressure_hpa / 1000.0, R / cp))
}

/// Get the vapor pressure of water.
///
/// Eq 5.1 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// * `dew_point_c` - the dew point in Celsius of the parcel. If saturation is assumed this is also
///                   the temperature of the parcel.
///
/// Returns: The vapor pressure of water vapor in hPa.
#[inline]
pub fn vapor_pressure_water(dew_point_c: f64) -> f64 {
    use std::f64;
    let kelvin = celsius_to_kelvin(dew_point_c);

    6.1078 * f64::exp(19.8 * dew_point_c / kelvin)
}

/// Get the dew point given the vapor pressure of water.
///
/// `vp_hpa` - The vapor pressure of water in hPa.
///
/// Returns: The dew point in Celsius.
fn dew_point_from_vapor_pressure(vp_hpa: f64) -> f64 {
    use std::f64;

    let kelvin = -273.15 / (f64::ln(vp_hpa / 6.1078) / 19.8 - 1.0);
    kelvin_to_celsius(kelvin)
}

/// Calculate the relative humidity.
///
/// * `temperature_c` - The temperature of the parcel in Celsius.
/// * `dew_point_c` - The dew point of the parcel in Celsius.
///
/// Returns: The relative humidity as a decimal, i.e. 0.95 instead of 95%.
#[inline]
pub fn rh(temperature_c: f64, dew_point_c: f64) -> f64 {
    let e = vapor_pressure_water(dew_point_c);
    let es = vapor_pressure_water(temperature_c);
    e / es
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
pub fn mixing_ratio(dew_point_c: f64, pressure_hpa: f64) -> f64 {
    let vp = vapor_pressure_water(dew_point_c);
    R / Rv * vp / (pressure_hpa - vp)
}

/// Given a mixing ratio and pressure, calculate the dew point temperature. If saturation is
/// assumed this is also the temperature.
///
/// * `pressure_hpa` - the pressure of the parcel in hPa
/// * `mw` - the mixing ratio in kg/kg.
///
/// Returns: The dew point in celsius
#[inline]
pub fn dew_point_from_p_and_mw(pressure_hpa: f64, mw: f64) -> f64 {
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
pub fn temperature_kelvin_at_lcl(temperature_c: f64, dew_point_c: f64) -> f64 {
    let celsius_lcl =
        dew_point_c - (0.001296 * dew_point_c + 0.1963) * (temperature_c - dew_point_c);
    celsius_to_kelvin(celsius_lcl)
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
) -> (f64, f64) {
    use std::f64;

    let t_lcl = temperature_kelvin_at_lcl(temperature_c, dew_point_c);
    let t_kelvin = celsius_to_kelvin(temperature_c);

    let p_lcl = pressure_hpa * f64::powf(t_lcl / t_kelvin, cp / R);

    (p_lcl, t_lcl)
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
pub fn pressure_hpa_at_lcl(temperature_c: f64, dew_point_c: f64, pressure_hpa: f64) -> f64 {
    pressure_and_temperature_at_lcl(temperature_c, dew_point_c, pressure_hpa).0
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
pub fn specific_humidity(dew_point_c: f64, pressure_hpa: f64) -> f64 {
    vapor_pressure_water(dew_point_c) / pressure_hpa * (R / Rv)
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
pub fn theta_e_kelvin(temperature_c: f64, dew_point_c: f64, pressure_hpa: f64) -> f64 {
    let theta = theta_kelvin(pressure_hpa, temperature_c);
    let t_lcl = temperature_kelvin_at_lcl(temperature_c, dew_point_c);
    let qs = specific_humidity(temperature_c, pressure_hpa);
    let lc = latent_heat_of_condensation(temperature_c);

    theta * (1.0 + lc * qs / (cp * t_lcl))
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
pub fn theta_e_saturated_kelvin(pressure_hpa: f64, temperature_c: f64) -> f64 {
    let theta = theta_kelvin(pressure_hpa, temperature_c);
    let qs = specific_humidity(temperature_c, pressure_hpa);
    let lc = latent_heat_of_condensation(temperature_c);

    theta * (1.0 + lc * qs / (cp * celsius_to_kelvin(temperature_c)))
}

/// Calculate the web bulb temperature.
///
/// * `temperature_c` - the temperature of the parcel in Celsius.
/// * `dew_point_c` - the dew point of the parcel in Celsius.
/// * `pressure_hpa` - the pressure of the parcel in hPa
///
/// Returns: The wet bulb temperature in Celsius.
#[inline]
pub fn wet_bulb_c(temperature_c: f64, dew_point_c: f64, pressure_hpa: f64) -> f64 {
    let (p_lcl, t_lcl) = pressure_and_temperature_at_lcl(temperature_c, dew_point_c, pressure_hpa);
    let theta_e = theta_e_saturated_kelvin(p_lcl, kelvin_to_celsius(t_lcl));

    find_root(
        &|celsius| theta_e_saturated_kelvin(pressure_hpa, celsius) - theta_e,
        dew_point_c,
        temperature_c,
    )
}

/// Latent heat of condensation for water.
///
/// Polynomial curve fits to Table 2.1. R. R. Rogers; M. K. Yau (1989). A Short Course in Cloud
/// Physics (3rd ed.). Pergamon Press. p. 16. ISBN 0-7506-3215-1.
///
/// `temperature_c` - temperature in Celsius.
///
/// Returns: the latent heat of condensation for water in J kg<sup>-1</sup>.
#[inline]
pub fn latent_heat_of_condensation(temperature_c: f64) -> f64 {
    let t = temperature_c;
    (2500.8 - 2.36 * t + 0.0016 * t * t - 0.00006 * t * t * t) * 1000.0
}

/// Bisection algorithm for finding the root of an equation given values bracketing a root. Used
/// when drawing moist adiabats.
fn find_root(f: &Fn(f64) -> f64, mut low_val: f64, mut high_val: f64) -> f64 {
    use std::f64;
    const MAX_IT: usize = 50;
    const EPS: f64 = 1.0e-10;

    if low_val > high_val {
        ::std::mem::swap(&mut low_val, &mut high_val);
    }

    let mut f_low = f(low_val);
    // let mut f_high = f(high_val);

    let mut mid_val = (high_val - low_val) / 2.0 + low_val;
    let mut f_mid = f(mid_val);
    for _ in 0..MAX_IT {
        if f_mid * f_low > 0.0 {
            low_val = mid_val;
            f_low = f_mid;
        } else {
            high_val = mid_val;
            // f_high = f_mid;
        }

        if (high_val - low_val).abs() < EPS {
            break;
        }

        mid_val = (high_val - low_val) / 2.0 + low_val;
        f_mid = f(mid_val);
    }

    mid_val
}

#[cfg(test)]
mod test {
    use super::*;

    fn approx_equal(left: f64, right: f64, tol: f64) -> bool {
        use std::f64;
        assert!(tol > 0.0);
        f64::abs(left - right) <= tol
    }

    #[test]
    fn test_find_root() {
        assert!(approx_equal(
            1.0,
            find_root(&|x| x * x - 1.0, 2.0, 0.0),
            1.0e-10
        ));
        assert!(approx_equal(
            -1.0,
            find_root(&|x| x * x - 1.0, -2.0, 0.0),
            1.0e-10
        ));
    }

    #[test]
    fn test_theta_kelvin() {
        for &t in [-20.0, -10.0, 0.0, 10.0, 20.0].iter() {
            let kelvin = celsius_to_kelvin(t);
            let theta = theta_kelvin(1000.0, t);

            assert!(approx_equal(kelvin, theta, 1.0e-9));
        }
    }

    #[test]
    fn test_vapor_pressure() {
        for &dp in [-20.0, -10.0, 0.0, 10.0, 20.0].iter() {
            let forward = vapor_pressure_water(dp);
            let back = dew_point_from_vapor_pressure(forward);

            println!("{} ~= {}", dp, back);
            assert!(approx_equal(dp, back, 1.0e-9), "Failed there and back!");
        }
    }

    #[test]
    fn test_mw() {
        for &press in [1000.0, 800.0, 700.0, 500.0, 300.0].iter() {
            for &dp in [-20.0, -10.0, 0.0, 10.0, 20.0].iter() {
                let mw = mixing_ratio(dp, press);
                let back = dew_point_from_p_and_mw(press, mw);

                println!("{} ~= {}", dp, back);
                assert!(approx_equal(dp, back, 1.0e-3), "Failed there and back!");
            }
        }
    }

    #[test]
    fn test_wet_bulb() {
        for &press in [1000.0, 906.4, 850.0, 700.0,314.6].iter() {
            for &t in [-45.46,-20.0, -10.0, -5.0, -1.46, 0.0, 5.0, 10.0, 20.0].iter() {
                for &dp in [-59.49,-21.0, -11.0, -5.1, -3.06, 0.0, 4.9, 9.9, 20.0].iter() {
                    if dp > t {
                        continue;
                    }
                    let wb = wet_bulb_c(t, dp, press);
                    println!("{} <= {} <= {} at {} hPa", dp, wb, t, press);
                    assert!(dp <= wb);
                    assert!(wb <= t);
                    assert!(dp <= t);
                }
            }
        }
    }
}
