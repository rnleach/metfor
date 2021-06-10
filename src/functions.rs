use crate::constants::*;
use crate::types::*;

macro_rules! debug_validate {
    ($($val:expr),+) => (
        #[cfg(debug_assertions)]{
            $(
                let _ = $val.unwrap();
            )+
        }
    );
}

/// Calculate potential temperature assuming a 1000 hPa reference level.
///
/// Returns: Potential temperature.
#[inline]
pub fn potential_temperature<P, T>(pressure: P, temperature: T) -> Kelvin
where
    P: Pressure,
    T: Temperature,
    HectoPascal: From<P>,
    Kelvin: From<T>,
{
    use std::f64;

    debug_validate!(pressure, temperature);

    let HectoPascal(p) = HectoPascal::from(pressure);
    let Kelvin(t) = Kelvin::from(temperature);
    let theta = t * f64::powf(1000.0 / p, Rd / cp);
    let theta = Kelvin(theta);

    debug_validate!(theta);

    theta
}

/// Given a potential temperature and pressure, calculate the temperature.
///
/// Returns: The temperature.
#[inline]
pub fn temperature_from_pot_temp<P, T>(theta: T, pressure: P) -> Kelvin
where
    P: Pressure,
    T: Temperature,
    HectoPascal: From<P>,
    Kelvin: From<T>,
{
    use std::f64;

    debug_validate!(pressure, theta);

    let HectoPascal(p) = HectoPascal::from(pressure);
    let Kelvin(theta) = Kelvin::from(theta);
    let t = theta * f64::powf(p / 1000.0, Rd / cp);
    let t = Kelvin(t);

    debug_validate!(t);

    t
}

// Constants used for the limits of applicability to the empirical relationships used for
// vapor pressure.
const MIN_T_VP_ICE: Celsius = Celsius(-80.0);
const MIN_T_VP_LIQUID: Celsius = Celsius(-80.0);
const MAX_T_VP_ICE: Celsius = Celsius(0.0);
const MAX_T_VP_LIQUID: Celsius = Celsius(50.0);

/// Get the vapor pressure over liquid water.
///
/// Alduchov, O.A., and Eskridge, R.E. Improved Magnus` form approximation of saturation vapor
/// pressure. United States: N. p., 1997. Web. doi:10.2172/548871.
///
/// Returns: The vapor pressure of water vapor.
#[inline]
pub fn vapor_pressure_water<T>(dew_point: T) -> Option<HectoPascal>
where
    T: Temperature + PartialOrd<Celsius>,
    Celsius: From<T>,
{
    debug_validate!(dew_point);

    if dew_point < MIN_T_VP_LIQUID || dew_point > MAX_T_VP_LIQUID {
        None
    } else {
        let Celsius(dp) = Celsius::from(dew_point);
        let vp = 6.1037 * f64::exp(17.641 * dp / (dp + 243.27));
        let vp = HectoPascal(vp);

        debug_validate!(vp);

        Some(vp)
    }
}

/// Get the dew point given the vapor pressure of water over liquid water. This function is the
/// inverse of `vapor_pressure_liquid_water`.
///
/// Returns: The dew point.
#[inline]
pub fn dew_point_from_vapor_pressure_water<P>(vp: P) -> Option<Celsius>
where
    P: Pressure,
    HectoPascal: From<P>,
{
    debug_validate!(vp);

    let HectoPascal(vp) = HectoPascal::from(vp);
    let a = f64::ln(vp / 6.1037) / 17.641;
    let dp = a * 243.27 / (1.0 - a);
    let dp = Celsius(dp);
    if dp < MIN_T_VP_LIQUID || dp > MAX_T_VP_LIQUID {
        None
    } else {
        Some(dp)
    }
}

/// Get the vapor pressure over ice.
///
/// Alduchov, O.A., and Eskridge, R.E. Improved Magnus` form approximation of saturation vapor
/// pressure. United States: N. p., 1997. Web. doi:10.2172/548871.
///
/// Returns: The vapor pressure of water vapor.
#[inline]
pub fn vapor_pressure_ice<T>(frost_point: T) -> Option<HectoPascal>
where
    T: Temperature + PartialOrd,
    Celsius: From<T>,
{
    debug_validate!(frost_point);

    let Celsius(fp) = Celsius::from(frost_point);
    // should not need unpack calls below, should be able to do this as temperatures but compiler
    // can't handle the type inference, even with help. rustc --version = 1.30.1
    if fp < MIN_T_VP_ICE.unpack() || fp > MAX_T_VP_ICE.unpack() {
        None
    } else {
        let vp = 6.1121 * f64::exp(22.587 * fp / (fp + 273.86));
        let vp = HectoPascal(vp);

        debug_validate!(vp);

        Some(vp)
    }
}

/// Get the frost point given the vapor pressure of water over ice. This function is the inverse of
/// `vapor_pressure_ice`.
///
/// Returns: The frost point.
#[inline]
pub fn frost_point_from_vapor_pressure_over_ice<P>(vp: P) -> Option<Celsius>
where
    P: Pressure,
    HectoPascal: From<P>,
{
    debug_validate!(vp);

    let HectoPascal(vp) = HectoPascal::from(vp);

    let a = f64::ln(vp / 6.1121) / 22.587;
    let fp = a * 273.86 / (1.0 - a);
    let fp = Celsius(fp);

    debug_validate!(fp);

    if fp < MIN_T_VP_ICE || fp > MAX_T_VP_ICE {
        None
    } else {
        Some(fp)
    }
}

/// Calculate the relative humidity with respect to liquid water.
///
/// Returns: The relative humidity as a decimal, i.e. 0.95 instead of 95%.
#[inline]
pub fn rh<T1, T2>(temperature: T1, dew_point: T2) -> Option<f64>
where
    T1: Temperature,
    T2: Temperature,
    Celsius: From<T1> + From<T2>,
{
    debug_validate!(temperature, dew_point);

    let t = Celsius::from(temperature);
    let dp = Celsius::from(dew_point);

    let (HectoPascal(es), HectoPascal(e)) = (vapor_pressure_water(t)?, vapor_pressure_water(dp)?);

    let rh_val = e / es;

    debug_assert!(rh_val >= 0.0);

    // Allow e > es for supersaturation.
    Some(rh_val)
}

/// Calculate the relative humidity with respect to ice.
///
/// Returns: The relative humidity as a decimal, i.e. 0.95 instead of 95%.
#[inline]
pub fn rh_ice<T1, T2>(temperature: T1, frost_point: T2) -> Option<f64>
where
    T1: Temperature,
    T2: Temperature,
    Celsius: From<T1> + From<T2>,
{
    debug_validate!(temperature, frost_point);

    let t = Celsius::from(temperature);
    let fp = Celsius::from(frost_point);

    let (HectoPascal(es), HectoPascal(e)) = (vapor_pressure_water(t)?, vapor_pressure_water(fp)?);

    let rh_val = e / es;

    debug_assert!(rh_val >= 0.0);

    // Allow e > es for supersaturation.
    Some(rh_val)
}

/// Calculate the mixing ratio of water.
///
/// Eq 5.9 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// Returns: The mixing ratio as a unitless value. Note this is often reported as g/kg, but this
///          function returns kg/kg.
#[inline]
pub fn mixing_ratio<T, P>(dew_point: T, pressure: P) -> Option<f64>
where
    T: Temperature,
    Celsius: From<T>,
    P: Pressure,
    HectoPascal: From<P>,
{
    debug_validate!(dew_point, pressure);

    let dp = Celsius::from(dew_point);
    let HectoPascal(p) = HectoPascal::from(pressure);

    let HectoPascal(vp) = vapor_pressure_water::<Celsius>(dp)?;
    if vp > p {
        None
    } else {
        Some(epsilon * vp / (p - vp))
    }
}

/// Given a mixing ratio and pressure, calculate the dew point temperature. If saturation is
/// assumed, this is also the temperature.
///
/// Returns: The dew point.
#[inline]
pub fn dew_point_from_p_and_mw<P>(pressure: P, mw: f64) -> Option<Celsius>
where
    P: Pressure,
    HectoPascal: From<P>,
{
    debug_validate!(pressure);
    debug_assert!(mw >= 0.0);

    let HectoPascal(p) = HectoPascal::from(pressure);

    let vp = HectoPascal(mw * p / (mw + epsilon));
    dew_point_from_vapor_pressure_water::<HectoPascal>(vp)
}

/// Calculate the specific humidity.
///
/// Eqs 5.11 and 5.12 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// * `dew_point` - the dew point, if this is the same as the temperature then this
///                 calculates the saturation specific humidity.
/// * `pressure` - the pressure in hPa.
///
/// Returns the specific humidity. (no units)
#[inline]
pub fn specific_humidity<DP, P>(dew_point: DP, pressure: P) -> Option<f64>
where
    DP: Temperature,
    P: Pressure,
    Celsius: From<DP>,
    HectoPascal: From<P>,
{
    debug_validate!(dew_point, pressure);

    let dp = Celsius::from(dew_point);
    let p = HectoPascal::from(pressure).into_option()?;
    let HectoPascal(vp) = vapor_pressure_water::<Celsius>(dp)?;

    if vp > p {
        None
    } else {
        Some(vp / p * epsilon)
    }
}

/// Convert specific humidity into mixing ratio.
pub fn mixing_ratio_from_specific_humidity(specific_humidity: f64) -> f64 {
    specific_humidity / (1.0 - specific_humidity)
}

/// Convert mixing ratio into specific humidity.
pub fn specific_humidity_from_mixing_ratio(mixing_ratio: f64) -> f64 {
    mixing_ratio / (1.0 + mixing_ratio)
}

/// Given a specific humidity and pressure, calculate the dew point temperature. If saturation is
/// assumed, this is also the temperature.
#[inline]
pub fn dew_point_from_p_and_specific_humidity<P>(pressure: P, q: f64) -> Option<Celsius>
where
    P: Pressure,
    HectoPascal: From<P>,
{
    debug_validate!(pressure);
    debug_assert!(q >= 0.0);

    let HectoPascal(p) = HectoPascal::from(pressure);

    let vp = HectoPascal(q * p / epsilon);
    dew_point_from_vapor_pressure_water::<HectoPascal>(vp)
}

/// Calculate the temperature and pressure at the lifting condensation level (LCL).
///
/// Eqs 5.17 and 5.18 from "Weather Analysis" by Dus&#780;an Dujric&#769;
///
/// * `temperature` - the initial temperature.
/// * `dew_point` - the initial dew point.
/// * `pressure` - the initial pressure of the parcel.
///
/// Returns: The pressure and the temperature at the LCL.
#[inline]
pub fn pressure_and_temperature_at_lcl<T, DP, P>(
    temperature: T,
    dew_point: DP,
    pressure: P,
) -> Option<(HectoPascal, Kelvin)>
where
    T: Temperature,
    DP: Temperature,
    Celsius: From<T> + From<DP>,
    P: Pressure,
    HectoPascal: From<P>,
{
    debug_validate!(temperature, dew_point, pressure);

    let t = Celsius::from(temperature);
    let dp = Celsius::from(dew_point);
    let p = HectoPascal::from(pressure);

    let plcl = pressure_at_lcl::<Celsius, Celsius, HectoPascal>(t, dp, p)?;
    let tlcl = temperature_from_pot_temp::<HectoPascal, Kelvin>(
        potential_temperature::<HectoPascal, Celsius>(p, t),
        plcl,
    );

    debug_validate!(plcl, tlcl);

    Some((plcl, tlcl))
}

/// Approximate pressure of the Lifting Condensation Level (LCL).
///
/// * `temperature` - the initial temperature.
/// * `dew_point` - the initial dew point.
/// * `pressure` - the initial pressure.
///
/// Returns: The pressure at the LCL.
#[inline]
pub fn pressure_at_lcl<T, DP, P>(temperature: T, dew_point: DP, pressure: P) -> Option<HectoPascal>
where
    T: Temperature,
    DP: Temperature,
    Celsius: From<T> + From<DP>,
    P: Pressure,
    HectoPascal: From<P>,
{
    debug_validate!(temperature, dew_point, pressure);

    let t = Celsius::from(temperature);
    let dp = Celsius::from(dew_point);
    let p = HectoPascal::from(pressure);

    if dp >= t {
        Some(p)
    } else {
        let theta = potential_temperature::<HectoPascal, Celsius>(p, t);
        let mw = mixing_ratio::<Celsius, HectoPascal>(dp, p)?;
        let theta_e = equiv_pot_temperature::<_, _, HectoPascal>(t, dp, p)?;

        let eq = |p: f64| -> Option<f64> {
            let p = HectoPascal(p);
            let Kelvin(t) = temperature_from_pot_temp::<HectoPascal, Kelvin>(theta, p);
            let Kelvin(dp) = Kelvin::from(dew_point_from_p_and_mw::<HectoPascal>(p, mw)?);

            Some(t - dp)
        };

        // Search between 1060 and 300 hPa. If it falls outside that range, give up!
        let first_guess = find_root(&eq, 300.0, 1060.0)?;

        let eq = |p: f64| -> Option<f64> {
            let p = HectoPascal(p);
            let Kelvin(t1) = temperature_from_pot_temp::<HectoPascal, Kelvin>(theta, p);
            let Kelvin(t2) =
                Kelvin::from(temperature_from_equiv_pot_temp_saturated_and_pressure::<
                    HectoPascal,
                    Kelvin,
                >(p, theta_e)?);

            Some(t1 - t2)
        };

        const DELTA: f64 = 10.0;
        let lclp = if let Some(lcl) = find_root(&eq, first_guess + DELTA, first_guess - DELTA) {
            lcl
        } else {
            first_guess
        };

        let lclp = HectoPascal(lclp);

        debug_validate!(lclp);

        Some(lclp)
    }
}

/// Calculate equivalent potential temperature.
///
/// Equation from ["The Glossary of Meteorology"]
/// (http://glossary.ametsoc.org/wiki/Equivalent_potential_temperature) online where the
/// approximation of ignoring the "total water mixing ratio" is used since most of the time we do
/// not have the necessary information to calculate that.
///
/// * `temperature` - the initial temperature.
/// * `dew_point` - the initial dew point.
/// * `pressure` - the initial pressure.
///
/// Returns: The equivalent potential temperature.
#[inline]
pub fn equiv_pot_temperature<T, DP, P>(temperature: T, dew_point: DP, pressure: P) -> Option<Kelvin>
where
    T: Temperature,
    DP: Temperature,
    P: Pressure,
    Celsius: From<T> + From<DP>,
    HectoPascal: From<P>,
{
    debug_validate!(temperature, dew_point, pressure);

    let t = Celsius::from(temperature);
    let Kelvin(naked_tk) = Kelvin::from(t);
    let dp = Celsius::from(dew_point);
    let p = HectoPascal::from(pressure);
    let HectoPascal(naked_p) = p;

    const P0: f64 = 1000.0; // Reference pressure for potential temperatures.
    let rv = mixing_ratio::<Celsius, HectoPascal>(dp, p)?;
    let pd = naked_p - vapor_pressure_water::<Celsius>(dp)?.into_option()?;

    if pd < 0.0 {
        return None;
    }

    let h = rh(t, dp)?;

    let JpKg(lv) = latent_heat_of_condensation(t)?;

    let theta_e_val = Kelvin(
        naked_tk
            * f64::powf(P0 / pd, Rd / cp)
            * f64::powf(h, -rv * (Rv / cp))
            * f64::exp(lv * rv / cp.unpack() / naked_tk),
    );

    debug_validate!(theta_e_val);

    Some(theta_e_val)
}

/// Given the pressure and equivalent potential temperature, assume saturation and calculate the
/// temperature.
///
/// This function is useful if you were trying to calculate the temperature for plotting a
/// moist adiabat on a skew-t log-p diagram.
///
/// Returns: The temperature.
#[inline]
pub fn temperature_from_equiv_pot_temp_saturated_and_pressure<P, T>(
    pressure: P,
    theta_e_val: T,
) -> Option<Celsius>
where
    P: Pressure,
    HectoPascal: From<P>,
    T: Temperature,
    Kelvin: From<T>,
{
    debug_validate!(pressure, theta_e_val);

    let p = HectoPascal::from(pressure);
    let Kelvin(theta_e_k) = Kelvin::from(theta_e_val);

    find_root(
        &|t_c| {
            let t = Celsius(t_c);
            let Kelvin(theta_e_calc) =
                equiv_pot_temperature::<Celsius, Celsius, HectoPascal>(t, t, p)?;
            Some(theta_e_calc - theta_e_k)
        },
        -80.0,
        50.0,
    )
    .map(Celsius)
}

/// Calculate the web bulb temperature.
///
/// Returns: The wet bulb temperature.
#[inline]
pub fn wet_bulb<T, DP, P>(temperature: T, dew_point: DP, pressure: P) -> Option<Celsius>
where
    T: Temperature,
    DP: Temperature,
    Celsius: From<T> + From<DP>,
    P: Pressure,
    HectoPascal: From<P>,
{
    debug_validate!(temperature, dew_point, pressure);

    let t = Celsius::from(temperature);
    let dp = Celsius::from(dew_point);
    let p = HectoPascal::from(pressure);

    let Celsius(temperature_c) = t;
    let Celsius(dew_point_c) = dp;

    let (p_lcl, t_lcl) = pressure_and_temperature_at_lcl::<_, _, HectoPascal>(t, dp, p)?;

    let Kelvin(theta_e_val) = equiv_pot_temperature::<_, _, HectoPascal>(t_lcl, t_lcl, p_lcl)?;

    find_root(
        &|t_c| {
            let t_c = Celsius(t_c);
            let Kelvin(theta_e_calc) = equiv_pot_temperature::<_, _, HectoPascal>(t_c, t_c, p)?;

            Some(theta_e_calc - theta_e_val)
        },
        dew_point_c,
        temperature_c,
    )
    .map(Celsius)
}

/// Latent heat of condensation for water.
///
/// Polynomial curve fit to Table 2.1. R. R. Rogers; M. K. Yau (1989). A Short Course in Cloud
/// Physics (3rd ed.). Pergamon Press. p. 16. ISBN 0-7506-3215-1.
///
/// Returns: the latent heat of condensation for water in J kg<sup>-1</sup>.
#[inline]
pub fn latent_heat_of_condensation<T>(temperature: T) -> Option<JpKg>
where
    T: Temperature,
    Celsius: From<T>,
{
    debug_validate!(temperature);

    let Celsius(t) = Celsius::from(temperature);

    // The table has values from -40.0 to 40.0. So from -100.0 to -40.0 is actually an exrapolation.
    // I graphed the values from the extrapolation, and the curve looks good, and is approaching the
    // latent heat of sublimation, but does not exceed it. This seems very reasonable to me,
    // especially considering that a common approximation is to just us a constant value.
    if t < -100.0 || t > 60.0 {
        None
    } else {
        Some(JpKg(
            (2500.8 - 2.36 * t + 0.0016 * t * t - 0.00006 * t * t * t) * 1000.0,
        ))
    }
}

/// Virtual temperature.
///
/// From the [Glossary of Meteorology].(http://glossary.ametsoc.org/wiki/Virtual_temperature)
///
/// Returns the virtual temperature.
#[inline]
pub fn virtual_temperature<T, DP, P>(temperature: T, dew_point: DP, pressure: P) -> Option<Kelvin>
where
    T: Temperature,
    DP: Temperature,
    P: Pressure,
    Celsius: From<T> + From<DP>,
    HectoPascal: From<P>,
{
    debug_validate!(temperature, dew_point, pressure);

    let t = Celsius::from(temperature);
    let dp = Celsius::from(dew_point);
    let p = HectoPascal::from(pressure);

    let rv = mixing_ratio::<Celsius, HectoPascal>(dp, p)?;
    let Kelvin(t_k) = potential_temperature::<HectoPascal, Celsius>(p, t);
    let vt_k = Kelvin(t_k * (1.0 + rv / epsilon) / (1.0 + rv));
    Some(temperature_from_pot_temp::<HectoPascal, _>(vt_k, p))
}

/// Find the root of an equation given values bracketing a root. Used when finding wet bulb
/// temperature among other functions.
fn find_root(f: &dyn Fn(f64) -> Option<f64>, mut a: f64, mut b: f64) -> Option<f64> {
    use std::f64;
    const MAX_IT: usize = 50;
    const EPS: f64 = 1.0e-9;

    let mut fa = f(a)?;
    let mut fb = f(b)?;

    // Check to make sure we are bracketed
    if fa * fb >= 0.0 {
        return None;
    }

    if fa.abs() < fb.abs() {
        ::std::mem::swap(&mut a, &mut b);
        ::std::mem::swap(&mut fa, &mut fb);
    }

    let mut c = a;
    let mut fc = fa;
    let mut d = c; // This value is never used, mflag ensures it will be initialized in the loop.
    let mut mflag = true;

    for _ in 0..MAX_IT {
        let mut s = if fa != fc && fb != fc {
            // Try inverse quadratic for next step
            (a * fb * fc) / ((fa - fb) * (fa - fc))
                + (b * fa * fc) / ((fb - fa) * (fb - fc))
                + (c * fa * fb) / ((fc - fa) * (fc - fb))
        } else {
            // Try secant method for next step
            b - fb * (b - a) / (fb - fa)
        };

        // Check to see if bisection would be a better idea
        let condition1 = if a < b {
            s > b || s < (3.0 * a + b) / 4.0
        } else {
            s < b || s > (3.0 * a + b) / 4.0
        };

        let condition2 = mflag && (s - b).abs() >= (b - c).abs() / 2.0;
        let condition3 = !mflag && (s - b).abs() >= (c - d).abs() / 2.0;
        let condition4 = mflag && (b - c).abs() < EPS;
        let condition5 = !mflag && (c - d).abs() < EPS;

        if condition1 || condition2 || condition3 || condition4 || condition5 {
            s = (a + b) / 2.0;
            mflag = true;
        } else {
            mflag = false;
        }

        let fs = f(s)?;
        d = c;
        c = b;
        fc = fb;

        if fa * fs < 0.0 {
            b = s;
            fb = fs;
        } else {
            a = s;
            fa = fs;
        }

        if fa.abs() < fb.abs() {
            std::mem::swap(&mut a, &mut b);
            std::mem::swap(&mut fa, &mut fb);
        }

        // Check for convergence and return
        if fb == 0.0 || (b - a).abs() < EPS {
            return Some(b);
        }
    }

    None
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::test_utils::{approx_equal, approx_lte, pressure_levels, temperatures};

    const TOL: f64 = 1.0e-9;

    #[test]
    fn test_theta_kelvin() {
        for celsius in temperatures() {
            let kelvin = Kelvin::from(celsius);
            let theta = potential_temperature(HectoPascal(1000.0), celsius);
            assert!(approx_equal(kelvin, theta, CelsiusDiff(TOL)));
        }
    }

    #[test]
    fn test_temperature_c_from_theta() {
        for celsius in temperatures() {
            let kelvin = Kelvin::from(celsius);
            let back_to_celsius = temperature_from_pot_temp(kelvin, HectoPascal(1000.0));

            assert!(approx_equal(celsius, back_to_celsius, CelsiusDiff(TOL)));
        }
    }

    #[test]
    fn test_temperature_c_from_theta_and_back_by_theta_kelvin() {
        for p in pressure_levels() {
            for celsius in temperatures() {
                let theta = potential_temperature(p, celsius);
                let back_to_celsius = Celsius::from(temperature_from_pot_temp(theta, p));
                assert!(approx_equal(celsius, back_to_celsius, CelsiusDiff(TOL)));
            }
        }
    }

    #[test]
    fn test_vapor_pressure_liquid_water_there_and_back() {
        for dp in [
            -80.0, -60.0, -40.0, -20.0, -10.0, 0.0, 10.0, 20.0, 40.0, 49.0,
        ]
        .iter()
        .map(|&dp| Celsius(dp))
        {
            let forward = vapor_pressure_water(dp).unwrap();
            let back = dew_point_from_vapor_pressure_water(forward).unwrap();
            assert!(approx_equal(dp, back, CelsiusDiff(TOL)));
        }
    }

    #[test]
    fn test_vapor_pressure_ice_there_and_back() {
        for dp in [-80.0, -60.0, -40.0, -20.0, -10.0, -5.0, 0.0]
            .iter()
            .map(|&dp| Celsius(dp))
        {
            let forward = vapor_pressure_ice(dp).unwrap();
            let back = frost_point_from_vapor_pressure_over_ice(forward).unwrap();
            assert!(approx_equal(dp, back, CelsiusDiff(TOL)));
        }
    }

    #[test]
    fn test_rh() {
        for t in temperatures() {
            for dp in temperatures().filter(|dp| *dp <= t) {
                if let Some(rh_result) = rh(t, dp) {
                    assert!(rh_result > 0.0, "Negative RH!");
                    assert!(rh_result <= 1.0);
                    if t == dp {
                        assert!(rh_result == 1.0);
                    }
                }
            }
        }
    }

    #[test]
    fn test_rh_ice() {
        for t in temperatures() {
            for dp in temperatures().filter(|dp| *dp <= t) {
                if let Some(rh_result) = rh_ice(t, dp) {
                    assert!(rh_result > 0.0, "Negative RH!");
                    assert!(rh_result <= 1.0);
                    if t == dp {
                        assert!(rh_result == 1.0);
                    }
                }
            }
        }
    }

    #[test]
    fn test_mixing_ratio() {
        for p in pressure_levels() {
            for dp in temperatures() {
                if let Some(vp) = vapor_pressure_water(dp) {
                    if vp >= p {
                        continue;
                    }

                    let mw = mixing_ratio(dp, p).unwrap();
                    assert!(mw > 0.0);
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
                    Some(mw) => {
                        if let Some(back) = dew_point_from_p_and_mw(press, mw) {
                            println!(
                                "{:?} == {:?} with tolerance <= {:?}",
                                dp,
                                back,
                                f64::abs(dp.unwrap() - back.unwrap())
                            );
                            assert!(approx_equal(dp, back, CelsiusDiff(1.0e-2)));
                        }
                    }
                    None => { /* Ignore this for now.*/ }
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
                        Some((p_lcl, t_lcl)) => {
                            let k = Kelvin::from(t);
                            assert!(p_lcl <= p);
                            if t == dp {
                                assert!(approx_equal(t_lcl, k, CelsiusDiff(0.001)));
                            }

                            // Check for valid values
                            p_lcl.unwrap();
                            t_lcl.unwrap();
                        }
                        None => { /* Ignore these for now. */ }
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
                    match pressure_at_lcl(t, dp, p) {
                        Some(p_lcl) => {
                            assert!(p_lcl <= p);
                            p_lcl.unwrap();
                        }
                        None => { /* Ignore */ }
                    }
                }
            }
        }
    }

    #[test]
    fn test_specific_humidity() {
        for p in pressure_levels() {
            for dp in temperatures() {
                if let Some(sh) = specific_humidity(dp, p) {
                    assert!(sh > 0.0 && sh <= 1.0);
                }
            }
        }
    }

    #[test]
    fn test_specific_humidity_and_back_by_dew_point_from_p_and_specific_humidity() {
        for press in pressure_levels() {
            for dp in temperatures() {
                let q = specific_humidity(dp, press);

                match q {
                    Some(q) => {
                        if let Some(back) = dew_point_from_p_and_specific_humidity(press, q) {
                            println!(
                                "{:?} == {:?} with tolerance <= {:?}",
                                dp,
                                back,
                                f64::abs(dp.unwrap() - back.unwrap())
                            );
                            assert!(approx_equal(dp, back, CelsiusDiff(1.0e-2)));
                        }
                    }
                    None => { /* Ignore this for now.*/ }
                }
            }
        }
    }

    #[test]
    fn test_theta_e() {
        for p in pressure_levels() {
            for t in temperatures() {
                for dp in temperatures() {
                    match equiv_pot_temperature(t, dp, p) {
                        Some(theta_e) => {
                            let theta = potential_temperature(p, t);
                            assert!(theta <= theta_e);
                        }
                        None => { /* Ignore forwarded errors. */ }
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
                    if let Some((theta, theta_e, theta_es)) = equiv_pot_temperature(t, dp, p)
                        .and_then(|theta_e_val| {
                            Some((theta_e_val, equiv_pot_temperature(t, t, p)?))
                        })
                        .and_then(|pair| Some((potential_temperature(p, t), pair.0, pair.1)))
                    {
                        println!("{:?} <= {:?} <= {:?}", theta, theta_e, theta_es);
                        assert!(
                            approx_lte(theta, theta_e, CelsiusDiff(TOL))
                                && approx_lte(theta_e, theta_es, CelsiusDiff(TOL))
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn test_specific_humidity_and_mixing_ratio_conversions() {
        for sh in (0..1000).map(|i| i as f64 / 10_000.0) {
            let mr = mixing_ratio_from_specific_humidity(sh);
            let sh_back = specific_humidity_from_mixing_ratio(mr);
            assert!(approx_equal(sh, sh_back, 1.0e-10));
        }
    }

    #[test]
    fn test_theta_e_saturated_there_and_back() {
        for p in pressure_levels() {
            for t in temperatures() {
                if let Some(t_back) = equiv_pot_temperature(t, t, p).and_then(|theta_es| {
                    temperature_from_equiv_pot_temp_saturated_and_pressure(p, theta_es)
                }) {
                    println!(
                        "{:?} {:?} {:?}",
                        t,
                        t_back,
                        (t.unwrap() - t_back.unwrap()).abs()
                    );
                    assert!(approx_equal(t, t_back, CelsiusDiff(1.0e-6)));
                }
            }
        }
    }

    #[test]
    fn test_wet_bulb_dp_temp_consistency() {
        for press in [1000.0, 906.4, 850.0, 700.0, 314.6]
            .iter()
            .map(|&p| HectoPascal(p))
        {
            for t in [-45.46, -20.0, -10.0, -5.0, -1.46, 0.0, 5.0, 10.0, 20.0]
                .iter()
                .map(|&t| Celsius(t))
            {
                for dp in [-59.49, -21.0, -11.0, -5.1, -3.06, 0.0, 4.9, 9.9, 20.0]
                    .iter()
                    .map(|&t| Celsius(t))
                {
                    if dp > t {
                        continue;
                    }
                    if let Some(wb) = wet_bulb(t, dp, press) {
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
        let t_c = [0.0, 10.0, 20.0, 25.0, 30.0].iter().map(|&t| Celsius(t));
        let dp_c = [-15.0, -15.0, 18.0, 19.0, 19.0].iter().map(|&t| Celsius(t));
        let p_hpa = [1000.0, 1000.0, 900.0, 875.0, 875.0]
            .iter()
            .map(|&p| HectoPascal(p));
        let vt = [0.2, 10.21, 22.57, 27.86, 32.91]
            .iter()
            .map(|&t| Celsius(t));

        for (t, (dp, (p, target))) in t_c.zip(dp_c.zip(p_hpa.zip(vt))) {
            println!(
                "target = {:?}, value = {:?}",
                target,
                virtual_temperature(t, dp, p).unwrap()
            );
            assert!(approx_equal(
                target,
                virtual_temperature(t, dp, p).unwrap(),
                CelsiusDiff(0.05)
            ));
        }
    }

    #[test]
    fn test_find_root() {
        assert!(approx_equal(
            1.0,
            find_root(&|x| Some(x * x - 1.0), 2.0, 0.0).unwrap(),
            1.0e-8
        ));
        assert!(approx_equal(
            -1.0,
            find_root(&|x| Some(x * x - 1.0), -2.0, 0.0).unwrap(),
            1.0e-8
        ));
    }
}
