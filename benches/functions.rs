//! Run these benches with `cargo bench --bench funtions -- --verbose`
//!
//! Run with `cargo bench --bench funtions -- --verbose vapor_pressure_over_ice` to select the 
//! single benchmark.

use criterion::{black_box, criterion_group, criterion_main, BenchmarkId, Criterion};
use metfor::*;

criterion_main!(potential_temperature, vapor_pressure_liquid, vapor_pressure_ice_benches);

/**************************************************************************************************
 *                                  Vapor Pressure (Ice) Group
 *************************************************************************************************/
criterion_group!(
    vapor_pressure_ice_benches,
    vapor_pressure_ice_bench,
    frost_point_from_vapor_pressure_over_ice_bench,
    rh_ice_bench
);

fn rh_ice_bench(c: &mut Criterion) {

    c.bench_function("rh_ice", |b| {
        b.iter(|| {
                rh_ice(black_box(Celsius(-22.0)), black_box(Celsius(-40.0)));
        })
    });
}

fn frost_point_from_vapor_pressure_over_ice_bench(c: &mut Criterion) {

    let vps: Vec<_> = (0..120).map(|i| HectoPascal(i as f64)).collect();

    c.bench_function("dew_point_from_vapor_pressure_over_water", |b| {
        b.iter(|| {
            for vp in &vps {
                frost_point_from_vapor_pressure_over_ice(black_box(*vp));
            }
        })
    });
}

fn vapor_pressure_ice_bench(c: &mut Criterion) {
    let temperatures: Vec<_> = (-90..10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("vapor_pressure_ice", |b| {
        b.iter(|| {
            for temperature in &temperatures {
                vapor_pressure_ice(black_box(*temperature));
            }
        })
    });
}

/**************************************************************************************************
 *                                  Vapor Pressure (Liquid) Group
 *************************************************************************************************/
criterion_group!(
    vapor_pressure_liquid,
    vapor_pressure_liquid_water_bench,
    dew_point_from_vapor_pressure_over_liquid_bench,
    rh_bench
);

fn rh_bench(c: &mut Criterion) {

    c.bench_function("rh", |b| {
        b.iter(|| {
                rh(black_box(Celsius(22.0)), black_box(Celsius(-1.0)));
        })
    });
}

fn dew_point_from_vapor_pressure_over_liquid_bench(c: &mut Criterion) {

    let vps: Vec<_> = (0..120).map(|i| HectoPascal(i as f64)).collect();

    c.bench_function("dew_point_from_vapor_pressure_over_water", |b| {
        b.iter(|| {
            for vp in &vps {
                dew_point_from_vapor_pressure_over_liquid(black_box(*vp));
            }
        })
    });
}

fn vapor_pressure_liquid_water_bench(c: &mut Criterion) {
    let dps: Vec<_> = (-90..60).map(|i| Celsius(i as f64)).collect();

    c.bench_function("vapor_pressure_liquid", |b| {
        b.iter(|| {
            for dew_point in &dps {
                vapor_pressure_liquid_water(black_box(*dew_point));
            }
        })
    });
}

/**************************************************************************************************
 *                                          Theta Group
 *************************************************************************************************/
criterion_group!(
    potential_temperature,
    theta_bench,
    temperature_from_theta_bench
);

fn theta_bench(c: &mut Criterion) {
    let p = HectoPascal(700.0);
    let t = Celsius(-2.8);

    let args = (p, t);

    c.bench_with_input(
        BenchmarkId::new("theta", "700hPa and -2.8C"),
        &args,
        |b, &(p, t)| b.iter(|| theta(black_box(p), black_box(t))),
    );
}

fn temperature_from_theta_bench(c: &mut Criterion) {
    let pressure = HectoPascal(700.0);
    let theta = Kelvin(298.0);
    let args = (theta, pressure);

    c.bench_with_input(
        BenchmarkId::new("temperature from theta", "700hPa and 298K"),
        &args,
        |b, &(theta, pressure)| {
            b.iter(|| temperature_from_theta(black_box(theta), black_box(pressure)))
        },
    );

    let pressure = HectoPascal(1000.0);
    let theta = Kelvin(312.0);
    let args = (theta, pressure);

    c.bench_with_input(
        BenchmarkId::new("temperature from theta", "1000hPa and 312"),
        &args,
        |b, &(theta, pressure)| {
            b.iter(|| temperature_from_theta(black_box(theta), black_box(pressure)))
        },
    );
}
