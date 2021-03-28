//! Run these benches with `cargo bench --bench funtions -- --verbose`
//!
//! Run with `cargo bench --bench funtions -- --verbose vapor_pressure_over_ice` to select the
//! single benchmark.

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use metfor::*;

criterion_main!(
    find_root_benches,
    potential_temperature,
    vapor_pressure_liquid,
    vapor_pressure_ice_benches,
    theta_e_benches
);

/**************************************************************************************************
 *                                     Find Root Group
 *************************************************************************************************/
// These test functions that depend on the internal, private function find_root.
criterion_group!(
    find_root_benches,
    temperature_from_theta_e_saturated_and_pressure_bench,
    pressure_at_lcl_bench,
    pressure_and_temperature_at_lcl_bench,
    wet_bulb_bench
);

fn temperature_from_theta_e_saturated_and_pressure_bench(c: &mut Criterion) {
    let pressures: Vec<_> = (300..1101)
        .step_by(100)
        .map(|i| HectoPascal(i as f64))
        .collect();
    let theta_es: Vec<_> = (260..321).step_by(10).map(|i| Kelvin(i as f64)).collect();

    c.bench_function("temperature_from_theta_e_saturated_and_pressure", |b| {
        b.iter(|| {
            for p in &pressures {
                for theta_e in &theta_es {
                    temperature_from_theta_e_saturated_and_pressure(
                        black_box(*p),
                        black_box(*theta_e),
                    );
                }
            }
        })
    });
}

fn pressure_at_lcl_bench(c: &mut Criterion) {
    let pressures: Vec<_> = (700..1101)
        .step_by(100)
        .map(|i| HectoPascal(i as f64))
        .collect();
    let temperatures: Vec<_> = (-40..41).step_by(10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("pressure_at_lcl", |b| {
        b.iter(|| {
            for p in &pressures {
                for t in &temperatures {
                    for dp in temperatures.iter().filter(|dp| *dp <= t) {
                        pressure_at_lcl(black_box(*t), black_box(*dp), black_box(*p));
                    }
                }
            }
        })
    });
}

fn pressure_and_temperature_at_lcl_bench(c: &mut Criterion) {
    let pressures: Vec<_> = (700..1101)
        .step_by(100)
        .map(|i| HectoPascal(i as f64))
        .collect();
    let temperatures: Vec<_> = (-40..41).step_by(10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("pressure_and_temperature_at_lcl", |b| {
        b.iter(|| {
            for p in &pressures {
                for t in &temperatures {
                    for dp in temperatures.iter().filter(|dp| *dp <= t) {
                        pressure_and_temperature_at_lcl(
                            black_box(*t),
                            black_box(*dp),
                            black_box(*p),
                        );
                    }
                }
            }
        })
    });
}

fn wet_bulb_bench(c: &mut Criterion) {
    let pressures: Vec<_> = (700..1101)
        .step_by(100)
        .map(|i| HectoPascal(i as f64))
        .collect();
    let temperatures: Vec<_> = (-40..41).step_by(10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("wet_bulb", |b| {
        b.iter(|| {
            for p in &pressures {
                for t in &temperatures {
                    for dp in temperatures.iter().filter(|dp| *dp <= t) {
                        wet_bulb(black_box(*t), black_box(*dp), black_box(*p));
                    }
                }
            }
        })
    });
}

/**************************************************************************************************
 *                                        Theta E Group
 *************************************************************************************************/
criterion_group!(
    theta_e_benches,
    theta_e_bench,
    latent_heat_of_condensation_bench,
    virtual_temperature_bench
);

fn theta_e_bench(c: &mut Criterion) {
    let pressures: Vec<_> = (700..1101)
        .step_by(100)
        .map(|i| HectoPascal(i as f64))
        .collect();
    let temperatures: Vec<_> = (-40..41).step_by(10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("theta_e", |b| {
        b.iter(|| {
            for p in &pressures {
                for t in &temperatures {
                    for dp in temperatures.iter().filter(|dp| *dp <= t) {
                        theta_e(black_box(*t), black_box(*dp), black_box(*p));
                    }
                }
            }
        })
    });
}

fn latent_heat_of_condensation_bench(c: &mut Criterion) {
    let temperatures: Vec<_> = (-40..41).step_by(10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("latent_heat_of_condensation", |b| {
        b.iter(|| {
            for t in &temperatures {
                latent_heat_of_condensation(black_box(*t));
            }
        })
    });
}

fn virtual_temperature_bench(c: &mut Criterion) {
    let pressures: Vec<_> = (700..1101)
        .step_by(100)
        .map(|i| HectoPascal(i as f64))
        .collect();
    let temperatures: Vec<_> = (-40..41).step_by(10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("virtual_temperature", |b| {
        b.iter(|| {
            for p in &pressures {
                for t in &temperatures {
                    for dp in temperatures.iter().filter(|dp| *dp <= t) {
                        virtual_temperature(black_box(*t), black_box(*dp), black_box(*p));
                    }
                }
            }
        })
    });
}

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
    let temperatures: Vec<_> = (-40..11).step_by(10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("rh_ice", |b| {
        b.iter(|| {
            for t in &temperatures {
                for dp in temperatures.iter().filter(|&dp| dp <= t) {
                    rh_ice(black_box(*t), black_box(*dp));
                }
            }
        })
    });
}

fn frost_point_from_vapor_pressure_over_ice_bench(c: &mut Criterion) {
    let vps: Vec<_> = (0..121)
        .step_by(10)
        .map(|i| HectoPascal(i as f64))
        .collect();

    c.bench_function("frost_point_from_vapor_pressure_over_ice", |b| {
        b.iter(|| {
            for vp in &vps {
                frost_point_from_vapor_pressure_over_ice(black_box(*vp));
            }
        })
    });
}

fn vapor_pressure_ice_bench(c: &mut Criterion) {
    let temperatures: Vec<_> = (-90..11).step_by(10).map(|i| Celsius(i as f64)).collect();

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
    rh_bench,
    mixing_ratio_bench,
    dew_point_from_p_and_mw_bench,
    specific_humidity_bench,
);

fn specific_humidity_bench(c: &mut Criterion) {
    let pressures: Vec<_> = (500..1101)
        .step_by(100)
        .map(|i| HectoPascal(i as f64))
        .collect();
    let temperatures: Vec<_> = (-40..41).step_by(10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("specific_humidity", |b| {
        b.iter(|| {
            for p in &pressures {
                for t in &temperatures {
                    specific_humidity(black_box(*t), black_box(*p));
                }
            }
        })
    });
}

fn dew_point_from_p_and_mw_bench(c: &mut Criterion) {
    let pressures: Vec<_> = (500..1050).map(|i| HectoPascal(i as f64)).collect();
    let mws: Vec<_> = (0..700).step_by(10).map(|i| (i as f64) / 100.0).collect();

    c.bench_function("dew_point_from_p_and_mw", |b| {
        b.iter(|| {
            for p in &pressures {
                for mw in &mws {
                    dew_point_from_p_and_mw(black_box(*p), black_box(*mw));
                }
            }
        })
    });
}

fn mixing_ratio_bench(c: &mut Criterion) {
    let pressures: Vec<_> = (500..1101)
        .step_by(100)
        .map(|i| HectoPascal(i as f64))
        .collect();
    let temperatures: Vec<_> = (-40..41).step_by(10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("mixing_ratio", |b| {
        b.iter(|| {
            for p in &pressures {
                for t in &temperatures {
                    mixing_ratio(black_box(*t), black_box(*p));
                }
            }
        })
    });
}

fn rh_bench(c: &mut Criterion) {
    let temperatures: Vec<_> = (-40..41).step_by(10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("rh", |b| {
        b.iter(|| {
            for t in &temperatures {
                for dp in temperatures.iter().filter(|&dp| dp <= t) {
                    rh_ice(black_box(*t), black_box(*dp));
                }
            }
        })
    });
}

fn dew_point_from_vapor_pressure_over_liquid_bench(c: &mut Criterion) {
    let vps: Vec<_> = (0..121)
        .step_by(10)
        .map(|i| HectoPascal(i as f64))
        .collect();

    c.bench_function("dew_point_from_vapor_pressure_over_liquid", |b| {
        b.iter(|| {
            for vp in &vps {
                dew_point_from_vapor_pressure_over_liquid(black_box(*vp));
            }
        })
    });
}

fn vapor_pressure_liquid_water_bench(c: &mut Criterion) {
    let dps: Vec<_> = (-90..61).step_by(10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("vapor_pressure_liquid_water", |b| {
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
    let pressures: Vec<_> = (500..1101)
        .step_by(100)
        .map(|i| HectoPascal(i as f64))
        .collect();
    let temperatures: Vec<_> = (-40..41).step_by(10).map(|i| Celsius(i as f64)).collect();

    c.bench_function("theta", |b| {
        b.iter(|| {
            for p in &pressures {
                for t in &temperatures {
                    theta(black_box(*p), black_box(*t));
                }
            }
        })
    });
}

fn temperature_from_theta_bench(c: &mut Criterion) {
    let pressures: Vec<_> = (500..1101)
        .step_by(100)
        .map(|i| HectoPascal(i as f64))
        .collect();
    let thetas: Vec<_> = (250..331).step_by(10).map(|i| Kelvin(i as f64)).collect();

    c.bench_function("temperature_from_theta", |b| {
        b.iter(|| {
            for p in &pressures {
                for t in &thetas {
                    theta(black_box(*p), black_box(*t));
                }
            }
        })
    });
}
