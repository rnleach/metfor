//! Utilities for running unit tests.
#![macro_use]
use crate::types::*;
use std::ops::Sub;

pub fn approx_equal<L, R, T>(left: L, right: R, tol: T) -> bool
where
    T: Quantity + From<L> + From<R> + Sub,
    <T as Sub>::Output: Quantity,
{
    use std::f64;
    let tol = tol.unpack();
    let left = T::from(left);
    let right = T::from(right);
    assert!(tol > 0.0);
    let passes = f64::abs((left - right).unpack()) <= tol;

    if !passes {
        println!("{} !~= {} within tolerance {}", left, right, tol);
    }

    passes
}

pub fn approx_lte<L, R, T>(a: L, b: R, tol: T) -> bool
where
    T: Quantity + From<L> + From<R> + Sub,
    <T as Sub>::Output: PartialOrd<T> + Quantity,
{
    let a = T::from(a);
    let b = T::from(b);
    let passes = (a - b) <= tol;

    if !passes {
        println!("{} !~<= {} within tolerance {}", a, b, tol);
    }

    passes
}

pub struct DRange {
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

pub fn pressure_levels() -> impl Iterator<Item = HectoPascal> {
    DRange {
        start: 1000.0,
        step: -10.0,
        stop: 100.0,
    }
    .map(HectoPascal)
}

pub fn temperatures() -> impl Iterator<Item = Celsius> {
    DRange {
        start: -100.0,
        step: 10.0,
        stop: 100.0,
    }
    .map(Celsius)
}
