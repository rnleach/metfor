[![Build Status](https://travis-ci.org/rnleach/metfor.svg?branch=master)](https://travis-ci.org/rnleach/metfor)
[![Build status](https://ci.appveyor.com/api/projects/status/aqsy0e5n90mw3lwh/branch/master?svg=true)](https://ci.appveyor.com/project/rnleach/metfor/branch/master)
[![Latest Version](https://img.shields.io/crates/v/metfor.svg)](https://crates.io/crates/metfor)
[![docs](https://docs.rs/metfor/badge.svg)](https://docs.rs/metfor)

# metfor


Meteorological constants and formulas.

This library includes wrapper types, or unit types, to help with using the proper units when doing
calculations, some common meteorological constants, and some functions for basic calculations. It
has thus far been developed to support libraries and an application used for displaying and
analyzing skew-t data. So you will find the functions and variables are mostly things that would
typically be used on a skew-t. Future versions may expand the intended use cases.

I investigated using some sort of dimensional analysis via types with a crate like [uom][uom] or
[dimensioned][dimensioned] instead of making my own _newtype_ types. However after experimentation,
I decided I did not want to make another crate part of the API since choosing to use one would force
that library on the users of this library.

I've found the [`optional`] [optional] crate to be very useful when dealing with `f64` types in
situations where there may be missing values, so I included a feature `use_optional` that enables
the newtypes in this crate to be used in the same manner as an `f64` is with `optional`.

[uom]: https://crates.io/crates/uom
[dimensioned]: https://crates.io/crates/dimensioned
[optional]: https://crates.io/crates/optional

## Examples
```rust
use metfor::{Kelvin, Celsius, HectoPascal, Millibar, theta};

let t1 = Kelvin(300.0);
let p1 = HectoPascal(1000.0);
let theta1 = theta(p1, t1);

let t2 = Celsius(0.0);
let p2 = Millibar(700.0);
let theta2 = theta(p2, t2);

println!("theta1 = {} and theta2 = {}", theta1, theta2);

```
