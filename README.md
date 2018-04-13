[![Build Status](https://travis-ci.org/rnleach/metfor.svg?branch=master)](https://travis-ci.org/rnleach/metfor)
[![Build status](https://ci.appveyor.com/api/projects/status/aqsy0e5n90mw3lwh/branch/master?svg=true)](https://ci.appveyor.com/project/rnleach/metfor/branch/master)
[![Latest Version](https://img.shields.io/crates/v/metfor.svg)](https://crates.io/crates/metfor)
[![docs](https://docs.rs/metfor/badge.svg)](https://docs.rs/metfor)

# metfor

Meteorological constants and formulas.

I investigated using some sort of dimensional analysis via types with a crate like [uom][uom]
or [dimensioned][dimensioned]. However after experimentation, neither of these work well for a
library. Choosing to use one would force that library on the users of this library. In the
future I may make a feature to use one or the other of these crates.

[uom]: https://crates.io/crates/uom
[dimensioned]: https://crates.io/crates/dimensioned
