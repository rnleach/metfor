use crate::types::Quantity;

impl Quantity for f64 {
    fn pack(val: f64) -> f64 {
        val
    }

    fn unpack(self) -> f64 {
        self
    }

    fn unwrap(self) -> f64 {
        self
    }

    fn into_option(self) -> Option<f64> {
        Some(self)
    }
}

// Added u8 so I could treat some simple integer indexes as quantities too.
impl Quantity for u8 {
    fn pack(val: f64) -> u8 {
        val as u8
    }

    fn unpack(self) -> f64 {
        self as f64
    }

    fn unwrap(self) -> f64 {
        self as f64
    }

    fn into_option(self) -> Option<f64> {
        Some(self as f64)
    }
}

// No need to use implOpsForQuantity macro here, all those ops should already be defined for
// f64 and u8!
