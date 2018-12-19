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

// No need to use implOpsForQuantity macro here, all those ops should already be defined for
// f64!
