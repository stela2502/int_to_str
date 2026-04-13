pub mod int_to_str;

#[cfg(feature = "alignment")]
pub mod alignement;

#[allow(unused_imports)]
pub use int_to_str::IntToStr as IntToStr;

#[cfg(feature = "alignment")]
pub use alignement::{AlignmentConfig, AlignmentResult, ExactOverlap};