//! # iq_terminal — U10a (peak features, find_peaks/peak_widths, XGBoost) and U10b
//! (per-block polyA/TSS accumulation → predictions)
//!
//! Design §7.9. Fixture runners: unit `U10a` (fns `detect_peaks_batch`, `xgb_predict`) and
//! `U10b` (fn `flush`), see [`register_fixtures`].

pub mod counter;
pub mod peaks;
pub mod xgb;

pub use counter::block_predictions;
pub use peaks::{detect_peaks, detect_peaks_detailed, Histogram};
pub use xgb::{bundled_models, XgbModel, XgbModels};

/// IsoQuant 4.1 `isoquant_lib/data/model_polya.json` (bundled, design §7.9).
pub const MODEL_POLYA_JSON: &[u8] = include_bytes!("../data/model_polya.json");
/// IsoQuant 4.1 `isoquant_lib/data/model_tss.json` (bundled, design §7.9).
pub const MODEL_TSS_JSON: &[u8] = include_bytes!("../data/model_tss.json");

/// Registers this crate's fixture runners.
pub fn register_fixtures(reg: &mut iq_types::fixture::Registry) {
    reg.add::<peaks::fixture::U10aDetectPeaks>();
    reg.add::<xgb::fixture::U10aXgbPredict>();
    counter::register_fixtures(reg);
}
