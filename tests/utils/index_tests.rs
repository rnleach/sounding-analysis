use super::*;
use sounding_analysis::Result;
use sounding_base::Sounding;

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
pub fn test_index<F: FnOnce(&Sounding) -> Result<f64>>(
    snd: &Sounding,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
    anal_func: F,
    index_key: &str,
    err_val: f64, // Value if result would return Err(_)
) {
    if let Some(target_vals) = tgt_float_vals.get(index_key) {
        assert_eq!(target_vals.len(), 1);
        let target_val = target_vals[0];

        let analysis = anal_func(&snd).unwrap_or(err_val);
        assert!(approx_equal(analysis, target_val, 0.5));
    }
}
