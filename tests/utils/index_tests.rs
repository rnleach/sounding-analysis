use super::*;
use metfor::Quantity;
use sounding_analysis::Result;
use sounding_base::Sounding;

#[allow(dead_code)] // False alarm - lint is done before macro expansion.
pub fn test_index<F, Q>(
    snd: &Sounding,
    tgt_float_vals: &HashMap<String, Vec<f64>>,
    anal_func: F,
    index_key: &str,
    tol: Q,     // tolerance
    err_val: Q, // Value if result would return Err(_)
) where
    F: FnOnce(&Sounding) -> Result<Q>,
    Q: Quantity,
{
    if let Some(target_vals) = tgt_float_vals.get(index_key) {
        assert_eq!(target_vals.len(), 1);
        let target_val = Q::pack(target_vals[0]);

        let analysis = anal_func(&snd);

        if let Err(anal_err) = analysis {
            println!("Analysis Error Value: {err:#?} => {err}", err = anal_err);
        }

        let analysis = analysis.unwrap_or(err_val);

        assert!(analysis.approx_eq(target_val, tol));
    }
}
