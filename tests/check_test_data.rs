use sounding_validate;

#[macro_use]
mod utils;

#[test]
fn test_load_test_csv_sounding() {
    let (_, ivals, fvals) = utils::load_test_file("standard.csv");

    assert!(Some(&1) == ivals.get("num dendritic zones"));
    assert!(Some(&0) == ivals.get("num warm dry bulb aloft"));
    assert!(Some(&0) == ivals.get("num warm wet bulb aloft"));
    assert!(Some(&0) == ivals.get("num inversions"));

    assert!(Some(&vec![604.2, 534.7]) == fvals.get("dendritic zone pressures"));
    assert!(Some(&vec![1080.0, 508.5]) == fvals.get("6km agl layer pressures"));
}

check_file_complete!(standard_file_complete, "standard.csv");
check_file_complete!(complex_dendritic_file_complete, "complex_dendritic.csv");
check_file_complete!(
    multiplie_warm_layers_aloft,
    "multiple_warm_layers_aloft.csv"
);
check_file_complete!(multiple_inversions_aloft, "multiple_inversions_aloft.csv");
