extern crate metfor;
extern crate optional;
extern crate sounding_analysis;
extern crate sounding_base;
extern crate sounding_validate;

#[macro_use]
mod utils;

test_file!(standard, "standard.csv");
test_file!(complex_dendritic, "complex_dendritic.csv");
test_file!(multiple_warm_layers, "multiple_warm_layers_aloft.csv");
test_file!(multiple_inversions_aloft, "multiple_inversions_aloft.csv");
