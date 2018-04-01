extern crate metfor;
extern crate sounding_analysis;
extern crate sounding_base;
extern crate sounding_validate;

#[macro_use]
mod utils;

test_file!(standard, "standard.csv");
test_file!(complex_dendritic, "complex_dendritic.csv");
