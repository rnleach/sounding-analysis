[![Build Status](https://travis-ci.org/rnleach/sounding-analysis.svg?branch=master)](https://travis-ci.org/rnleach/sounding-analysis)
[![Build status](https://ci.appveyor.com/api/projects/status/jb5joubn8bendk7s?svg=true)](https://ci.appveyor.com/project/rnleach/sounding-analysis)
[![Latest Version](https://img.shields.io/crates/v/sounding-analysis.svg)](https://crates.io/crates/sounding-analysis)
[![docs](https://docs.rs/sounding-analysis/badge.svg)](https://docs.rs/sounding-analysis)

# sounding-analysis

Functions and data types for analyzing soundings from the
[sounding-base](https://github.com/rnleach/sounding-base.git) crate.

### Purpose
Provides analysis capabilities for the [sounding-base](https://github.com/rnleach/sounding-base.git)
crate.

### Error handling strategy
Right now error handling is mixed up. Some functions/methods return `Option` while others return
`Result`, and all missing values are indicated by `Option::None`. I need to clean this up so that
funciton/methods return `Result`s. Then when building a data structure up, if a function returns an
error result, I will set the value to `Option::None`. If the user wants to know why it is
`Option::None`, they will have to call the function/method and inspect the error.

