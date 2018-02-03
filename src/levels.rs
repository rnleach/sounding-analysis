//! This module finds significant levels such as the freezing level and wet bulb zero level. It also
//! has functions for finding critical values at a single level, such as the maximum wet bulb
//! temperature aloft.  It does not include functions for finding levels related to parcel analysis
//! and convection, those are found in the `parcel` module.

use error::*;
use smallvec::SmallVec;

use sounding_base::Sounding;
use sounding_base::Profile::*;

// TODO: Wet bulb zero height Return multiple if needed.
// TODO: Freezing. Return multiple if needed.
// TODO: Max Tw aloft.
