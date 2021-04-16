#!/usr/bin/env zsh -l

setopt verbose

BASELINE_TAG = "opt-start"
COMPARE_TAG = "master"

BENCHES = (indexes layers levels parcel parcel_profile profile wind fire)


# Check out the baseline tag.
git checkout ${BASELINE_TAG}

# Force generic compile target.
export RUSTFLAGS="-C target-cpu=generic"

# Do a warm up run.
cargo bench

# Do the real benches
BASELINE="${BASELINE_TAG}-generic"
for bench_name in ${BENCHES}; do
    cargo bench --bench ${bench_name} -- --save-baseline ${BASELINE}
done

# Force to compile to target CPU
export RUSTFLAGS="-C target-cpu=native"

# Do a warm up run.
cargo bench

# Do the real benches
BASELINE="${BASELINE_TAG}-native"
for bench_name in ${BENCHES}; do
    cargo bench --bench ${bench_name} -- --save-baseline ${BASELINE}
done

# Check out the comparison tag
git checkout ${COMPARE_TAG}

# Force generic compile target.
export RUSTFLAGS="-C target-cpu=generic"

# Do a warm up run.
cargo bench

# Do the real benches
BASELINE="${COMPARE_TAG}-generic"
for bench_name in ${BENCHES}; do
    cargo bench --bench ${bench_name} -- --save-baseline ${BASELINE}
done

# Force to compile to target CPU
export RUSTFLAGS="-C target-cpu=native"

# Do a warm up run.
cargo bench

# Do the real benches
BASELINE="${COMPARE_TAG}-native"
for bench_name in ${BENCHES}; do
    cargo bench --bench ${bench_name} -- --save-baseline ${BASELINE}
done


