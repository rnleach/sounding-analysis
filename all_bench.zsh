#!/usr/bin/env zsh

BASELINE_TAG="opt-start"
COMPARE_TAG="master"

BENCHES=("indexes" "layers" "levels" "parcel" "parcel_profile" "profile" "wind" "fire")

echo "Cleaning up...."
cargo clean

function do_bench {
    
    # $1 git tag
    # $2 target-cpu, either "generic" or "native"
    
    RUSTFLAGS_VAL="-C target-cpu=$2"
    BASELINE="$1-$2"

    echo "Testing tag $1 with target-cpu=$2"
    echo "RUSTFLAGS=${RUSTFLAGS_VAL} and BASELINE=${BASELINE}"

    git checkout $1
    if [ $? -ne 0 ]; then
        echo "Git failed, aborting."
        exit 1
    fi

    if [ $# -gt 2 ]; then
            echo ""
            echo "Doing bench $3"
            echo ""

            RUSTFLAGS="${RUSTFLAGS_VAL}" cargo bench --bench $3 -- --save-baseline ${BASELINE}
    else
        for BNAME in ${BENCHES}; do
            echo ""
            echo "Doing bench ${BNAME}"
            echo ""

            RUSTFLAGS="${RUSTFLAGS_VAL}" cargo bench --bench ${BNAME} -- --save-baseline ${BASELINE}
        done
    fi

    echo "Completed baseline ${BASELINE}"

    return 0
}

if [ $# -eq 0 ]; then
    do_bench ${BASELINE_TAG}  "generic"
    do_bench ${BASELINE_TAG}  "native"
    do_bench ${COMPARE_TAG}  "generic"
    do_bench ${COMPARE_TAG}  "native"
else
    do_bench ${BASELINE_TAG}  "generic" $1
    do_bench ${BASELINE_TAG}  "native" $1
    do_bench ${COMPARE_TAG}  "generic" $1
    do_bench ${COMPARE_TAG}  "native" $1
fi

