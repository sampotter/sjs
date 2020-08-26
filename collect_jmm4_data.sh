#!/usr/bin/env sh

BUILD_DIR=./build
OUTPUT_DIR=./data

collect_cell_marching_data () {
    pmin=$1
    pmax=$2
    for slo_fun in 1 p v m g
    do
        echo "- slo_fun: ${slo_fun}"
        mkdir -p $OUTPUT_DIR/$slo_fun
        for p in {4..12}
        do
            echo "  * p = ${p}"
            N=$((2**$p + 1))
            $BUILD_DIR/cell_marching $slo_fun $N || exit 1
            for bin_path in $(ls *_$N.bin)
            do
                mv $bin_path $OUTPUT_DIR/$slo_fun
            done
        done
    done
}

collect_cell_marching_data
