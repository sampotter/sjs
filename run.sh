#!/usr/bin/env sh

PYTHON=python

BUILD_DIR=./build
OUTPUT_DIR=./data
PLOTS_DIR=./plots

mkdir -p $BUILD_DIR
cd $BUILD_DIR
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd ..

mkdir -p $OUTPUT_DIR
$BUILD_DIR/JMMmeasurements $OUTPUT_DIR
$BUILD_DIR/olim8mp0 $OUTPUT_DIR
$BUILD_DIR/FMM $OUTPUT_DIR

mkdir -p $PLOTS_DIR
$PYTHON make_plots.py
$PYTHON make_tables.py > $PLOTS_DIR/table.tex
