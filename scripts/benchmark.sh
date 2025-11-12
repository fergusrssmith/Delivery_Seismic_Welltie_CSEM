#!/usr/bin/env bash

# Benchmark of colt Linear Algebra library routines (dgemm = double general matrix multiply)

# All CAPITALISED variables imported from this file
source environment_setup.sh

echo "COLT BENCHMARKS"
java -classpath ${SEISMIC_INVERTER_HOME}/jars/colt.jar:${SEISMIC_INVERTER_HOME}/jars/concurrent.jar \
     cern.colt.matrix.bench.BenchmarkMatrix dgemm dense 1 2 0.99 false true 5 20 50 100 300 500

echo "JAVA STANDARD MATRICES (A[][] style) BENCHMARKS"

java -classpath ${SEISMIC_INVERTER_HOME}:${SEISMIC_INVERTER_HOME}/jars/ejml-0.25.jar au.csiro.JamesGunning.lib.Numerical.NativeArraySpeedTest
