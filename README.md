# METAL

[![Build Status](https://travis-ci.org/statgen/METAL.svg?branch=master)](https://travis-ci.org/statgen/METAL)

METAL is a tool for the meta-analysis of genome-wide association studies. More information about METAL is available at https://genome.sph.umich.edu/wiki/METAL.

## Requirements

- CMake cross-platform make system

## Compilation

Download and unarchive latest stable version of METAL source code or clone this GitHub repository to obtain development version.
Then execute the following commands inside the directory with METAL source code:

```
mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
make test
make install
```

The compiled `metal` executable will be installed into `build/bin` directory.

## New features

- 2020-05-05
  - Implemented new options `EFFECT_PRINT_PRECISION` and `STDERR_PRINT_PRECISION` to control print precision of `Effect` and `StdErr` output columns, correspondingly. For example, setting `EFFECT_PRINT_PRECISION 12` will force METAL to print 12 digits after the decimal point in the `Effect` column. The default is to print 4 digits after the decimal point for both `Effect` and `StdErr`. Increasing precision is useful when you plan to use `Effect` and `StdErr` values computed by METAL in other downstream analyses.

- 2018-08-28
  - Implemented new option `TRACKPOSITIONS`. If enabled, METAL checks if chromosome and position of a variant match across studies. The `Chromosome` and `Position` fields are propagated to the meta-analysis result file (useful when generating Manhattan or LocusZoom plots).

- 2017-12-21
  - Correction for sample overlap in sample size weighted meta-analysis (developed by Sebanti Sengupta). First, METAL estimates the number of individuals that are common among two or more studies based on Z-statistics from each study. Then, METAL adjusts for sample overlap when calculating overall Z-statistics by correcting the weights with the estimated number of individuals in common. 

    To enable correction for sample overlap in your sample size weighted meta-analysis, use `OVERLAP ON` command (valid only with `SCHEME SAMPLESIZE`) before any `PROCESS` commands. By default, METAL uses Z-statistics <1 for esimating the number of individuals that are common among studies. To change this threshold, use `ZCUTOFF [number]` command. The `N` column in the result file will store sample size corrected for overlap between studies.

    More information on this specific method can be found on our wiki: https://genome.sph.umich.edu/wiki/METAL_Documentation#Sample_Overlap_Correction

## Documentation

Complete documentation is available at https://genome.sph.umich.edu/wiki/METAL_Documentation.
