# PYCA #

PYCA is a SAGE based toolkit for the study of cellular automata.

The aim of the project is to be an extensible toolkit that covers practical uses of CA and theoretical analysis.

## Implemented features ##

This are the features that are already implemented, note that currently the processing speed is unacceptable for any practical problem size. The plan is to keep the the suboptimal implementation for verification of faster implementations.

Basic classes:
  * an 1D CA rule
  * an 1D CA lattice (for cyclic and bounded lattices)

Practical:
  * forward transition step

Theoretical:
  * computation of preimages
  * checking if a configuration is a "Garden of Eden"

## Planed features ##

Basic classes:
  * optimized binary 1D CA
  * optimized multi value 1D CA
  * 2D CA
  * nD CA

Practical:
  * compressed raster images (PNG) of 1D CA transitions
  * compressed video of 2D CA transitions

Theoretical:
  * visualization of attraction basins

SAGE integration:
  * packages for downloading
  * documentation inside source files
  * automatic regression test suite