# VASP Scripts
Various Python files for post-processing of VASP calculations.

## xdatcar2xyz.py
Script for converting an XDATCAR file to an xyz file, accounting for changes in the lattice vectors (e.g., ISIF=3). 

## gofr.py
Class for computing the pair distribution function g(r) for any trajectory readable by ASE. 

## polyhedra.py
Class for characterizing the polyhedra network for any trajectory file readable by ASE.

## orbital_pdos.py
Script for extracting the orbital-projected pdos from the DOSCAR file (LORBIT = 11). Requires input file called "projection" - see example. 
