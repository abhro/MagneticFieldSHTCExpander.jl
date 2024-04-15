# MagneticFieldSHTCExpander

[![Build Status](https://github.com/abhro/MagneticFieldSHTCExpander.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/abhro/MagneticFieldSHTCExpander.jl/actions/workflows/CI.yml?query=branch%3Amain)

Create grids of the solar magnetic field based on provided Spherical Harmonic Transform Coefficients (SHTC)

## Physical motivation
See [Physical-theory.md](Physical-theory.md).

## To-do
- [ ] Implement caching for `assoc_legendre_func_table`
- [ ] Provide a `collect`-like function for magnetic field at many `(r, θ, φ)`
- [ ] Write tests to cover real-world SHTC files
