# Changelog

## 0.25.01

### Message

- Initial release

### Added

- Implementation of the core library.
- Module parse_cmdline provides a class `ArgObj` and related methods to process commandline arguments.  It works in ways similar to Python's `argparse` with less features.
- New module for electronic-structure calculations providing tools to compute 1-electron integrals and convert orbitals.
- New module to process basis-sets, allowing conversions between pure and Cartesian forms.
- Numerical-precision-related constants and procedures are stored in the module numeric.
- New module `physics` containing physical constants and conversion formulas.
- Atomic information are provided by the module `atominfo`.
- Module `datatypes` provides several derived types to act as container of simulation-related data (molecule specifications, basis sets, orbitals, electronic excitations, vibrations...).
- New parser module to parse Gaussian formatted checkpoint files.
- New drivers modules providing interfaces to common libraries like BLAS and LAPACK
- Implementation of the SOS formalism (P. Bour, Chem. Phys. Lett. 1997, 265, 65-70) to compute excited-excited transition moments.
- New program `mcd_tensor` to compute the MCD G tensor using the SOS formalism.
