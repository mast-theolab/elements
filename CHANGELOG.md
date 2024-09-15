# Changelog

## [Unreleased]

### Added

- New exception `QuantityError` and the associated procedure `RaiseQuantityError` for errors related to the processing of data/quantities.
- New module `vibdata` to store vibrational data extracted from input files.
- New procedure `build_vibdata` to extract vibrational data.
- The library can now support data sets for multiple molecules or states.

### Changed

- Module `transdata` is renamed `excdata` for consistency.
- `build_moldata` has optional arguments to control which data have to be loaded, grouped by categories: molecular specification, basis set specifications, molecular orbitals.
- Data effectively loaded in the `moldata` module can be probed with logical variables of the form `dat_XXXX_loaded`, with XXXX=spec, morb, bset.
- Modules `moldata`, `vibdata` and `excdata` have been dropped.  Data are now stored in derived types, contained in `workdata`.
- Improved error messages when the number of values given to a command-line argument does not match requirements.