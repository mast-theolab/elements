# Changelog

## [Unreleased]

### Added

- New exception `QuantityError` and the associated procedure `RaiseQuantityError` for errors related to the processing of data/quantities.

### Changed

- Module `transdata` is renamed `excdata` for consistency.
- `build_moldata` has optional arguments to control which data have to be loaded, grouped by categories: molecular specification, basis set specifications, molecular orbitals.
- Data effectively loaded in the `moldata` module can be probed with logical variables of the form `dat_XXXX_loaded`, with XXXX=spec, morb, bset.
