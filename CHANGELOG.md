# Changelog

## [Unreleased]

### Added

- Added support of non-adiabatic couplings in `Datafile.getdata`.
- New character function `timestamp` to return a formatted string of the current date and time.
- New "real work precision" (`realwp`) provided by module `numeric`.  The precision can be set through preprocessor directives (`-DUSE_R4`, `-DUSE_R8`).
- Core and data modules that do not provide versions for different real precisions use `realwp` internally.
- New method `get_data` for `DataFile` provides a simple way to extract data selectively from data files.
- Added dependency on xmake `openblas` package to support calls to BLAS/LAPACK routines.
- New exception `QuantityError` and the associated procedure `RaiseQuantityError` for errors related to the processing of data/quantities.
- New module `vibdata` to store vibrational data extracted from input files.
- New procedure `build_vibdata` to extract vibrational data.
- The library can now support data sets for multiple molecules or states.
- Extraction of electronic energies of current states (`moleculeDB` block) and pure excited-state energies (`ExcitationDB` block) from Gaussian fchk file.

### Changed

- `DataFile.getdata` extracts explicitly the number of atoms internally for derivatives, since the transition scalar with the total number of derivatives may be missing in some older chk/fchk.
- Software-specific parsing features are moved to new directory `parsers/`; this includes the parsing tools and the interfaces for `input`.
- `input_extract` is renamed `input_data` for better clarity.
- Module `transdata` is renamed `excdata` for consistency.
- `build_moldata` has optional arguments to control which data have to be loaded, grouped by categories: molecular specification, basis set specifications, molecular orbitals.
- Data effectively loaded in the `moldata` module can be probed with logical variables of the form `dat_XXXX_loaded`, with XXXX=spec, morb, bset.
- Modules `moldata`, `vibdata` and `excdata` have been dropped.  Data are now stored in derived types, contained in `workdata`.
- Improved error messages when the number of values given to a command-line argument does not match requirements.
