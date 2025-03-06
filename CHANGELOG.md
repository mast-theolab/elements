# Changelog

## 0.25.03

### Added
* An instance of `PhysFact`, `phys_conv`, is now provided by the `physics` module as well.
* New conversion method from atomic unit of period to wavenumbers provided as `PhysFact % au2cm1` in `physics`.
* `Eckart_orient` in `geometry` was rewritten to be more robust regarding the overwriting of geometry arrays.  The default functions are now supposed to avoid any side effects, with new and old geometries stored in separate structures, and new routines were added to consent overriding the old data set.
* New interface to LAPACK sorting routines `xLASRT` for quick sorting of lists of real numbers.
* New interfaces to LAPACK routines dedicated to QR factorization, `xGEQRF` and `xORGQR`.
* Added support of potential energy in `PropertyDB` (`datatypes`).
* Added extraction of potential energy from Gaussian formatted checkpoint file in `get_data` (`input`), as electronic transition moments or state-specific quantities, including derivatives.
* Added operator `.x.` for cross product in `math`.
* Cross product can now be applied to list(s) of vectors (`math`).
* Derived type `ErrorHandler` (`exception`) provides a structure to manage runtime errors or warnings.  It can also manage automatic termination in case of error.
* Module `exception` provides a variable, `runstat` to keep track of the run status between modules and program units.
* Added conversion procedure bound to `PhysFact` (`physics`) to convert derivatives of energies with respect to mass-weighted nuclear coordinates in atomic units to wavenumbers scale, `dE_au2cm`.  The derivation order can be chosen.
* New `vibrational` module for operations related to vibrational modes and spectroscopy.
* New interface `build_modes` (`vibrational`) can construct the vibrational energies and normal coordinates from the force constants matrix.
* Most conversion functions in `PhysFact` (`physics`) support the reverse operation, with the optional keyword `reverse`.
* New conversion function from atomic unit of mass to unified atomic mass in `PhysFact` (`physics`).
* New procedure `set_orientation` (`vibrational`) to fix the orientation of normal modes, setting the largest component to positive.

### Fixed
* `prt_coord` (`output`) failed if masses were not provided.

### Changed
* `write_err` (`output`) has been improved to be more consistent in the output between the different kind of errors.  The names of the arguments have been changed to be (hopefully) clearer.
* Module `geometry` now uses the new `runstat` from `exception` to manage exceptions.

### Removed
* The conversion method from angstroms to Bohr in `PhysFact` (`physics`) is superseded by `bohr2Ang(reverse=.true.)` and has been removed.


## 0.25.02

### Added
* Operator `.iscloseto.` (`numeric`) tests closeness of 2 values or arrays with sensible thresholds for numeric simulations.
* Function `is_close_to` (`numeric`) tests the closeness, with absolute and relative thresholds that can be set by developers.
* Module `numeric` now has internal "sensible" thresholds for numeric simulations, to be used as defaults in "closeness" tests.
* New methods `copy_to` and `copy_from` added to `MoleculeDB` facilitate copies between databases.  Test on overlaps between source and destination must be done a priori to avoid an error.
* New module `geometry` related to operations on molecular structures.
* New function interface `center_of_mass` (`geometry`) to compute the center of mass for a given geometry, given as array or though a `MoleculeDB` type.
* New function interface `inertia moment` (`geometry`) to compute the inertia moments tensor related to a given geometry, given as an array or through a `MoleculeDB` type.
* New subroutine interface `Eckart_orient` (`geometry`) to transform a given geometry, given as an array or as a `MoleculeDB` type, into an orientation satisfying Eckart's conditions.  The transformation may not be unique and depends on the original orientation.
* New subroutine interface `superpose` (`geometry`) to superpose a given structure, given as an array or as a `MoleculeDB` type, onto an existing one, given as an array.  The subroutine supports different weight models and masks to perform the superposition on a subset of atoms.

### Fixed
* Mass-weighted normal coordinates vectors was not properly normalized after extraction from fchk files.
* Calculation of the number of Cartesian atomic orbitals/basis functions (`basisset: num_cart_AOs`).


## 0.25.01

### Message

* Initial release

### Added

* Implementation of the core library.
* Module parse_cmdline provides a class `ArgObj` and related methods to process commandline arguments.  It works in ways similar to Python's `argparse` with less features.
* New module for electronic-structure calculations providing tools to compute 1-electron integrals and convert orbitals.
* New module to process basis-sets, allowing conversions between pure and Cartesian forms.
* Numerical-precision-related constants and procedures are stored in the module numeric.
* New module `physics` containing physical constants and conversion formulas.
* Atomic information are provided by the module `atominfo`.
* Module `datatypes` provides several derived types to act as container of simulation-related data (molecule specifications, basis sets, orbitals, electronic excitations, vibrations...).
* New parser module to parse Gaussian formatted checkpoint files.
* New drivers modules providing interfaces to common libraries like BLAS and LAPACK
* Implementation of the SOS formalism (P. Bour, Chem. Phys. Lett. 1997, 265, 65-70) to compute excited-excited transition moments.
* New program `mcd_tensor` to compute the MCD G tensor using the SOS formalism.
