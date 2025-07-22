# Changelog

## 0.25.07

### Added
* ***arrays***: new interface `indexes_sym2lt` to convert the coordinates (2-4) in a nD symmetric tensor to the linear storage index of the lower-triangular block.
* ***basisset***: new function to compute the derivative of the angular momentum for a shell made of Cartesian atomic orbitals at a given point: `get_cart_L_der_sh_at`.
* ***basisset***: new function to compute the normalization factor for the angular momentum component for a shell made of Cartesian atomic orbitals at a given point: `get_cart_L_norms_sh`.
* ***basisset***: new function to compute the angular momentum for a shell made of Cartesian atomic orbitals at a given point: `get_cart_L_sh_at`.
* ***basisset***: new function to compute the derivative of the radial component for a shell made of Cartesian atomic orbitals at a given point: `get_cart_r_der_sh_at`.
* ***basisset***: new function to compute the radial component for a shell made of Cartesian atomic orbitals at a given point: `get_cart_r_sh_at`.
* ***basisset***: new function to compute the overlap between primitives in a given shell made of Cartesian atomic orbitals: `get_primitives_overlap_sh`.
* ***basisset***: new function to return the sequence of powers for a given angular momentum: `list_L_powers`.
* ***basisset***: function `build_bset_DB` now returns the highest angular moment found while building the database.
* ***basisset***: `convert_pure2cart` can now produce a full basis set database in output.
* ***calcites***: new program to compute the thermal contributions to vibrational spectra within the harmonic approximation for a given type of spectroscopy.
* ***datatypes***: new component `L_max` to derived type `BasisSetDB` storing the highest component in the whole basis set.
* ***datatypes***: added methods `is_pure` and `is_cart` to `BasisSetDB` to check the type of basis set.  Note that mixed basis sets could lead to both functions yielding `.false.`.
* ***math***: new function `double_factorial` to compute the double factorial "n!!".
* ***orbital***: new module dedicated to operation on atomic/molecular orbitals.
* ***orbital***: new subroutine `eval_AOs_chi_at` to evaluate atomic orbitals at a chosen position.
* ***orbital***: new subroutine `eval_AOs_nabla_chi_at` to evaluate atomic orbitals and their first derivatives at a chosen position.
* ***vertex***: new program to compute and print the extrapolated geometry for vibronic vertical models.
* ***vibronic***: new `extrapolate_geom` to compute the extrapolated minimum geometry for vertical vibronic models (VG, VH).

### Fixed
* ***input***: the fchk parser now handles better missing quantities related to electronic excited states.
* ***parse_cmdline***: `ArgObj%is_set` was not properly set if no default value was provided and the user had not provided any value, resulting in a risk of wrong interpretation if the value was randomly set to 1 or 2.
* ***vibronic***: Extrapolated shift vector (VH, VG) now supports imaginary frequencies given as negative numbers.

### Changed
* ***basisset***: removed unused argument `n_ao` from `num_cart_AOs`.
* ***build***: OpenMP is now enabled by default in `exc_sos`.
* ***electronic***: `convert_AO2MO` is transferred to module `orbital` and thus removed from module `electronic`.
* ***math***: Definition of linear-algebra routines has been moved to submodules based on the type and kind (real32, real64, complex32, complex64).
* ***math***: Procedures now ordered by alphabetical order to facilitate search of procedures.


## 0.25.06

### Added
* ***autoclave***: New program to compute VPT2 vibrational energies from a (Gaussian) logfile.  The program handles file parsing itself.
* ***datatypes***: New method in `PropertyDB` to easily clear the stored data and free memory.
* ***numeric***: Precision-specific thresholds and value of pi are now publicly available for modules that do not use the _math_ library otherwise.
* ***physics***: New derived type `SpectroConv` and new instance `spec_conv` to store common conversion operations related to computational spectroscopy.
* ***physics***: New type-bound procedure `mwq2q` to facilitate conversions of quantity expressed with respect to normal coordinates.
* ***vibrational***: New function interface `full_boltz_pop` to compute the analytic total Boltzmann population for a given set of frequencies at any temperature.
* ***vibrational***: New module `vibrational_pt2` dedicated to vibrational perturbation theory at second order.
* ***vibrational***: New function interface `calc_en_vib` in `vibrational_pt2` to compute anharmonic vibrational energies without variational correction.

### Fixed
* ***build_boltz_pop***: The program now computes the analytical total Boltzmann population and compares it to the computed value, giving some measure of the convergence.
* ***input***: Fixed parsing of Gaussian version in fchk files generated from checkpoint files obtained with the `c86dv` (form: CDVRev-X.XX).
* ***vibrational***: Fixed missing initialization of `nqi` in `boltz_pop_max_quanta` when checking duplicate modes, which could lead to the false assertion that the list of mode indexes contained duplicates.
* ***vibrational***: Fixed logic associated to optional argument `is_weighted`, which was internally interpreted with the opposite meaning.

### Changed
* ***input***: To speed up the file type identification process, by default the file is only analyzed if the extension is not conclusive.  The older behavior can be restored with `soft_check` set to `.false.`.
* ***build***: Build file `xmake.lua` has been split into 4 files, with three sub-files storing recipes for the library (`xmake_libs.lua`), the internal tests (`xmake_tests.lua`) and the stand-alone programs (`xmake_progs.lua`).


## 0.25.04

### Note
* Change in format: now the program or module name on which the change applies is listed first: `name: change`.  Special names are:
    * `build`: for the general build tool chain
    * `lib`: the ELEMENTS library in general
    * `misc`: unclassified change
* Modifications are sorted by category, not anymore purely historically (last change last).

### Added
* ***lib***: Renaming of some private module procedures to use `_db` for the version using databases (derived types) instead of pure arrays.
* ***numeric***: added keywords for work-precision powers of 10, as `f10xy` with `x`= p(lus) or m(inus), and `y` an integer. Ex.: `f10p2` = `1.0e2_realwp`.
* ***output***: it is possible to set the number of leading blank characters before a header with optional keyword `lead_spaces` in `sec_header`.
* ***string***: `num_chars_int` calculates the number of characters needed to store an integer.
* ***vibrational***: new subroutine `boltz_pop_max_quanta` to compute the populated vibrational states based on a given temperature and minimum population with respect to the ground state.  The procedure can run in several modes, only reporting the maximum number of quanta for each mode considering at least one quanta per mode, or explicitly listing all possible states and their respective populations as well.
* ***build_boltz_pop***: new program to build the list of populated states above a given threshold based on the input temperature.
* ***vibronic***: new module for operations related to vibrationally-resolved electronic spectroscopy.
* ***vibronic***: new function `Duschinsky_matrix` to compute the Duschinsky matrix.
* ***vibronic***: new function `Duschinsky_shift` to compute the shift vector in the Duschinsky transformation

### Changed
* ***output***: error messages are now preceded by a blank line, so they are more visible in the output.
* ***output***: headers do not have a leading blank character anymore by default.
* ***string***: the default for module procedures is private, as will be enforced for every module.


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
* New procedure `prt_vec` (`output`) to print vectors.
* New component `red_freq` in `VibrationsDB` to store reduced harmonic frequencies in atomic units.
* Added conversion function `to_int` and `to_real` (`numeric`) to convert an arbitrary object to the internal default kind (as defined in `numeric`).

### Fixed
* `prt_coord` (`output`) failed if masses were not provided.
* `xgemm` properly loaded in `basisset_purecart`.
* Fixed race condition when using `inertia_tensor` (`geometry`) in parallel, giving garbage tensors.
* Copy methods of `MoleculeDB` objects (`datatypes`) now properly check that an array is allocated before trying to copy them in another object.
* Datafiles were opened twice at the same time because of a missing closure.
* The number of Cartesian AOs was incorrectly computed in `num_cart_AOs_bsetBF` (`basisset`).
* Fixed errors in the definition of the kind in `to_int64` and `to_real64` (`numeric`).
* Fixed value of 0! in `factorial` (`math`).

### Changed
* `write_err` (`output`) has been improved to be more consistent in the output between the different kind of errors.  The names of the arguments have been changed to be (hopefully) clearer.
* Module `geometry` now uses the new `runstat` from `exception` to manage exceptions.
* Improved error messages in `mcd_tensor`.

### Improved
* Better efficiency of `int_xn_e2ax2` (`math`).
* Better performance in the building of basis set coefficients, especially beyond f (`basisset`).

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
