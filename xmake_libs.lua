target("corelib")
    -- Core library: contains basic routines/constants
    set_kind("static")
    add_packages("openmp")
    add_files("src/core/numeric.F90")
    add_files("src/core/string.f90")
    add_files("src/core/arrays.f90")
    add_files("src/core/physics.f90")
    add_files("src/core/output.f90")
    add_files("src/core/exception.f90")
    add_files("src/core/parse_cmdline.f90")
    add_files("src/core/parse_cmdline_*.f90")
    add_files("src/data/datatypes.f90")
    add_files("src/data/atominfo.f90")


target("mathlib")
    -- Math and lineary algebra library
    set_kind("static")
    add_packages("openmp")
    add_packages("openblas")
    add_files("src/drivers/blas.f90")
    add_files("src/drivers/lapack.f90")
    add_files("src/core/math.f90")
    add_files("src/core/math_*.f90")


target("datalib")
    -- Library to process and store input data.
    set_kind("static")
    add_packages("openmp")
    add_deps("corelib")
    add_deps("mathlib")
    add_files("src/parsers/parse_fchk.f90")
    add_files("src/data/propinfo.f90")
    add_files("src/core/basisset.f90")
    add_files("src/core/basisset_*.f90")
    add_files("src/core/input.f90")
    add_files("src/core/input_*.f90")
    add_files("src/parsers/input_data_*.f90")


target("molelib")
    -- Molecular structure-related resources
    set_kind("static")
    -- add_packages("openmp")
    add_deps("corelib")
    add_deps("datalib")
    add_deps("mathlib")
    add_files("src/core/geometry.f90")


target("eleclib")
    -- Electronic structure-related resources
    set_kind("static")
    add_deps("corelib")
    add_deps("datalib")
    add_deps("mathlib")
    add_files("src/core/orbital.f90")
    add_files("src/core/electronic.f90")


target("speclib")
    -- Spectroscopy-related resources
    set_default(false)
    set_kind("static")
    add_deps("corelib")
    add_deps("molelib")
    add_deps("mathlib")
    add_packages("openmp")
    add_files("src/spectro/vibrational.f90")
    add_files("src/spectro/vibrational_*.f90")
    add_files("src/spectro/vibronic.f90")
    add_files("src/spectro/vibronic_*.f90")


target("elements")
    -- Full ELEMENTS library
    set_kind("static")
    add_deps("corelib")
    add_deps("mathlib")
    add_deps("datalib")
    add_deps("molelib")
    add_deps("eleclib")
    add_deps("speclib")
    add_packages("openmp")
    add_files("src/core/exc_sos.f90")
