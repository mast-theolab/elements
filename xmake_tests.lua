
target("test_symm_array")
    set_default(false)
    set_rundir("$(projectdir)/tests")
    add_deps("corelib")
    add_files("src/tests/symm_array.f90")
    add_tests("default")

target("test_read_file")
    set_default(false)
    set_rundir("$(projectdir)/tests")
    add_deps("datalib")
    add_files("src/tests/read_file.f90")
    add_tests("default")

target("test_read_vib")
    set_default(false)
    set_rundir("$(projectdir)/tests")
    add_deps("datalib")
    add_deps("corelib")
    add_files("src/tests/read_vibdat.f90")
    add_tests("default",
              {runargs = {"-f", "H2CO_S0_frq.fchk", "-o", "test_vib_def.txt"}})
    add_tests("twofiles",
              {runargs = {"-f", "H2CO_S0_frq.fchk", "-o", "H2CO_S2_frq.fchk",
               "-o", "test_vib_two.txt"}})

target("test_blas_ops")
    set_default(false)
    set_rundir("$(projectdir)/tests")
    add_deps("corelib")
    add_deps("mathlib")
    add_files("src/tests/blas_ops.f90")
    add_tests("default")

target("test_get_data")
    set_default(false)
    add_packages("openmp")
    set_rundir("$(projectdir)/tests")
    add_deps("datalib")
    add_deps("corelib")
    add_files("src/tests/test_getdata.f90")
    add_tests("h2co",
              {runargs = {"-f", "H2CO_S2_frq.fchk",
                          "-o", "test_getdata_H2CO.txt"}})
    add_tests("biphenyl",
              {runargs = {"-f", "biphenyl_S1_frq.fchk",
                          "-o", "test_getdata_biphenyl.txt"}})

target("test_timestamp")
    set_default(false)
    set_rundir("$(projectdir)/tests")
    add_deps("corelib")
    add_files("src/tests/test_timestamp.f90")
    add_tests("default")

target("test_geom_ops")
    set_default(false)
    -- add_packages("openmp")
    set_rundir("$(projectdir)/tests")
    add_deps("molelib")
    add_files("src/tests/geom_ops.f90")
    add_tests("default")

target("test_vibronic")
    set_default(false)
    set_rundir("$(projectdir)/tests")
    add_deps("speclib")
    add_files("src/tests/test_vibronic.f90")
    add_tests("default")
