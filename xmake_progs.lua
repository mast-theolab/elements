target("autoclave")
    set_default(false)
    add_packages("openmp")
    set_rundir("$(projectdir)/tests")
    add_deps("elements")
    add_files("src/progs/autoclave.f90")

target("build_boltz_pop")
    set_default(false)
    add_packages("openmp")
    set_rundir("$(projectdir)/tests")
    add_deps("elements")
    add_files("src/progs/build_boltz_pop.f90")
    add_tests("h2co",
              {runargs = {"H2CO_S0_frq.fchk",
                          "-t", "1000",
                          "-p", "0.001",
                          "-o", "bzpop_h2co_T1000_P001.txt"}})
    add_tests("meox",
              {runargs = {"meox.S0.vac.B3PW91.junTZ.frq-ROA.fchk",
                          "-t", "400",
                          "-p", "0.005",
                          "-o", "bzpop_meox_T400_P005.txt"}})

target("calcites")
    set_default(false)
    add_packages("openmp")
    set_rundir("$(projectdir)/tests")
    add_deps("elements")
    add_files("src/progs/calcites.f90")
    -- add_tests("h2co",
    --           {runargs = {"H2CO_S0_frq.fchk",
    --                       "-t", "1000",
    --                       "-p", "0.001",
    --                       "-o", "bzpop_h2co_T1000_P001.txt"}})
    -- add_tests("meox",
    --           {runargs = {"meox.S0.vac.B3PW91.junTZ.frq-ROA.fchk",
    --                       "-t", "400",
    --                       "-p", "0.005",
    --                       "-o", "bzpop_meox_T400_P005.txt"}})

target("density_cube")
    set_default(false)
    set_kind("binary")
    add_packages("openmp")
    add_deps("elements")
    add_files("src/progs/density_cube.f90")
-- TO DO: ADD TESTS


target("gen_py_atomDB")
    set_default(false)
    set_rundir("$(projectdir)/tests")
    add_deps("corelib")
    add_deps("datalib")
    add_files("src/progs/gen_py_atomdb.f90")
    add_tests("to_ang", {runargs = {"atom_ang.py"}})
    add_tests("to_au", {runargs = {"atom_bohr.py"}})

target("mcd_tensor")
    -- Program to compute the MCD tensor.
    set_default(false)
    set_kind("binary")
    add_packages("openmp")
    add_deps("elements")
    add_files("src/extras/mcd/gmcd_output.f90")
    add_files("src/extras/mcd/gmcd_legacy.f90")
    add_files("src/progs/mcd_tensor.f90")
    set_rundir("$(projectdir)/tests")
    -- Check that input data are consistent with formchk from G16/GDV
    add_tests("HOF:default with G16 fchk",
              {runargs = {"HOF.vac.B3LYP.321G.TD.G16.fchk", "--no-giao",
                          "-o", "mcd_hof_default_g16.txt",
                          "--no-timestamp"}})
    add_tests("HOF:default with GDV fchk",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.fchk", "--no-giao",
                          "-o", "mcd_hof_default_gdv.txt",
                          "--no-timestamp"}})
    -- Check that unrestricted and closed-shell 50-50 singlet-triplet are
    --   consistent, GIAO correction deactivated (transition S0 -> T1).
    add_tests("HOF:openshell, final=T1, no GIAO",
              {runargs = {"HOF.vac.UB3LYP.321G.TD.GDV.fchk", "--no-giao",
                          "-o", "mcd_hof_openshell_no-GIAO.txt",
                          "--no-timestamp"}})
    add_tests("HOF:closed 50-50, final=T1, no GIAO",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.5050.fchk", "--no-giao",
                          "-o", "mcd_hof_closed_50-50_no-GIAO.txt",
                          "--no-timestamp"}})
    -- Check that unrestricted and closed-shell 50-50 are consistent (S0->S1)
    add_tests("HOF:openshell, final=S1, no GIAO",
              {runargs = {"HOF.vac.UB3LYP.321G.TD.GDV.fchk", "--no-giao",
                          "-o", "mcd_hof_openshell_no-GIAO_S1.txt",
                          "--final=2", "--debug=ijaa",
                          "--no-timestamp"}})
    add_tests("HOF:closed 50-50, final=S1, no GIAO",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.5050.fchk", "--no-giao",
                          "-o", "mcd_hof_closed_50-50_no-GIAO_S1.txt",
                          "--final=2", "--debug=ijaa",
                          "--no-timestamp"}})
    -- Same as before but deactivates the component <ia|ja> in the prefactor
    --   as done in GUVCDE for closed-shells but not open-shells calculations.
    add_tests("HOF:openshell, final=S1, no GIAO, no <ia|ja>",
              {runargs = {"HOF.vac.UB3LYP.321G.TD.GDV.fchk", "--no-giao",
                          "-o", "mcd_hof_openshell_no-GIAO_no-ijaa_S1.txt",
                          "--final=2", "--debug=noijaa",
                          "--no-timestamp"}})
    add_tests("HOF:closed 50-50, final=S1, no GIAO, no <ia|ja>",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.5050.fchk", "--no-giao",
                          "-o", "mcd_hof_closed_50-50_no-GIAO_no-ijaa_S1.txt",
                          "--final=2", "--debug=noijaa",
                          "--no-timestamp"}})
    -- Check that unrestricted and closed-shell 50-50 are consistent with GIAO
    add_tests("HOF:openshell, final=S1, GIAO",
              {runargs = {"HOF.vac.UB3LYP.321G.TD.GDV.fchk",
                          "-o", "mcd_hof_openshell_S1.txt",
                          "--final=2",
                          "--no-timestamp"}})
    add_tests("HOF:closed 50-50, final=S1, GIAO",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.5050.fchk",
                          "-o", "mcd_hof_closed_50-50_S1.txt",
                          "--final=2",
                          "--no-timestamp"}})
    -- Test pure Slater model for excitations
    add_tests("HOF:openshell, final=S1, Slater",
              {runargs = {"HOF.vac.UB3LYP.321G.TD.GDV.fchk",
                          "--exc-model=slater",
                          "-o", "mcd_hof_openshell_S1_Slater.txt",
                          "--final=2",
                          "--no-timestamp"}})
    add_tests("HOF:closed 50-50, final=S1, Slater",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.5050.fchk",
                          "--exc-model=slater",
                          "-o", "mcd_hof_closed_50-50_S1_Slater.txt",
                          "--final=2",
                          "--no-timestamp"}})
    -- Test model relying purely on TD amplitudes for excitations
    add_tests("HOF:openshell, final=S1, TD amplitudes",
              {runargs = {"HOF.vac.UB3LYP.321G.TD.GDV.fchk",
                          "--exc-model=amplitudes",
                          "-o", "mcd_hof_openshell_S1_TDampl.txt",
                          "--final=2",
                          "--no-timestamp"}})
    add_tests("HOF:closed 50-50, final=S1, TD amplitudes",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.5050.fchk",
                          "--exc-model=amplitudes",
                          "-o", "mcd_hof_closed_50-50_S1_TDampl.txt",
                          "--final=2",
                          "--no-timestamp"}})
    -- Test using hybrid model for excitations
    add_tests("HOF:openshell, final=S1, TD hybrid model",
              {runargs = {"HOF.vac.UB3LYP.321G.TD.GDV.fchk",
                          "--exc-model=hybrid",
                          "-o", "mcd_hof_openshell_S1_hybrid-exc.txt",
                          "--final=2",
                          "--no-timestamp"}})
    add_tests("HOF:closed 50-50, final=S1, TD hybrid model",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.5050.fchk",
                          "--exc-model=hybrid",
                          "-o", "mcd_hof_closed_50-50_S1_hybrid-exc.txt",
                          "--final=2",
                          "--no-timestamp"}})

target("tcd_cube")
    set_default(false)
    set_kind("binary")
    add_links("blas")
    add_packages("openmp")
    add_deps("molelib")
    add_deps("mathlib")
    add_deps("cubelib")
    add_files("src/progs/tcd_cube.f90")
    set_rundir("$(projectdir)/tests")
    add_tests("meaz_ECD",
              {runargs = {"methylaziridine+.vac.B3LYP.631Gd.TD.fchk",
                          "-s", "2", "-e", "2", "-d", "scarce",
                          "-o", "tcd_cube_meax_scarce.txt"
                          }})

target("vertex")
    set_default(false)
    add_packages("openmp")
    set_rundir("$(projectdir)/tests")
    add_deps("elements")
    add_files("src/progs/vertex.f90")

target("vtcd_cube")
    set_default(false)
    set_kind("binary")
    add_packages("openmp")
    add_links("blas")
    add_deps("corelib")
    add_deps("molelib")
    add_deps("mathlib")
    add_deps("cubelib")
    add_files("src/progs/vtcd_cube.f90")
    set_rundir("$(projectdir)/tests")

    before_test(function (target, opt)
        if opt.name == "vtcd_cube/meox_v05_moldata" then
            import("lib.detect.find_file")
            local moldatfile = "vtcd_meox_mol.dat"
            local cubdatfile = "vtcd_meox_cube.dat"
            local workdir = path.join(os.projectdir(), "tests")
            os.cd(workdir)
            if find_file(moldatfile, ".") then
                os.rm(moldatfile)
            end
            if find_file(cubdatfile, ".") then
                os.rm(cubdatfile)
            end
            os.runv(target:targetfile(),
                    {"vtcd_meox_TD.fchk",
                     "vtcd_meox_nac_###-mod.fchk",
                     "vtcd_meox_freq.fchk",
                     "-q", "6",
                     "-s", "1",
                     "-e", "3",
                     "-m", moldatfile,
                     "-b", cubdatfile,
                     "-d", "scarce"})
            os.rm("vtcd_meox_TD_v06_s1-3.cube")
            os.cd(os.projectdir())
            return
        end
    end)

    add_tests("meox_v05_moldata",
              {runargs = {
                  "-u", "vtcd_meox_mol.dat",
                  "-c", "vtcd_meox_cube.dat",
                  "-q", "5",
                  "-s", "1",
                  "-e", "3",
                  "-o", "vtcd_cube_meox_v05_data.txt"
            }})

    add_tests("meox_v04_scarce",
              {runargs = {
                    "vtcd_meox_TD.fchk",
                    "vtcd_meox_nac_###-mod.fchk",
                    "vtcd_meox_freq.fchk",
                    "-q", "4",
                    "-s", "1",
                    "-e", "3",
                    "-d", "scarce",
                    "-o", "vtcd_cube_meox_v04_scarce.txt"
            }})

    add_tests("meox_v18_vlow",
              {runargs = {
                  "vtcd_meox_TD.fchk",
                  "vtcd_meox_nac_###-mod.fchk",
                  "vtcd_meox_freq.fchk",
                  "-q", "18",
                  "-s", "1",
                  "-e", "3",
                  "-d", "scarce",
                  "-o", "vtcd_cube_meox_v18_vlow.txt"
            }})
