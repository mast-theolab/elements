add_rules("mode.debug", "mode.release")

add_repositories("local-repo xmake")
add_requires("openmp")


rule("ford")
    set_extensions(".md", ".markdown")
    on_build_file(function (target, sourcefile, opt)
        import("core.project.depend")
        import("utils.progress")
        -- make sure build directory exists
        os.mkdir(target:targetdir())
        -- systematically rebuild the API system
        progress.show(opt.progress, "${color.build.object}reading %s", sourcefile)
        os.execv("ford", {"-g", sourcefile, "-o", target:targetdir()})
    end)

rule("adoc")
    set_extensions(".adoc")
    on_build_file(function (target, sourcefile, opt)
        import("core.project.depend")
        import("utils.progress")
        local docdir = target:targetdir()
        local libdir = path.join(".", docdir, 'lib')
        -- replace .adoc with .html
        local targetfile = path.join(docdir, path.basename(sourcefile) .. ".html")
        -- only rebuild the file if its changed since last run
        depend.on_changed(function ()
            -- call ford to generate the documentation
            progress.show(opt.progress, "${color.build.object}reading %s", sourcefile)
            os.execv("asciidoctor",
                     {"-r", path.join(libdir, "custom-admonition-block.rb"),
                      "-a", "stylesheet=mastdoc.css",
                      "-o", targetfile,
                      sourcefile})
            -- os.execv("asciidoctor", {"-o", targetfile, sourcefile})
        end, {files = sourcefile})
    end)


target("apidoc")
    on_clean(function (target)
        os.rm(target:targetdir())
    end)

    set_kind("object")
    -- deactive from default build --
    set_default(false)
    set_targetdir("doc/api")
    -- make the test target support the construction rules of the markdown file
    add_rules("ford")
    -- adding a markdown file to build
    add_files("elements.md")


target("userdoc")
    on_clean(function (target)
        for _, file in pairs(target:sourcefiles()) do
            os.rm(path.join(target:targetdir(), path.basename(file) .. '.html'))
        end
        -- print(path.basename(target:sourcefile()) .. '.html')
        -- print(path.basename(target:sourcebatches()) .. '.html')
    end)

    set_kind("object")
    -- deactive from default build --
    set_targetdir("doc")
    -- make the test target support the construction rules of the markdown file
    add_rules("adoc")
    -- adding ASCIIDOC files to build
    add_files("doc/*.adoc")


target("corelib")
    set_kind("static")
    add_packages("openmp")
    add_files("src/lib/numeric.f90")
    add_files("src/lib/arrays.f90")
    add_files("src/lib/atomic.f90")
    add_files("src/lib/exception.f90")
    add_files("src/lib/math.f90")
    add_files("src/lib/physics.f90")
    add_files("src/lib/string.f90")
    add_files("src/lib/output.f90")
    add_files("src/lib/basisset.f90")
    add_files("src/lib/electronic.f90")
    add_files("src/lib/parsefchk.f90")


target("parselib")
    set_kind("static")
    add_packages("openmp")
    add_deps("corelib")
    add_files("src/lib/parse_cmdline.f90")
    add_files("src/lib/parse_cmdline_*.f90")


target("datalib")
    set_kind("static")
    add_packages("openmp")
    add_deps("corelib")
    add_files("src/lib/moldata.f90")
    add_files("src/lib/transdata.f90")
    add_files("src/lib/input.f90")
    add_files("src/lib/input_*.f90")


target("mcd_tensor")
    set_kind("binary")
    add_packages("openmp")
    add_deps("corelib")
    add_deps("parselib")
    add_deps("datalib")
    add_files("src/lib/exc_sos.f90")
    add_files("src/prog/gmcd_output.f90")
    add_files("src/prog/gmcd_legacy.f90")
    add_files("src/prog/mcd_tensor.f90")
    set_rundir("$(projectdir)/tests")
    -- Check that input data are consistent with formchk from G16/GDV
    add_tests("HOF:default with G16 fchk",
              {runargs = {"HOF.vac.B3LYP.321G.TD.G16.fchk", "--no-giao",
                          "-o", "mcd_hof_default_g16.txt"}})
    add_tests("HOF:default with GDV fchk",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.fchk", "--no-giao",
                          "-o", "mcd_hof_default_gdv.txt"}})
    -- Check that unrestricted and closed-shell 50-50 singlet-triplet are
    --   consistent, GIAO correction deactivated (transition S0 -> T1).
    add_tests("HOF:openshell, final=T1, no GIAO",
              {runargs = {"HOF.vac.UB3LYP.321G.TD.GDV.fchk", "--no-giao",
                          "-o", "mcd_hof_openshell_no-GIAO.txt"}})
    add_tests("HOF:closed 50-50, final=T1, no GIAO",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.5050.fchk", "--no-giao",
                          "-o", "mcd_hof_closed_50-50_no-GIAO.txt"}})
    -- Check that unrestricted and closed-shell 50-50 are consistent (S0->S1)
    add_tests("HOF:openshell, final=S1, no GIAO",
              {runargs = {"HOF.vac.UB3LYP.321G.TD.GDV.fchk", "--no-giao",
                          "-o", "mcd_hof_openshell_no-GIAO_S1.txt",
                          "--final=2", "--debug=ijaa"}})
    add_tests("HOF:closed 50-50, final=S1, no GIAO",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.5050.fchk", "--no-giao",
                          "-o", "mcd_hof_closed_50-50_no-GIAO_S1.txt",
                          "--final=2", "--debug=ijaa"}})
    -- Same as before but deactivates the component <ia|ja> in the prefactor
    --   as done in GUVCDE for closed-shells but not open-shells calculations.
    add_tests("HOF:openshell, final=S1, no GIAO, no <ia|ja>",
              {runargs = {"HOF.vac.UB3LYP.321G.TD.GDV.fchk", "--no-giao",
                          "-o", "mcd_hof_openshell_no-GIAO_no-ijaa_S1.txt",
                          "--final=2", "--debug=noijaa"}})
    add_tests("HOF:closed 50-50, final=S1, no GIAO, no <ia|ja>",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.5050.fchk", "--no-giao",
                          "-o", "mcd_hof_closed_50-50_no-GIAO_no-ijaa_S1.txt",
                          "--final=2", "--debug=noijaa"}})
    -- Check that unrestricted and closed-shell 50-50 are consistent with GIAO
    add_tests("HOF:openshell, final=S1, GIAO",
              {runargs = {"HOF.vac.UB3LYP.321G.TD.GDV.fchk",
                          "-o", "mcd_hof_openshell_S1.txt",
                          "--final=2"}})
    add_tests("HOF:closed 50-50, final=S1, GIAO",
              {runargs = {"HOF.vac.B3LYP.321G.TD.GDV.5050.fchk",
                          "-o", "mcd_hof_closed_50-50_S1.txt",
                          "--final=2"}})


target("tcd_cube")
    set_default(false)
    set_kind("binary")
    add_packages("openmp")
    add_deps("corelib")
    add_deps("parselib")
    add_deps("datalib")
    add_files("src/prog/tcd_cube.f90")

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