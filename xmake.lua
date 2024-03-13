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
    set_targetdir("doc")
    -- make the test target support the construction rules of the markdown file
    add_rules("adoc")
    -- adding ASCIIDOC files to build
    add_files("doc/*.adoc")


target("corelib")
    set_kind("static")
    add_packages("openmp")
    add_files("src/lib/numeric.f90")
    add_files("src/lib/atomic.f90")
    add_files("src/lib/exception.f90")
    add_files("src/lib/math.f90")
    add_files("src/lib/physic.f90")
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
    add_files("src/lib/parse_cmdline_parser.f90")
    add_files("src/lib/parse_cmdline_addarg.f90")
    add_files("src/lib/parse_cmdline_setval.f90")
    add_files("src/lib/parse_cmdline_getval.f90")


target("datalib")
    set_kind("static")
    add_packages("openmp")
    add_deps("corelib")
    add_files("src/lib/moldata.f90")
    add_files("src/lib/transdata.f90")
    add_files("src/lib/input.f90")


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

