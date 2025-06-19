add_rules("mode.debug", "mode.release")

add_requires("openmp")
add_requires("openblas")

rule("ford")
    set_extensions(".md", ".markdown")
    on_build_file(function (target, sourcefile, opt)
        import("core.project.depend")
        import("lib.detect.find_tool")
        import("utils.progress")
        -- make sure build directory exists
        os.mkdir(target:targetdir())
        -- systematically rebuild the API system
        progress.show(opt.progress, "${color.build.object}reading %s",
                      sourcefile)
        local extra_args = {}
        local comp = find_tool("gfortran")
        if comp then
            table.insert(extra_args, "preprocessor='gfortran -E'")
        elseif not find_tool("pcpp") then
            table.insert(extra_args, "preprocess=false")
        end
        if extra_args then
            config = string.format(" --config=\"%s\"",
                                   table.concat(extra_args, "; "))
        else
            config = ''
        end
        local base_args = {"-g", sourcefile, "-o", target:targetdir()}
        cmd = "ford " .. table.concat(base_args, ' ') .. config
        os.exec(cmd)
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
    set_default(false)
    -- deactive from default build --
    set_targetdir("doc")
    -- make the test target support the construction rules of the markdown file
    add_rules("adoc")
    -- adding ASCIIDOC files to build
    add_files("doc/*.adoc")

-- Rules for specific compilations
-- * library system
includes("xmake_libs.lua")
-- * tests
includes("xmake_tests.lua")
-- * programs
includes("xmake_progs.lua")
