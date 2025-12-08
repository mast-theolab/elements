-- Compilation rules/targets for GauOpen.
-- GauOpen supports different versions internally (32 and 64 bits labels)
-- The configuration below ensures they are properly built and linked for the
-- general interface.



task("check_gauopen_files")
    on_run(function (target, version)
        import("lib.detect.find_file")
        import("lib.detect.find_program")

        local libpath = "src/extlibs/gauopen"
        local files = {"qcmatrix", "qcmatrixio"}
        local ext = '.F'
        -- Check basic file
        print("Checking if 'gauopen' library is installed")
        for _, file in ipairs(files) do
            if not find_file(file .. ext, libpath) then
                raise(string.format(
                    "Missing gauopen file '%s' in '%s'. Stopping",
                    file .. ext, libpath))
            end
        end
        print("'gauopen' library found.")
        print("Checking if 'gauopen' derived files have been generated.")
        local absent = false
        for _, file in ipairs(files) do
            absent = absent or not find_file(file .. version .. ext, libpath)
        end
        if absent then
            print("Derived files absent, attempting to build them.")
            if find_program("python3") then
                python_exe = "python3"
            elseif find_program("python") then
                python_exe = "python"
            else
                raise("Missing Python executable to generate derived files.")
            end
            os.cd(libpath)
            os.runv(python_exe, {"process_qcmat.py"})
            os.cd("$(projectdir)")
        else
            print("Derived files found. Proceeding with compilation.")
        end
    end)

-- Basic rules
rule("gauopen_F")
    -- set_extensions(".F")
    on_load(function (target)
        target:add("fcflags", "gfortran::-ffixed-form", {force = true})
        target:add("fcflags", "gfortran::-std=legacy")
        target:add("fcflags", "-O0")
        target:add("fcflags", "-g")
    end)

rule("gauopen64_F")
    -- set_extensions(".F")
    on_load(function (target)
        -- target:add("fcflags", "-cpp")
        -- target:add("fcflags", "-DUSE_I8")
        -- target:add("defines", "USE_I8")
        target:add("fcflags", "gfortran::-fdefault-integer-8")
    end)

-- Build 32 and 64-bits versions of GauOpen files.
target("gauopen_I4")
    set_default(false)
    set_kind("object")
    add_rules("gauopen_F")
    add_files("src/extlibs/gauopen/qcmatrix*4.F")

    before_build(function (target)
        import("core.project.task")
        task.run("check_gauopen_files", {}, target, "4")
    end)

target("gauopen_I8")
    set_default(false)
    set_kind("object")
    add_rules("gauopen_F")
    add_rules("gauopen64_F")
    add_files("src/extlibs/gauopen/qcmatrix*8.F")

    before_build(function (target)
        import("core.project.task")
        task.run("check_gauopen_files", {}, target, "8")
    end)

target("gauopen")
    -- Binding and library for gauopen
    set_default(false)
    set_kind("static")

    add_deps("corelib")
    add_deps("gauopen_I4")
    add_deps("gauopen_I8")
    add_files("src/drivers/gauopen.F90")
    add_files("src/parsers/gfaf_io.f90")

    before_build(function (target)
        -- import("core.project.task")
        -- task.run("check_gauopen_files", {}, target, "4")
        -- check_gauopen_files (target, "4")
        qcm_ver = io.readfile("src/extlibs/gauopen/gauopen.version")
        target:add("fcflags", "-DGAUOPEN=" .. qcm_ver)
    end)

-- target("gauopen64")
--     -- Binding and library for gauopen
--     set_default(false)
--     set_kind("static")

--     add_deps("qcmatrix8")
--     add_deps("qcmatrixio8")
--     add_files("src/drivers/gauopen_io64.f90")

--     on_load(function (target)
--         import("lib.detect.find_file")
--         if not (find_file("qcmatrix.F", "src/extlibs/gauopen")
--                 and find_file("qcmatrixio.F", "src/extlibs/gauopen")) then
--             raise("Missing gauopen files in 'src/extlibs/gauopen'. Stopping")
--         end
--     end)

-- target("gauopen")
--     set_default(false)
--     -- set_kind("static")
--     add_deps("gauopen32")
--     add_deps("gauopen64")
--     add_files("src/drivers/gauopen.f90")
