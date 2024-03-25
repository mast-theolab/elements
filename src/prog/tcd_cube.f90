program tcd_mesh
    use iso_fortran_env, only: real64, output_unit
    use parse_cmdline, only: CmdArgDB
    use output, only: iu_out, sec_header, len_int, prt_coord, prt_mat, &
      write_err
    use exception, only: BaseException, Error, AllocateError, ArgumentError, &
        FileError, ValueError
    use input, only: build_moldata, build_transdata
    use moldata
    use transdata

    implicit none

    integer :: max_state = -1
    integer, dimension(3) :: &
        cube_npoints = [0, 0, 0]
    real(real64), dimension(3) :: &
        cube_min = [0.0_real64, 0.0_real64, 0.0_real64], &
        cube_delta = [0.0_real64, 0.0_real64, 0.0_real64]
    logical :: DEBUG=.false., TIMEIT=.false.
    character(len=256) :: fchk_file
    character(len=*), parameter :: &
        PROGTITLE = 'Transition Current Density Cube Generator', &
        PROGNAME = 'tcd_cube'

    class(BaseException), allocatable :: err

    interface write_param
        procedure write_param_bool, write_param_int, write_param_real, &
            write_param_char
    end interface write_param

    ! Write title
    call sec_header(-1, PROGTITLE)

    ! Read user options
    call parse_argopts(fchk_file, cube_npoints, cube_min, cube_delta, &
                       max_state)

    ! Read molecular specifications
    call sec_header(1, 'Data on Molecular System')
    if (TIMEIT) call write_time('Molecular data')
    call build_moldata(fchk_file, err)
    if (err%raised()) then
        select type(err)
            class is (FileError)
                call write_err('std', 'Error found with input file', err%msg())
                stop
            class is (Error)
                call write_err('std', &
                    'Error found while parsing basic molecular data', &
                    err%msg() &
                )
                stop
            class default
                call write_err('gen', &
                    'Something went wrong while parsing molecular data.' &
                )
                stop
        end select
    end if
    if (openshell) then
        call write_err('gen', &
            'Sorry, open-shell systems not yet supported.  Working on it.' &
            )
        stop
    end if
    call write_moldata


contains
    function fmt_param(label, sub) result(res)
        !! Formats the parameter label for output
        implicit none

        integer, parameter :: maxchar = 44

        character(len=*), intent(in) :: label
        !! Label of the parameter
        logical, intent(in), optional :: sub
        !! Build parameter label as suboption.
        character(len=maxchar) :: res
        !! Result of the function

        logical :: main

        if (present(sub)) then
            main = .not.sub
        else
            main = .true.
        end if

        if (main) then
            res = ' - ' // trim(label) // ' ' // repeat('-', maxchar)
            res(maxchar:maxchar) = '>'
        else
            res = '   > ' // trim(label)
            res(maxchar:maxchar) = ':'
        end if
    end function fmt_param

    subroutine parse_argopts(infile, cube_npts, xyz_min, cube_step, hi_state)
        !! Parses options from commandline arguments.
        !!
        !! Builds options parser and parses user-options
        implicit none

        character(len=*), intent(out) :: infile
        integer, dimension(3), intent(out) :: cube_npts
        real(real64), dimension(3), intent(out) :: xyz_min, cube_step
        integer, intent(inout) :: hi_state
        
        integer :: ios, iu
        integer, dimension(:), allocatable :: ivals
        real(real64), dimension(:), allocatable :: rvals
        logical :: exists, rd_minpos, rd_npoints, rd_delta
        character(len=80) :: str
        character(len=512) :: fname
        class(CmdArgDB), allocatable :: opts
        data rd_minpos, rd_npoints, rd_delta/3*.false./

        ! Build parser and check arguments
        opts = CmdArgDB(progname='tcd_mesh')
        if (opts%error%raised()) then
            write(*, '(a)') 'Failed to initialize the command-line parser'
            select type (err => opts%error)
                class is (AllocateError)
                    write(*, '(a)') trim(err%msg())
                class is (ArgumentError)
                    write(*, '(a)') trim(err%msg())
                class is (Error)
                    write(*, '(a)') trim(err%msg())
                class default
                    write(*, '(a)') 'Unknown error.'
            end select
            stop
        end if
        call opts%add_arg_char( &
            'string', err, label='filename', &
            help='Gaussian formatted checkpoint file')
        call opts%add_arg_real( &
            'list', err, longname='--minpos', &
            help='Coordinates of the lowest point of the cube',&
            min_nvals=3, max_nvals=3)
        call opts%add_arg_int( &
            'list', err, longname='--npoints', &
            help='Number of points along each coordinate of the cube', &
            min_nvals=3, max_nvals=3)
        call opts%add_arg_int( &
            'list', err, longname='--delta', &
            help='Number of points along each coordinate of the cube', &
            min_nvals=3, max_nvals=3)
        call opts%add_arg_char( &
            'string', err, longname='--cubespec', &
            help='specifies a file with cube specifications')
        call opts%add_arg_int( &
            'scalar', err, shortname='-m', longname='--max-state', &
            help='Highest electronic state to include (upper bound of summation)', &
            min_value=1)
        call opts%add_arg_char( &
            'string', err, label='output', shortname='-o', longname='--output', &
            help='Name of the file to store the output. Existing content will &
            &be overwritten.')
        
        call opts%parse_args(err)
        if (err%raised()) then
            select type (err)
                class is (ValueError)
                    write(*, '(a)') trim(err%msg())
                class is (Error)
                    write(*, '(a)') trim(err%msg())
                class default
                    write(*, '(a)') 'Unknown error while reading options'
            end select
            stop
        end if

        if (opts%is_user_set('output', err)) then
            call opts%get_value('output', fname, err)
            open(newunit=iu_out, file=fname, action='write')
            call sec_header(-1, PROGTITLE)
        end if

        call sec_header(1, 'Parameters for the TCD Meshing')
        write(iu_out, '(1x)')
        if (DEBUG) call write_param("Debugging mode", .true.)
        if (TIMEIT) call write_param("Timer", .true.)

        if (opts%is_user_set('cubespec', err)) then
            call opts%get_value('cubespec', fname, err)
            inquire(file=fname, exist=exists)
            if (.not.exists) then
                print '("ERROR: File ",a," does not exist.")', trim(fname)
                stop
            end if
            open(newunit=iu, file=fname, action='read')
            call write_param('Reading basic cube parameters from', fname)
            ! Reading position of "lowest" point
            read(iu, *, iostat=ios) xyz_min
            if (ios /= 0) then
                print *, 'Error met while reading position of lowest point.'
                stop
            end if
            rd_minpos = .true.
            read(iu, *, iostat=ios) cube_npts
            if (ios /= 0) then
                print *, 'Error met while reading number of points in cube.'
                stop
            end if
            rd_npoints = .true.
            read(iu, *, iostat=ios) cube_step
            if (ios /= 0) then
                print *, 'Error met while reading cube delta parameters.'
                stop
            end if
            rd_delta = .true.
        end if

        if (opts%is_user_set('minpos', err)) then
            call opts%get_value('minpos', rvals, err)
            if (rd_minpos) then
                print *, &
                    '  > NOTE: Lowest-position of cube overridden by user.'
            else
                rd_minpos = .true.
            end if
            xyz_min = rvals
        end if

        if (opts%is_user_set('npoints', err)) then
            call opts%get_value('npoints', ivals, err)
            if (rd_npoints) then
                print *, &
                    '  > NOTE: Number of points in cube overridden by user.'
            else
                rd_npoints = .true.
            end if
            cube_npts = ivals
        end if

        if (opts%is_user_set('delta', err)) then
            call opts%get_value('delta', rvals, err)
            if (rd_delta) then
                print *, &
                    '  > NOTE: Cube delta R parameters overridden by user.'
            else
                rd_delta = .true.
            end if
            cube_step = rvals
        end if

        if (rd_minpos) then
            write(str, '("X=",f0.6,", Y=",f0.6,", Z=",f0.6)') xyz_min
            call write_param('User-defined lowest position of cube', str)
        end if

        if (rd_npoints) then
            write(str, '("Nx=",i0,", Ny=",i0,", Nz=",i0)') cube_npts
            call write_param('User-defined size of the cube', str)
        end if

        if (rd_delta) then
            write(str, '("X=",f0.6,", Y=",f0.6,", Z=",f0.6)') cube_step
            call write_param('User-defined spacing for the cube', str)
        end if

        if (opts%is_user_set('max-state', err)) then
            call opts%get_value('max-state', hi_state, err)
            call write_param('Highest excited state', hi_state)
        else
            hi_state = -1
            call write_param('Highest excited state', 'include all')
        end if

        call opts%get_value('filename', infile, err)
        inquire(file=infile, exist=exists)
        if (.not.exists) then
            write(*, '("Error: File ",a," does not exist.")') trim(infile)
            stop
        end if
        call write_param('Input filename', infile)

    end subroutine parse_argopts

    subroutine write_moldata
        !! Writes molecular data.
        !!
        !! Write some of the molecular data stored in the moldata module
        !! for tcd_mesh.
        implicit none

        character(len=:), allocatable :: tag_bfunD, tag_bfunF, tag_open

        1100 format(/, &
            ' > Number of atoms           : ',i0,/, &
            ' > Charge                    : ',i0,/, &
            ' > Multiplicity              : ',i0,/, &
            ' > Open-/Closed-shell        : ',a,//, &
            ' > Number of basis functions : ',i0,/, &
            ' > Type of D functions       : ',a,/, &
            ' > Type of F functions       : ',a,/, &
            ' > Number of atomic orbitals : ',i0)

        if (openshell) then
            tag_open = 'open'
        else
            tag_open = 'closed'
        end if

        if (pureD) then
            tag_bfunD = 'pure'
        else
            tag_bfunD = 'Cartesian'
        end if
        if (pureF) then
            tag_bfunF = 'pure'
        else
            tag_bfunF = 'Cartesian'
        end if

        write(iu_out, 1100) n_at, charge, multip, tag_open, n_basis, &
            tag_bfunD, tag_bfunF, n_ao

        call sec_header(2, 'Atomic coordinates')
        call prt_coord(n_at, at_lab, at_crd, at_mas)
    end subroutine write_moldata

    subroutine write_param_bool(param, status, subparam)
        !! Writes status for logical-type parameter.
        !!
        !! Formats and writes one parameter of simulation with the
        !! specified status.
        implicit none

        character(len=*), intent(in) :: param
        !! Parameter name (with info like unit).
        logical, intent(in) :: status
        !! Status
        logical, intent(in), optional :: subparam
        !! Parameter is a sub-parameter type

        1000 format(a,1x,"enabled")
        1010 format(a,1x,"disabled")

        if (status) then
            write(iu_out, 1000) fmt_param(param, subparam)
        else
            write(iu_out, 1010) fmt_param(param, subparam)
        end if
    end subroutine write_param_bool

    subroutine write_param_int(param, value, subparam)
        !! Writes value of integer-type parameter.
        !!
        !! Formats and writes one parameter of simulation with the
        !! indicated value.
        implicit none

        character(len=*), intent(in) :: param
        !! Parameter name (with info like unit).
        integer, intent(in) :: value
        !! Value associated to parameter
        logical, intent(in), optional :: subparam
        !! Parameter is a sub-parameter type

        1000 format(a,1x,i0)

        write(iu_out, 1000) fmt_param(param, subparam), value
    end subroutine write_param_int

    subroutine write_param_real(param, value, subparam)
        !! Writes value of real-type parameter.
        !!
        !! Formats and writes one parameter of simulation with the
        !! indicated value.
        implicit none

        character(len=*), intent(in) :: param
        !! Parameter name (with info like unit).
        real(real64), intent(in) :: value
        !! Value associated to parameter
        logical, intent(in), optional :: subparam
        !! Parameter is a sub-parameter type

        1000 format(a,1x,f0.4)
        1010 format(a,es13.6)

        if (abs(value) < 1.0e-2_real64 .or. abs(value) > 1.0e4_real64) then
            write(iu_out, 1010) fmt_param(param, subparam), value
        else
            write(iu_out, 1000) fmt_param(param, subparam), value
        end if
    end subroutine write_param_real

    subroutine write_param_char(param, value, subparam)
        !! Writes value of character-type parameter.
        !!
        !! Formats and writes one parameter of simulation with the
        !! indicated value.
        implicit none

        character(len=*), intent(in) :: param
        !! Parameter name (with info like unit).
        character(len=*), intent(in) :: value
        !! Value associated to parameter
        logical, intent(in), optional :: subparam
        !! Parameter is a sub-parameter type

        1000 format(a,1x,a)

        write(iu_out, 1000) fmt_param(param, subparam), trim(value)
    end subroutine write_param_char

    subroutine write_time(label)
        !! Writes time information.
        !!
        !! Writes time data at call time, displaying the label as well.
        implicit none

        character(len=*), intent(in) :: label
        !! Label to display in printing information

        integer, dimension(8) :: ia_dtime

        1000 format('Entering: ',a,' - Date: ',i4,'/',i2.2,'/',i2.2,' at ', &
            i2.2,':',i2.2,':',i2.2)
        call date_and_time(values=ia_dtime)
        write(iu_out, 1000) trim(label), ia_dtime(1:3), ia_dtime(5:7)
    end subroutine write_time

end program tcd_mesh
