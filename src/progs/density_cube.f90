program density_cube
    use basisset, only: convert_pure2cart, fix_norm_AOs
    use cubegen, only: CubeDensDB, CubeGridDB
    use datatypes, only: BasisSetDB, MoleculeDB, OrbitalsDB
    use input, only: DataFile
    use numeric, only: realwp, f10m2, f10p4
    use output, only: iu_out, sec_header, prt_coord, prt_mat
    use parse_cmdline, only: CmdLineArgsDB
    use run_env, only: ErrorHandle, run
    ! use exception, only: BaseException, Error, AllocateError, ArgumentError, &
    !                            FileError, ValueError, RaiseFileError, InitError

    ! use datatypes, only: MoleculeDB, BasisSetDB, OrbitalsDB, ExcitationDB, PropertyDB

    ! use arrays,          only: symm_tri_arr_r64

    ! use cubegen,         only: moldb, orbdb, bsetdb, generate_cube, initialize_grid_parameters, generate_partioning_cube

    implicit none

    type Params
        character(len=:), allocatable :: file_fchk
        character(len=:), allocatable :: file_grid
        character(len=:), allocatable :: grid_sparsity
        character(len=:), allocatable :: print_level
        logical :: partition = .false.
    end type Params

    logical, parameter :: TIMEIT = .False.
    logical, parameter :: DEBUG = .False.

    integer :: ios, iu
    real(realwp), dimension(:,:), allocatable :: density
    real(realwp):: density_sum
    character(len=:), allocatable :: cube_fname
    character(len=*), parameter :: &
        PROGTITLE = 'Electonic Density Cube Generator'
    type(ErrorHandle), allocatable :: err
    type(BasisSetDB), target :: bsetDB_orig, bsetDB_cart
    type(BasisSetDB), pointer :: bsetDB
    type(CubeDensDB) :: cubeDB
    type(DataFile) :: dfile
    type(MoleculeDB) :: molDB
    type(OrbitalsDB) :: orbDB
    type(Params) :: opts

    interface write_param
        procedure write_param_bool, write_param_int, write_param_real, &
            write_param_char
    end interface write_param

    call sec_header(-1, PROGTITLE)

    call parse_opts(opts)

    dfile = DataFile(opts%file_fchk)
    call run%check(dfile%error, 'Error found while initializing data file')

    call sec_header(1, 'Data on Molecular System')

    molDB = dfile%get_mol_data(get_dens=.true.)
    call run%check(dfile%error, &
        'Error found while parsing molecular data in file')
    if (.not.molDB%loaded) then
        call run%error%raise_error('data', 'missing', &
            'Failure to extract molecular data')
    end if

    bsetDB_orig = dfile%get_bset_data()
    call run%check(dfile%error, 'Parsing basis set data in file failed')
    if (.not.bsetDB_orig%loaded) then
        call run%error%raise_error('data', 'missing', &
            'Failure to extract basis-set data')
    end if

    orbDB = dfile%get_orb_data()
    call run%check(dfile%error, &
        'Parsing molecular orbitals data in file failed')
    if (.not.orbDB%loaded) then
        call run%error%raise_error('data', 'missing', &
            'Failure to extract orbitals-related data')
    end if

    call write_moldata(molDB, bsetDB_orig, orbDB)

    if (orbDB%openshell) &
        call run%error%raise_deverror('nyi', &
            'Sorry, open-shell systems not yet supported.  Working on it.')

    if (TIMEIT) call write_time('Check AO normalization')
    call fix_norm_AOs(iu_out, molDB%n_at, orbDB%n_ao, molDB%at_crd, &
                      bsetDB_orig, err, DEBUG)
    call run%check(err)

    ! Check if basis set Cartesian.
    if (.not.bsetDB_orig%is_cart()) then
        call convert_pure2cart(bsetDB_orig, bsetDB_cart)
        bsetDB => bsetDB_cart
        allocate(density(bsetDB_cart%n_basis, bsetDB_cart%n_basis))
        call convert_pure2cart(bsetDB_orig, molDB%el_dens(:,:,1), density)
    else
        bsetDB => bsetDB_orig
        allocate(density(bsetDB_orig%n_basis, bsetDB_orig%n_basis))
        density = molDB%el_dens(:,:,1)
    end if

    ! Initialize cube
    call cubeDB%init(molDB)

    call sec_header(1, 'Grid Parameters')
    if (allocated(opts%file_grid)) then
        call cubeDB%grid%read(opts%file_grid)
        call run%check(cubeDB%grid%error, &
            'Error found while parsing grid input file')
    else
        call cubeDB%grid%init('density', opts%grid_sparsity, molDB=molDB, &
                              bsetDB=bsetDB, e_dens=density)
    end if
    call write_griddata(iu_out, cubeDB%grid)

    call cubeDB%gendata(bsetDB, density, density_sum, .true.)
    if (opts%partition) then
        call run%error%raise_deverror('nyi', 'Partition not yet available')
    end if

    cube_fname = 'density.cube'
    open(newunit=iu, file=cube_fname, status='replace', action='write', &
         iostat=ios)
    if (ios /= 0) &
        call run%error%raise_error('file', 'open', &
            'Failed to open cube file for writing')
    call cubeDB%write(iu)
    close(iu)

    write(iu_out,'("File ",a," has been generated!")') cube_fname

    if (opts%print_level == 'debug') then
        write(iu_out, '(/A, E14.6)') 'Integrated electronic density: ', &
            density_sum
    endif

contains

! ======================================================================

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

! ======================================================================

subroutine parse_opts(optsDB)
    !! Parses command line options and updates information.
    !!
    !! Builds options parser and parse user-options, setting
    !! default values where necessary.
    type(Params), intent(out) :: optsDB

    logical :: exists
    character(len=1024) :: argval
    class(CmdLineArgsDB), allocatable :: parser

    ! Build options parser and check arguments
    parser = CmdLineArgsDB(progname='density_cube')
    call run%check(parser%error, 'Failed to initialize the command-line parser')

    call parser%add_arg_char('string', label='grid_filename', shortname='-g', &
                             longname='--grid', &
                             help='Input file to set grid parameters', &
                             required=.false.)

    call parser%add_arg_char('string', label='checkpointfile_filename', &
                             help='Gaussian formatted checkpoint file')

    call parser%add_arg_char('string', label='grid_sparsity', shortname='-d', &
                             longname='--sparsity', &
                             help='Set grid sparsity', required=.false.)

    call parser%add_arg_char('string', label='print_level', shortname='-p', &
                             longname='--print-level', &
                             help='Set the print level', required=.false.)

    call parser%add_arg_char('string', label='partitioning', shortname='-p', &
                             longname='--partitioning', &
                             help='Generate partitioning cube file', &
                             required=.false.)

    if (TIMEIT) call write_time('Option parser')
    call parser%parse_args()
    call run%check(parser%error, 'Failed to parse the command-line arguments')

    call sec_header(1, 'Simulation Parameters')
    write(iu_out, '(1x)')

    if (parser%is_user_set('grid_filename')) then
        call parser%get_value('grid_filename', argval)
        optsDB%file_grid = trim(argval)
        call write_param('Grid input filename', optsDB%file_grid)
    end if

    if (parser%is_user_set('print_level')) then
        call parser%get_value('print_level', argval)
        if (trim(argval) == 'debug') then
            optsDB%print_level = 'debug'
            call write_param('Print level', optsDB%print_level)
        end if
    endif

    if (parser%is_user_set('partitioning')) then
        call parser%get_value('partitioning', argval)
        if (trim(argval) == 'true') then
            optsDB%partition = .true.
            call write_param('Partitioning .cube', optsDB%partition)
        end if
    end if

    if (parser%is_user_set('grid_sparsity')) then
        call parser%get_value('grid_sparsity', argval)
        optsDB%grid_sparsity = trim(argval)
        call write_param('Grid sparsity', optsDB%grid_sparsity)
    else
        optsDB%grid_sparsity = 'default'
    end if

    call parser%get_value('checkpointfile_filename', argval)
    optsDB%file_fchk = trim(argval)
    inquire(file=optsDB%file_fchk, exist=exists)
    if (.not.exists) then
        argval = ' '
        write(argval, '("File ",a," does not exist.")') &
            trim(optsDB%file_fchk)
        call run%error%raise_error('file', 'not found', trim(argval))
    end if
    call write_param('Input filename', optsDB%file_fchk)

end subroutine parse_opts

! ======================================================================

subroutine write_moldata(mol_db, bset_db, orb_db)
    !! Write molecular data.
    !!
    !! Writes some of the molecular data stored in the moldata module
    !! for mcd_tensor.
    class(MoleculeDB), intent(in) :: mol_db
    !! Molecule specifications database.
    class(BasisSetDB), intent(in) :: bset_db
    !! Basis set specifications database.
    class(OrbitalsDB), intent(in) :: orb_db
    !! Molecular orbitals database.

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

    if (orb_db%openshell) then
        tag_open = 'open'
    else
        tag_open = 'closed'
    end if

    if (bset_db%pureD) then
        tag_bfunD = 'pure'
    else
        tag_bfunD = 'Cartesian'
    end if
    if (bset_db%pureF) then
        tag_bfunF = 'pure'
    else
        tag_bfunF = 'Cartesian'
    end if

    write(iu_out, 1100) mol_db%n_at, mol_db%charge, mol_db%multip, &
        tag_open, bset_db%n_basis, tag_bfunD, tag_bfunF, orb_db%n_ao

    call sec_header(2, 'Atomic coordinates')
    call prt_coord(mol_db%n_at, mol_db%at_lab, mol_db%at_crd, &
                    mol_db%at_mas)

end subroutine write_moldata

! ======================================================================

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

! ======================================================================

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

! ======================================================================

    subroutine write_param_real(param, value, subparam)
        !! Writes value of real-type parameter.
        !!
        !! Formats and writes one parameter of simulation with the
        !! indicated value.
        implicit none

        character(len=*), intent(in) :: param
        !! Parameter name (with info like unit).
        real(realwp), intent(in) :: value
        !! Value associated to parameter
        logical, intent(in), optional :: subparam
        !! Parameter is a sub-parameter type

        1000 format(a,1x,f0.4)
        1010 format(a,es13.6)

        if (abs(value) < f10m2 .or. abs(value) > f10p4) then
            write(iu_out, 1010) fmt_param(param, subparam), value
        else
            write(iu_out, 1000) fmt_param(param, subparam), value
        end if
    end subroutine write_param_real

! ======================================================================

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

! ======================================================================

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

! ======================================================================

subroutine write_griddata(iunit, gridDB)
    !! Write grid data.
    !!
    !! Writes grid parameters stored in a grid specifications
    !!database.
    integer, intent(in) :: iunit
        !! Identifier unit to the file connected for output.
    type(CubeGridDB), intent(in) :: gridDB
        !! Grid specifications database.

    1100 format(/, &
            'Grid type:                       ',a,/, &
            'Grid origin x:                  ',f10.6,/, &
            'Grid origin y:                  ',f10.6,/, &
            'Grid origin z:                  ',f10.6,/, &
            'Number of points along x-axis:  ',i4,/,    &
            'Number of points along y-axis:  ',i4,/,    &
            'Number of points along z-axis:  ',i4,/,    &
            'Step size along x-axis:         ',f10.6,/, &
            'Step size along y-axis:         ',f10.6,/, &
            'Step size along z-axis:         ',f10.6,/)

    write(iunit, 1100) gridDB%shape, gridDB%min, gridDB%n_points, &
        gridDB%step_size

end subroutine write_griddata

! ======================================================================

end program density_cube
