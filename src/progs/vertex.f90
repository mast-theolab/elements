program vertex
    !! Vertical Excitation-Related Transformation to EXtrapolated Geometry
    !!
    !! Simple program to compute the extrapolated geometry from a
    !! vertical position in a vertical excitation.

    use datatypes, only: MoleculeDB, PropertyDB, VibrationsDB
    use exception, only: BaseException, runstat
    use geometry, only: Eckart_orient, superpose
    use input, only: DataFile
    use numeric, only: f0, f1, realwp
    use output, only: iu_out, prt_coord, prt_mat, sec_header
    use parse_cmdline, only: CmdArgDB
    use string, only: locase, upcase
    use vibrational, only: build_modes, set_orientation
    use vibronic, only: Duschinsky_matrix, Duschinsky_shift, extrapolate_geom

    implicit none

    type Params
        character(len=:), allocatable :: file_eq
        character(len=:), allocatable :: file_ve
        character(len=:), allocatable :: file_out
        logical :: do_VG = .false.
    end type Params

    type(DataFile) :: dfile
    type(MoleculeDB) :: mol_eq, mol_ve
    type(Params) :: opts
    type(PropertyDB) :: prop
    type(VibrationsDB) :: vib_eq, vib_ve

    integer :: i, n_atoms, n_vib, ia
    real(realwp), dimension(3,3) :: rot_eq, rot_ve
    real(realwp), dimension(:), allocatable :: kvec
    real(realwp), dimension(:,:), allocatable :: en_grad, ex_geom, jmat
    character(len=*), parameter :: PROGNAME = 'VERTEX'

    ! Start program
    call sec_header(-1, PROGNAME)
1000 format(/, &
        'This program constructs and prints extrapolated atomic coordinates &
        &based on a',/, &
        'vertical, reference geometry and the curvature of the potential &
        &energy surface',/, &
        'of a different electronic state at the same geometry.')
    write(iu_out, 1000)

    call sec_header(1, 'Parsing user options')

    call parse_options(opts)

    call sec_header(1, 'Reading input data')

    call sec_header(2, 'State at equilibrium')

    dfile = DataFile(opts%file_eq)
    if (dfile%has_error()) then
        call runstat%raise_error( &
            'Failed to initialize data file', dfile%get_error())
    end if

    write(iu_out, '(/,"1. Reading molecular data")')
    mol_eq = dfile%get_mol_data()
    if (dfile%has_error()) then
        call runstat%raise_error( &
            'Unable to parse molecular data, check file.', &
            dfile%get_error())
    end if

    write(iu_out, '("2. Reading vibrational data")')
    vib_eq = dfile%get_vib_data()
    if (dfile%has_error()) then
        call runstat%raise_error( &
            'Unable to parse vibrational data, check file.', &
            dfile%get_error())
    end if

    write(iu_out, '("3. Construction of normal modes")')
    prop = dfile%get_data(1, derorder=2)
    call build_modes(prop%data, mol_eq, vib_eq)

    call sec_header(2, 'State in vertical region')

    dfile = DataFile(opts%file_ve)
    if (dfile%has_error()) then
        call runstat%raise_error( &
            'Failed to initialize data file', dfile%get_error())
    end if

    write(iu_out, '(/,"1. Reading molecular data")')
    mol_ve = dfile%get_mol_data()
    if (dfile%has_error()) then
        call runstat%raise_error( &
            'Unable to parse molecular data, check file.', &
            dfile%get_error())
    end if

    write(iu_out, '("2. Reading vibrational data")')
    vib_ve = dfile%get_vib_data()
    if (dfile%has_error()) then
        call runstat%raise_error( &
            'Unable to parse vibrational data, check file.', &
            dfile%get_error())
    end if

    write(iu_out, '("3. Construction of normal modes")')
    prop = dfile%get_data(1, derorder=2)
    call build_modes(prop%data, mol_ve, vib_ve)
    print *, vib_ve%freq

    ! Extract gradient
    write(iu_out, '("4. Extraction of energy gradient")')
    prop = dfile%get_data(1, derorder=1)

    ! Finalization
    n_atoms = mol_eq%n_at
    n_vib = vib_eq%n_vib

    call sec_header(1, 'Superposition of structures')
    ! Set mol_eq to Eckart orientation
    write(iu_out, '(/,"1. Setting equilibrium state to Eckart")')
    call Eckart_orient(mol_eq, .true., rot_mat=rot_eq)
    call prt_coord(mol_eq%n_at, mol_eq%at_lab, mol_eq%at_crd)

    ! Superpose second structure to it
    write(iu_out, '("2. Superposing second structure")')
    call superpose(mol_ve, mol_eq%at_crd, rot_mat=rot_ve)
    call prt_coord(mol_ve%n_at, mol_ve%at_lab, mol_ve%at_crd)

    write(iu_out, '("3. Reorienting relevant quantities")')

    call set_orientation(vib_eq)
    call set_orientation(vib_ve)
    do i = 1, n_vib
        do ia = 1, 3*n_atoms, 3
            vib_eq%L_mat(ia:ia+2,i) = matmul(vib_eq%L_mat(ia:ia+2,i), rot_eq)
            vib_ve%L_mat(ia:ia+2,i) = matmul(vib_ve%L_mat(ia:ia+2,i), rot_ve)
        end do
    end do
    ! Rotate the gradient
    en_grad = matmul(rot_ve, reshape(prop%data, [3, n_atoms]))

    call sec_header(1, 'Computation of Duschinsky matrix and shift vector')

    if (opts%do_VG) then
        jmat = Duschinsky_matrix(vib_eq, vib_ve, is_identity=.true.)
        kvec = Duschinsky_shift(vib_eq, mol_eq, mode='VG', grad2=en_grad)
    else
        jmat = Duschinsky_matrix(vib_eq, vib_ve)
        kvec = Duschinsky_shift(vib_eq, mol_eq, mode='VH', vib2=vib_ve, &
                                Jmat=jmat, grad2=en_grad)
    end if

    call sec_header(1, 'Computation of extrapolated geometry')
    call prt_mat(jmat, n_vib, n_vib)

    ex_geom = extrapolate_geom(mol_eq, vib_eq, kvec)
    call prt_coord(mol_ve%n_at, mol_ve%at_lab, ex_geom)

contains

! ======================================================================

subroutine parse_options(opts_db)
    type(Params), intent(out) :: opts_db

    character(len=1024) :: argval, msg
    class(BaseException), allocatable :: err
    type(CmdArgDB) :: parser

    ! Build option parser for commandline
    parser = CmdArgDB(progname=locase(PROGNAME))
    if (parser%has_error()) then
        err = parser%exception()
        call runstat%raise_error( &
            'Unable to initialize the commandline parser', &
            err%msg())
    end if
    call parser%add_arg_char( &
        'string', label='file_eq', &
        help='Gaussian formatted checkpoint file containing the description &
            &of the state at equilibrium (geom+vib).')
    call parser%add_arg_char( &
        'string', label='file_ve', &
        help='Gaussian formatted checkpoint file containing data of &
            &other state out-of-equilibrium.')
    call parser%add_arg_char( &
        'string', shortname='-m', longname='--method', &
        def_value='VH', &
        help='Model to use: VG (vertical gradient), VH (vertical Hessian, &
             &default)')
    call parser%add_arg_char( &
        'string', shortname='-o', longname='--output', &
        help='Output filename.')

    call parser%parse_args()
    if (parser%has_error()) then
        err = parser%exception()
        call runstat%raise_error( &
            'Failure to parse commandline options', &
            err%msg())
    end if

    ! Check commandline arguments and set information
    call parser%get_value('file_eq', argval)
    opts_db%file_eq = trim(argval)
    call parser%get_value('file_ve', argval)
    opts_db%file_ve = trim(argval)
    if (parser%is_user_set('output')) then
        call parser%get_value('output', argval)
        opts_db%file_out = trim(argval)
        open(newunit=iu_out, file=opts_db%file_out, action='write')
    end if

    if (parser%is_user_set('method')) then
        call parser%get_value('method', argval)
        select case(upcase(trim(argval)))
            case('VG')
                opts_db%do_VG = .true.
            case('VH')
                opts_db%do_VG = .false.
            case default
                write(msg, '("Unrecognized method: ",a)') trim(argval)
                call runstat%raise_error(trim(msg))
        end select
    end if

end subroutine parse_options

! ======================================================================

end program vertex
