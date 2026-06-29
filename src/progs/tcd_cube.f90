program tcd_cube
    use arrays, only: antisymm_sum
    use basisset, only: convert_pure2cart, fix_norm_AOs, get_cart_L_norms_sh
    use cubegen, only: CubeTCDDB, CubeGridDB
    use datatypes, only: BasisSetDB, ExcitationDB, MoleculeDB, OrbitalsDB, &
        PrimitiveFunction
    use input, only: DataFile
    use math, only: cross
    use numeric, only: realwp, f0, f1, f2, f10m2, f10p4, small
    use output, only: iu_out, prt_coord, prt_mat, sec_header
    use parse_cmdline, only: CmdLineArgsDB
    use run_env, only: ErrorHandle, run
    use string, only: num_chars_int

    use datatypes, only: PropertyDB

    implicit none

    type Params
        character(len=:), allocatable :: file_fchk
        character(len=:), allocatable :: file_grid
        character(len=:), allocatable :: grid_sparsity
        character(len=:), allocatable :: print_level
        character(len=:), allocatable :: outfile
            !! Optional output filename.
        logical :: get_axial = .false.
        logical :: length_gauge = .false.
        integer :: min_state = -1
        integer :: max_state = -1
        logical :: grid_all = .false.
        integer :: grid_state = -1
        logical :: align_edtm = .false.
        real(realwp), dimension(3) :: edtm_ref = f0
        integer :: edtm_axis = 0
    end type Params


    logical, parameter :: TIMEIT = .False.
    logical, parameter :: DEBUG = .False.
    integer :: i, ios, istate, iu, lstate
    character(len=256) :: header_base, header_with_integral
    real(realwp), dimension(:,:), allocatable :: tmp_dens, trans_dens, ao_rot
    real(realwp), dimension(3) :: dtm_sum, edtm, mdtm, cube_integral
    logical :: convert
    character(len=:), allocatable :: cube_fname, base_name
    character(len=*), parameter :: &
        PROGTITLE = 'Transition Current Density Cube Generator'
    character(len=80) :: fmt
    character(len=512) :: msg
    type(ErrorHandle) :: err
    type(BasisSetDB), target :: bsetDB_orig, bsetDB_cart
    type(BasisSetDB), pointer :: bsetDB
    class(CubeGridDB), allocatable :: gridDB, grid_tmp
    type(CubeTCDDB) :: cubeDB
    type(DataFile) :: dfile
    type(ExcitationDB) :: excDB
    type(MoleculeDB) :: molDB
    type(OrbitalsDB) :: orbDB
    type(Params) :: opts
    real(realwp), dimension(3,3) :: edtm_rotmat

    class(PropertyDB), allocatable :: propDB

    interface write_param
        procedure write_param_bool, write_param_int, write_param_real, &
            write_param_char
    end interface write_param

    1200 format(/, &
    ' > Number of excited states : ',i0,/, &
    ' > Reference excited state  : ',i0)

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

    excDB = dfile%get_exc_data(get_dens=.true.)
    call run%check(dfile%error, &
        'Parsing electronic excitation data in file failed')
    if (.not.excdb%prop_loaded) then
        call run%error%raise_error('data', 'missing', &
            'Failure to extract excited-state data')
    else if (.not.excdb%dens_loaded) then
        call run%error%raise_error('data', 'missing', &
            'Failure to extract transition density data')
    end if

    ! Let us now correct the final state if not set
    if (opts%max_state == -1) opts%max_state = excDB%n_states

    if (opts%grid_state > excDB%n_states) then
        write(msg, '("Requested grid state ",i0, &
            " exceeds number of excited states")') opts%grid_state
        call run%error%raise_error('value', 'wrong', trim(msg))
    end if

    if (opts%align_edtm) then
        if (allocated(opts%file_grid)) then
            call run%error%raise_error('opt', 'conflict', &
                'The -t/--align-edtm option cannot be used together with &
                &-g/--grid because explicit grid parameters are not rotated.')
        end if
        call build_alignment_rotation(opts%edtm_ref, opts%edtm_axis, edtm_rotmat)
        call rotate_molecule(molDB, edtm_rotmat)
        call write_param('EDTM alignment axis', axis_label(opts%edtm_axis))
        call write_param('EDTM alignment vector', fmt_real3(opts%edtm_ref))
    end if

    call write_moldata(molDB, bsetDB_orig, orbDB)

    call sec_header(1, 'Transition Data')
    write(iu_out, 1200) excDB%n_states, excDB%id_state

    if (TIMEIT) call write_time('Check AO normalization')
    call fix_norm_AOs(iu_out, molDB%n_at, orbDB%n_ao, molDB%at_crd, &
                      bsetDB_orig, err, DEBUG)
    call run%check(err)

    ! Check if basis set Cartesian.
    if (.not.bsetDB_orig%is_cart()) then
        convert = .true.
        call convert_pure2cart(bsetDB_orig, bsetDB_cart)
        bsetDB => bsetDB_cart
    else
        convert = .false.
        bsetDB => bsetDB_orig
    end if

    if (opts%align_edtm) then
        call build_ao_rotation_matrix(bsetDB, edtm_rotmat, ao_rot)
    end if

    allocate(tmp_dens(bsetDB%n_basis, bsetDB%n_basis), &
             trans_dens(bsetDB%n_basis, bsetDB%n_basis))

    ! Store grid data if given in filename
    if (allocated(opts%file_grid)) then
        allocate(gridDB)
        call gridDB%read(opts%file_grid)
        call run%check(gridDB%error, &
            'Error found while parsing grid input file')
    end if

    call sec_header(1, 'Cube Generation')

    ! Construct base name for cube filenames
    i = index(opts%file_fchk, '.', back=.true.)
    base_name = opts%file_fchk(:i-1)
    ! find length to print highest state
    lstate = num_chars_int(opts%max_state)
    ! Format: <base>_AT/PT_I<lstate>.<lstate>.cube or <base>_TDD_<state>.cube
    ! Allocate enough space to store partitioning file name as well
    allocate(character(len=len(base_name)+20+2*lstate) :: cube_fname)
    if (opts%length_gauge) then
        write(fmt, '(''(a,"_TDD_",i0,".cube")'')')
    else if (opts%get_axial) then
        write(fmt, '(''(a,"_AT_",i'',i0,''.'',i0,'',".cube")'')') &
            lstate, lstate
    else
        write(fmt, '(''(a,"_PT_",i'',i0,''.'',i0,'',".cube")'')') &
            lstate, lstate
    end if

    if (opts%length_gauge) then
        call cubeDB%init(molDB, 'ETDD')
    else
        call cubeDB%init(molDB, 'ETCD')
    end if
    if (.not.allocated(gridDB)) then
        allocate(gridDB)

    end if

    if (.not.allocated(opts%file_grid)) then
        if (opts%grid_all) then
            istate = opts%min_state
            if (orbDB%openshell) then
                tmp_dens = (excDB%g2e_dens(:,:,1,istate) &
                            + excDB%g2e_dens(:,:,2,istate)) / f2
            else
                tmp_dens = excDB%g2e_dens(:,:,1,istate)
            end if
            if (convert) then
                call convert_pure2cart(bsetDB_orig, tmp_dens, trans_dens)
            else
                trans_dens = tmp_dens
            end if
            if (opts%align_edtm) call rotate_ao_density(ao_rot, trans_dens)
            call antisymm_sum(bsetDB%n_basis, trans_dens)
            call gridDB%init(cubeDB%type, opts%grid_sparsity, &
                             molDB=molDB, bsetDB=bsetDB, &
                             e_trans_dens=trans_dens)
            if (opts%max_state > opts%min_state) then
                allocate(grid_tmp)
                do istate = opts%min_state+1, opts%max_state
                    if (orbDB%openshell) then
                        tmp_dens = (excDB%g2e_dens(:,:,1,istate) &
                                    + excDB%g2e_dens(:,:,2,istate)) / f2
                    else
                        tmp_dens = excDB%g2e_dens(:,:,1,istate)
                    end if
                    if (convert) then
                        call convert_pure2cart(bsetDB_orig, tmp_dens, trans_dens)
                    else
                        trans_dens = tmp_dens
                    end if
                    if (opts%align_edtm) call rotate_ao_density(ao_rot, trans_dens)
                    call antisymm_sum(bsetDB%n_basis, trans_dens)
                    call grid_tmp%init(cubeDB%type, opts%grid_sparsity, &
                                       molDB=molDB, bsetDB=bsetDB, &
                                       e_trans_dens=trans_dens)
                    if (grid_tmp%min(1) < gridDB%min(1)) gridDB%min(1) = grid_tmp%min(1)
                    if (grid_tmp%min(2) < gridDB%min(2)) gridDB%min(2) = grid_tmp%min(2)
                    if (grid_tmp%min(3) < gridDB%min(3)) gridDB%min(3) = grid_tmp%min(3)
                    if (grid_tmp%max(1) > gridDB%max(1)) gridDB%max(1) = grid_tmp%max(1)
                    if (grid_tmp%max(2) > gridDB%max(2)) gridDB%max(2) = grid_tmp%max(2)
                    if (grid_tmp%max(3) > gridDB%max(3)) gridDB%max(3) = grid_tmp%max(3)
                end do
                gridDB%step_size = (gridDB%max - gridDB%min) / &
                    real(gridDB%n_points - 1, realwp)
            end if
            call write_griddata(iu_out, gridDB)
        else if (opts%grid_state > 0) then
            istate = opts%grid_state
            if (orbDB%openshell) then
                tmp_dens = (excDB%g2e_dens(:,:,1,istate) &
                            + excDB%g2e_dens(:,:,2,istate)) / f2
            else
                tmp_dens = excDB%g2e_dens(:,:,1,istate)
            end if
            if (convert) then
                call convert_pure2cart(bsetDB_orig, tmp_dens, trans_dens)
            else
                trans_dens = tmp_dens
            end if
            if (opts%align_edtm) call rotate_ao_density(ao_rot, trans_dens)
            call antisymm_sum(bsetDB%n_basis, trans_dens)
            call gridDB%init(cubeDB%type, opts%grid_sparsity, &
                             molDB=molDB, bsetDB=bsetDB, &
                             e_trans_dens=trans_dens)
            call write_griddata(iu_out, gridDB)
        end if
    else
        call write_griddata(iu_out, gridDB)
    end if

    do istate = opts%min_state, opts%max_state
        ! Initialize cube
        if (opts%length_gauge) then
            call cubeDB%init(molDB, 'ETDD', [istate])
        else
            call cubeDB%init(molDB, 'ETCD', [istate])
        end if

        write(iu_out, '(/,"Doing state num. ",i0)') istate

        if (orbDB%openshell) then
            ! sum alpha and beta densities, 2.0 comes from sqrt(2.0)
            tmp_dens = (excDB%g2e_dens(:,:,1,istate) &
                        + excDB%g2e_dens(:,:,2,istate)) / f2
        else
            tmp_dens = excDB%g2e_dens(:,:,1,istate)
        end if
        if (convert) then
            call convert_pure2cart(bsetDB_orig, tmp_dens, trans_dens)
        else
            trans_dens = tmp_dens
        end if
        if (opts%align_edtm) call rotate_ao_density(ao_rot, trans_dens)
        call antisymm_sum(bsetDB%n_basis, trans_dens)

        if (.not.allocated(opts%file_grid)) then
            if (.not.opts%grid_all .and. opts%grid_state <= 0) then
                call gridDB%init(cubeDB%type, opts%grid_sparsity, &
                                 molDB=molDB, bsetDB=bsetDB, &
                                 e_trans_dens=trans_dens)
                call write_griddata(iu_out, gridDB)
            end if
        end if

        cubeDB%grid = gridDB
        call cubeDB%gendata(bsetDB, trans_dens, dtm_sum, &
                            gen_axial=opts%get_axial, &
                            length_gauge=opts%length_gauge)
        cube_integral = [sum(cubeDB%data(1,:,:,:)), &
                         sum(cubeDB%data(2,:,:,:)), &
                         sum(cubeDB%data(3,:,:,:))]
        cube_integral = cube_integral * &
            product((cubeDB%grid%max + cubeDB%grid%step_size) - cubeDB%grid%min) &
            / product(cubeDB%grid%n_points)

        write(cube_fname, fmt) base_name, istate
        open(newunit=iu, file=cube_fname, status='replace', action='write', &
             iostat=ios)
        if (ios /= 0) &
            call run%error%raise_error('file', 'open', &
                'Failed to open cube file for writing')

        header_base = cubeDB%header
        write(header_with_integral, '(a," | Integral: ",3es20.12)') &
            trim(header_base), cube_integral
        cubeDB%header = header_with_integral
        call cubeDB%write(iu)
        cubeDB%header = header_base
        close(iu)

        write(iu_out, '("File ",a," has been generated!")') cube_fname

        if (opts%print_level == 'debug') then

            if (opts%get_axial) then

                write(iu_out, '(/A, I0, A, E14.6)') 'Excitation energy (au) for the state ', &
                      istate, ' is:', excDB%g2e_energy(istate)

                mdtm = dtm_sum

                write(iu_out, '(/A, I0, A, E14.6, E14.6, E14.6)') 'MDTM (au) for the state ', &
                      istate, ' is:', mdtm

            else

                write(iu_out, '(/A, I0, A, E14.6)') 'Excitation energy (au) for the state ', &
                      istate, ' is:', excDB%g2e_energy(istate)

                if (opts%length_gauge) then
                    edtm = dtm_sum
                    write(iu_out, '(/A, I0, A, E14.6, E14.6, E14.6)') 'Length gauge EDTM (au) for the state ', &
                          istate, ' is:', edtm
                else
                    edtm = dtm_sum / excDB%g2e_energy(istate)
                    write(iu_out, '(/A, I0, A, E14.6, E14.6, E14.6)') 'Velocity gauge EDTM (au) for the state ', &
                          istate, ' is:', edtm
                end if

            endif
        endif

    end do

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
    character(len=:), dimension(:), allocatable :: argvals
    class(CmdLineArgsDB), allocatable :: parser
    integer :: ios

    ! Build options parser and check arguments
    parser = CmdLineArgsDB(progname='tcd')
    call run%check(parser%error, 'Failed to initialize the command-line parser')

    call parser%add_arg_char('string', label='checkpointfile_filename', &
                             help='Gaussian formatted checkpoint file')

    call parser%add_arg_char('string', label='axial_tensor', &
                             shortname='-a', longname='--axial-tensor', &
                             required=.false., &
                             help='Generate axial tensor .cube file')

    call parser%add_arg_char('string', label='grid_sparsity', &
                             shortname='-d', longname='--sparsity', &
                             required=.false., &
                             help='Set grid sparsity')

    call parser%add_arg_int('scalar', &
                            min_value=1, &
                            shortname='-e', longname='--max-state', &
                            help='Highest electronic state to include (upper &
                                &bound of summation)')

    call parser%add_arg_char('string', label='grid_filename', &
                             shortname='-g', longname='--grid', &
                             required=.false., &
                             help='Input file to set grid parameters')

    call parser%add_arg_bool('store_true', &
                             shortname='-j', longname='--grid_all', &
                             required=.false., &
                             help='Use a grid encompassing all selected states')

    call parser%add_arg_int('scalar', &
                            shortname='-k', longname='--grid_specific', &
                            min_value=1, required=.false., &
                            help='Use grid parameters generated for the &
                                &specified state')

    call parser%add_arg_bool('store_true', label='length_gauge', &
                             shortname='-L', longname='--length-gauge', &
                             required=.false., &
                             help='Use the length-gauge electric-dipole &
                                 &representation. This generates transition &
                                 &dipole-density maps instead of transition &
                                 &current-density maps.')

    call parser%add_arg_char('string', label='output', &
                             shortname='-o', longname='--output', &
                             help='Name of the file to store the output. &
                                 &Existing content will be overwritten.')

    call parser%add_arg_char('string', label='print_level', &
                             shortname='-p', longname='--print-level', &
                             required=.false., &
                             help='Set the print level')

    call parser%add_arg_int('scalar', &
                            shortname='-s', longname='--min-state', &
                            min_value=1, &
                            help='Lowest electronic state to include (upper &
                                &bound of summation)')

    call parser%add_arg_char('list', label='x y z axis', &
                             shortname='-t', longname='--align-edtm', &
                             min_nvals=4, max_nvals=4, required=.false., &
                             help='Rotate system so EDTM (x y z) aligns to &
                                 &axis int (1=x, 2=y, 3=z)')

    if (TIMEIT) call write_time('Option parser')
    call parser%parse_args()
    call run%check(parser%error, 'Failed to parse command-line arguments')

    ! Output file
    if (parser%is_user_set('output')) then
        call parser%get_value('output', argval)
        optsDB%outfile = trim(argval)
        open(newunit=iu_out, file=optsDB%outfile, action='write')
        call sec_header(-1, PROGTITLE)
    end if

    call sec_header(1, 'Simulation Parameters')
    write(iu_out, '(1x)')

    call parser%get_value('checkpointfile_filename', argval)
    optsDB%file_fchk = trim(argval)
    inquire(file=optsDB%file_fchk, exist=exists)
    if (.not.exists) then
        argval = ' '
        write(argval, '("Error: File ",a," does not exist.")') &
            trim(optsDB%file_fchk)
        call run%error%raise_error('file', 'not found', trim(argval))
    end if
    call write_param('TD input filename', optsDB%file_fchk)

    if (parser%is_user_set('grid_filename')) then
        call parser%get_value('grid_filename', argval)
        optsDB%file_grid = trim(argval)
        call write_param('Grid input filename', optsDB%file_grid)
    endif

    if (parser%is_user_set('grid_sparsity')) then
        call parser%get_value('grid_sparsity', argval)
        optsDB%grid_sparsity = trim(argval)
        call write_param('Grid sparsity', optsDB%grid_sparsity)
    else
        optsDB%grid_sparsity = 'default'
    end if

    if (parser%is_user_set('print_level')) then
        call parser%get_value('print_level', argval)
        if (trim(argval) == 'debug') then
            optsDB%print_level = 'debug'
            call write_param('Print level', optsDB%print_level)
        end if
    endif

    if (parser%is_user_set('axial_tensor')) then
        call parser%get_value('axial_tensor', argval)
        if (trim(argval) == 'true') then
            optsDB%get_axial = .true.
            call write_param('Axial tensor .cube', optsDB%get_axial)
        end if
    endif

    if (parser%is_user_set('length_gauge')) optsDB%length_gauge = .true.

    if (parser%is_user_set('length_gauge') .and. &
            parser%is_user_set('axial_tensor')) then
        call run%error%raise_error('opt', 'conflict', &
            'Options --length-gauge and --axial-tensor cannot be used together')
    end if

    if (optsDB%length_gauge) then
        call write_param('Dipole gauge', 'length')
    else
        call write_param('Dipole gauge', 'velocity')
    end if

    if (parser%is_user_set('min-state')) then
        call parser%get_value('min-state', optsDB%min_state)
    end if

    if (parser%is_user_set('max-state')) then
        call parser%get_value('max-state', optsDB%max_state)
    end if

    if (parser%is_user_set('grid_all') .and. &
            parser%is_user_set('grid_specific')) then
        call run%error%raise_error('opt', 'conflict', &
            'Options --grid_all and --grid_specific cannot be used together')
    end if

    if (parser%is_user_set('grid_all')) then
        optsDB%grid_all = .true.
        call write_param('Common grid for all states', optsDB%grid_all)
    end if

    if (parser%is_user_set('grid_specific')) then
        call parser%get_value('grid_specific', optsDB%grid_state)
        call write_param('Grid taken from state', optsDB%grid_state)
    end if

    if (parser%is_user_set('min-state') &
            .and. parser%is_user_set('max-state')) then
        if (optsDB%min_state == optsDB%max_state) then
            call write_param('Requested excited state', optsDB%min_state)
        else
            call write_param('Lowest excited state', optsDB%min_state)
            call write_param('Highest excited state', optsDB%max_state)
        endif

    elseif (parser%is_user_set('max-state')) then
        optsDB%min_state = 1
        call write_param('Lowest excited state', 'default to first state')
        call write_param('Highest excited state', optsDB%max_state)

    elseif (parser%is_user_set('min-state')) then
        call write_param('Lowest excited state', optsDB%min_state)
        call write_param('Highest excited state', 'include all')

    else
        optsDB%min_state = 1
        call write_param('Lowest excited state', 'default to first state')
        call write_param('Highest excited state', 'include all')

    endif

    if (parser%is_user_set('x y z axis')) then
        call parser%get_value('x y z axis', argvals)
        read(argvals(1), *, iostat=ios) optsDB%edtm_ref(1)
        if (ios /= 0) call run%error%raise_error('opt', 'value', &
            'Invalid x component for -t/--align-edtm.')
        read(argvals(2), *, iostat=ios) optsDB%edtm_ref(2)
        if (ios /= 0) call run%error%raise_error('opt', 'value', &
            'Invalid y component for -t/--align-edtm.')
        read(argvals(3), *, iostat=ios) optsDB%edtm_ref(3)
        if (ios /= 0) call run%error%raise_error('opt', 'value', &
            'Invalid z component for -t/--align-edtm.')
        read(argvals(4), *, iostat=ios) optsDB%edtm_axis
        if (ios /= 0) call run%error%raise_error('opt', 'value', &
            'Invalid axis for -t/--align-edtm. Expected 1, 2 or 3.')
        if (optsDB%edtm_axis < 1 .or. optsDB%edtm_axis > 3) then
            call run%error%raise_error('opt', 'value', &
                'The alignment axis must be 1, 2 or 3.')
        end if
        if (sum(abs(optsDB%edtm_ref)) <= 10*small) then
            call run%error%raise_error('opt', 'value', &
                'The EDTM vector provided with -t/--align-edtm must be non-zero.')
        end if
    end if

end subroutine parse_opts

! ======================================================================

subroutine build_alignment_rotation(edtm_vec, axis_id, rotmat)
    implicit none

    real(realwp), dimension(3), intent(in) :: edtm_vec
    integer, intent(in) :: axis_id
    real(realwp), dimension(3,3), intent(out) :: rotmat

    real(realwp) :: cval, sval
    real(realwp), dimension(3) :: source_vec, target_vec, axis_vec
    real(realwp), dimension(3,3) :: kmat, kmat2

    source_vec = edtm_vec / norm2(edtm_vec)
    target_vec = f0
    target_vec(axis_id) = f1
    cval = dot_product(source_vec, target_vec)

    if (abs(cval - f1) <= 10*small) then
        rotmat = f0
        rotmat(1,1) = f1
        rotmat(2,2) = f1
        rotmat(3,3) = f1
        return
    end if

    if (abs(cval + f1) <= 10*small) then
        axis_vec = orthogonal_unit_vector(source_vec)
        call cross_matrix(axis_vec, kmat)
        kmat2 = matmul(kmat, kmat)
        rotmat = -reshape([f1, f0, f0, f0, f1, f0, f0, f0, f1], [3,3]) + &
                 f2*matmul(reshape(axis_vec, [3,1]), reshape(axis_vec, [1,3]))
        return
    end if

    axis_vec = cross(source_vec, target_vec)
    sval = norm2(axis_vec)
    axis_vec = axis_vec / sval
    call cross_matrix(axis_vec, kmat)
    kmat2 = matmul(kmat, kmat)
    rotmat = reshape([f1, f0, f0, f0, f1, f0, f0, f0, f1], [3,3]) + &
             sval*kmat + (f1 - cval)*kmat2
end subroutine build_alignment_rotation

! ======================================================================

subroutine rotate_molecule(mol_db, rotmat)
    implicit none

    class(MoleculeDB), intent(inout) :: mol_db
    real(realwp), dimension(3,3), intent(in) :: rotmat

    mol_db%at_crd = matmul(rotmat, mol_db%at_crd)
end subroutine rotate_molecule

! ======================================================================

subroutine rotate_ao_density(ao_rotmat, density)
    implicit none

    real(realwp), dimension(:,:), intent(in) :: ao_rotmat
    real(realwp), dimension(:,:), intent(inout) :: density

    real(realwp), dimension(:,:), allocatable :: tmp

    allocate(tmp(size(density, 1), size(density, 2)))
    tmp = matmul(ao_rotmat, density)
    density = matmul(tmp, transpose(ao_rotmat))
    deallocate(tmp)
end subroutine rotate_ao_density

! ======================================================================

subroutine build_ao_rotation_matrix(bset_db, rotmat, ao_rotmat)
    implicit none

    class(BasisSetDB), intent(in) :: bset_db
    real(realwp), dimension(3,3), intent(in) :: rotmat
    real(realwp), dimension(:,:), allocatable, intent(out) :: ao_rotmat

    integer :: ia, iprim, beg_ao, end_ao, nprim_sh
    integer :: ndim_sh

    allocate(ao_rotmat(bset_db%n_basis, bset_db%n_basis))
    ao_rotmat = f0

    beg_ao = 1
    do ia = 1, size(bset_db%nprim_per_at)
        iprim = 1
        do while (iprim <= bset_db%nprim_per_at(ia))
            nprim_sh = 1
            do while (.not.bset_db%info(ia, iprim + nprim_sh - 1)%shell_last)
                nprim_sh = nprim_sh + 1
            end do
            ndim_sh = bset_db%info(ia, iprim)%ndim
            end_ao = beg_ao + ndim_sh - 1
            call build_shell_rotation_matrix( &
                bset_db%info(ia, iprim:iprim+nprim_sh-1), rotmat, &
                ao_rotmat(beg_ao:end_ao, beg_ao:end_ao))
            beg_ao = end_ao + 1
            iprim = iprim + nprim_sh
        end do
    end do
end subroutine build_ao_rotation_matrix

! ======================================================================

subroutine build_shell_rotation_matrix(shell_prims, rotmat, shell_rot)
    implicit none

    type(PrimitiveFunction), dimension(:), intent(in) :: shell_prims
    real(realwp), dimension(3,3), intent(in) :: rotmat
    real(realwp), dimension(:,:), intent(out) :: shell_rot

    integer :: l_ang

    shell_rot = f0
    if (shell_prims(1)%shelltype == 'SP') then
        shell_rot(1,1) = f1
        call build_cart_rotation_block(shell_prims, 1, rotmat, &
            shell_rot(2:4, 2:4))
    else
        l_ang = shell_prims(1)%L
        call build_cart_rotation_block(shell_prims, l_ang, rotmat, shell_rot)
    end if
end subroutine build_shell_rotation_matrix

! ======================================================================

subroutine build_cart_rotation_block(shell_prims, l_ang, rotmat, dblock)
    implicit none

    type(PrimitiveFunction), dimension(:), intent(in) :: shell_prims
    integer, intent(in) :: l_ang
    real(realwp), dimension(3,3), intent(in) :: rotmat
    real(realwp), dimension(:,:), intent(out) :: dblock

    integer :: iold, inew, npow
    integer, dimension(:,:), allocatable :: powers
    real(realwp), dimension(:), allocatable :: norms
    real(realwp), dimension(3,3) :: invrot
    real(realwp), dimension(:), allocatable :: coeffs

    npow = (l_ang + 1) * (l_ang + 2) / 2
    allocate(powers(3, npow), norms(npow), coeffs(npow))
    powers = list_cart_powers(shell_prims, l_ang)
    norms = get_cart_L_norms_sh(l_ang, npow)
    invrot = transpose(rotmat)

    dblock = f0
    do iold = 1, npow
        call expand_rotated_monomial(invrot, powers(:, iold), powers, coeffs)
        do inew = 1, npow
            dblock(inew, iold) = coeffs(inew) * norms(inew) / norms(iold)
        end do
    end do

    deallocate(powers, norms, coeffs)
end subroutine build_cart_rotation_block

! ======================================================================

subroutine expand_rotated_monomial(rotmat, power_old, powers_all, coeffs)
    implicit none

    real(realwp), dimension(3,3), intent(in) :: rotmat
    integer, dimension(3), intent(in) :: power_old
    integer, dimension(:,:), intent(in) :: powers_all
    real(realwp), dimension(:), intent(out) :: coeffs

    integer :: iterm, ipow, nterm
    integer, dimension(:,:), allocatable :: cur_pows, new_pows
    real(realwp), dimension(:), allocatable :: cur_coef, new_coef

    nterm = 1
    allocate(cur_pows(3, 1), cur_coef(1))
    cur_pows(:, 1) = 0
    cur_coef(1) = f1

    do ipow = 1, 3
        if (power_old(ipow) <= 0) cycle
        call multiply_linear_power(cur_pows, cur_coef, nterm, rotmat(ipow, :), &
                                   power_old(ipow), new_pows, new_coef)
        deallocate(cur_pows, cur_coef)
        call move_alloc(new_pows, cur_pows)
        call move_alloc(new_coef, cur_coef)
        nterm = size(cur_coef)
    end do

    coeffs = f0
    do iterm = 1, nterm
        ipow = find_power_index(cur_pows(:, iterm), powers_all)
        coeffs(ipow) = coeffs(ipow) + cur_coef(iterm)
    end do

    deallocate(cur_pows, cur_coef)
end subroutine expand_rotated_monomial

! ======================================================================

subroutine multiply_linear_power(cur_pows, cur_coef, ncur, linear_coef, npower, &
                                 out_pows, out_coef)
    implicit none

    integer, dimension(:,:), intent(in) :: cur_pows
    real(realwp), dimension(:), intent(in) :: cur_coef
    integer, intent(in) :: ncur
    real(realwp), dimension(3), intent(in) :: linear_coef
    integer, intent(in) :: npower
    integer, dimension(:,:), allocatable, intent(out) :: out_pows
    real(realwp), dimension(:), allocatable, intent(out) :: out_coef

    integer :: ix, iy, iz, iterm, nout
    real(realwp) :: pref

    nout = ncur * ((npower + 1) * (npower + 2) / 2)
    allocate(out_pows(3, nout), out_coef(nout))
    nout = 0
    do iterm = 1, ncur
        do ix = 0, npower
            do iy = 0, npower - ix
                iz = npower - ix - iy
                nout = nout + 1
                pref = real(multinomial3(npower, ix, iy, iz), realwp) * &
                       linear_coef(1)**ix * linear_coef(2)**iy * linear_coef(3)**iz
                out_pows(:, nout) = cur_pows(:, iterm) + [ix, iy, iz]
                out_coef(nout) = cur_coef(iterm) * pref
            end do
        end do
    end do
end subroutine multiply_linear_power

! ======================================================================

function multinomial3(n, nx, ny, nz) result(res)
    implicit none

    integer, intent(in) :: n, nx, ny, nz
    integer :: res

    res = factorial_int(n) &
        / (factorial_int(nx) * factorial_int(ny) * factorial_int(nz))
end function multinomial3

! ======================================================================

function factorial_int(n) result(res)
    implicit none

    integer, intent(in) :: n
    integer :: res, i

    res = 1
    do i = 2, n
        res = res * i
    end do
end function factorial_int

! ======================================================================

function find_power_index(power, powers_all) result(idx)
    implicit none

    integer, dimension(3), intent(in) :: power
    integer, dimension(:,:), intent(in) :: powers_all
    integer :: idx

    do idx = 1, size(powers_all, 2)
        if (all(powers_all(:, idx) == power)) return
    end do
    call run%error%raise_deverror('gen', &
        'Internal error while building AO rotation matrix.')
end function find_power_index

! ======================================================================

function list_cart_powers(shell_prims, l_ang) result(powers)
    implicit none

    type(PrimitiveFunction), dimension(:), intent(in) :: shell_prims
    integer, intent(in) :: l_ang
    integer, dimension(3, (l_ang + 1) * (l_ang + 2) / 2) :: powers

    select case (l_ang)
    case (0)
        powers(:, 1) = [0, 0, 0]
    case (1)
        powers(:, 1) = [1, 0, 0]
        powers(:, 2) = [0, 1, 0]
        powers(:, 3) = [0, 0, 1]
    case (2)
        powers(:, 1) = [2, 0, 0]
        powers(:, 2) = [0, 2, 0]
        powers(:, 3) = [0, 0, 2]
        powers(:, 4) = [1, 1, 0]
        powers(:, 5) = [1, 0, 1]
        powers(:, 6) = [0, 1, 1]
    case (3)
        powers(:, 1) = [3, 0, 0]
        powers(:, 2) = [0, 3, 0]
        powers(:, 3) = [0, 0, 3]
        powers(:, 4) = [2, 1, 0]
        powers(:, 5) = [2, 0, 1]
        powers(:, 6) = [1, 2, 0]
        powers(:, 7) = [0, 2, 1]
        powers(:, 8) = [1, 0, 2]
        powers(:, 9) = [0, 1, 2]
        powers(:,10) = [1, 1, 1]
    case default
        if (allocated(shell_prims(1)%lxyz)) then
            powers = shell_prims(1)%lxyz
        else
            call run%error%raise_deverror('case', &
                'Missing Cartesian powers for high-L shell rotation.')
        end if
    end select
end function list_cart_powers

! ======================================================================

function axis_label(axis_id) result(label)
    implicit none

    integer, intent(in) :: axis_id
    character(len=1) :: label

    select case (axis_id)
    case (1)
        label = 'x'
    case (2)
        label = 'y'
    case (3)
        label = 'z'
    case default
        label = '?'
    end select
end function axis_label

! ======================================================================

function fmt_real3(vec) result(label)
    implicit none

    real(realwp), dimension(3), intent(in) :: vec
    character(len=96) :: label

    write(label, '("(",es13.6,", ",es13.6,", ",es13.6,")")') vec
end function fmt_real3

! ======================================================================

subroutine cross_matrix(axis_vec, kmat)
    implicit none

    real(realwp), dimension(3), intent(in) :: axis_vec
    real(realwp), dimension(3,3), intent(out) :: kmat

    kmat = reshape([ &
        f0, -axis_vec(3), axis_vec(2), &
        axis_vec(3), f0, -axis_vec(1), &
        -axis_vec(2), axis_vec(1), f0 &
    ], [3,3])
end subroutine cross_matrix

! ======================================================================

function orthogonal_unit_vector(vec) result(res)
    implicit none

    real(realwp), dimension(3), intent(in) :: vec
    real(realwp), dimension(3) :: res, trial

    if (abs(vec(1)) <= abs(vec(2)) .and. abs(vec(1)) <= abs(vec(3))) then
        trial = [f1, f0, f0]
    else if (abs(vec(2)) <= abs(vec(3))) then
        trial = [f0, f1, f0]
    else
        trial = [f0, f0, f1]
    end if

    res = cross(vec, trial)
    res = res / norm2(res)
end function orthogonal_unit_vector

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

        if (abs(value) < 1.0e-2_realwp .or. abs(value) > 1.0e4_realwp) then
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

subroutine write_griddata(iunit, grid_db)
    !! Write grid data.
    !!
    !! Writes grid parameters stored in a grid specifications
    !!database.
    integer, intent(in) :: iunit
        !! Identifier unit to the file connected for output.
    type(CubeGridDB), intent(in) :: grid_db
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

    write(iunit, 1100) grid_db%shape, grid_db%min, grid_db%n_points, &
        grid_db%step_size

end subroutine write_griddata

! ======================================================================

end program tcd_cube
