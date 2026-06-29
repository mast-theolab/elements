program vtcd_cube
    use iso_fortran_env, only: int32
    use arrays, only: antisymm_sum
    use basisset, only: convert_pure2cart, fix_norm_AOs, get_cart_L_norms_sh
    use cubegen, only: CubeTCDDB, CubeGridDB, FileTCDCube
    use datatypes, only: BasisSetDB, ExcitationDB, MoleculeDB, OrbitalsDB, &
        PrimitiveFunction, PropertyDB, VibrationsDB
    use geometry, only: Eckart_orient
    use input, only: DataFile
    use math, only: cross
    use numeric, only: realwp, f0, f1, f2, fhalf, small
    use output, only: iu_out, prt_coord, sec_header
    use parse_cmdline, only: CmdLineArgsDB
    use physics, only: bohr_radius, phys_conv
    use run_env, only: ErrorHandle, run
    use string, only: locase, num_chars_int

    implicit none

    type BaseParams
        character(len=:), allocatable :: file_td
            !! Name of file containing basic TD data.
        character(len=:), allocatable :: file_frq
            !! Name of file containing frequency data.
        character(len=:), allocatable :: template_nac
            !! Template for the NAC files.
        character(len=:), allocatable :: file_grid
            !! Name of the file containing the grid parameters.
        character(len=:), allocatable :: file_cubes
            !! Name of the file containing the cube specifications.
        character(len=:), allocatable :: moldata_file
            !! Name of the file containing moldata output.
        character(len=:), allocatable :: moldata_in_file
            !! Name of the file containing moldata input.
        character(len=:), allocatable :: outfile
            !! Optional output filename.
        integer :: mode = -1
            !! Vibrational mode of interest.
        integer :: which_grid = 0
            !! Mode used to build grid
            !! -1: read data from file
            !!  0: build automatically using all chosen states
            !! >0: use state `which_grid` to build grid.
        integer :: min_state = -1
            !! Lowest state in the summation (-1: automatic).
        integer :: max_state = -1
            !! Highest state in the summation (-1: automatic).
        logical :: min_state_set = .false.
            !! True if min_state was set by user.
        logical :: max_state_set = .false.
            !! True if max_state was set by user.
        character(len=512) :: grid_sparsity = 'medium'
            !! Sparsity of the cube grid.
            !! TODO: transform into integer: >0 level, -1 read params.
        logical :: eckart = .false.
            !! Rotate to Eckart orientation.
        logical :: read_cubes = .false.
            !! Read cube data from `file_cube`.
        logical :: save_cubes = .false.
            !! Save cube data in `file_cube`.
        character(len=:), allocatable :: print_level
            !! Set the print level (default, debug).
        logical :: get_axial = .false.
            !! Get axial tensor.
        logical :: length_gauge = .false.
            !! Use the length-gauge electric-dipole representation.
        logical :: write_moldata = .false.
            !! Write moldata output.
        logical :: read_moldata = .false.
            !! Read moldata input.
        logical :: align_edtm = .false.
            !! Rotate the system so a reference EDTM points along a Cartesian axis.
        real(realwp), dimension(3) :: edtm_ref = f0
            !! Reference EDTM components used to define the alignment rotation.
        integer :: edtm_axis = 0
            !! Target Cartesian axis: 1=x, 2=y, 3=z.
    end type BaseParams

    logical, parameter :: TIMEIT = .false.
    logical, parameter :: DEBUG = .false.

    integer :: i, ios, istate, iu, lmode, lstate, iat, offset
    character(len=256) :: header_base, header_with_integral
    real(realwp), dimension(3) :: dtm_sum, dtm_1der, cube_integral, &
        edtm_nuc_1der, mdtm_nuc_1der
    real(realwp), dimension(:,:), allocatable :: apt_nuc, aat_nuc
    real(realwp), dimension(:), allocatable :: nac_dQ, rvec
    real(realwp), dimension(:,:), allocatable :: nac_data_all
    real(realwp), dimension(:), allocatable :: nac_vec
    real(realwp), dimension(:,:), allocatable :: tmp_dens, trans_dens, ao_rot
    logical :: convert
    character(len=*), parameter :: &
        PROGTITLE = 'Vibrational Transition Current Density Cube Generator'
    character(len=8) :: nac_ftype
    character(len=80) :: fmt
    character(len=512) :: msg
    character(len=:), allocatable :: base_name, cube_fname
    character(len=:), dimension(:), allocatable :: nac_fnames
    type(ErrorHandle) :: err
    class(CubeGridDB), allocatable :: gridDB
    class(FileTCDCube), allocatable :: tcdDB
    class(PropertyDB), allocatable :: nacDB
    type(BaseParams) :: opts
    type(BasisSetDB), target :: bsetDB_orig, bsetDB_cart
    type(BasisSetDB), pointer :: bsetDB
    type(CubeTCDDB) :: cube_tcd, cube_vtcd
    type(DataFile) :: dfile, nac_file
    type(ExcitationDB) :: excDB
    type(MoleculeDB) :: molDB
    type(OrbitalsDB) :: orbDB
    type(VibrationsDB) :: vibDB
    integer :: nac_min_state, nac_max_state
    logical :: nac_data_is_dq
    real(realwp), dimension(3,3) :: edtm_rotmat

    interface write_param
        procedure write_param_bool, write_param_int, write_param_real, &
            write_param_char
    end interface write_param

    1200 format(/, &
    ' > Number of excited states : ',i0,/, &
    ' > Reference excited state  : ',i0)

    call sec_header(-1, PROGTITLE)

    call parse_options(opts)

    if (opts%read_moldata) then
        call read_moldata_binary(opts%moldata_in_file, excDB, molDB, orbDB, &
                                 vibDB, bsetDB_orig, nac_data_all, &
                                 nac_min_state, nac_max_state, nac_data_is_dq)
    else
        ! == Extract density-related data from TD file
        dfile = DataFile(opts%file_td)
        call run%check(dfile%error, 'Error found while initializing data file')

        call sec_header(1, 'Data on Molecular System')

        molDB = dfile%get_mol_data()
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

        ! Let us now correct the orientation
        if (opts%eckart) then
           call Eckart_orient(moldb, .true.)
           write(iu_out, '(a)') 'The molecular system is now set in Eckart orientation'
        endif

        ! == Now extract vibrationa-specific data
        dfile = DataFile(opts%file_frq)
        call run%check(dfile%error, &
            'Unable to initialize file containing vibrational data')
        vibDB = dfile%get_vib_data()
        call run%check(dfile%error, &
            'Parsing vibrational data in file failed')
        if (.not.vibDB%loaded) then
            call run%error%raise_error('data', 'missing', &
                'Failure to extract vibrational data')
        end if
    end if

    ! Let us now correct the final state if not set
    if (opts%max_state == -1) opts%max_state = excDB%n_states
    if (opts%read_moldata) then
        if (.not.opts%min_state_set .and. nac_min_state > 0) then
            opts%min_state = nac_min_state
        end if
        if (.not.opts%max_state_set .and. nac_max_state > 0) then
            opts%max_state = min(opts%max_state, nac_max_state)
        end if
    end if

    if (opts%read_moldata) then
        call sec_header(1, 'Data on Molecular System')
        if (opts%eckart) then
            call Eckart_orient(moldb, .true.)
            write(iu_out, '(a)') 'The molecular system is now set in Eckart orientation'
        end if
    end if

    if (opts%align_edtm) then
        if (opts%read_cubes) then
            call run%error%raise_error('opt', 'conflict', &
                'The -t/--align-edtm option cannot be used together with ' // &
                '-c/--cubefname because precomputed cubes cannot be rotated ' // &
                'consistently.')
        end if
        call build_alignment_rotation(opts%edtm_ref, opts%edtm_axis, edtm_rotmat)
        call rotate_molecule(molDB, edtm_rotmat)
        call rotate_vibrations(vibDB, edtm_rotmat)
        if (allocated(nac_data_all) .and. .not.nac_data_is_dq) then
            call rotate_3n_matrix(nac_data_all, edtm_rotmat)
        end if
        call write_param('EDTM alignment axis', axis_label(opts%edtm_axis))
        call write_param('EDTM alignment vector', fmt_real3(opts%edtm_ref))
    end if

    ! == Loading finished
    call write_moldata(molDB, bsetDB_orig, orbDB)

    call sec_header(1, 'Transition Data')
    write(iu_out, 1200) excDB%n_states, excDB%id_state

    if (TIMEIT) call write_time('Check AO normalization')
    call fix_norm_AOs(iu_out, molDB%n_at, orbDB%n_ao, molDB%at_crd, &
                      bsetDB_orig, err, DEBUG)
    call run%check(err)

    if (.not.opts%read_moldata) then
        ! Extract list of NACs files
        call list_nac_files(opts%template_nac, opts%min_state, opts%max_state, &
                            nac_ftype, nac_fnames)
    end if

    if (opts%read_moldata) then
        if (.not.allocated(nac_data_all)) then
            call run%error%raise_error('data', 'missing', &
                'Missing NAC data in moldata input file.')
        end if
        if (opts%min_state < nac_min_state .or. &
            opts%max_state > nac_max_state) then
            call run%error%raise_error('data', 'missing', &
                'Moldata input NAC range does not cover requested states.')
        end if
    else if (opts%write_moldata) then
        allocate(nac_data_all(3*molDB%n_at, opts%max_state - opts%min_state + 1))
        nac_data_all = f0
    end if

    call sec_header(1, 'Cube Generation')

    ! Initialize cubes
    if (opts%length_gauge) then
        call cube_tcd%init(molDB, 'ETDD')
        call cube_vtcd%init(molDB, 'VTDD', &
                            [opts%min_state, opts%max_state, opts%mode])
    else
        call cube_tcd%init(molDB, 'ETCD')
        call cube_vtcd%init(molDB, 'VTCD', &
                            [opts%min_state, opts%max_state, opts%mode])
    end if

    dtm_1der = f0

    if (opts%read_cubes) then
        ! Read Cube data stored in data file.
        tcdDB = FileTCDCube(opts%file_cubes)
        call run%check(tcdDB%error, 'Unable to initialized TCD data cube file')
        if (opts%length_gauge) then
            write(iu_out, '(a)') &
                'Interpreting cube data file as length-gauge TDD data.'
        end if
        cube_tcd%grid%min = tcdDB%grid_min
        cube_tcd%grid%max = tcdDB%grid_max
        cube_tcd%grid%step_size = tcdDB%grid_step
        cube_tcd%grid%n_points = tcdDB%n_points
        cube_vtcd%grid = cube_tcd%grid
        cube_tcd%vector_size = tcdDB%len_point
        if (allocated(cube_tcd%data)) deallocate(cube_tcd%data)
        allocate(cube_tcd%data(cube_tcd%vector_size, cube_tcd%grid%n_points(3), &
                               cube_tcd%grid%n_points(2), &
                               cube_tcd%grid%n_points(1)))
        cube_tcd%data = f0
        cube_vtcd%vector_size = cube_tcd%vector_size
        allocate(cube_vtcd%data(3,cube_vtcd%grid%n_points(3), &
                                cube_vtcd%grid%n_points(2), &
                                cube_vtcd%grid%n_points(1)))
        cube_vtcd%data = f0

        if (opts%read_moldata) then
            allocate(nac_dQ(vibDB%n_vib))
        else
            if (nac_ftype == 'fchk') allocate(rvec(3*molDB%n_at))
            allocate(nac_vec(3*molDB%n_at))
            allocate(nac_dQ(vibDB%n_vib))
        end if
        do istate = opts%min_state, opts%max_state
            call tcdDB%read(istate, cube_tcd%data)
            write(msg, '("Unable to read TCD cube data for state ",i0)') &
                istate
            call run%check(tcdDB%error, msg)
            if (opts%read_moldata) then
                if (nac_data_is_dq) then
                    nac_dQ = nac_data_all(:, istate-nac_min_state+1)
                else
                    nac_dQ = matmul(nac_data_all(:, istate-nac_min_state+1), vibDB%L_mwg)
                end if
            else
                if (nac_ftype == 'fchk') then
                    nac_file = DataFile(nac_fnames(istate-opts%min_state+1))
                    nacDB = nac_file%get_data(50, 0, -1)
                    nac_vec = nacDB%data
                else
                    open(newunit=iu, file=nac_fnames(istate-opts%min_state+1), &
                         status='old', action='read')
                    if (ios /= 0) then
                        write(msg, '("Unable to open file ",a,".")') &
                            trim(nac_fnames(istate-opts%min_state+1))
                        call run%error%raise_error('file', 'open', msg)
                    end if
                    read(iu, *) rvec(:3*molDB%n_at)
                    nac_vec = rvec
                end if
                if (opts%align_edtm) call rotate_3n_vector(nac_vec, edtm_rotmat)
                nac_dQ = matmul(nac_vec, vibDB%L_mwg)
                if (opts%write_moldata) then
                    nac_data_all(:, istate-opts%min_state+1) = nac_vec
                end if
            end if

            dtm_sum = [sum(cube_tcd%data(1,:,:,:)), &
                       sum(cube_tcd%data(2,:,:,:)), &
                       sum(cube_tcd%data(3,:,:,:))]
            dtm_sum = dtm_sum * &
                product((cube_tcd%grid%max + cube_tcd%grid%step_size) - cube_tcd%grid%min) &
                / product(cube_tcd%grid%n_points)

            if (opts%length_gauge .or. opts%get_axial) then
                cube_vtcd%data = cube_vtcd%data + cube_tcd%data * nac_dQ(opts%mode)
                dtm_1der = dtm_1der + dtm_sum * nac_dQ(opts%mode)
            else
                cube_vtcd%data = cube_vtcd%data + cube_tcd%data * nac_dQ(opts%mode) &
                               / excDB%g2e_energy(istate)
                dtm_1der = dtm_1der + dtm_sum * nac_dQ(opts%mode) &
                         / excDB%g2e_energy(istate)
            end if

        end do

        call tcdDB%close()
    else
        ! Construct cube
        if (opts%read_moldata) then
            call run%error%raise_error('data', 'missing', &
                'Moldata input does not include transition densities; &
                &use -c/--cubefname to read cubes.')
        end if
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
        ! Build grid specifications
        if (opts%which_grid == -1) then
            call cube_tcd%grid%read(opts%file_grid, ignore_unknown=.true.)
            call run%check(cube_tcd%grid%error, &
                'Error found while parsing grid input file')
        else if (opts%which_grid > 0) then
            if (opts%which_grid < opts%min_state &
                    .or. opts%which_grid > opts%max_state) &
                call run%error%raise_error('opt', 'value', &
                    'State chosen as reference for grid construction outside &
                    &range of states considered for VTCD')
            istate = opts%which_grid
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
            call cube_tcd%grid%init(cube_tcd%type, opts%grid_sparsity, &
                                    molDB=molDB, bsetDB=bsetDB, &
                                    e_trans_dens=trans_dens)
        else
            if (orbDB%openshell) then
                ! sum alpha and beta densities, 2.0 comes from sqrt(2.0)
                tmp_dens = (excDB%g2e_dens(:,:,1,opts%min_state) &
                           + excDB%g2e_dens(:,:,2,opts%min_state)) / f2
            else
                tmp_dens = excDB%g2e_dens(:,:,1,opts%min_state)
            end if
            if (convert) then
                call convert_pure2cart(bsetDB_orig, tmp_dens, trans_dens)
            else
                trans_dens = tmp_dens
            end if
            if (opts%align_edtm) call rotate_ao_density(ao_rot, trans_dens)
            call antisymm_sum(bsetDB%n_basis, trans_dens)
            call cube_tcd%grid%init(cube_tcd%type, opts%grid_sparsity, &
                                    molDB=molDB, bsetDB=bsetDB, &
                                    e_trans_dens=trans_dens)
            allocate(gridDB)
            do istate = opts%min_state+1, opts%max_state
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
                call gridDB%init(cube_tcd%type, opts%grid_sparsity, &
                                 molDB=molDB, bsetDB=bsetDB, &
                                 e_trans_dens=trans_dens)
                if (gridDB%min(1) < cube_tcd%grid%min(1)) &
                    cube_tcd%grid%min(1) = gridDB%min(1)
                if (gridDB%min(2) < cube_tcd%grid%min(2)) &
                    cube_tcd%grid%min(2) = gridDB%min(2)
                if (gridDB%min(3) < cube_tcd%grid%min(3)) &
                    cube_tcd%grid%min(3) = gridDB%min(3)
                if (gridDB%max(1) > cube_tcd%grid%max(1)) &
                    cube_tcd%grid%max(1) = gridDB%max(1)
                if (gridDB%max(2) > cube_tcd%grid%max(2)) &
                    cube_tcd%grid%max(2) = gridDB%max(2)
                if (gridDB%max(3) > cube_tcd%grid%max(3)) &
                    cube_tcd%grid%max(3) = gridDB%max(3)
            end do
            ! Now correct steps based on the grid specifications
            cube_tcd%grid%step_size = (cube_tcd%grid%max - cube_tcd%grid%min) &
                / real(cube_tcd%grid%n_points-1, realwp)
        end if
        call write_griddata(iu_out, cube_tcd%grid)
        cube_vtcd%grid = cube_tcd%grid
        allocate(cube_vtcd%data(3,cube_vtcd%grid%n_points(3), &
                                cube_vtcd%grid%n_points(2), &
                                cube_vtcd%grid%n_points(1)))
        cube_vtcd%data = f0
        ! Now constructs cubes
        if (opts%save_cubes) then
            tcdDB = FileTCDCube(cube_vtcd%grid, 3, min_state=opts%min_state)
            call tcdDB%init_file(opts%file_cubes)
            call run%check(tcdDB%error, &
                'Unable to initialize TCD cube data file')
        end if
        if (nac_ftype == 'fchk') allocate(rvec(3*molDB%n_at))
        allocate(nac_vec(3*molDB%n_at))
        allocate(nac_dQ(vibDB%n_vib))
        do istate = opts%min_state, opts%max_state
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
            call cube_tcd%gendata(bsetDB, trans_dens, dtm_sum, &
                                  gen_axial=opts%get_axial, &
                                  length_gauge=opts%length_gauge)
            if (opts%save_cubes) then
                call tcdDB%write(istate, cube_tcd%data)
            end if
            if (nac_ftype == 'fchk') then
                nac_file = DataFile(nac_fnames(istate-opts%min_state+1))
                nacDB = nac_file%get_data(50, 0, -1)
                nac_vec = nacDB%data
            else
                open(newunit=iu, file=nac_fnames(istate-opts%min_state+1), &
                     status='old', action='read')
                if (ios /= 0) then
                    write(msg, '("Unable to open file ",a,".")') &
                        trim(nac_fnames(istate-opts%min_state+1))
                    call run%error%raise_error('file', 'open', msg)
                end if
                read(iu, *) rvec(:3*molDB%n_at)
                nac_vec = rvec
            end if
            if (opts%align_edtm) call rotate_3n_vector(nac_vec, edtm_rotmat)
            nac_dQ = matmul(nac_vec, vibDB%L_mwg)
            if (opts%write_moldata) then
                nac_data_all(:, istate-opts%min_state+1) = nac_vec
            end if

            if (opts%length_gauge .or. opts%get_axial) then
                cube_vtcd%data = cube_vtcd%data + cube_tcd%data * nac_dQ(opts%mode)
                dtm_1der = dtm_1der + dtm_sum * nac_dQ(opts%mode)
            else
                cube_vtcd%data = cube_vtcd%data + cube_tcd%data * nac_dQ(opts%mode) &
                               / excDB%g2e_energy(istate)
                dtm_1der = dtm_1der + dtm_sum * nac_dQ(opts%mode) &
                         / excDB%g2e_energy(istate)
            end if

        end do
        if (opts%save_cubes) call tcdDB%close()
    end if

    if (opts%write_moldata) then
        call write_moldata_binary(opts%moldata_file, excDB, molDB, orbDB, &
                                  vibDB, bsetDB_orig, opts%min_state, &
                                  opts%max_state, nac_data_all)
    end if

    ! Construct base name for cube filenames
    if (opts%read_moldata) then
        i = index(opts%moldata_in_file, '.', back=.true.)
        if (i > 1) then
            base_name = opts%moldata_in_file(:i-1)
        else
            base_name = opts%moldata_in_file
        end if
    else
        i = index(opts%file_td, '.', back=.true.)
        if (i > 1) then
            base_name = opts%file_td(:i-1)
        else
            base_name = opts%file_td
        end if
    end if
    ! find length to print highest state
    lmode = num_chars_int(vibDB%n_vib)
    lstate = num_chars_int(opts%max_state)
    ! Format: <base>_v<lmode>.<lmode>_s<lstate>.<lstate>-<lstate>.<lstate>.cube
    allocate(character(len=len(base_name)+26+2*lmode+4*lstate) :: cube_fname)
    if (opts%length_gauge) then
        write(fmt, &
              '(''(a,"_VTDD_v",i'',i0,''.'',i0,'',"_s",i'',i0,''.'',i0,'',"-",i'',i0,&
              &''.'',i0,'',".cube")'')') &
            lmode, lmode, lstate, lstate, lstate, lstate
    else
        write(fmt, &
              '(''(a,"_v",i'',i0,''.'',i0,'',"_s",i'',i0,''.'',i0,'',"-",i'',i0,&
              &''.'',i0,'',".cube")'')') &
            lmode, lmode, lstate, lstate, lstate, lstate
    end if

    write(cube_fname, fmt) base_name, opts%mode, opts%min_state, opts%max_state

    open(newunit=iu, file=cube_fname, status='replace', action='write', &
         iostat=ios)
    if (ios /= 0) &
        call run%error%raise_error('file', 'open', &
            'Failed to open cube file for writing')

    cube_integral = [sum(cube_vtcd%data(1,:,:,:)), &
                     sum(cube_vtcd%data(2,:,:,:)), &
                     sum(cube_vtcd%data(3,:,:,:))]
    cube_integral = cube_integral * &
        product((cube_vtcd%grid%max + cube_vtcd%grid%step_size) - cube_vtcd%grid%min) &
        / product(cube_vtcd%grid%n_points)

    header_base = cube_vtcd%header
    write(header_with_integral, '(a," | Integral: ",3es20.12)') &
        trim(header_base), cube_integral
    cube_vtcd%header = header_with_integral
    call cube_vtcd%write(iu)
    cube_vtcd%header = header_base
    close(iu)

    write(iu_out, '("File ",a," has been generated!")') cube_fname

    if (allocated(tmp_dens)) deallocate(tmp_dens)
    if (allocated(trans_dens)) deallocate(trans_dens)
    if (opts%print_level == 'debug') then
        if (opts%get_axial) then

            write(iu_out, '(/A, I0, A, E14.6, E14.6, E14.6/)') &
                  'Electronic MDTM (au) for the mode ', opts%mode, ' is:', - dtm_1der &
                                                                     / sqrt(vibDB%red_mass(opts%mode)) &
                                                                     * sqrt(phys_conv%au2cm1(vibDB%freq(opts%mode),.true.)/f2)

        else if (opts%length_gauge) then

            write(iu_out, '(/A, I0, A, E14.6, E14.6, E14.6/)') &
                  'Length gauge electronic EDTM (au) for the mode ', opts%mode, ' is:', dtm_1der &
                                                                     / sqrt(vibDB%red_mass(opts%mode)) &
                                                                     * sqrt(phys_conv%au2cm1(vibDB%freq(opts%mode),.true.)/f2)

        else

            write(iu_out, '(/A, I0, A, E14.6, E14.6, E14.6/)') &
                  'Velocity gauge electronic EDTM (au) for the mode ', opts%mode, ' is:', dtm_1der &
                                                                     / sqrt(vibDB%red_mass(opts%mode)) &
                                                                     * sqrt(phys_conv%au2cm1(vibDB%freq(opts%mode),.true.)/f2)

        end if

    end if

    allocate(apt_nuc(3*molDB%n_at,3))
    allocate(aat_nuc(3*molDB%n_at,3))

    apt_nuc = f0
    aat_nuc = f0

    edtm_nuc_1der = f0
    mdtm_nuc_1der = f0

    do iat = 1, molDB%n_at

        offset = (iat-1)*3

        apt_nuc(offset+1,1) = molDB%at_num(iat)
        apt_nuc(offset+2,2) = molDB%at_num(iat)
        apt_nuc(offset+3,3) = molDB%at_num(iat)

        aat_nuc(offset+2,1) = -molDB%at_crd(3, iat) * bohr_radius
        aat_nuc(offset+3,1) =  molDB%at_crd(2, iat) * bohr_radius
        aat_nuc(offset+1,2) =  molDB%at_crd(3, iat) * bohr_radius
        aat_nuc(offset+3,2) = -molDB%at_crd(1, iat) * bohr_radius
        aat_nuc(offset+1,3) = -molDB%at_crd(2, iat) * bohr_radius
        aat_nuc(offset+2,3) =  molDB%at_crd(1, iat) * bohr_radius

        aat_nuc(offset+1:offset+3,:) = aat_nuc(offset+1:offset+3,:) * molDB%at_num(iat) * fhalf**2

    end do

    edtm_nuc_1der = matmul(transpose(apt_nuc), vibDB%L_mwg(:,opts%mode)) / sqrt(vibDB%red_mass(opts%mode))
    mdtm_nuc_1der = matmul(transpose(aat_nuc), vibDB%L_mwg(:,opts%mode)) / sqrt(vibDB%red_mass(opts%mode))

    if (opts%print_level == 'debug') then

        if (opts%get_axial) then

            write(iu_out, '(/A, I0, A, E14.6, E14.6, E14.6/)') &
                  'Nuclear MDTM (au) for the mode ', opts%mode, ' is:', mdtm_nuc_1der
        else

            write(iu_out, '(A, I0, A, E14.6, E14.6, E14.6)') &
                  'Nuclear EDTM (au) for the mode ', opts%mode, ' is:', edtm_nuc_1der

        end if

    end if

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

subroutine parse_options(optsDB)
    !! Parses command line options and updates information.
    !!
    !! Builds options parser and parse user-options, setting
    !! default values where necessary.
    type(BaseParams), intent(out) :: optsDB

    logical :: exists
    character(len=1024) :: argval
    character(len=:), dimension(:), allocatable :: argvals
    class(CmdLineArgsDB), allocatable :: parser

    ! Build options parser and check arguments
    parser = CmdLineArgsDB(progname='vtcd_cube')
    call run%check(parser%error, 'Failed to initialize the command-line parser')

    call parser%add_arg_char( &
        'list', label='data_files', &
        min_nvals=0, max_nvals=3, &
        help='TD_filename nac_template freq_filename:\n &
            &- Gaussian fchk file with electronic excitation data\n&
            &- Template for NAC files (use "#" for the placeholder of the &
            &indexes, matching the necessary length: 001 -> ###)\n&
            &- Gaussian fchk file storing ground-state vibrational data')

    call parser%add_arg_char('string', label='axial_tensor', &
                             shortname='-a', longname='--axial-tensor', &
                             required=.false., &
                             help='Generate axial tensor .cube file')

    call parser%add_arg_char('string', label='save_cubes_in_file', &
                             shortname='-b', longname='--savecube', &
                             required=.false., &
                             help='Save .cube in binary file')

    call parser%add_arg_char('string', label='cubes_from_file', &
                             shortname='-c', longname='--cubefname', &
                             required=.false., &
                             help='Read .cube from file')

    call parser%add_arg_char('string', label='grid_sparsity', &
                             shortname='-d', longname='--sparsity', &
                             required=.false., &
                             help='Set grid sparsity')

    call parser%add_arg_int('scalar', &
                            shortname='-e', longname='--max-state', &
                            min_value=1, &
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
                            help='Use grid parameters generated for the specified state')

    call parser%add_arg_char('string', label='print_level', &
                             shortname='-l', longname='--print-level', &
                             help='Set the print level', required=.false.)

    call parser%add_arg_bool('store_true', label='length_gauge', &
                             shortname='-L', longname='--length-gauge', &
                             required=.false., &
                             help='Use the length-gauge electric-dipole &
                             &representation. This generates transition &
                             &dipole-density maps instead of transition &
                             &current-density maps.')

    call parser%add_arg_char('string', label='moldata_output', &
                             shortname='-m', longname='--moldata', &
                             required=.false., &
                             help='Write moldata output file')

    call parser%add_arg_char('string', label='output', &
                             shortname='-o', longname='--output', &
                             help='Name of the file to store the output. &
                                 &Existing content will be overwritten.')

    call parser%add_arg_int('scalar', label='normal_mode', &
                            shortname='-q', longname='--normal-mode', &
                            min_value=1, required=.true., &
                            help='Normal mode to include')

    call parser%add_arg_char('string', label='eckart_orient', &
                             shortname='-r', longname='--eckart', &
                             required=.false., &
                             help='Rotate the molecule to Eckart orientation')

    call parser%add_arg_int('scalar', &
                            shortname='-s', longname='--min-state', &
                            min_value=1, &
                            help='Lowest electronic state to include (upper &
                                 &bound of summation)')

    call parser%add_arg_char('list', label='x y z axis', &
                             shortname='-t', longname='--align-edtm', &
                             required=.false., &
                             help='Rotate system so EDTM (x y z) aligns to &
                                  &axis int (1=x, 2=y, 3=z)')

    call parser%add_arg_char('string', label='moldata_input', &
                             shortname='-u', longname='--moldata-input', &
                             required=.false., &
                             help='Read moldata input file')

    call parser%parse_args()
    call run%check(parser%error, 'Parsing of command-line arguments failed')

    ! Output file
    if (parser%is_user_set('output')) then
        call parser%get_value('output', argval)
        optsDB%outfile = trim(argval)
        open(newunit=iu_out, file=optsDB%outfile, action='write')
        call sec_header(-1, PROGTITLE)
    end if

    call sec_header(1, 'Simulation Parameters')
    write(iu_out, '(1x)')

    if (parser%is_user_set('moldata_input')) then
        if (.not.parser%is_user_set('cubes_from_file')) then
            call run%error%raise_error('opt', 'missing', &
                '`-u` (moldata input) requires `-c` (TCD cube data input).')
        end if
        optsDB%read_moldata = .true.
        call parser%get_value('moldata_input', argval)
        optsDB%moldata_in_file = trim(argval)
        inquire(file=optsDB%moldata_in_file, exist=exists)
        if (.not.exists) then
            write(argval, '("Moldata file ''",a,"'' does not exist")') &
                optsDB%moldata_in_file
            call run%error%raise_error('file', 'not found', trim(argval))
        end if
        call write_param('Reading moldata from file', optsDB%moldata_in_file)
    else if (parser%is_user_set('data_files')) then
        call parser%get_value('data_files', argvals)
        if (size(argvals) /= 3) then
            call run%error%raise_error('opt', 'missing', &
                'Exactly 3 parameters expected: TD_filename nac_template &
                &freq_filename')
        end if
        optsDB%file_td = trim(argvals(1))
        inquire(file=optsDB%file_td, exist=exists)
        if (.not.exists) then
            write(argval, '("File with TD data ''",a,"'' does not exist")') &
                optsDB%file_td
            call run%error%raise_error('file', 'not found', trim(argval))
        end if
        call write_param('TD input filename', optsDB%file_td)
        optsDB%template_nac = trim(argvals(2))
        call write_param('NAC filenames template', optsDB%template_nac)
        optsDB%file_frq = trim(argvals(3))
        inquire(file=optsDB%file_frq, exist=exists)
        if (.not.exists) then
            write(argval, &
                  '("File with frequency data ''",a,"'' does not exist")') &
                optsDB%file_frq
            call run%error%raise_error('file', 'not found', trim(argval))
        end if
        call write_param('Frequencies input filename', optsDB%file_frq)
    else
        call run%error%raise_error('opt', 'missing', &
            'Missing data files')
    endif

    if (parser%is_user_set('moldata_output')) then
        if (.not.parser%is_user_set('save_cubes_in_file')) then
            call run%error%raise_error('opt', 'missing', &
                '`-m` (moldata output) requires `-b` (TCD cube data output).')
        end if
        optsDB%write_moldata = .true.
        call parser%get_value('moldata_output', argval)
        if (len_trim(argval) == 0) then
            write(*, '(a)') 'Missing filename for -m option'
            stop
        end if
        optsDB%moldata_file = trim(argval)
        call write_param('Writing moldata in file', optsDB%moldata_file)
    end if

    call parser%get_value('normal_mode', optsDB%mode)

    if (parser%is_user_set('grid_filename')) then
        call parser%get_value('grid_filename', argval)
        optsDB%file_grid = trim(argval)
        call write_param('Grid input filename', optsDB%file_grid)
        optsDB%which_grid = -1
    endif

    if (parser%is_user_set('grid_sparsity')) then
        if (optsDB%which_grid == -1) then
            optsDB%grid_sparsity = 'ignored'
        else
            call parser%get_value('grid_sparsity', optsDB%grid_sparsity)
        end if
        call write_param('Grid sparsity', optsDB%grid_sparsity)
    endif

    if (parser%is_user_set('grid_all') .and. &
        parser%is_user_set('grid_specific')) then
        call run%error%raise_error('opt', 'conflict', &
            'Options --grid_all and --grid_specific cannot be used together')
    end if

    if (parser%is_user_set('moldata_output') .and. &
        parser%is_user_set('moldata_input')) then
        call run%error%raise_error('opt', 'conflict', &
            '`-m` (write) and `-u` (read) cannot be used together.')
    end if

    if (parser%is_user_set('grid_all')) then
        if (optsDB%which_grid == -1) then
            call run%error%raise_error('opt', 'conflict', &
                'Conflicting options in definition of grid.')
        end if
        optsDB%which_grid = 0
        call write_param('Common grid for all states', .true.)
    end if

    if (parser%is_user_set('grid_specific')) then
        if (optsDB%which_grid == -1) then
            call run%error%raise_error('opt', 'conflict', &
                'Conflicting options in definition of grid.')
        end if
        call parser%get_value('grid_specific', optsDB%which_grid)
        if (optsDB%which_grid <= 0) &
            call run%error%raise_error('opt', 'value', &
                'Strictly positive value expected for `--grid_specific`')
        call write_param('Grid taken from state', optsDB%which_grid)
    end if

    if (parser%is_user_set('min-state')) then
        call parser%get_value('min-state', optsDB%min_state)
        optsDB%min_state_set = .true.
    end if

    if (parser%is_user_set('max-state')) then
        call parser%get_value('max-state', optsDB%max_state)
        optsDB%max_state_set = .true.
    end if

    if (parser%is_user_set('min-state') .and. &
        parser%is_user_set('max-state')) then
        if (optsDB%min_state == optsDB%max_state) then
            call write_param('Requested excited state', optsDB%min_state)
        else
            call write_param('Lowest excited state', optsDB%min_state)
            call write_param('Highest excited state', optsDB%max_state)
        endif
    else if (parser%is_user_set('max-state')) then
        optsDB%min_state = 1
        call write_param('Lowest excited state', 'default to first state')
        call write_param('Highest excited state', optsDB%max_state)
    else if (parser%is_user_set('min-state')) then
        call write_param('Lowest excited state', optsDB%min_state)
        call write_param('Highest excited state', 'include all')
    else
        optsDB%min_state = 1
        call write_param('Lowest excited state', 'default to first state')
        call write_param('Highest excited state', 'include all')
    endif

    if (parser%is_user_set('cubes_from_file')) then
        optsDB%read_cubes = .true.
        call parser%get_value('cubes_from_file', argval)
        optsDB%file_cubes = trim(argval)
        inquire(file=optsDB%file_cubes, exist=exists)
        if (.not.exists) then
            write(argval, &
                  '("File with cubes data ''",a,"'' does not exist")') &
                optsDB%file_cubes
            call run%error%raise_error('file', 'not found', trim(argval))
        end if
        call write_param('Reading cubes from file', optsDB%file_cubes)
    endif

    if (parser%is_user_set('eckart_orient')) then
        call parser%get_value('eckart_orient', argval)
        if (trim(argval) == 'true') then
            optsDB%eckart = .true.
            call write_param('Eckart orientation', optsDB%eckart)
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

    if (parser%is_user_set('save_cubes_in_file')) then
        if (parser%is_user_set('cubes_from_file')) then
            call run%error%raise_error('opt', 'conflict', &
                'Cannot save and read cube files at the same time')
        end if
        call parser%get_value('save_cubes_in_file', argval)
        optsDB%save_cubes = .true.
        if (trim(argval) == 'true') then
            optsDB%file_cubes = 'TCD_cubes.dat'
        else
            optsDB%file_cubes = trim(argval)
        end if
        call write_param('Saving cubes data in file', &
            optsDB%save_cubes)
    endif

    if (parser%is_user_set('print_level')) then
        call parser%get_value('print_level', argval)
        if (trim(argval) == 'debug') then
            optsDB%print_level = 'debug'
            call write_param('Print level', optsDB%print_level)
        end if
    endif

    ! if (optsDB%align_edtm) then
    !     if (.not.align_vec_set) then
    !         call run%error%raise_error('opt', 'missing', &
    !             'Missing EDTM alignment data. Use -t/--align-edtm x y z axis.')
    !     end if
    !     if (.not.align_axis_set) then
    !         call run%error%raise_error('opt', 'missing', &
    !             'Missing alignment axis. Use -t/--align-edtm x y z axis.')
    !     end if
    !     if (optsDB%edtm_axis < 1 .or. optsDB%edtm_axis > 3) then
    !         call run%error%raise_error('opt', 'missing', &
    !             'The alignment axis must be 1, 2 or 3.')
    !     end if
    !     if (sum(abs(optsDB%edtm_ref)) <= 10*small) then
    !         call run%error%raise_error('opt', 'value', &
    !             'The EDTM vector provided with -t/--align-edtm must be non-zero.')
    !     end if
    ! end if

end subroutine parse_options

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

subroutine rotate_vibrations(vib_db, rotmat)
    implicit none

    class(VibrationsDB), intent(inout) :: vib_db
    real(realwp), dimension(3,3), intent(in) :: rotmat

    call rotate_3n_matrix(vib_db%L_mwg, rotmat)
    call rotate_3n_matrix(vib_db%L_mat, rotmat)
end subroutine rotate_vibrations

! ======================================================================

subroutine rotate_3n_matrix(data, rotmat)
    implicit none

    real(realwp), dimension(:,:), intent(inout) :: data
    real(realwp), dimension(3,3), intent(in) :: rotmat

    integer :: irow

    do irow = 1, size(data, 1), 3
        data(irow:irow+2, :) = matmul(rotmat, data(irow:irow+2, :))
    end do
end subroutine rotate_3n_matrix

! ======================================================================

subroutine rotate_3n_vector(data, rotmat)
    implicit none

    real(realwp), dimension(:), intent(inout) :: data
    real(realwp), dimension(3,3), intent(in) :: rotmat

    integer :: irow

    do irow = 1, size(data), 3
        data(irow:irow+2) = matmul(rotmat, data(irow:irow+2))
    end do
end subroutine rotate_3n_vector

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
        call build_cart_rotation_block(shell_prims, 1, rotmat, shell_rot(2:4, 2:4))
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

    res = factorial_int(n) / (factorial_int(nx) * factorial_int(ny) * factorial_int(nz))
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
    implicit none

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

    write(iu_out, 1100) mol_db%n_at, mol_db%charge, mol_db%multip, tag_open, &
        bset_db%n_basis, tag_bfunD, tag_bfunF, orb_db%n_ao

    call sec_header(2, 'Atomic coordinates')
    call prt_coord(mol_db%n_at, mol_db%at_lab, mol_db%at_crd, mol_db%at_mas)

end subroutine write_moldata

! ======================================================================

subroutine write_moldata_binary(filename, exc_db, mol_db, orb_db, vib_db, &
                                bset_db, min_state, max_state, nac_data_all)
    implicit none

    character(len=*), intent(in) :: filename
    class(ExcitationDB), intent(in) :: exc_db
    class(MoleculeDB), intent(in) :: mol_db
    class(OrbitalsDB), intent(in) :: orb_db
    class(VibrationsDB), intent(in) :: vib_db
    class(BasisSetDB), intent(in) :: bset_db
    integer, intent(in) :: min_state
    integer, intent(in) :: max_state
    real(realwp), dimension(:,:), allocatable, intent(in) :: nac_data_all

    integer :: iunit, ios
    character(len=8) :: magic
    integer(int32) :: version

    magic = 'MOLDATA'
    version = 2_int32

    open(newunit=iunit, file=filename, status='replace', action='write', &
         form='unformatted', access='stream', iostat=ios)
    if (ios /= 0) then
        write(msg, '("Unable to open moldata output file ''",a,"''")') &
            filename
        call run%error%raise_error('file', 'open', trim(msg))
    end if

    write(iunit) magic
    write(iunit) version

    write(iunit) int(min_state, int32)
    write(iunit) int(max_state, int32)

    call write_excitation_db(iunit, exc_db)
    call write_molecule_db(iunit, mol_db)
    call write_orbitals_db(iunit, orb_db)
    call write_vibrations_db(iunit, vib_db)
    call write_basisset_db(iunit, bset_db)
    call write_real_2d(iunit, nac_data_all)

    close(iunit)
end subroutine write_moldata_binary

! ======================================================================

subroutine read_moldata_binary(filename, exc_db, mol_db, orb_db, vib_db, &
                               bset_db, nac_data_all, min_state, max_state, &
                               nac_data_is_dq)
    implicit none

    character(len=*), intent(in) :: filename
    class(ExcitationDB), intent(out) :: exc_db
    class(MoleculeDB), intent(out) :: mol_db
    class(OrbitalsDB), intent(out) :: orb_db
    class(VibrationsDB), intent(out) :: vib_db
    class(BasisSetDB), intent(out) :: bset_db
    real(realwp), dimension(:,:), allocatable, intent(out) :: nac_data_all
    integer, intent(out) :: min_state
    integer, intent(out) :: max_state
    logical, intent(out) :: nac_data_is_dq

    integer :: iunit, ios
    character(len=8) :: magic
    integer(int32) :: version

    open(newunit=iunit, file=filename, status='old', action='read', &
         form='unformatted', access='stream', iostat=ios)
    if (ios /= 0) then
        write(msg, '("Unable to open moldata input file ''",a,"''")') &
            filename
        call run%error%raise_error('file', 'open', trim(msg))
    end if

    read(iunit) magic
    read(iunit) version
    if (trim(magic) /= 'MOLDATA') then
        call run%error%raise_error('key', 'not found', &
            'Invalid moldata file header')
    end if
    select case (version)
    case (1_int32)
        min_state = 1
        max_state = -1
        nac_data_is_dq = .true.
    case (2_int32)
        read(iunit) min_state
        read(iunit) max_state
        nac_data_is_dq = .false.
    case default
        call run%error%raise_error('file', 'version', &
            'Unsupported moldata file version.')
    end select

    call read_excitation_db(iunit, exc_db)
    call read_molecule_db(iunit, mol_db)
    call read_orbitals_db(iunit, orb_db)
    call read_vibrations_db(iunit, vib_db)
    call read_basisset_db(iunit, bset_db)
    call read_real_2d(iunit, nac_data_all)

    if (version == 1_int32) then
        max_state = min_state + size(nac_data_all, 2) - 1
    end if

    close(iunit)
end subroutine read_moldata_binary

! ======================================================================

subroutine write_excitation_db(iunit, exc_db)
    implicit none

    integer, intent(in) :: iunit
    class(ExcitationDB), intent(in) :: exc_db
    integer :: nvals, n1, n2, n3, n4

    call write_int32_val(iunit, exc_db%n_states)
    call write_int32_val(iunit, exc_db%id_state)
    write(iunit) exc_db%gs_energy

    if (allocated(exc_db%ispin_exc)) then
        nvals = size(exc_db%ispin_exc)
        call write_int32_val(iunit, nvals)
        call write_int32_array(iunit, exc_db%ispin_exc)
    else
        call write_int32_val(iunit, 0)
    end if

    if (allocated(exc_db%exc_dens)) then
        n1 = size(exc_db%exc_dens, 1)
        n2 = size(exc_db%exc_dens, 2)
        n3 = size(exc_db%exc_dens, 3)
        n4 = size(exc_db%exc_dens, 4)
        call write_int32_val(iunit, n1)
        call write_int32_val(iunit, n2)
        call write_int32_val(iunit, n3)
        call write_int32_val(iunit, n4)
        write(iunit) exc_db%exc_dens
    else
        call write_int32_val(iunit, 0)
        call write_int32_val(iunit, 0)
        call write_int32_val(iunit, 0)
        call write_int32_val(iunit, 0)
    end if

    call write_real_2d(iunit, exc_db%g2e_eldip)
    call write_real_2d(iunit, exc_db%g2e_magdip)
    call write_real_1d(iunit, exc_db%exc_energy)
    call write_real_1d(iunit, exc_db%g2e_energy)
    call write_logical_val(iunit, exc_db%dens_loaded)
    call write_logical_val(iunit, exc_db%prop_loaded)
end subroutine write_excitation_db

! ======================================================================

subroutine read_excitation_db(iunit, exc_db)
    implicit none

    integer, intent(in) :: iunit
    class(ExcitationDB), intent(out) :: exc_db
    integer :: nvals, n1, n2, n3, n4

    call read_int32_val(iunit, exc_db%n_states)
    call read_int32_val(iunit, exc_db%id_state)
    read(iunit) exc_db%gs_energy

    call read_int32_val(iunit, nvals)
    if (nvals > 0) then
        allocate(exc_db%ispin_exc(nvals))
        call read_int32_array(iunit, exc_db%ispin_exc)
    end if

    call read_int32_val(iunit, n1)
    call read_int32_val(iunit, n2)
    call read_int32_val(iunit, n3)
    call read_int32_val(iunit, n4)
    if (n1 > 0 .and. n2 > 0 .and. n3 > 0 .and. n4 > 0) then
        allocate(exc_db%exc_dens(n1, n2, n3, n4))
        read(iunit) exc_db%exc_dens
    end if

    call read_real_2d(iunit, exc_db%g2e_eldip)
    call read_real_2d(iunit, exc_db%g2e_magdip)
    call read_real_1d(iunit, exc_db%exc_energy)
    call read_real_1d(iunit, exc_db%g2e_energy)
    call read_logical_val(iunit, exc_db%dens_loaded)
    call read_logical_val(iunit, exc_db%prop_loaded)
end subroutine read_excitation_db

! ======================================================================

subroutine write_molecule_db(iunit, mol_db)
    implicit none

    integer, intent(in) :: iunit
    class(MoleculeDB), intent(in) :: mol_db
    integer :: nvals, n1, n2, n3

    call write_int32_val(iunit, mol_db%charge)
    call write_int32_val(iunit, mol_db%multip)
    call write_int32_val(iunit, mol_db%n_at)
    call write_int32_val(iunit, mol_db%n_el)
    write(iunit) mol_db%energy

    if (allocated(mol_db%at_num)) then
        nvals = size(mol_db%at_num)
        call write_int32_val(iunit, nvals)
        call write_int32_array(iunit, mol_db%at_num)
    else
        call write_int32_val(iunit, 0)
    end if

    call write_real_1d(iunit, mol_db%at_chg)
    call write_real_1d(iunit, mol_db%at_mas)
    call write_real_2d(iunit, mol_db%at_crd)

    if (allocated(mol_db%el_dens)) then
        n1 = size(mol_db%el_dens, 1)
        n2 = size(mol_db%el_dens, 2)
        n3 = size(mol_db%el_dens, 3)
        call write_int32_val(iunit, n1)
        call write_int32_val(iunit, n2)
        call write_int32_val(iunit, n3)
        write(iunit) mol_db%el_dens
    else
        call write_int32_val(iunit, 0)
        call write_int32_val(iunit, 0)
        call write_int32_val(iunit, 0)
    end if

    call write_char_array(iunit, mol_db%at_lab)
    call write_logical_val(iunit, mol_db%loaded)
    call write_logical_val(iunit, mol_db%dens_loaded)
end subroutine write_molecule_db

! ======================================================================

subroutine read_molecule_db(iunit, mol_db)
    implicit none

    integer, intent(in) :: iunit
    class(MoleculeDB), intent(out) :: mol_db
    integer :: nvals, n1, n2, n3

    call read_int32_val(iunit, mol_db%charge)
    call read_int32_val(iunit, mol_db%multip)
    call read_int32_val(iunit, mol_db%n_at)
    call read_int32_val(iunit, mol_db%n_el)
    read(iunit) mol_db%energy

    call read_int32_val(iunit, nvals)
    if (nvals > 0) then
        allocate(mol_db%at_num(nvals))
        call read_int32_array(iunit, mol_db%at_num)
    end if

    call read_real_1d(iunit, mol_db%at_chg)
    call read_real_1d(iunit, mol_db%at_mas)
    call read_real_2d(iunit, mol_db%at_crd)

    call read_int32_val(iunit, n1)
    call read_int32_val(iunit, n2)
    call read_int32_val(iunit, n3)
    if (n1 > 0 .and. n2 > 0 .and. n3 > 0) then
        allocate(mol_db%el_dens(n1, n2, n3))
        read(iunit) mol_db%el_dens
    end if

    call read_char_array(iunit, mol_db%at_lab)
    call read_logical_val(iunit, mol_db%loaded)
    call read_logical_val(iunit, mol_db%dens_loaded)
end subroutine read_molecule_db

! ======================================================================

subroutine write_orbitals_db(iunit, orb_db)
    implicit none

    integer, intent(in) :: iunit
    class(OrbitalsDB), intent(in) :: orb_db

    call write_int32_val(iunit, orb_db%n_ab)
    call write_int32_val(iunit, orb_db%n_ao)
    call write_int32_val(iunit, orb_db%n_mo)
    call write_int32_val(iunit, orb_db%n_ao_cart)
    call write_int32_array(iunit, orb_db%n_mos)
    call write_int32_array(iunit, orb_db%n_els)
    call write_real_2d(iunit, orb_db%en_mos)
    call write_real_3d(iunit, orb_db%coef_mos)
    call write_logical_val(iunit, orb_db%openshell)
    call write_logical_val(iunit, orb_db%loaded)
end subroutine write_orbitals_db

! ======================================================================

subroutine read_orbitals_db(iunit, orb_db)
    implicit none

    integer, intent(in) :: iunit
    class(OrbitalsDB), intent(out) :: orb_db

    call read_int32_val(iunit, orb_db%n_ab)
    call read_int32_val(iunit, orb_db%n_ao)
    call read_int32_val(iunit, orb_db%n_mo)
    call read_int32_val(iunit, orb_db%n_ao_cart)
    call read_int32_array(iunit, orb_db%n_mos)
    call read_int32_array(iunit, orb_db%n_els)
    call read_real_2d(iunit, orb_db%en_mos)
    call read_real_3d(iunit, orb_db%coef_mos)
    call read_logical_val(iunit, orb_db%openshell)
    call read_logical_val(iunit, orb_db%loaded)
end subroutine read_orbitals_db

! ======================================================================

subroutine write_vibrations_db(iunit, vib_db)
    implicit none

    integer, intent(in) :: iunit
    class(VibrationsDB), intent(in) :: vib_db

    call write_int32_val(iunit, vib_db%n_vib)
    call write_real_1d(iunit, vib_db%freq)
    call write_real_1d(iunit, vib_db%red_freq)
    call write_real_1d(iunit, vib_db%red_mass)
    call write_real_2d(iunit, vib_db%L_mwg)
    call write_real_2d(iunit, vib_db%L_mat)
    call write_logical_val(iunit, vib_db%loaded)
end subroutine write_vibrations_db

! ======================================================================

subroutine read_vibrations_db(iunit, vib_db)
    implicit none

    integer, intent(in) :: iunit
    class(VibrationsDB), intent(out) :: vib_db

    call read_int32_val(iunit, vib_db%n_vib)
    call read_real_1d(iunit, vib_db%freq)
    call read_real_1d(iunit, vib_db%red_freq)
    call read_real_1d(iunit, vib_db%red_mass)
    call read_real_2d(iunit, vib_db%L_mwg)
    call read_real_2d(iunit, vib_db%L_mat)
    call read_logical_val(iunit, vib_db%loaded)
end subroutine read_vibrations_db

! ======================================================================

subroutine write_basisset_db(iunit, bset_db)
    implicit none

    integer, intent(in) :: iunit
    class(BasisSetDB), intent(in) :: bset_db
    integer :: n1, n2, i, j, nvals

    call write_int32_val(iunit, bset_db%n_basis)
    call write_int32_val(iunit, bset_db%n_basok)
    call write_int32_val(iunit, bset_db%n_shells)
    call write_int32_val(iunit, bset_db%L_max)

    if (allocated(bset_db%nprim_per_at)) then
        nvals = size(bset_db%nprim_per_at)
        call write_int32_val(iunit, nvals)
        call write_int32_array(iunit, bset_db%nprim_per_at)
    else
        call write_int32_val(iunit, 0)
    end if

    call write_logical_val(iunit, bset_db%pureD)
    call write_logical_val(iunit, bset_db%pureF)

    if (allocated(bset_db%info)) then
        n1 = size(bset_db%info, 1)
        n2 = size(bset_db%info, 2)
        call write_int32_val(iunit, n1)
        call write_int32_val(iunit, n2)
        do j = 1, n2
            do i = 1, n1
                call write_primitive(iunit, bset_db%info(i, j))
            end do
        end do
    else
        call write_int32_val(iunit, 0)
        call write_int32_val(iunit, 0)
    end if

    call write_logical_val(iunit, bset_db%loaded)
end subroutine write_basisset_db

! ======================================================================

subroutine read_basisset_db(iunit, bset_db)
    implicit none

    integer, intent(in) :: iunit
    class(BasisSetDB), intent(out) :: bset_db
    integer :: n1, n2, i, j, nvals

    call read_int32_val(iunit, bset_db%n_basis)
    call read_int32_val(iunit, bset_db%n_basok)
    call read_int32_val(iunit, bset_db%n_shells)
    call read_int32_val(iunit, bset_db%L_max)

    call read_int32_val(iunit, nvals)
    if (nvals > 0) then
        allocate(bset_db%nprim_per_at(nvals))
        call read_int32_array(iunit, bset_db%nprim_per_at)
    end if

    call read_logical_val(iunit, bset_db%pureD)
    call read_logical_val(iunit, bset_db%pureF)

    call read_int32_val(iunit, n1)
    call read_int32_val(iunit, n2)
    if (n1 > 0 .and. n2 > 0) then
        allocate(bset_db%info(n1, n2))
        do j = 1, n2
            do i = 1, n1
                call read_primitive(iunit, bset_db%info(i, j))
            end do
        end do
    end if

    call read_logical_val(iunit, bset_db%loaded)
end subroutine read_basisset_db

! ======================================================================

subroutine write_primitive(iunit, prim)
    implicit none

    integer, intent(in) :: iunit
    type(PrimitiveFunction), intent(in) :: prim
    integer :: nvals, n1, n2

    call write_fixed_char(iunit, prim%shelltype)
    call write_int32_val(iunit, prim%L)
    call write_int32_val(iunit, prim%shellid)
    call write_int32_val(iunit, prim%ndim)
    call write_logical_val(iunit, prim%pure)
    call write_logical_val(iunit, prim%shell_first)
    call write_logical_val(iunit, prim%shell_last)
    write(iunit) prim%alpha

    if (allocated(prim%coeff)) then
        nvals = size(prim%coeff)
        call write_int32_val(iunit, nvals)
        write(iunit) prim%coeff
    else
        call write_int32_val(iunit, 0)
    end if

    if (allocated(prim%lxyz)) then
        n1 = size(prim%lxyz, 1)
        n2 = size(prim%lxyz, 2)
        call write_int32_val(iunit, n1)
        call write_int32_val(iunit, n2)
        call write_int32_matrix(iunit, prim%lxyz)
    else
        call write_int32_val(iunit, 0)
        call write_int32_val(iunit, 0)
    end if
end subroutine write_primitive

! ======================================================================

subroutine read_primitive(iunit, prim)
    implicit none

    integer, intent(in) :: iunit
    type(PrimitiveFunction), intent(out) :: prim
    integer :: nvals, n1, n2

    call read_fixed_char(iunit, prim%shelltype)
    call read_int32_val(iunit, prim%L)
    call read_int32_val(iunit, prim%shellid)
    call read_int32_val(iunit, prim%ndim)
    call read_logical_val(iunit, prim%pure)
    call read_logical_val(iunit, prim%shell_first)
    call read_logical_val(iunit, prim%shell_last)
    read(iunit) prim%alpha

    call read_int32_val(iunit, nvals)
    if (nvals > 0) then
        allocate(prim%coeff(nvals))
        read(iunit) prim%coeff
    end if

    call read_int32_val(iunit, n1)
    call read_int32_val(iunit, n2)
    if (n1 > 0 .and. n2 > 0) then
        allocate(prim%lxyz(n1, n2))
        call read_int32_matrix(iunit, prim%lxyz)
    end if
end subroutine read_primitive

! ======================================================================

subroutine write_int32_val(iunit, val)
    implicit none

    integer, intent(in) :: iunit
    integer, intent(in) :: val
    integer(int32) :: tmp

    tmp = int(val, int32)
    write(iunit) tmp
end subroutine write_int32_val

! ======================================================================

subroutine read_int32_val(iunit, val)
    implicit none

    integer, intent(in) :: iunit
    integer, intent(out) :: val
    integer(int32) :: tmp

    read(iunit) tmp
    val = int(tmp)
end subroutine read_int32_val

! ======================================================================

subroutine write_logical_val(iunit, val)
    implicit none

    integer, intent(in) :: iunit
    logical, intent(in) :: val
    integer :: ival

    ival = 0
    if (val) ival = 1
    call write_int32_val(iunit, ival)
end subroutine write_logical_val

! ======================================================================

subroutine read_logical_val(iunit, val)
    implicit none

    integer, intent(in) :: iunit
    logical, intent(out) :: val
    integer :: ival

    call read_int32_val(iunit, ival)
    val = (ival /= 0)
end subroutine read_logical_val

! ======================================================================

subroutine write_int32_array(iunit, arr)
    implicit none

    integer, intent(in) :: iunit
    integer, dimension(:), intent(in) :: arr
    integer(int32), allocatable :: tmp(:)

    allocate(tmp(size(arr)))
    tmp = int(arr, int32)
    write(iunit) tmp
    deallocate(tmp)
end subroutine write_int32_array

! ======================================================================

subroutine read_int32_array(iunit, arr)
    implicit none

    integer, intent(in) :: iunit
    integer, dimension(:), intent(out) :: arr
    integer(int32), allocatable :: tmp(:)

    allocate(tmp(size(arr)))
    read(iunit) tmp
    arr = int(tmp)
    deallocate(tmp)
end subroutine read_int32_array

! ======================================================================

subroutine write_int32_matrix(iunit, arr)
    implicit none

    integer, intent(in) :: iunit
    integer, dimension(:,:), intent(in) :: arr
    integer(int32), allocatable :: tmp(:,:)

    allocate(tmp(size(arr, 1), size(arr, 2)))
    tmp = int(arr, int32)
    write(iunit) tmp
    deallocate(tmp)
end subroutine write_int32_matrix

! ======================================================================

subroutine read_int32_matrix(iunit, arr)
    implicit none

    integer, intent(in) :: iunit
    integer, dimension(:,:), intent(out) :: arr
    integer(int32), allocatable :: tmp(:,:)

    allocate(tmp(size(arr, 1), size(arr, 2)))
    read(iunit) tmp
    arr = int(tmp)
    deallocate(tmp)
end subroutine read_int32_matrix

! ======================================================================

subroutine write_real_1d(iunit, arr)
    implicit none

    integer, intent(in) :: iunit
    real(realwp), dimension(:), intent(in), allocatable :: arr
    integer :: nvals

    if (allocated(arr)) then
        nvals = size(arr)
        call write_int32_val(iunit, nvals)
        write(iunit) arr
    else
        call write_int32_val(iunit, 0)
    end if
end subroutine write_real_1d

! ======================================================================

subroutine read_real_1d(iunit, arr)
    implicit none

    integer, intent(in) :: iunit
    real(realwp), dimension(:), allocatable, intent(out) :: arr
    integer :: nvals

    call read_int32_val(iunit, nvals)
    if (nvals > 0) then
        allocate(arr(nvals))
        read(iunit) arr
    end if
end subroutine read_real_1d

! ======================================================================

subroutine write_real_2d(iunit, arr)
    implicit none

    integer, intent(in) :: iunit
    real(realwp), dimension(:,:), intent(in), allocatable :: arr
    integer :: n1, n2

    if (allocated(arr)) then
        n1 = size(arr, 1)
        n2 = size(arr, 2)
        call write_int32_val(iunit, n1)
        call write_int32_val(iunit, n2)
        write(iunit) arr
    else
        call write_int32_val(iunit, 0)
        call write_int32_val(iunit, 0)
    end if
end subroutine write_real_2d

! ======================================================================

subroutine read_real_2d(iunit, arr)
    implicit none

    integer, intent(in) :: iunit
    real(realwp), dimension(:,:), allocatable, intent(out) :: arr
    integer :: n1, n2

    call read_int32_val(iunit, n1)
    call read_int32_val(iunit, n2)
    if (n1 > 0 .and. n2 > 0) then
        allocate(arr(n1, n2))
        read(iunit) arr
    end if
end subroutine read_real_2d

! ======================================================================

subroutine write_real_3d(iunit, arr)
    implicit none

    integer, intent(in) :: iunit
    real(realwp), dimension(:,:,:), intent(in), allocatable :: arr
    integer :: n1, n2, n3

    if (allocated(arr)) then
        n1 = size(arr, 1)
        n2 = size(arr, 2)
        n3 = size(arr, 3)
        call write_int32_val(iunit, n1)
        call write_int32_val(iunit, n2)
        call write_int32_val(iunit, n3)
        write(iunit) arr
    else
        call write_int32_val(iunit, 0)
        call write_int32_val(iunit, 0)
        call write_int32_val(iunit, 0)
    end if
end subroutine write_real_3d

! ======================================================================

subroutine read_real_3d(iunit, arr)
    implicit none

    integer, intent(in) :: iunit
    real(realwp), dimension(:,:,:), allocatable, intent(out) :: arr
    integer :: n1, n2, n3

    call read_int32_val(iunit, n1)
    call read_int32_val(iunit, n2)
    call read_int32_val(iunit, n3)
    if (n1 > 0 .and. n2 > 0 .and. n3 > 0) then
        allocate(arr(n1, n2, n3))
        read(iunit) arr
    end if
end subroutine read_real_3d

! ======================================================================

subroutine write_fixed_char(iunit, val)
    implicit none

    integer, intent(in) :: iunit
    character(len=*), intent(in) :: val
    integer :: nlen

    nlen = len(val)
    call write_int32_val(iunit, nlen)
    write(iunit) val
end subroutine write_fixed_char

! ======================================================================

subroutine read_fixed_char(iunit, val)
    implicit none

    integer, intent(in) :: iunit
    character(len=*), intent(out) :: val
    integer :: nlen

    call read_int32_val(iunit, nlen)
    read(iunit) val
end subroutine read_fixed_char

! ======================================================================

subroutine write_char_array(iunit, arr)
    implicit none

    integer, intent(in) :: iunit
    character(len=*), dimension(:), allocatable, intent(in) :: arr
    integer :: nvals, nlen

    if (allocated(arr)) then
        nvals = size(arr)
        nlen = len(arr)
        call write_int32_val(iunit, nlen)
        call write_int32_val(iunit, nvals)
        write(iunit) arr
    else
        call write_int32_val(iunit, 0)
        call write_int32_val(iunit, 0)
    end if
end subroutine write_char_array

! ======================================================================

subroutine read_char_array(iunit, arr)
    implicit none

    integer, intent(in) :: iunit
    character(len=*), dimension(:), allocatable, intent(out) :: arr

    integer :: nvals, nlen
    character(len=:), dimension(:), allocatable :: tmp

    call read_int32_val(iunit, nlen)
    call read_int32_val(iunit, nvals)
    if (nvals > 0 .and. nlen > 0) then
        allocate(character(len=nlen) :: tmp(nvals))
        read(iunit) tmp
    end if
    arr = tmp
    if (nlen > len(arr(1))) &
        call run%error%raise_warning('arg', 'size', &
            'Extracted character strings may be truncated.')
end subroutine read_char_array

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

subroutine list_nac_files(template, num_start, num_end, ftype, fnames, &
                          ignore_missing)
    !! Search for NAC files.
    !!
    !! Looks for all files matching the input template between chosen
    !! indexes.
    !! The general filetype is also reported, between: text, fchk
    !! If `ignore_missing` is true, missing files are tolerated
    !! (default: no)
    character(len=*), intent(in) :: template
        !! Template
    integer, intent(in) :: num_start
        !! Starting index.  If `-1`, use lowest index from file list.
    integer, intent(in) :: num_end
        !! End index. If `-1`, use highest index from file list.
    character(len=*), intent(out) :: ftype
        !! Generic file type.
    character(len=:), dimension(:), allocatable, intent(out) :: fnames
        !! list of files to return.
    logical, intent(in), optional :: ignore_missing
        !! Ignore missing files, simply printing a warning.

    integer :: i, i0, i1, lnum, nmax
    logical :: exists, missing_ok
    character(len=512) :: fmt, msg
    character(len=:), allocatable :: fname

    if (len_trim(template) == 0) &
        call run%error%raise_deverror('gen', &
            'Missing template.  Cannot proceed')

    if (num_end < num_start) &
        call run%error%raise_argerror('gen', &
            'Inverted indexes in list of states to read: `num_end < num_start`')

    ! Look for extension
    i0 = index(template, '.', back=.true.)
    if (i0 > len(template)) then
        ftype = 'text'
    else
        select case (locase(template(i0+1:)))
        case ('fchk', 'fch', 'fck')
            ftype = 'fchk'
        case default
            ftype = 'text'
        end select
    end if

    i0 = index(template, '#')
    if (i0 > 0) then
        i1 = i0
        do while (template(i1:i1) == '#')
            i1 = i1 + 1
        end do
        i1 = i1 - 1
        if (index(template(i1+1:), '#') > 0) &
            call run%error%raise_error('val', 'wrong', &
                'Only one index block allowed in template name.')
    else
        call run%error%raise_error('val', 'wrong', &
            'Missing index block "#" in template name.')
    end if

    lnum = i1 - i0 + 1
    1000 format('("',a,'",i',i0,'.',i0,',"',a,'")')
    write(fmt, 1000) template(:i0-1), lnum, lnum, trim(template(i1+1:))

    if (present(ignore_missing)) then
        missing_ok = ignore_missing
    else
        missing_ok = .false.
    end if

    if (num_start == -1 .or. num_end == -1) then
        nmax = 10**lnum - 1
    end if
    allocate(character(len=len_trim(template)) :: fname)

    ! Check if num_start set or not
    if (num_start == -1) then
        i0 = 0
        exists = .false.
        do while (.not.exists)
            i0 = i0 + 1
            if (i0 > nmax) &
                call run%error%raise_error('file', 'not found', &
                    'Could not find any file fitting the given pattern')
            write(fname, fmt) i0
            inquire(file=fname, exist=exists)
        end do
    else
        i0 = num_start
    end if

    ! Check if num_end set or not
    if (num_end == -1) then
        i1 = nmax + 1
        exists = .false.
        do while (.not.exists)
            i1 = i1 - 1
            if (i1 == 0) &
                call run%error%raise_error('val', 'missing', &
                    'Could not find any file fitting the given pattern')
            write(fname, fmt) i1
            inquire(file=fname, exist=exists)
        end do
    else
        i1 = num_end
    end if

    allocate(character(len=len_trim(template)) :: fnames(i1-i0+1))

    nmax = 0
    do i = i0, i1
        write(fname, fmt) i
        inquire(file=fname, exist=exists)
        if (.not.exists) then
            if (missing_ok) then
                write(msg, '("File ",a," does not exist. Ignored.")') fname
                call run%error%raise_warning('file', 'not found', trim(msg))
            else
                write(msg, '("File ",a," does not exist. Aborting.")') fname
                call run%error%raise_error('file', 'not found', trim(msg))
            end if
            fnames(i-i0+1) = ' '
        else
            fnames(i-i0+1) = fname
            nmax = nmax + 1
        end if
    end do
    if (nmax == 0) &
        call run%error%raise_error('file', 'not found', &
            'Could not find any file fitting the given pattern')

end subroutine list_nac_files

! ======================================================================

end program vtcd_cube
