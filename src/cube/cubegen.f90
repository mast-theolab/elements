module cubegen
    !! This module provides the necessary functionality to generate cubes.
    !! It includes the necessary subroutines and functions
    !! to create and manipulate cubes in a Fortran program.
    use blas_drv, only: xdot, xgemv, xgemm
    use datatypes, only: BasisSetDB, MoleculeDB
    use math, only: cross
    use numeric, only: realwp, f0, f1, f2
    use orbital, only: eval_AOs_chi_at, eval_AOs_nabla_chi_at
    use run_env, only: CoreExecObject, run
    use string, only: locase

    implicit none

    character(len=*), parameter, private :: &
        cube_type_dens = 'density', &
        cube_type_trdens = 'transition density', &
        cube_type_etcd = 'ETCD', &
        cube_type_vtcd = 'VTCD', &
        cube_type_etdd = 'ETDD', &
        cube_type_vtdd = 'VTDD'

    type, public, extends(CoreExecObject) :: CubegridDB
        !! Cube grid specification database
        character(len=:), allocatable :: shape
            !! Shape of the grid: cubic.
        real(realwp), dimension(3) :: min
            !! Minimum bounds of the grid.
        real(realwp), dimension(3) :: max
            !! Maximum bounds of the grid.
        integer, dimension(3) :: n_points
            !! Number of points along each direction.
        real(realwp), dimension(3) :: step_size
            !! Step size along each direction.
    contains
        procedure, private :: init_gridDB_params_lab, init_gridDB_params_int
        generic :: init => init_gridDB_params_lab, init_gridDB_params_int
        procedure :: read => read_gridDB_params
    end type CubegridDB

    type, private, abstract, extends(CoreExecObject) :: CubeBaseDB
        !! Basis cube database specification
        integer :: itype
            !! Type of data contained in cube, as integer.
            !! Supported forms:
            !! 1: density
            !! 2: electronic transition density
            !! 3: electronic transition current density
            !! 4: vibrational transition current density
            !! 5: electronic transition dipole density
            !! 6: vibrational transition dipole density
        character(len=:), allocatable :: type
            !! Type of data contained in cube.
        character(len=256) :: header
            !! Header line, with basic parameters.
        character(len=:), allocatable :: extra_head
            !! Extra header line, below header line.
        type(CubegridDB) :: grid
            !! Grid specification.
        type(MoleculeDB), pointer :: mol => null()
            !! Database for the molecule associated to cube.
    contains
        procedure(write_cube_data), deferred :: write
        procedure(init_cube_dbase), deferred :: init
    end type CubeBaseDB

    abstract interface
        subroutine write_cube_data(cubeDB, iunit, molDB)
            use datatypes, only: MoleculeDB
            import :: CubeBaseDB
            class(CubeBaseDB), intent(inout) :: cubeDB
                !! Cube database.
            integer, intent(in) :: iunit
                !! Unit identifier to file opened in writing.
            class(MoleculeDB), intent(in), target, optional :: molDB
                !! Molecule database.
        end subroutine write_cube_data
        subroutine init_cube_dbase(cubeDB, molDB, cube_type, iparams, &
                                   rparams, extra_head)
            use datatypes, only : MoleculeDB
            use numeric, only: realwp
            import :: CubeBaseDB
            class(CubeBaseDB), intent(inout) :: cubeDB
                !! Cube database.
            class(MoleculeDB), intent(in), target, optional :: molDB
                !! Molecule database.
            character(len=*), intent(in), optional :: cube_type
                !! Type of cube.
            integer, dimension(:), intent(in), optional :: iparams
                !! Integer-type parameters for header, cube-dependent.
            real(realwp), dimension(:), intent(in), optional :: rparams
                !! Real-type parameters for header, cube-dependent.
            character(len=*), intent(in), optional :: extra_head
                !! Extra header, printed below generated header line.
        end subroutine init_cube_dbase
    end interface

    type, private, extends(CubeBaseDB) :: CubeScalarDB
        !! Cube database for scalar quantities.
        real(realwp), dimension(:,:,:), allocatable :: data
            !! Data stored in cube.
    contains
        procedure :: write => write_cube_scalar_data
        procedure :: init => init_cube_base_scalar
    end type CubeScalarDB

    type, private, extends(CubeBaseDB) :: CubeVectorDB
        !! Cube database for vectorial quantities.
        integer :: vector_size
            !! Size of the vectorial quantity.
        real(realwp), dimension(:,:,:,:), allocatable :: data
            !! Data stored in cube.
    contains
        procedure :: write => write_cube_vector_data
        procedure :: init => init_cube_base_vector
    end type CubeVectorDB

    type, public, extends(CubeScalarDB) :: CubeDensDB
    contains
        procedure :: gendata => gen_cube_data_density
    end type CubeDensDB

    type, public, extends(CubeVectorDB) :: CubeTCDDB
    contains
        procedure :: gendata => gen_cube_data_tcd
    end type CubeTCDDB

    type, public, extends(CoreExecObject) :: FileTCDCube
        integer :: version
            !! Version of the file object.
        integer :: len_header
            !! Number of records dedicated to header in file.
        integer :: min_state = 1
            !! Lowest state stored in file.
        integer :: max_state = 0
            !! Highest state stored in file.
        integer :: len_record
            !! Record length.
        integer, dimension(2), private :: len_data
            !! Length of integer (leading) and real (actual) data in record.
        character(len=:), allocatable :: file
            !! File name associated to object.
        integer :: len_point
            !! Data length at a grid point.
        integer, dimension(3) :: n_points
            !! Number of grid points in each direction.
        integer :: grid_shape_id
            !! Grid shape identifier:
            !! 1. cubic
        real(realwp), dimension(3) :: grid_min
            !! Grid minimum point coordinate.
        real(realwp), dimension(3) :: grid_max
            !! Grid maximum point coordinate.
        real(realwp), dimension(3) :: grid_step
            !! Grid step along each direction.
        integer, private :: iunit = 0
            !! Unit to connected file.
    contains
        procedure, private :: file_tcdcube_read_data_single, &
            file_tcdcube_read_data_single_scal, &
            file_tcdcube_read_data_single_vec, &
            file_tcdcube_read_data_multi, &
            file_tcdcube_read_data_multi_scal, &
            file_tcdcube_read_data_multi_vec, &
            file_tcdcube_write_data_single, &
            file_tcdcube_write_data_single_scal, &
            file_tcdcube_write_data_single_vec, &
            file_tcdcube_write_data_multi, &
            file_tcdcube_write_data_multi_scal, &
            file_tcdcube_write_data_multi_vec
        generic :: read => file_tcdcube_read_data_single, &
            file_tcdcube_read_data_single_scal, &
            file_tcdcube_read_data_single_vec, &
            file_tcdcube_read_data_multi, &
            file_tcdcube_read_data_multi_scal, &
            file_tcdcube_read_data_multi_vec
        generic :: write => file_tcdcube_write_data_single, &
            file_tcdcube_write_data_single_scal, &
            file_tcdcube_write_data_single_vec, &
            file_tcdcube_write_data_multi, &
            file_tcdcube_write_data_multi_scal, &
            file_tcdcube_write_data_multi_vec
        procedure :: init_file => file_tcdcube_init_file
        procedure :: close => file_tcdcube_close_file
    end type FileTCDCube

    interface FileTCDCube
        module procedure file_tcdcube_load, file_tcdcube_init
    end interface FileTCDCube

interface

! ----------------------------------------------------------------------

module subroutine init_gridDB_params_int(gridDB, type, sparsity, n_points, &
                                         molDB, bsetDB, e_dens, e_trans_dens)
    !! Initialize a CubeGridDB object, computing its parameters.
    !!
    !! Compute the parameters for a Cube grid.
    !! Several types of cubes are supported to define the grid
    !!
    !! * 'dens': electronic density
    !! * 'tcd_dens': electronic transition density for TCD calcs.
    !!
    !! @warning
    !! A Cartesian basis set is necessary.
    !! @endwarning
    !!
    !! @note "version"
    !! `Sparsity` is given as an integer.
    !! @endnote
    class(CubegridDB), intent(inout) :: gridDB
        !! Cube grid database.
    character(len=*), intent(in) :: type
        !! Type of cube for the grid definition.
    integer, intent(in) :: sparsity
        !! Grid sparsity, as integer
    integer, dimension(:), intent(in), optional :: n_points
        !! Number of grid points in each direction, only if `sparsity=-1`.
    class(MoleculeDB), intent(in), optional :: molDB
        !! Molecule database.
    class(BasisSetDB), intent(in), optional :: bsetDB
        !! Basis set database.
    real(realwp), dimension(:,:), intent(in), optional :: e_dens
        !! Electronic density.
    real(realwp), dimension(:,:), intent(in), optional :: e_trans_dens
        !! Electronic transition density.

end subroutine init_gridDB_params_int

! ----------------------------------------------------------------------

module subroutine init_gridDB_params_lab(gridDB, type, sparsity, n_points, &
                                         molDB, bsetDB, e_dens, e_trans_dens)
    !! Initialize a CubeGridDB object, computing its parameters.
    !!
    !! Compute the parameters for a Cube grid.
    !! Several types of cubes are supported to define the grid
    !!
    !! * 'dens': electronic density
    !! * 'tcd_dens': electronic transition density for TCD calcs.
    !!
    !! @warning
    !! A Cartesian basis set is necessary.
    !! @endwarning
    !!
    !! @note "version"
    !! `Sparsity` is given as a label.
    !! It acts as a wrapper to [init_gridDB_params_lab], which contains
    !! the actual code.
    !! @endnote
    class(CubegridDB), intent(inout) :: gridDB
        !! Cube grid database.
    character(len=*), intent(in) :: type
        !! Type of cube for the grid definition.
    character(len=*), intent(in) :: sparsity
        !! Grid sparsity (see [sparsity_lab2int] for supported keys).
    integer, dimension(:), intent(in), optional :: n_points
        !! Number of grid points in each direction, only if `sparsity=read`.
    class(MoleculeDB), intent(in) :: molDB
        !! Molecule database.
    class(BasisSetDB), intent(in) :: bsetDB
        !! Basis set database.
    real(realwp), dimension(:,:), intent(in), optional :: e_dens
        !! Electronic density.
    real(realwp), dimension(:,:), intent(in), optional :: e_trans_dens
        !! Electronic transition density.

end subroutine init_gridDB_params_lab

! ----------------------------------------------------------------------

module subroutine read_gridDB_params(gridDB, file_grid, ignore_unknown)
    !! Initialize a CubeGridDB object, by reading its parameters
    class(CubegridDB), intent(inout) :: gridDB
        !! Cube grid database.
    character(len=*), intent(in) :: file_grid
        !! Grid definition file.
    logical, intent(in), optional :: ignore_unknown
        !! Ignore unknown keywords, without raising errors.

end subroutine read_gridDB_params

! ----------------------------------------------------------------------

end interface

contains

! ======================================================================

subroutine file_tcdcube_close_file(fileDB)
    !! Close file after finalizing writing, cutting connection.
    !!
    !! Finalizes data in header, closes file and resets connection
    !! information in file database.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.

    integer, dimension(:), allocatable :: ivec
    real(realwp), dimension(:), allocatable :: rvec

    if (fileDB%iunit == 0) then
        ! Nothing to do, resetting file name if set
        if (allocated(fileDB%file)) deallocate(fileDB%file)
        return
    end if

    ! Read first record and update it
    allocate(ivec(fileDB%len_data(1)), rvec(fileDB%len_data(2)))
    read(fileDB%iunit, rec=1) ivec, rvec
    ! Now update rvec with updated data
    ! max_state can be finalized.
    rvec(6) = fileDB%max_state
    write(fileDB%iunit, rec=1) ivec, rvec
    close(fileDB%iunit)
    fileDB%iunit = 0
    deallocate(fileDB%file)

end subroutine file_tcdcube_close_file

! ======================================================================

function file_tcdcube_init(gridDB, len_data, min_state) result(fileDB)
    !! Initialize a TCD cube file object from grid specifications.
    !!
    !! Initializes a TCD cube file object from a grid specifications
    !! database.
    !! The lowest state can be provided and cannot be modified
    !! afterwards since it would require moving all records each time.
    class(CubegridDB), intent(in) :: gridDB
        !! Grid specifications database.
    integer, intent(in) :: len_data
        !! Length of dataset at each point (1 for scalar)
    integer, intent(in), optional :: min_state
        !! Lowest state index.
    type(FileTCDCube) :: fileDB
        !! TCD Cube File database.

    integer, dimension(:), allocatable :: ivec
    real(realwp), dimension(:), allocatable :: rvec

    fileDB%version = 20250731
    fileDB%len_header = 1
    if (present(min_state)) then
        if (min_state < 0) then
            call fileDB%error%raise_argerror('value', &
                'Negative value not accepted for min_state')
            return
        end if
        fileDB%min_state = min_state
    end if
    fileDB%len_point = len_data
    fileDB%grid_min = gridDB%min
    fileDB%grid_max = gridDB%max
    fileDB%grid_step = gridDB%step_size
    fileDB%n_points = gridDB%n_points
    ! Grid shape is converted to integer value for easier storage
    select case (gridDB%shape)
    case ('cube', 'cubic')
        fileDB%grid_shape_id = 1
    case default
        fileDB%grid_shape_id = 1
    end select

    fileDB%len_data(1) = 1
    fileDB%len_data(2) = fileDB%len_point*product(fileDB%n_points)
    if (fileDB%len_data(2) <= 0) then
        call fileDB%error%raise_error('calc', 'invalid', &
            'Could not set the size of the cube datasets')
        return
    end if
    allocate(ivec(fileDB%len_data(1)), rvec(fileDB%len_data(2)))
    inquire(iolength=fileDB%len_record) ivec, rvec
    deallocate(ivec, rvec)

end function file_tcdcube_init

! ======================================================================

subroutine file_tcdcube_init_file(fileDB, fname, overwrite_OK)
    !! Connect and write header to target file.
    !!
    !! Connects `fname` with the `fileDB` instance and writes the header
    !! content.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    character(len=*), intent(in) :: fname
        !! Name of the file to connect.
    logical, intent(in), optional :: overwrite_OK
        !! If `fname` exists, the routine will overwrite it.

    integer :: ios
    integer, dimension(:), allocatable :: ivec
    real(realwp), dimension(:), allocatable :: rvec
    logical :: exists, overwrite
    character(len=512) :: msg

    if (present(overwrite_OK)) then
        overwrite = overwrite_OK
    else
        ! For safety, we do not consent overwriting an existing file by default
        overwrite = .false.
    end if

    ! Check if database already connected to a file and close.
    if (fileDB%iunit /= 0) call fileDB%close()

    ! Check if given file exists and overwriting OK.
    inquire(file=fname, exist=exists)
    if (exists .and..not. overwrite) then
        write(msg, '("File ",a," already exists.")') trim(fname)
        call fileDB%error%raise_error('file', 'write', trim(msg))
        return
    end if

    fileDB%file = trim(fname)
    open(newunit=fileDB%iunit, file=fileDB%file, access='direct', &
         action='readwrite', recl=fileDB%len_record)
    allocate(ivec(fileDB%len_data(1)), rvec(fileDB%len_data(2)))
    ! Build arrays for data
    ivec(1) = fileDB%version
    rvec(1) = fileDB%len_header
    rvec(2) = fileDB%len_record
    rvec(3) = fileDB%len_data(1)
    rvec(4) = fileDB%len_data(2)
    rvec(5) = fileDB%min_state
    rvec(6) = 0
    rvec(7) = fileDB%grid_shape_id
    rvec(8) = fileDB%len_point
    rvec(9:11) = fileDB%n_points
    rvec(12:14) = fileDB%grid_min
    rvec(15:17) = fileDB%grid_max
    rvec(18:20) = fileDB%grid_step
    write(fileDB%iunit, rec=1, iostat=ios) ivec, rvec
    if (ios /= 0) then
        call fileDB%error%raise_error('file', 'write', &
            'Unable to write into TCD cube data file.')
        return
    end if

end subroutine file_tcdcube_init_file

! ======================================================================

function file_tcdcube_load(fname, prt_warn) result(fileDB)
    !! Initialize a TCD cube file object by reading data from a file.
    !!
    !! Initializes a TCD cube file object by reading header data from
    !! an existing file
    character(len=*), intent(in) :: fname
        !! Filename with TCD cube data.
    logical, intent(in), optional :: prt_warn
        !! Print warning messages.
    type(FileTCDCube) :: fileDB
        !! TCD Cube File database.

    integer, parameter :: LENDIM = 4
    integer :: iu, n
    integer, dimension(:), allocatable :: ivec
    real(realwp), dimension(LENDIM) :: vhead
    real(realwp), dimension(:), allocatable :: rvec
    logical :: do_warn, exists

    inquire(file=fname, exist=exists)
    if (.not.exists) then
        call fileDB%error%raise_error('file', 'not found', 'File not found')
        return
    end if

    if (present(prt_warn)) then
        do_warn = prt_warn
    else
        do_warn = .false.
    end if

    fileDB%file = trim(fname)
    inquire(iolength=n) iu, vhead
    open(newunit=iu, file=fileDB%file, access='direct', status='old', recl=n)
    read(iu, rec=1) fileDB%version, vhead
    close(iu)
    ! Parsing of the lead head, which contains information to properly open and
    ! parse the file.
    select case (fileDB%version)
    case (20250731)
        fileDB%len_header = int(vhead(1))
        fileDB%len_record = int(vhead(2))
        fileDB%len_data(1) = int(vhead(3)) ! Number of leading integers elements
        fileDB%len_data(2) = int(vhead(4)) ! Number of data elements (real)
        allocate(ivec(fileDB%len_data(1)), rvec(fileDB%len_data(2)))
        ! Now check record length is consistent (may be system-dependent)
        inquire(iolength=n) ivec, rvec
        if (n /= fileDB%len_record) then
            if (do_warn) &
                call run%error%raise_warning('data', 'consistency', &
                    'Inconsistency between computed and stored records &
                    &lengths for data file.')
            fileDB%len_record = n
        end if
        open(newunit=fileDB%iunit, file=fileDB%file, access='direct', &
             status='old', recl=fileDB%len_record)
        read(fileDB%iunit, rec=1) ivec, rvec
        fileDB%min_state = int(rvec(LENDIM+1))
        fileDB%max_state = int(rvec(LENDIM+2))
        if (fileDB%max_state == 0 .and. do_warn) then
            call run%error%raise_warning('data', 'missing', &
                'Maximum state unset. The file may not have been finalized &
                &properly!')
        end if
        fileDB%grid_shape_id = int(rvec(LENDIM+3))
        fileDB%len_point = int(rvec(LENDIM+4))
        fileDB%n_points = int(rvec(LENDIM+5:LENDIM+7))
        fileDB%grid_min = rvec(LENDIM+8:LENDIM+10)
        fileDB%grid_max = rvec(LENDIM+11:LENDIM+13)
        fileDB%grid_step = rvec(LENDIM+14:LENDIM+16)
        deallocate(ivec, rvec)
    case default
        call fileDB%error%raise_error('data', 'value', &
            'Unsupported version of the TCD cube file.')
    end select

end function file_tcdcube_load

! ======================================================================

subroutine file_tcdcube_read_data_multi(fileDB, states, tcd_data)
    !! Read multiple TCD cube data entries from fileDB.
    !!
    !! Reads a list of TCD cube data from an existing fileDB.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    integer, dimension(:), intent(in) :: states
        !! List of states for the data to extract
    real(realwp), dimension(:,:), intent(out) :: tcd_data
        !! TCD data extracted from file.

    integer :: i, ios, istat
    character(len=256) :: msg

    if (fileDB%iunit == 0) then
        call fileDB%error%raise_error('file', 'missing', 'No file connected')
        return
    end if

    if (size(tcd_data, 2) /= size(states) &
            .or. size(tcd_data, 1) /= fileDB%len_data(2)) then
        call fileDB%error%raise_error('data', 'inconsistent', &
            'Size inconsistency for data extraction.')
        return
    end if

    do i = 1, size(states)
        read(fileDB%iunit, &
             rec=fileDB%len_header+states(i)+1-fileDB%min_state, iostat=ios) &
            istat, tcd_data(:,i)
        if (ios /= 0) then
            write(msg, '("Unable to read data for state ",i0)') states(i)
            call fileDB%error%raise_error('file', 'read', trim(msg))
            return
        else if (istat == 0) then
            write(msg, '("Missing data for state ",i0)') states(i)
            call fileDB%error%raise_error('data', 'missing', trim(msg))
            return
        end if
    end do

end subroutine file_tcdcube_read_data_multi

! ======================================================================

subroutine file_tcdcube_read_data_multi_scal(fileDB, states, tcd_data)
    !! Read multiple TCD scalar cube data entries from fileDB.
    !!
    !! Reads a list of TCD scalar cube data with true dimensions from
    !! an existing fileDB.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    integer, dimension(:), intent(in) :: states
        !! List of states for the data to extract
    real(realwp), dimension(:,:,:,:), intent(out) :: tcd_data
        !! Cube scalar data extracted from file, with full dimension.

    integer :: i, ios, istat, ncubes, ndata
    character(len=256) :: msg

    if (fileDB%iunit == 0) then
        call fileDB%error%raise_error('file', 'missing', 'No file connected')
        return
    end if

    ncubes = size(tcd_data, 4)
    ndata = size(tcd_data) / ncubes
    if (ncubes /= size(states) .or. ndata /= fileDB%len_data(2)) then
        call fileDB%error%raise_error('data', 'inconsistent', &
            'Size inconsistency for data extraction.')
        return
    end if

    do i = 1, size(states)
        read(fileDB%iunit, &
             rec=fileDB%len_header+states(i)+1-fileDB%min_state, iostat=ios) &
            istat, tcd_data(:,:,:,i)
        if (ios /= 0) then
            write(msg, '("Unable to read data for state ",i0)') states(i)
            call fileDB%error%raise_error('file', 'read', trim(msg))
            return
        else if (istat == 0) then
            write(msg, '("Missing data for state ",i0)') states(i)
            call fileDB%error%raise_error('data', 'missing', trim(msg))
            return
        end if
    end do

end subroutine file_tcdcube_read_data_multi_scal

! ======================================================================

subroutine file_tcdcube_read_data_multi_vec(fileDB, states, tcd_data)
    !! Read multiple TCD vectorial cube data entries from fileDB.
    !!
    !! Reads a list of TCD vectorial cube data with true dimensions from
    !! an existing fileDB.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    integer, dimension(:), intent(in) :: states
        !! List of states for the data to extract
    real(realwp), dimension(:,:,:,:,:), intent(out) :: tcd_data
        !! Cube vectorial data extracted from file, with full dimension.

    integer :: i, ios, istat, ncubes, ndata
    character(len=256) :: msg

    if (fileDB%iunit == 0) then
        call fileDB%error%raise_error('file', 'missing', 'No file connected')
        return
    end if

    ncubes = size(tcd_data, 5)
    ndata = size(tcd_data) / ncubes
    if (ncubes /= size(states) .or. ndata /= fileDB%len_data(2)) then
        call fileDB%error%raise_error('data', 'inconsistent', &
            'Size inconsistency for data extraction.')
        return
    end if

    do i = 1, size(states)
        read(fileDB%iunit, &
             rec=fileDB%len_header+states(i)+1-fileDB%min_state, iostat=ios) &
            istat, tcd_data(:,:,:,:,i)
        if (ios /= 0) then
            write(msg, '("Unable to read data for state ",i0)') states(i)
            call fileDB%error%raise_error('file', 'read', trim(msg))
            return
        else if (istat == 0) then
            write(msg, '("Missing data for state ",i0)') states(i)
            call fileDB%error%raise_error('data', 'missing', trim(msg))
            return
        end if
    end do

end subroutine file_tcdcube_read_data_multi_vec

! ======================================================================

subroutine file_tcdcube_read_data_single(fileDB, state, tcd_data)
    !! Read a single TCD cube data entry from fileDB.
    !!
    !! Reads a single entry of TCD cube data from an existing fileDB.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    integer, intent(in) :: state
        !! State of reference for data to extract.
    real(realwp), dimension(:), intent(out) :: tcd_data
        !! TCD data extracted from file.

    integer :: ios, istat
    character(len=256) :: msg

    if (fileDB%iunit == 0) then
        call fileDB%error%raise_error('file', 'missing', 'No file connected')
        return
    end if

    if (size(tcd_data) /= fileDB%len_data(2)) then
        call fileDB%error%raise_error('data', 'inconsistent', &
            'Size inconsistency for data extraction.')
        return
    end if

    read(fileDB%iunit, rec=fileDB%len_header+state+1-fileDB%min_state, &
         iostat=ios) istat, tcd_data
    if (ios /= 0) then
        write(msg, '("Unable to read data for state ",i0)') state
        call fileDB%error%raise_error('file', 'read', trim(msg))
        return
    else if (istat == 0) then
        write(msg, '("Missing data for state ",i0)') state
        call fileDB%error%raise_error('data', 'missing', trim(msg))
        return
    end if

end subroutine file_tcdcube_read_data_single

! ======================================================================

subroutine file_tcdcube_read_data_single_scal(fileDB, state, tcd_data)
    !! Read a single TCD scalar cube data entry from fileDB.
    !!
    !! Reads a single entry of scalar cube data with true dimensions
    !! from an existing fileDB.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    integer, intent(in) :: state
        !! State of reference for data to extract.
    real(realwp), dimension(:,:,:), intent(out) :: tcd_data
        !! Cube scalar data extracted from file, with full dimension.

    integer :: ios, istat
    character(len=256) :: msg

    if (fileDB%iunit == 0) then
        call fileDB%error%raise_error('file', 'missing', 'No file connected')
        return
    end if

    if (size(tcd_data) /= fileDB%len_data(2)) then
        call fileDB%error%raise_error('data', 'inconsistent', &
            'Size inconsistency for data extraction.')
        return
    end if

    read(fileDB%iunit, rec=fileDB%len_header+state+1-fileDB%min_state, &
         iostat=ios) istat, tcd_data
    if (ios /= 0) then
        write(msg, '("Unable to read data for state ",i0)') state
        call fileDB%error%raise_error('file', 'read', trim(msg))
        return
    else if (istat == 0) then
        write(msg, '("Missing data for state ",i0)') state
        call fileDB%error%raise_error('data', 'missing', trim(msg))
        return
    end if

end subroutine file_tcdcube_read_data_single_scal

! ======================================================================

subroutine file_tcdcube_read_data_single_vec(fileDB, state, tcd_data)
    !! Read a single TCD vectorial cube data entry from fileDB.
    !!
    !! Reads a single entry of vectorial cube data with true dimensions
    !! from an existing fileDB.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    integer, intent(in) :: state
        !! State of reference for data to extract.
    real(realwp), dimension(:,:,:,:), intent(out) :: tcd_data
        !! Cube vectorial data extracted from file, with full dimension.

    integer :: ios, istat
    character(len=256) :: msg

    if (fileDB%iunit == 0) then
        call fileDB%error%raise_error('file', 'missing', 'No file connected')
        return
    end if

    if (size(tcd_data) /= fileDB%len_data(2)) then
        call fileDB%error%raise_error('data', 'inconsistent', &
            'Size inconsistency for data extraction.')
        return
    end if

    read(fileDB%iunit, rec=fileDB%len_header+state+1-fileDB%min_state, &
         iostat=ios) istat, tcd_data
    if (ios /= 0) then
        write(msg, '("Unable to read data for state ",i0)') state
        call fileDB%error%raise_error('file', 'read', trim(msg))
        return
    else if (istat == 0) then
        write(msg, '("Missing data for state ",i0)') state
        call fileDB%error%raise_error('data', 'missing', trim(msg))
        return
    end if

end subroutine file_tcdcube_read_data_single_vec

! ======================================================================

subroutine file_tcdcube_write_data_multi(fileDB, states, tcd_data, fname)
    !! Write multiple TCD cube data entries into fileDB.
    !!
    !! Writes a list of TCD cube data into a fileDB.
    !! The routine can initialize the file if needed.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    integer, dimension(:), intent(in) :: states
        !! List of states for the data to write.
    real(realwp), dimension(:,:), intent(in) :: tcd_data
        !! TCD data to save in file.
    character(len=*), intent(in), optional :: fname
        !! Optional fname if file has not been yet initialized.

    integer :: i, ios
    character(len=256) :: msg

    if (fileDB%iunit == 0) then
        if (present(fname)) then
            call fileDB%init_file(fname)
        else
            call fileDB%error%raise_error('file', 'missing', &
                'No connected file to store TCD cube data.')
            return
        end if
    end if

    if (size(tcd_data, 2) /= size(states) &
        .or. size(tcd_data, 1) /= fileDB%len_data(2) ) then
        call fileDB%error%raise_error('data', 'inconsistent', &
            'Size inconsistency for data storage.')
        return
    end if

    do i = 1, size(states)
        write(fileDB%iunit, &
              rec=fileDB%len_header+states(i)+1-fileDB%min_state, iostat=ios) &
            1, tcd_data(:,i)
        if (ios /= 0) then
            write(msg, '("Unable to write data for state ",i0)') states(i)
            call fileDB%error%raise_error('file', 'write', trim(msg))
            return
        end if
        if (states(i) > fileDB%max_state) fileDB%max_state = states(i)
    end do

end subroutine file_tcdcube_write_data_multi

! ======================================================================

subroutine file_tcdcube_write_data_multi_scal(fileDB, states, tcd_data, fname)
    !! Write multiple TCD scalar cube data entries into fileDB.
    !!
    !! Write a list of TCD scalar cube data with true dimensions into
    !! an existing fileDB.
    !! The routine can initialize the file if needed.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    integer, dimension(:), intent(in) :: states
        !! List of states for the data to write.
    real(realwp), dimension(:,:,:,:), intent(in) :: tcd_data
        !! TCD data to save in file.
    character(len=*), intent(in), optional :: fname
        !! Optional fname if file has not been yet initialized.

    integer :: i, ios, ncubes, ndata
    character(len=256) :: msg

    if (fileDB%iunit == 0) then
        if (present(fname)) then
            call fileDB%init_file(fname)
        else
            call fileDB%error%raise_error('file', 'missing', &
                'No connected file to store TCD cube data.')
            return
        end if
    end if

    ncubes = size(tcd_data, 4)
    ndata = size(tcd_data) / ncubes
    if (ncubes /= size(states) .or. ndata /= fileDB%len_data(2)) then
        call fileDB%error%raise_error('data', 'inconsistent', &
            'Size inconsistency for data storage.')
        return
    end if

    do i = 1, size(states)
        write(fileDB%iunit, &
              rec=fileDB%len_header+states(i)+1-fileDB%min_state, iostat=ios) &
            1, tcd_data(:,:,:,i)
        if (ios /= 0) then
            write(msg, '("Unable to write data for state ",i0)') states(i)
            call fileDB%error%raise_error('file', 'write', trim(msg))
            return
        end if
        if (states(i) > fileDB%max_state) fileDB%max_state = states(i)
    end do

end subroutine file_tcdcube_write_data_multi_scal

! ======================================================================

subroutine file_tcdcube_write_data_multi_vec(fileDB, states, tcd_data, fname)
    !! Write multiple TCD vectorial cube data entries into fileDB.
    !!
    !! Write a list of TCD vectorial cube data with true dimensions into
    !! an existing fileDB.
    !! The routine can initialize the file if needed.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    integer, dimension(:), intent(in) :: states
        !! List of states for the data to write.
    real(realwp), dimension(:,:,:,:,:), intent(in) :: tcd_data
        !! TCD data to save in file.
    character(len=*), intent(in), optional :: fname
        !! Optional fname if file has not been yet initialized.

    integer :: i, ios, ncubes, ndata
    character(len=256) :: msg

    if (fileDB%iunit == 0) then
        if (present(fname)) then
            call fileDB%init_file(fname)
        else
            call fileDB%error%raise_error('file', 'missing', &
                'No connected file to store TCD cube data.')
            return
        end if
    end if

    ncubes = size(tcd_data, 5)
    ndata = size(tcd_data) / ncubes
    if (ncubes /= size(states) .or. ndata /= fileDB%len_data(2)) then
        call fileDB%error%raise_error('data', 'inconsistent', &
            'Size inconsistency for data storage.')
        return
    end if

    do i = 1, size(states)
        write(fileDB%iunit, &
              rec=fileDB%len_header+states(i)+1-fileDB%min_state, iostat=ios) &
            1, tcd_data(:,:,:,:,i)
        if (ios /= 0) then
            write(msg, '("Unable to write data for state ",i0)') states(i)
            call fileDB%error%raise_error('file', 'write', trim(msg))
            return
        end if
        if (states(i) > fileDB%max_state) fileDB%max_state = states(i)
    end do

end subroutine file_tcdcube_write_data_multi_vec

! ======================================================================

subroutine file_tcdcube_write_data_single(fileDB, state, tcd_data, fname)
    !! Write a single TCD cube data entry into fileDB.
    !!
    !! Writes a single entry of TCD cube data into an existing fileDB.
    !! The routine can initialize the file if needed.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    integer, intent(in) :: state
        !! State of reference for data to store.
    real(realwp), dimension(:), intent(in) :: tcd_data
        !! TCD data to save in file.
    character(len=*), intent(in), optional :: fname
        !! Optional fname if file has not been yet initialized.

    integer :: ios
    character(len=256) :: msg

    if (fileDB%iunit == 0) then
        if (present(fname)) then
            call fileDB%init_file(fname)
        else
            call fileDB%error%raise_error('file', 'missing', &
                'No connected file to store TCD cube data.')
            return
        end if
    end if

    if (size(tcd_data) /= fileDB%len_data(2)) then
        call fileDB%error%raise_error('data', 'inconsistent', &
            'Size inconsistency for data storage.')
        return
    end if

    write(fileDB%iunit, rec=fileDB%len_header+state+1-fileDB%min_state, &
          iostat=ios) 1, tcd_data
    if (ios /= 0) then
        write(msg, '("Unable to write data for state ",i0)') state
        call fileDB%error%raise_error('file', 'write', trim(msg))
        return
    end if
    if (state > fileDB%max_state) fileDB%max_state = state

end subroutine file_tcdcube_write_data_single

! ======================================================================

subroutine file_tcdcube_write_data_single_scal(fileDB, state, tcd_data, fname)
    !! Write a single TCD scalar cube data entry into fileDB.
    !!
    !! Writes a single entry of scalar cube data with true dimensions
    !! into an existing fileDB.
    !! The routine can initialize the file if needed.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    integer, intent(in) :: state
        !! State of reference for data to store.
    real(realwp), dimension(:,:,:), intent(in) :: tcd_data
        !! TCD data to save in file.
    character(len=*), intent(in), optional :: fname
        !! Optional fname if file has not been yet initialized.

    integer :: ios
    character(len=256) :: msg

    if (fileDB%iunit == 0) then
        if (present(fname)) then
            call fileDB%init_file(fname)
        else
            call fileDB%error%raise_error('file', 'missing', &
                'No connected file to store TCD cube data.')
            return
        end if
    end if

    if (size(tcd_data) /= fileDB%len_data(2)) then
        call fileDB%error%raise_error('data', 'inconsistent', &
            'Size inconsistency for data storage.')
        return
    end if

    write(fileDB%iunit, rec=fileDB%len_header+state+1-fileDB%min_state, &
          iostat=ios) 1, tcd_data
    if (ios /= 0) then
        write(msg, '("Unable to write data for state ",i0)') state
        call fileDB%error%raise_error('file', 'write', trim(msg))
        return
    end if
    if (state > fileDB%max_state) fileDB%max_state = state

end subroutine file_tcdcube_write_data_single_scal

! ======================================================================

subroutine file_tcdcube_write_data_single_vec(fileDB, state, tcd_data, fname)
    !! Write a single TCD vectorial cube data entry into fileDB.
    !!
    !! Writes a single entry of vectorial cube data with true dimensions
    !! into an existing fileDB.
    !! The routine can initialize the file if needed.
    class(FileTCDCube), intent(inout) :: fileDB
        !! TCD Cube File database.
    integer, intent(in) :: state
        !! State of reference for data to store.
    real(realwp), dimension(:,:,:,:), intent(in) :: tcd_data
        !! TCD data to save in file.
    character(len=*), intent(in), optional :: fname
        !! Optional fname if file has not been yet initialized.

    integer :: ios
    character(len=256) :: msg

    if (fileDB%iunit == 0) then
        if (present(fname)) then
            call fileDB%init_file(fname)
        else
            call fileDB%error%raise_error('file', 'missing', &
                'No connected file to store TCD cube data.')
            return
        end if
    end if

    if (size(tcd_data) /= fileDB%len_data(2)) then
        call fileDB%error%raise_error('data', 'inconsistent', &
            'Size inconsistency for data storage.')
        return
    end if

    write(fileDB%iunit, rec=fileDB%len_header+state+1-fileDB%min_state, &
          iostat=ios) 1, tcd_data
    if (ios /= 0) then
        write(msg, '("Unable to write data for state ",i0)') state
        call fileDB%error%raise_error('file', 'write', trim(msg))
        return
    end if
    if (state > fileDB%max_state) fileDB%max_state = state

end subroutine file_tcdcube_write_data_single_vec

! ======================================================================

subroutine gen_cube_data_density(cubeDB, bsetDB, e_dens, e_dens_sum, prt_warn, molDB)
    !! Generate cube data with electronic density.
    !!
    !! Generates the cube data with electronic density.
    !! If `cubefile` is provided, the data are printed in output.
    class(CubeDensDB), intent(inout) :: cubeDB
        !! Density cube database.
    class(BasisSetDB), intent(in) :: bsetDB
        !! Basis set database.
    real(realwp), dimension(:,:), intent(in), optional :: e_dens
        !! Electronic density.
    real(realwp), intent(inout) :: e_dens_sum
        !! Integrated dlectronic density.
    logical, intent(in), optional :: prt_warn
        !! Print warning messages if risk of divergences.
    class(MoleculeDB), intent(in), target, optional :: molDB
        !! Molecule database.

    integer :: ix, iy, iz
    real(realwp) :: x, y, z
    real(realwp), dimension(:), allocatable :: chi_at, tmp

    if (present(molDB)) then
        if(associated(cubeDB%mol)) then
            call run%error%raise_argerror('conflict', &
                'Molecular specification specified twice', &
                source='gen_cube_data_density')
            return
        else
            cubeDB%mol => molDB
        end if
    else if(.not.associated(cubeDB%mol)) then
        call run%error%raise_argerror('missing', &
            'Missing molecular specification', source='gen_cube_data_density')
        return
    end if

    if (allocated(cubeDB%data)) deallocate(cubeDB%data)
    allocate(cubeDB%data(cubeDB%grid%n_points(3),cubeDB%grid%n_points(2), &
                         cubeDB%grid%n_points(1)))
    allocate(chi_at(bsetDB%n_basis), tmp(bsetDB%n_basis))

    chi_at = f0
    tmp = f0
    e_dens_sum = f0

    !$omp parallel do collapse(3) private(x, y, z, chi_at, tmp) reduction(+:e_dens_sum)
    do ix = 1, cubeDB%grid%n_points(1)
        do iy = 1, cubeDB%grid%n_points(2)
            do iz = 1, cubeDB%grid%n_points(3)
                x = cubeDB%grid%min(1) + (ix-1)*cubeDB%grid%step_size(1)
                y = cubeDB%grid%min(2) + (iy-1)*cubeDB%grid%step_size(2)
                z = cubeDB%grid%min(3) + (iz-1)*cubeDB%grid%step_size(3)

                call eval_AOs_chi_at(cubeDB%mol, bsetDB, x, y, z, chi_at, &
                                     prt_warn)
                call xgemv('N', bsetDB%n_basis, bsetDB%n_basis, f1, e_dens, &
                           bsetDB%n_basis, chi_at, 1, f0, tmp, 1)
                cubeDB%data(iz,iy,ix) = xdot(bsetDB%n_basis, chi_at, 1, tmp, 1)
                e_dens_sum = e_dens_sum + cubeDB%data(iz,iy,ix)
            end do
        end do
    end do
    !$omp end parallel do

    ! NOTE: Add half the step size to both min and max
    ! to accurately compute the integration volume.

    e_dens_sum = e_dens_sum * product((cubeDB%grid%max + cubeDB%grid%step_size) - cubeDB%grid%min) &
        / product(cubeDB%grid%n_points)

    deallocate(chi_at, tmp)

end subroutine gen_cube_data_density

! ======================================================================

subroutine gen_cube_data_tcd(cubeDB, bsetDB, e_trans_dens,  dip, gen_axial, &
                             length_gauge, prt_warn, molDB)
    !! Generate cube data with transition current density.
    !!
    !! Generates the cube data with transition current density.
    !! If `cubefile` is provided, the data are printed in output.
    class(CubeTCDDB), intent(inout) :: cubeDB
        !! Density cube database.
    class(BasisSetDB), intent(in) :: bsetDB
        !! Basis set database.
    real(realwp), dimension(:,:), intent(in), optional :: e_trans_dens
        !! Electronic transition density.
    logical, intent(in), optional :: gen_axial
        !! Generate axial tensor instead of polar tensor.
    logical, intent(in), optional :: length_gauge
        !! Generate length-gauge transition dipole density.
    real(realwp), dimension(3), intent(inout) :: dip
        !! Integrated dipole moment tensor.
    logical, intent(in), optional :: prt_warn
        !! Print warning messages if risk of divergences.
    class(MoleculeDB), intent(in), target, optional :: molDB

    integer :: i, ix, iy, iz
    real(realwp) :: x, y, z, rho_at
    real(realwp), dimension(3) :: r_at, tcd_at
    real(realwp), dimension(:), allocatable :: chi_at, tmp
    real(realwp), dimension(:,:), allocatable :: d1_chi_at
    logical :: do_axial, do_length

    if (present(gen_axial)) then
        do_axial = gen_axial
    else
        do_axial = .false.
    end if

    if (present(length_gauge)) then
        do_length = length_gauge
    else
        do_length = .false.
    end if

    if (do_length .and. do_axial) then
        call run%error%raise_argerror('conflict', &
            'Length-gauge transition dipole density cannot be combined with &
            &axial tensor.')
        return
    end if

    if (present(molDB)) then
        if(associated(cubeDB%mol)) then
            call run%error%raise_argerror('conflict', &
                'Molecular specification specified twice', &
                source='gen_cube_data_density')
            return
        else
            cubeDB%mol => molDB
        end if
    else if(.not.associated(cubeDB%mol)) then
        call run%error%raise_argerror('missing', &
            'Missing molecular specification', source='gen_cube_data_density')
        return
    end if

    if (allocated(cubeDB%data)) deallocate(cubeDB%data)
    cubeDB%vector_size = 3

    allocate(cubeDB%data(cubeDB%vector_size,cubeDB%grid%n_points(3), &
                         cubeDB%grid%n_points(2),cubeDB%grid%n_points(1)))
    allocate(chi_at(bsetDB%n_basis), tmp(bsetDB%n_basis), &
             d1_chi_at(3,bsetDB%n_basis))

    chi_at = f0
    d1_chi_at = f0
    tmp = f0
    dip = f0
    tcd_at = f0

    !$omp parallel do collapse(3) private(x, y, z, &
    !$omp chi_at, d1_chi_at, tmp, rho_at, r_at, tcd_at) reduction(+:dip)
    do ix = 1, cubeDB%grid%n_points(1)
        do iy = 1, cubeDB%grid%n_points(2)
            do iz = 1, cubeDB%grid%n_points(3)
                x = cubeDB%grid%min(1) + (ix-1)*cubeDB%grid%step_size(1)
                y = cubeDB%grid%min(2) + (iy-1)*cubeDB%grid%step_size(2)
                z = cubeDB%grid%min(3) + (iz-1)*cubeDB%grid%step_size(3)
                r_at = [x, y, z]

                if (do_length) then
                    call eval_AOs_chi_at(cubeDB%mol, bsetDB, x, y, z, &
                                         chi_at, prt_warn)
                else
                    call eval_AOs_nabla_chi_at(cubeDB%mol, bsetDB, x, y, z, &
                                               chi_at, d1_chi_at, prt_warn)
                end if

                call xgemv('N', bsetDB%n_basis, bsetDB%n_basis, f1, &
                           e_trans_dens, bsetDB%n_basis, chi_at, 1, f0, tmp, 1)

                if (do_length) then
                    rho_at = xdot(bsetDB%n_basis, chi_at, 1, tmp, 1)
                    tcd_at = -r_at*rho_at
                else
                    do i = 1, 3
                        tcd_at(i) = - xdot(bsetDB%n_basis, d1_chi_at(i,:), 1, &
                                           tmp, 1)
                    end do
                    if (do_axial) tcd_at = cross(r_at, tcd_at) / f2
                end if

                cubeDB%data(:,iz,iy,ix) = tcd_at
                dip = dip + tcd_at
            end do
        end do
    end do
    !$omp end parallel do

    ! NOTE: Add half the step size to both min and max
    ! to accurately compute the integration volume.

    dip = dip * product((cubeDB%grid%max + cubeDB%grid%step_size) - cubeDB%grid%min) &
        / product(cubeDB%grid%n_points)

    deallocate(chi_at, d1_chi_at, tmp)

end subroutine gen_cube_data_tcd

subroutine init_cube_base_scalar(cubeDB, molDB, cube_type, iparams, rparams, &
                                 extra_head)
    !! Initialize basic content of cube database to store scalar data.
    class(CubeScalarDB), intent(inout) :: cubeDB
        !! Cube database.
    class(MoleculeDB), intent(in), target, optional :: molDB
        !! Molecule database.
    character(len=*), intent(in), optional :: cube_type
        !! Type of cube.
    integer, dimension(:), intent(in), optional :: iparams
        !! Integer-type parameters for header, cube-dependent.
    real(realwp), dimension(:), intent(in), optional :: rparams
        !! Real-type parameters for header, cube-dependent.
    character(len=*), intent(in), optional :: extra_head
        !! Extra header, printed below generated header line.

    integer :: key
    if (present(cube_type)) then
        key = test_cubetype_key(cube_type)
    else
        key = 0
    end if

    select type (cube => cubeDB)
    type is (CubeDensDB)
        select case(key)
        case (0, 1)
            cubeDB%type = cube_type_dens
            cubeDB%itype = 1
            cubeDB%header = '! Density cube file from ELEMENTS'
        case (2)
            cubeDB%type = cube_type_trdens
            cubeDB%itype = 2
            cubeDB%header = '! Transition density cube file from ELEMENTS'
        case default
            call run%error%raise_deverror('case', 'Unsupported type of cube', &
                                          source='init_cube_base')
            return
        end select
    class default
        if (present(cube_type)) then
            cubeDB%type = trim(cube_type)
            cubeDB%itype = -1
        else
            cubeDB%type = 'unknown'
            cubeDB%itype = -1
        end if
        cubeDB%header = '! Generic scalar cube file from ELEMENTS'
    end select
    if (present(extra_head)) cubeDB%extra_head = trim(extra_head)

    if (present(molDB)) cubeDB%mol => molDB

end subroutine init_cube_base_scalar

! ======================================================================

subroutine init_cube_base_vector(cubeDB, molDB, cube_type, iparams, rparams, &
                                 extra_head)
    !! Initialize basic content of cube database to store vectorial data.
    class(CubeVectorDB), intent(inout) :: cubeDB
        !! Cube database.
    class(MoleculeDB), intent(in), target, optional :: molDB
        !! Molecule database.
    character(len=*), intent(in), optional :: cube_type
        !! Type of cube.
    integer, dimension(:), intent(in), optional :: iparams
        !! Integer-type parameters for header, cube-dependent.
    real(realwp), dimension(:), intent(in), optional :: rparams
        !! Real-type parameters for header, cube-dependent.
    character(len=*), intent(in), optional :: extra_head
        !! Extra header, printed below generated header line.

    integer :: key

    if (present(cube_type)) then
        key = test_cubetype_key(cube_type)
    else
        key = 0
    end if

    select type (cube => cubeDB)
    type is (CubeTCDDB)
        select case (key)
        case (0, 3)
            cubeDB%type = cube_type_etcd
            cubeDB%itype = 3
            cubeDB%vector_size = 3
            if (present(iparams)) then
                if (size(iparams) < 1) then
                    call cubeDB%error%raise_argerror('size', &
                        'Not enough parameters for cube header')
                    return
                end if
                cubeDB%header = ' '
                write(cubeDB%header, &
                      '("! ETCD cube file from ELEMENTS &
                      &| Excited state: ",i0)') &
                    iparams(1)
            else
                cubeDB%header = '! ETCD cube file from ELEMENTS'
            end if
        case (4)
            cubeDB%type = cube_type_vtcd
            cubeDB%itype = 4
            cubeDB%vector_size = 3
            if (present(iparams)) then
                if (size(iparams) < 3) then
                    call cubeDB%error%raise_argerror('size', &
                        'Not enough parameters for cube header')
                    return
                end if
                cubeDB%header = ' '
                write(cubeDB%header, &
                      '("! VTCD cube file from ELEMENTS | Initial state: ",i0,&
                      &" | Final state: ",i0," | Mode: ",i0)') iparams(:3)
            else
                cubeDB%header = '! VTCD cube file from ELEMENTS'
            end if
        case (5)
            cubeDB%type = cube_type_etdd
            cubeDB%itype = 5
            cubeDB%vector_size = 3
            if (present(iparams)) then
                if (size(iparams) < 1) then
                    call cubeDB%error%raise_argerror('size', &
                        'Not enough parameters for cube header')
                    return
                end if
                cubeDB%header = ' '
                write(cubeDB%header, &
                      '("! ETDD cube file from ELEMENTS &
                      &| Excited state: ",i0)') &
                    iparams(1)
            else
                cubeDB%header = '! ETDD cube file from ELEMENTS'
            end if
        case (6)
            cubeDB%type = cube_type_vtdd
            cubeDB%itype = 6
            cubeDB%vector_size = 3
            if (present(iparams)) then
                if (size(iparams) < 3) then
                    call cubeDB%error%raise_argerror('size', &
                        'Not enough parameters for cube header')
                    return
                end if
                cubeDB%header = ' '
                write(cubeDB%header, &
                      '("! VTDD cube file from ELEMENTS | Initial state: ",i0,&
                      &" | Final state: ",i0," | Mode: ",i0)') iparams(:3)
            else
                cubeDB%header = '! VTDD cube file from ELEMENTS'
            end if
        case default
            call run%error%raise_deverror('case', 'Unsupported type of cube', &
                                          source='init_cube_base')
            return
        end select
    class default
        if (present(cube_type)) then
            cubeDB%type = trim(cube_type)
            cubeDB%itype = -1
        else
            cubeDB%type = 'unknown'
            cubeDB%itype = -1
        end if
        cubeDB%header = '! Generic vectorial cube file from ELEMENTS'
    end select
    if (present(extra_head)) cubeDB%extra_head = trim(extra_head)

    if (present(molDB)) cubeDB%mol => molDB

end subroutine init_cube_base_vector

! ======================================================================

function InitCube(cube_type, molDB) result(cubeDB)
    !! Initialize a cube database.
    !!
    !! Initializes a cube database, associated a molecular database if
    !! provided.
    character(len=*), intent(in) :: cube_type
        !! Type of cube.
    type(MoleculeDB), target, intent(in), optional :: molDB
        !! Molecule database.
    class(CubeBaseDB), allocatable :: cubeDB
        !! Cube database.

    integer :: key
    key = test_cubetype_key(cube_type)

    select case (key)
    case (1)
        allocate(CubeDensDB :: cubeDB)
        cubeDB%type = cube_type_dens
        cubeDB%itype = 1
    ! case ('trans_dens', 'transition_density', 'e_trans_dens', &
    !       'e_trans_density')
    !     allocate(CubeScalarDB :: cubeDB)
    !     cubeDB%type = 'transition density'
    !     cubeDB%itype = 2
    case (3)
        allocate(CubeTCDDB :: cubeDB)
        cubeDB%type = cube_type_etcd
        cubeDB%itype = 3
    case (4)
        allocate(CubeTCDDB :: cubeDB)
        cubeDB%type = cube_type_vtcd
        cubeDB%itype = 4
    case (5)
        allocate(CubeTCDDB :: cubeDB)
        cubeDB%type = cube_type_etdd
        cubeDB%itype = 5
    case (6)
        allocate(CubeTCDDB :: cubeDB)
        cubeDB%type = cube_type_vtdd
        cubeDB%itype = 6
    case default
        call run%error%raise_deverror('case', 'Unsupported type of cube', &
                                      source='InitCube')
        return
    end select

    if (present(molDB)) then
        cubeDB%mol => molDB
    else if (associated(cubeDB%mol)) then
        cubeDB%mol => null()
    end if

end function InitCube

! ======================================================================

function test_cubetype_key(cube_type, check_dens, check_tcd) result(ikey)
    !! Check cubetype label and return identifier key.
    !!
    !! The procedure tests different aliases and should be preferred for
    !! standardized tests of keywords.
    !! The test is always performed case insensitive.
    !! The function returns -1 if the keyword is not supported.
    !!
    !! Supported values:
    !! * `-1`: unsupported type
    !! * `1`: electronic density
    !! * `2`: electronic transition density
    !! * `3`: electronic transition current density
    !! * `4`: vibrational transition current density
    !! * `5`: electronic transition dipole density
    !! * `6`: vibrational transition dipole density
    character(len=*), intent(in) :: cube_type
        !! Type of cube.
    logical, intent(in), optional :: check_dens
        !! Check against density-related keywords.
    logical, intent(in), optional :: check_tcd
        !! Check against transition current density-related keywords.
    integer :: ikey
        !! Identifier key.

    logical :: test_dens, test_tcd
    character(len=:), allocatable :: keyword

    if (present(check_dens)) then
        test_dens = check_dens
    else
        test_dens = .true.
    end if

    if (present(check_tcd)) then
        test_tcd = check_tcd
    else
        test_tcd = .true.
    end if

    keyword = locase(trim(cube_type))

    ikey = -1
    if (test_dens) then
        select case(keyword)
        case ('dens', 'density', 'e_dens', 'e_density')
            ikey = 1
        case ('trans_dens', 'transition_density', 'e_trans_dens', &
              'e_trans_density')
            ikey = 2
        end select
    end if
    if (test_tcd .and. ikey == -1) then
        select case(keyword)
        case ('tcd', 'trans_curr_dens', 'etcd')
            ikey = 3
        case ('vtcd', 'vib_trans_curr_dens')
            ikey = 4
        case ('etdd', 'tdd', 'transition dipole density')
            ikey = 5
        case ('vtdd', 'vibrational transition dipole density')
            ikey = 6
        end select
    end if

end function test_cubetype_key

! ======================================================================

subroutine write_cube_scalar_data(cubeDB, iunit, molDB)
    !! Generate a .cube file from scalar data contained in `CubeDB`.
    class(CubeScalarDB), intent(inout) :: cubeDB
        !! Cube database.
    integer, intent(in) :: iunit
        !! Unit identifier to file opened in writing.
    class(MoleculeDB), intent(in), target, optional :: molDB
        !! Molecule database.

    integer :: ia, ix, iy, iz
    logical :: opened
    character(len=8) :: writable
    class(MoleculeDB), pointer :: mol => null()

    ! Check unit is properly connected to writable file
    inquire(iunit, opened=opened, write=writable)
    if (.not.opened .or. writable == 'NO') then
        call cubeDB%error%raise_error('file', 'write', 'Cannot write on file')
        return
    end if

    if (present(molDB)) then
        mol => molDB
    else if (associated(cubeDB%mol)) then
        mol => cubeDB%mol
    else
        call cubeDB%error%raise_argerror('missing', &
            'Missing molecular specifications')
        return
    end if

    1100 format(2x, i3, 2x, f10.6, 2x, f10.6, 2x, f10.6, 3x, i2)
    1010 format(2x, i3, 2x, f10.6, 2x, f10.6, 2x, f10.6)
    1110 format(3x, i2, 2x, f10.6, 2x, f10.6, 2x, f10.6, 2x, f10.6)

    write(iunit, '(a)') trim(cubeDB%header)
    if (allocated(cubeDB%extra_head)) &
        write(iunit, '(a)') trim(cubeDB%extra_head)
    write(iunit, '("! Grid parameters (Bohr):")')

    write(iunit, 1100) mol%n_at, cubeDB%grid%min, 1
    write(iunit, 1010) &
        cubeDB%grid%n_points(1), cubeDB%grid%step_size(1), f0, f0
    write(iunit, 1010) &
        cubeDB%grid%n_points(2), f0, cubeDB%grid%step_size(2), f0
    write(iunit, 1010) &
        cubeDB%grid%n_points(3), f0, f0, cubeDB%grid%step_size(3)

    do ia = 1, mol%n_at
        write(iunit, 1110) mol%at_num(ia), mol%at_chg(ia), &
                mol%at_crd(:,ia)
    end do

    do ix = 1, cubeDB%grid%n_points(1)
        do iy = 1, cubeDB%grid%n_points(2)
            write(iunit, '(6(e14.6,:))') &
                (cubeDB%data(iz,iy,ix), iz = 1, cubeDB%grid%n_points(3))

        end do
    end do

end subroutine write_cube_scalar_data

! ======================================================================

subroutine write_cube_vector_data(cubeDB, iunit, molDB)
    !! Generate a .cube file from scalar data contained in `CubeDB`.
    class(CubeVectorDB), intent(inout) :: cubeDB
        !! Cube database.
    integer, intent(in) :: iunit
        !! Unit identifier to file opened in writing.
    class(MoleculeDB), intent(in), target, optional :: molDB
        !! Molecule database.

    integer :: i, ix, iy, iz, ia
    logical :: opened
    character(len=8) :: writable
    class(MoleculeDB), pointer :: mol => null()

    ! Check unit is properly connected to writable file
    inquire(iunit, opened=opened, write=writable)
    if (.not.opened .or. writable == 'NO') then
        call cubeDB%error%raise_error('file', 'write', 'Cannot write on file')
        return
    end if

    if (present(molDB)) then
        mol => molDB
    else if (associated(cubeDB%mol)) then
        mol => cubeDB%mol
    else
        call cubeDB%error%raise_argerror('missing', &
            'Missing molecular specifications')
        return
    end if

    1100 format(2x, i3, 2x, f10.6, 2x, f10.6, 2x, f10.6, 3x, i2)
    1010 format(2x, i3, 2x, f10.6, 2x, f10.6, 2x, f10.6)
    1110 format(3x, i2, 2x, f10.6, 2x, f10.6, 2x, f10.6, 2x, f10.6)

    write(iunit, '(a)') trim(cubeDB%header)
    if (allocated(cubeDB%extra_head)) &
        write(iunit, '(a)') trim(cubeDB%extra_head)
    write(iunit, '("! Grid parameters (Bohr):")')

    write(iunit, 1100) mol%n_at, cubeDB%grid%min, cubeDB%vector_size
    write(iunit, 1010) &
        cubeDB%grid%n_points(1), cubeDB%grid%step_size(1), f0, f0
    write(iunit, 1010) &
        cubeDB%grid%n_points(2), f0, cubeDB%grid%step_size(2), f0
    write(iunit, 1010) &
        cubeDB%grid%n_points(3), f0, f0, cubeDB%grid%step_size(3)

    do ia = 1, mol%n_at
        write(iunit, 1110) mol%at_num(ia), mol%at_chg(ia), &
                mol%at_crd(:,ia)
    end do

    do ix = 1, cubeDB%grid%n_points(1)
        do iy = 1, cubeDB%grid%n_points(2)
            write(iunit, '(6(e14.6,:))') &
                ((cubeDB%data(i,iz,iy,ix), i = 1, cubeDB%vector_size), &
                 iz = 1, cubeDB%grid%n_points(3))
        end do
    end do

end subroutine write_cube_vector_data

! ======================================================================

end module cubegen
