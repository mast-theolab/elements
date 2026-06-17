module gfaf_io
    !! Gaussian Fortran Array Files Input/Ouput module.
    !!
    !! The module provides objects (derived types) to operate on Gaussian's
    !! Fortran array files (FAF).

    use iso_fortran_env, only: int32, int64, real64
    use run_env, only: CoreExecObject
    use numeric, only: intwp
    use string, only: upcase
    use gauopen_drv, only: AOInts, DAOInts, LenArr, Rd_2EN, Rd_CBuf, Rd_ChBuf, &
        Rd_IBuf, Rd_RBuf, Rd_RInd, Rd_SpA, Rd_SpAC, Rd_Labl

    private
    public :: gfaf_data, gfaf_parser

    integer, parameter :: MAXLENLAB = 64

    type :: gfaf_data
        character(len=:), allocatable :: key
        character(len=1) :: dtype
            !! Data type, as string.  Possible values are:
            !! * '0': unset (data was not extracted)
            !! * 'I': integer
            !! * 'R': real
            !! * 'C': complex
            !! * 'L': logical
            !! * 'S': character string
        logical :: packed
            !! Matrix is packed, index array sets
        integer :: ni = 0, nr = 0
            !! Number of integer and real used for each element
        integer :: size = 0
            !! Total number of elements.
        integer :: n1 = 0, n2 = 0, n3 = 0, n4 = 0, n5 = 0
            !! Number of elements along each dimension (0 if unused).
        integer :: isym = 0
            !! standard (0), symmetric (1) or antisymm (-1) matrix.
        integer, dimension(:), allocatable :: idata
        real(real64), dimension(:), allocatable :: rdata
        complex(real64), dimension(:), allocatable :: cdata
        logical, dimension(:), allocatable :: ldata
        character(len=:), allocatable :: sdata
        integer, dimension(:), allocatable :: map
    end type gfaf_data

    type, extends (CoreExecObject) :: gfaf_obj
        !! Basic class to handle Gaussian Fortran Array Files.
        !!
        !! @note
        !! While the object supports different labels for 1D/2D matrices and
        !! 4D matrices, like the current `GauOpen` library, it expects them to
        !! be consistent.
        !! Internally, `int_size` is in bits, `len12D` and `len4D` in bytes.
        !! @endnote
        private
        character(len=:), allocatable :: name
            !! Filename.
        integer :: unit = 0
            !! Fortran unit.
        integer :: int_size = 0
            !! Integer size used in FAF.
        integer :: rec = 0
            !! Current record read/written.
        integer :: len12D = 0
            !! Integer size used in FAF for labels of sparse 1D/2D matrices.
        integer :: len4D = 0
            !! Integer size used in FAF for labels of sparse 4D matrices.
        integer :: n_base_recs = 0
            !! Number of basic records (labels) stored in FAF.
        integer :: version = 0
            !! File version, as stored in FAF.
        character(len=:), allocatable :: title
            !! Route title, as stored in FAF.
        character(len=:), allocatable :: ftype
            !! File type, as stored in FAF.
        character(len=:), allocatable :: gaussian
            !! Gaussian versione, as stored in FAF.
    contains
        procedure :: filename => gfaf_get_name
        procedure :: close => gfaf_close

    end type gfaf_obj

    type, extends(gfaf_obj) :: gfaf_parser
        private
        integer :: n_at = -1
            !! Number of atoms.
        integer :: n_el = -1
            !! Number of electrons
        integer :: n_basis = -1
            !! Number of basis functions.
        integer :: n_basok = -1
            !! Number of linearly independent basis functions.
        integer :: charge = 0
            !! Molecular charge
        integer :: multip = 0
            !! Spin multiplicity.
        integer :: flag_opcl = 0
            !! Gaussian internal open/closed shell flag.
        integer :: flag_cgu = -1
            !! Gaussian internal complex/general/unsigned flag.
        logical :: preloaded = .false.
            !! If true, keys have been preloaded with their record positions.
        character(len=MAXLENLAB), dimension(:), allocatable :: label_keys
            !! Store the label keys contained in FAF.
        integer, dimension(:), allocatable :: label_recs
            !! Store the label positions in FAF for faster retrieval.
    contains
        procedure, private :: gfaf_read_item, gfaf_read_items
        procedure, private :: read_data => gfaf_read_data
        procedure :: head_pars => gfaf_get_head
        procedure :: head_data => gfaf_get_hdata
        procedure :: skip => gfaf_skip_records
        procedure :: keys => gfaf_read_keys
        generic :: get => gfaf_read_item, gfaf_read_items
    end type gfaf_parser

    !type, public :: gfafbuilder
    !end type gfafbuilder

    interface gfaf_parser
        module procedure init_gfaf_parser
    end interface gfaf_parser

contains

! ======================================================================

function init_gfaf_parser(fname, int_size, fix_int_size, preload, &
                          exit_on_error, silent) result(gfaf)
    !! Constructor-like function to create a FAF parser instance.
    !!
    !! `int_size` should correspond to the labels size used to generate the
    !! FAF.  If unknown, the routine tries to guess it, first by using the
    !! default integer kind, and then trying to query the beginning of the
    !! file content.
    character(len=*), intent(in) :: fname
        !! Name of the Fortran Array File.
    integer, intent(in), optional :: int_size
        !! Size of integers stored in FAF file.
        !! If not provided, the program tries to guess it.
    logical, intent(in), optional :: fix_int_size
        !! Correct `int_size` if inconsistent with internal size.
        !! Default: true.
    logical, intent(in), optional :: preload
        !! Preload keys and their positions for faster searches.
    logical, intent(in), optional :: exit_on_error
        !! Exit if an error is encountered.  By default, `.true.`
    logical, intent(in), optional :: silent
        !! Do not print messages.  By default, `.false.`.
    type(gfaf_parser) :: gfaf
        !! Gaussian FAF parser instance.

    integer, parameter :: len_lab = 64, max_reclen = 10, ikind = kind(1)
    integer :: ios, ival, len_rec
    integer(int32) :: idata32(max_reclen), ivers32, nlab32
    integer(int64) :: idata64(max_reclen), ivers64, nlab64
    logical :: change_ok, exit_ok, no_print, preload_keys, exists
    character(len=len_lab) :: label
    character(len=256) :: msg
    
    ! Set error handling policy.
    if (present(exit_on_error)) then
        exit_ok = exit_on_error
    else
        exit_ok = .true.
    end if

    if (present(silent)) then
        no_print = silent
    else
        no_print = .false.
    end if

    call gfaf%error%init(exit_on_error=exit_ok, no_printing=no_print)

    ! File operations
    if (present(fix_int_size)) then
        change_ok = fix_int_size
    else
        change_ok = .true.
    end if

    if (present(preload)) then
        preload_keys = preload
    else
        preload_keys = .false.
    end if

    inquire(file=fname, exist=exists)
    if (.not.exists) then
        call gfaf%error%raise_error('file', 'missing', 'file not found')
        return
    else
        allocate(character(len=len_trim(fname)) :: gfaf%name)
        gfaf%name = trim(fname)
        open(newunit=gfaf%unit, file=gfaf%name, iostat=ios, &
             form='unformatted', action='read', status='old')
        if (ios /= 0) then
            call gfaf%error%raise_error('file', 'open', 'operation failed')
            return
        end if
    end if

    ! File appears valid, now check int_size
    ! We first try to guess it and will then check if the user provided a
    ! value and if this value is acceptable.
    ! We use size misalignment to exclude directly 32 bits specification
    ! with 64-bits storage
    block
        character(len=len_lab) :: gversion
        read(gfaf%unit, iostat=ios) label, ivers32, nlab32, gversion
        if (ios /= 0) then
            call gfaf%error%raise_error('file', 'read', &
                                        'failed to read 1st record of FAF')
            return
        end if
        if (ivers32 <= 0 .or. nlab32 <= 0 .or. nlab32 > 10000) then
            gfaf%int_size = 64
            rewind gfaf%unit
            read(gfaf%unit, iostat=ios) label, ivers64, nlab64, gversion
            gfaf%n_base_recs = int(nlab64, kind=ikind)
            gfaf%version = int(ivers64, kind=ikind)
        else
            gfaf%int_size = 32
            gfaf%n_base_recs = int(nlab32, kind=ikind)
            gfaf%version = int(ivers64, kind=ikind)
        end if
        gfaf%ftype = trim(label)
        gfaf%gaussian = trim(gversion)
    end block
    gfaf%rec = 1

    ! Still possible to have mismatch in size storage, so try to extract next
    ! record with storage length info.
    if (gfaf%version == 1) then
        len_rec = 8
    else
        len_rec = 10
    end if
    if (gfaf%int_size == 64) then
        read(gfaf%unit, iostat=ios) label, idata64(:len_rec)
        if (ios /= 0) then
            call gfaf%error%raise_error( &
                'file', 'read', 'failed to read 2nd record of FAF')
            return
        end if
        gfaf%title = trim(label)
        gfaf%n_at = int(idata64(1), kind=ikind)
        gfaf%n_basis = int(idata64(2), kind=ikind)
        gfaf%n_basok = int(idata64(3), kind=ikind)
        gfaf%charge = int(idata64(4), kind=ikind)
        gfaf%multip = int(idata64(5), kind=ikind)
        gfaf%n_el = int(idata64(6), kind=ikind)
        gfaf%len12D = int(idata64(7), kind=ikind)
        gfaf%len4D = int(idata64(8), kind=ikind)
        if (len_rec > 8) then
            gfaf%flag_opcl = int(idata64(9), kind=ikind)
            gfaf%flag_cgu = int(idata64(10), kind=ikind)
        end if
    else
        read(gfaf%unit, iostat=ios) label, idata32(:len_rec)
        if (ios /= 0) then
            call gfaf%error%raise_error( &
                'file', 'read', 'failed to read 2nd record of FAF')
            return
        end if
        gfaf%title = trim(label)
        gfaf%n_at = int(idata32(1), kind=ikind)
        gfaf%n_basis = int(idata32(2), kind=ikind)
        gfaf%n_basok = int(idata32(3), kind=ikind)
        gfaf%charge = int(idata32(4), kind=ikind)
        gfaf%multip = int(idata32(5), kind=ikind)
        gfaf%n_el = int(idata32(6), kind=ikind)
        gfaf%len12D = int(idata32(7), kind=ikind)
        gfaf%len4D = int(idata32(8), kind=ikind)
        if (len_rec > 8) then
            gfaf%flag_opcl = int(idata32(9), kind=ikind)
            gfaf%flag_cgu = int(idata32(10), kind=ikind)
        end if
    end if
    gfaf%rec = 2

    if (gfaf%len12D /= gfaf%len4D) then
        call gfaf%error%raise_error( &
            'file', 'read', 'Inconsistency in storage specification')
        return
    end if
    if (gfaf%len12D == 4) then
        gfaf%int_size = 32
    else if (gfaf%len12D == 8) then
        gfaf%int_size = 64
    else
        call gfaf%error%raise_error( &
            'file', 'read', 'Unknown storage specification')
        return
    end if

    if (present(int_size)) then
        if (int_size == 4 .or. int_size == 32) then
            ival = 32
        else if (int_size == 8 .or. int_size == 64) then
            ival = 64
        else
            call gfaf%error%raise_deverror('argval', &
                'wrong value in int_size', source='init_gfaf_parser')
            return
        end if
        if (ival /= gfaf%int_size) then
            write(msg,'("labels stored in file on ",i0," bits, not ",i0)') &
                  gfaf%int_size, ival
            if (change_ok) then
                call gfaf%error%raise_warning('data', 'inconsistency', msg)
            else
                call gfaf%error%raise_error('data', 'inconsistency', msg)
                return
            end if
        end if
    end if

    if (preload_keys) then
        call gfaf%keys()
    end if

end function init_gfaf_parser

! ======================================================================

function gfaf_close(gfaf) result(status)
    !! Close connection to Gaussian FAF.
    class(gfaf_obj), intent(inout), target :: gfaf
        !! Gaussian generic FAF handler instance.
    logical :: status
        !! Return if file closed properly.

    integer :: ios

    close(gfaf%unit, iostat=ios)
    status = ios == 0
    if (.not.status) &
        call gfaf%error%raise_error('file', 'close', &
                                    'could not close properly')

end function gfaf_close

! ======================================================================

function gfaf_get_name(gfaf) result(fname)
    !! Query name of the Gaussian FAF.
    class(gfaf_obj), intent(in), target :: gfaf
        !! Gaussian generic FAF handler instance.
    character(len=:), pointer :: fname
        !! Filename.
    fname => gfaf%name
end function gfaf_get_name

! ======================================================================

subroutine gfaf_get_head(gfaf, ftype, fversion, nlabels, gversion, title, &
                         natoms, nbasis, nbsuse, charge, multip, nel, &
                         flag_opcl, flag_cgu)
    !! Extract header information from FAF file.
    !!
    !! Returns header information from the file connected to the `gfaf_parser`
    !! instance.
    class(gfaf_parser), intent(inout) :: gfaf
        !! Gaussian FAF parser instance.
    character(len=*), intent(out), optional :: ftype
        !! Label containing the type of file.
    integer, intent(out), optional :: fversion
        !! Version of the generated binary file.
    integer, intent(out), optional :: nlabels
        !! Number of labels stored in the file.
    character(len=*), intent(out), optional :: gversion
        !! Gaussian version.
    character(len=*), intent(out), optional :: title
        !! First 64 characters of the route tile.
    integer, intent(out), optional :: natoms
        !! Number of atoms.
    integer, intent(out), optional :: nbasis
        !! Number of basis functions.
    integer, intent(out), optional :: nbsuse
        !! Number of linearly-independent basis functions
    integer, intent(out), optional :: charge
        !! Total charge.
    integer, intent(out), optional :: multip
        !! Total multiplicity.
    integer, intent(out), optional :: nel
        !! Number of electrons.
    integer, intent(out), optional :: flag_opcl
        !! Close/open-shell flags.
    integer, intent(out), optional :: flag_cgu
        !! Complex/generalized/unrestricted flag.

    if(present(ftype))     ftype     = gfaf%ftype
    if(present(fversion))  fversion  = gfaf%version
    if(present(nlabels))   nlabels   = gfaf%n_base_recs
    if(present(gversion))  gversion  = gfaf%gaussian
    if(present(title))     title     = gfaf%title
    if(present(natoms))    natoms    = gfaf%n_at
    if(present(nbasis))    nbasis    = gfaf%n_basis
    if(present(nbsuse))    nbsuse    = gfaf%n_basok
    if(present(charge))    charge    = gfaf%charge
    if(present(multip))    multip    = gfaf%multip
    if(present(nel))       nel       = gfaf%n_el
    if(present(flag_opcl)) flag_opcl = gfaf%flag_opcl
    if(present(flag_cgu))  flag_cgu  = gfaf%flag_cgu

end subroutine gfaf_get_head

! ======================================================================

subroutine gfaf_get_hdata(gfaf, at_num, at_type, at_crd, at_chg, at_wgt, &
                          bfunc_atm_map, bfunc_type, n_froz_core, &
                          n_froz_virt, flag_transf, n_ao_shells, n_ao_prim, &
                          n_dens_fit_shells, n_dens_fit_prim, n_tot_bonds)
    !! Read data stored in head records.
    !!
    !! Reads and returns data contained in header, related to molecular
    !! specifications and basis functions.
    !! It is possible to choose which data to extract since all are optional.
    class(gfaf_parser), intent(inout) :: gfaf
        !! Gaussian FAF parser instance.
    integer, dimension(:), allocatable, optional, intent(out) :: at_num
        !! Atomic numbers.
    integer, dimension(:), allocatable, optional, intent(out) :: at_type
        !! Atomic types, intended for multi-scale or with frozen atoms.
    real(real64), dimension(:,:), allocatable, optional, intent(out) :: at_crd
        !! Atomic coordinates.
    real(real64), dimension(:), allocatable, optional, intent(out) :: at_chg
        !! Atomic charges.
    real(real64), dimension(:), allocatable, optional, intent(out) :: at_wgt
        !! Atomic weights.
    integer, dimension(:), allocatable, optional, intent(out) :: bfunc_atm_map
        !! Mapping information on basis functions -> atoms.
    integer, dimension(:), allocatable, optional, intent(out) :: bfunc_type
        !! Information on the type of basis functions, as "lllmmm".
        !! * "lll" is the angular momentum
        !! * "mmm" is the component number
        !! * negative values for pure functions, positive for Cartesian.
    integer, optional, intent(out) :: n_froz_core
        !! Number of frozen core orbitals for 2e molecular orbitals.
    integer, optional, intent(out) :: n_froz_virt
        !! Number of frozen virtual orbitals for 2e molecular orbitals.
    integer, optional, intent(out) :: flag_transf
        !! Flag on transformation rule for stored 2e molecular orbitals.
        !! * 0: no MO integrals stored.
        !! * 4: MOs involving at least one occupied orbital are stored
        !! * 5: full transformation was done before storage.
    integer, optional, intent(out) :: n_ao_shells
        !! Number of contracted shells of AO basis functions
    integer, optional, intent(out) :: n_ao_prim
        !! Number of primitive AO shells.
    integer, optional, intent(out) :: n_dens_fit_shells
        !! Number of contracted shells of density fitting functions.
    integer, optional, intent(out) :: n_dens_fit_prim
        !! Number of primitive density fitting shells.
    integer, optional, intent(out) :: n_tot_bonds
        !! Total number of bonds in connectivity data.

    integer, parameter :: ikind = kind(1)
    integer :: i, nrecs
    integer(int32), dimension(:), allocatable :: iarr32, lrecs32
    integer(int64), dimension(:), allocatable :: iarr64, lrecs64

    ! Check that we are about to read 3rd record
    if (gfaf%rec /= 2) call gfaf%skip(-2)

    ! Read record num. 03
    if (present(at_num)) then
        if (gfaf%n_base_recs < 3) then
            call gfaf%error%raise_error( &
                'data', 'missing', 'missing atomic numbers in FAF head data')
            return
        end if
        if (gfaf%int_size == 32) then
            allocate(iarr32(gfaf%n_at))
            read(gfaf%unit) iarr32
            at_num = int(iarr32, kind=ikind)
            deallocate(iarr32)
        else
            allocate(iarr64(gfaf%n_at))
            read(gfaf%unit) iarr64
            at_num = int(iarr64, kind=ikind)
            deallocate(iarr64)
        end if
    else
        read(gfaf%unit)
    end if
    gfaf%rec = gfaf%rec + 1

    ! Read record num. 04
    if (present(at_type)) then
        if (gfaf%n_base_recs < 4) then
            call gfaf%error%raise_error( &
                'data', 'missing', 'missing atomic types in FAF head data')
            return
        end if
        if (gfaf%int_size == 32) then
            allocate(iarr32(gfaf%n_at))
            read(gfaf%unit) iarr32
            at_type = int(iarr32, kind=ikind)
            deallocate(iarr32)
        else
            allocate(iarr64(gfaf%n_at))
            read(gfaf%unit) iarr64
            at_type = int(iarr64, kind=ikind)
            deallocate(iarr64)
        end if
    else
        read(gfaf%unit)
    end if
    gfaf%rec = gfaf%rec + 1

    ! Read record num. 05
    if (present(at_chg)) then
        if (gfaf%n_base_recs < 5) then
            call gfaf%error%raise_error( &
                'data', 'missing', 'missing atomic charges in FAF head data')
            return
        end if
        allocate(at_chg(gfaf%n_at))
        read(gfaf%unit) at_chg
    else
        read(gfaf%unit)
    end if
    gfaf%rec = gfaf%rec + 1

    ! Read record num. 06
    if (present(at_crd)) then
        if (gfaf%n_base_recs < 6) then
            call gfaf%error%raise_error( &
                'data', 'missing', &
                'missing atomic coordinates in FAF head data')
            return
        end if
        allocate(at_crd(3,gfaf%n_at))
        read(gfaf%unit) at_crd
    else
        read(gfaf%unit)
    end if
    gfaf%rec = gfaf%rec + 1

    ! Read record num. 07
    if (present(bfunc_atm_map) .or. present(bfunc_type)) then
        if (gfaf%n_base_recs < 7) then
            call gfaf%error%raise_error( &
                'data', 'missing', &
                'missing data on basis functions in FAF head data')
            return
        end if
        if (gfaf%int_size == 32) then
            allocate(iarr32(2*gfaf%n_basis))
            read(gfaf%unit) iarr32
            if (present(bfunc_atm_map)) &
                bfunc_atm_map = int(iarr32(1:n_basis), kind=ikind)
            if (present(bfunc_type)) &
                bfunc_type = int(iarr32(n_basis+1:2*n_basis), kind=ikind)
            deallocate(iarr32)
        else
            allocate(iarr64(2*gfaf%n_basis))
            read(gfaf%unit) iarr64
            if (present(bfunc_atm_map)) &
                bfunc_atm_map = int(iarr64(1:n_basis), kind=ikind)
            if (present(bfunc_type)) &
                bfunc_type = int(iarr64(n_basis+1:2*n_basis), kind=ikind)
            deallocate(iarr64)
        end if
    else
        read(gfaf%unit)
    end if
    gfaf%rec = gfaf%rec + 1

    ! Read record num. 08
    if (present(at_wgt)) then
        if (gfaf%n_base_recs < 8) then
            call gfaf%error%raise_error( &
                'data', 'missing', 'missing atomic masses in FAF head data')
            return
        end if
        allocate(at_wgt(gfaf%n_at))
        read(gfaf%unit) at_wgt
    else
        read(gfaf%unit)
    end if
    gfaf%rec = gfaf%rec + 1

    ! Read record num. 09
    if (present(n_froz_core) .or. present(n_froz_virt) &
            .or. present(flag_transf)) then
        if (gfaf%n_base_recs < 9) then
            call gfaf%error%raise_error( &
                'data', 'missing', &
                'missing 2e MO integrals information missing in FAF head data')
            return
        end if
        if (gfaf%int_size == 32) then
            allocate(iarr32(4))
            read(gfaf%unit) iarr32
            if (present(n_froz_core)) n_froz_core = int(iarr32(1), kind=ikind)
            if (present(n_froz_virt)) n_froz_virt = int(iarr32(2), kind=ikind)
            if (present(flag_transf)) flag_transf = int(iarr32(3), kind=ikind)
            deallocate(iarr32)
        else
            allocate(iarr64(4))
            read(gfaf%unit) iarr64
            if (present(n_froz_core)) n_froz_core = int(iarr64(1), kind=ikind)
            if (present(n_froz_virt)) n_froz_virt = int(iarr64(2), kind=ikind)
            if (present(flag_transf)) flag_transf = int(iarr64(3), kind=ikind)
            deallocate(iarr64)
        end if
    else
        read(gfaf%unit)
    end if
    gfaf%rec = gfaf%rec + 1

    ! Read record num. 10
    ! The record contains the lengths of all records beyond it.
    if (gfaf%n_base_recs > 10) then
        nrecs = gfaf%n_base_recs - 10
        if (gfaf%int_size == 32) then
            allocate(lrecs32(nrecs))
            read(gfaf%unit) lrecs32
        else
            allocate(lrecs64(nrecs))
            read(gfaf%unit) lrecs64
        end if
    end if
    gfaf%rec = gfaf%rec + 1

    ! Read record num. 11
    if (present(n_ao_shells) &
            .or. present(n_ao_prim) &
            .or. present(n_dens_fit_shells) &
            .or. present(n_dens_fit_prim) &
            .or. present(n_tot_bonds)) then
        if (gfaf%n_base_recs < 11) then
            call gfaf%error%raise_error( &
                'data', 'missing', 'record 11 in FAF is missing')
            return
        end if
        if (gfaf%int_size == 32) then
            allocate(iarr32(lrecs32(1)))
            read(gfaf%unit) iarr32
            if (present(n_ao_shells)) &
                n_ao_shells = int(iarr32(1), kind=ikind)
            if (present(n_ao_prim)) &
                n_ao_prim = int(iarr32(2), kind=ikind)
            if (present(n_dens_fit_shells)) &
                n_dens_fit_shells = int(iarr32(3), kind=ikind)
            if (present(n_dens_fit_prim)) &
                n_dens_fit_prim = int(iarr32(4), kind=ikind)
            if (present(n_tot_bonds)) &
                n_tot_bonds = int(iarr32(5), kind=ikind)
            deallocate(iarr32)
        else
            allocate(iarr64(lrecs64(1)))
            read(gfaf%unit) iarr64
            if (present(n_ao_shells)) &
                n_ao_shells = int(iarr64(1), kind=ikind)
            if (present(n_ao_prim)) &
                n_ao_prim = int(iarr64(2), kind=ikind)
            if (present(n_dens_fit_shells)) &
                n_dens_fit_shells = int(iarr64(3), kind=ikind)
            if (present(n_dens_fit_prim)) &
                n_dens_fit_prim = int(iarr64(4), kind=ikind)
            if (present(n_tot_bonds)) &
                n_tot_bonds = int(iarr64(5), kind=ikind)
            deallocate(iarr64)
        end if
    else
        read(gfaf%unit)
    end if
    gfaf%rec = gfaf%rec + 1

    ! Records 12 and later (currently not supported)
    do i = 2, gfaf%n_base_recs - 10
        read(gfaf%unit)
        gfaf%rec = gfaf%rec + 1
    end do

end subroutine gfaf_get_hdata

! ======================================================================

function gfaf_read_data(gfaf, label, NI, NR, NTot, LenBuf, N1, N2, N3, N4, &
                        N5, IType, NRI) result(dbase)
    !! Read 1 dataset from a Gaussian Fortran Array Files.
    !!
    !! Reads data associated to a label based on the size specifications from
    !! the label record (parsed before).
    class(gfaf_parser), intent(inout) :: gfaf
        !! Gaussian FAF parser instance.
    character(len=*), intent(in) :: label
        !! Label title, as stored in FAF.
    integer, intent(in) :: NI
        !! Number of integers used for each element.
    integer, intent(in) :: NR
        !! Number of reals used for each element. `<0` for complex numbers.
    integer, intent(in) :: NTot
        !! Total number of elements.
    integer, intent(in) :: LenBuf
        !! Number of elements per record.
    integer, intent(in) :: N1, N2, N3, N4, N5
        !! Dimensions of the matrix. `0` if unused, `<0` for triangular form.
    integer, intent(in) :: IType
        !! Type of matrix:
        !! * `0`: standard/symmmetric matrices
        !! * `-1`: anti-symmetric/anti-Hermetian matrices
        !! * `>0` array contains character data of length `TypeA`.
    integer, intent(in) :: NRI
        !! Number of real/imaginary components.
        !! * `1`: real numbers are stored.
        !! * `2`: complex numbers are stored.
    type(gfaf_data) :: dbase
        !! Extracted data.

    integer :: LR, LenBC, nrecs, NTotC
    character(len=256) :: msg

    dbase%key = trim(label)

    LR = LenArr(N1, N2, N3, N4, N5)
    nrecs = (NTot+LenBuf-1) / LenBuf
    ! Cases we do no now how to proceed
    ! 1. Inconsistency between size from dimensions and total num. elements
    if(((NI*NR) == 0 .and. LR /= NTot) &  ! 1
            ) then
        call gfaf%skip(nrecs)
        dbase%dtype = '0'
        write(msg, '("Inconsistency in array specifications for label &
                    &""",a,""".")') trim(label)
        call gfaf%error%raise_error('data', 'inconsistency', trim(msg))
        return
    end if

    dbase%size = NTot
    dbase%ni = NI
    dbase%nr = NR
    dbase%n1 = N1
    dbase%n2 = N2
    dbase%n3 = N3
    dbase%n4 = N4
    dbase%n5 = N5
    if (N1 < 0 .or. N2 < 0 .or. N3 < 0 .or. N4 < 0 .or. N5 < 0) then
        if (TypeA >= 0) then
            dbase%isym = 1
        else
            dbase%isym = -1
        end if
    else
        dbase%isym = 0
    end if

    ! Processing of specific data
    if (AOInts(label)) then
        dbase%dtype = 'R'
        allocate(dbase%rdata(LR*NR))
        call Rd_2EN(gfaf%unit, NR, LR, NTot, LenBuf, dbase%rdata, &
                    gfaf%int_size)
    else if (DAOInts(label)) then
        dbase%dtype = 'R'
        allocate(dbase%rdata(NR*NTot), dbase%map(NI*NTot))
        call Rd_SpA(gfaf%unit, NI, NR, NTot, LenBuf, dbase%idata, &
                    dbase%rdata, gfaf%int_size)
    ! Generic processing
    else if (NI == 1 .and. NR == 0) then
        if (IType > 0) then
            dbase%dtype = 'S'
            NTotC = gfaf%len4D*NTot
            LenBC = gfaf%len4D*LenBuf
            allocate(character(len=NTotC) :: dbase%sdata)
            call Rd_ChBuf(gfaf%unit, NTotC, LenBC, dbase%sdata)
        else
            dbase%dtype = 'I'
            allocate(dbase%idata(NTot))
            call Rd_IBuf(gfaf%unit, NTot, LenBuf, dbase%idata, gfaf%int_size)
        end if
    else if (NI == 0 .and. NR == 1 .and. NRI == 1) then
        dbase%dtype = 'R'
        allocate(dbase%rdata(NTot))
        call Rd_RBuf(gfaf%unit, NTot, LenBuf, dbase%rdata)
    else if (NI == 0 .and. NR == 1 .and. NRI == 2) then
        dbase%dtype = 'C'
        allocate(dbase%cdata(NTot))
        call Rd_CBuf(gfaf%unit, NTot, LenBuf, dbase%cdata)
    else if (NI == 1) then
        if (NRI == 1) then
            dbase%dtype = 'R'
            allocate(dbase%rdata(NR*LR))
            call Rd_RInd(gfaf%unit, NR, LR, NR*LR, NTot, LenBuf, dbase%rdata, &
                         gfaf%int_size)
        else
            dbase%dtype = 'C'
            ! Rd_RInd is intended for real arrays, so we cannot pass a complex
            ! arrays in principle.
            ! To bypass the problem, we allocate rdata as well to act as a temp
            ! array and then copy back the real and imaginary parts into cdata.
            allocate(dbase%cdata(NR*LR), dbase%rdata(2*NR*LR))
            call Rd_RInd(gfaf%unit, NRI*NR, LR, NRI*NR*LR, NTot, LenBuf, &
                         dbase%rdata, gfaf%int_size)
            dbase%cdata%re = dbase%rdata(1::2)
            dbase%cdata%im = dbase%rdata(2::2)
            deallocate(dbase%rdata)
        end if
    else if (NI > 0 .and. NR > 0) then
        if (NRI == 1) then
            dbase%dtype = 'R'
            allocate(dbase%rdata(NR*NTot), dbase%map(NI*NTot))
            call Rd_SpA(gfaf%unit, NI, NR, NTot, LenBuf, dbase%idata, &
                        dbase%rdata, gfaf%int_size)
        else
            dbase%dtype = 'C'
            allocate(dbase%cdata(NR*NTot), dbase%map(NI*NTot))
            call Rd_SpAC(gfaf%unit, NI, NR, NTot, LenBuf, dbase%idata, &
                         dbase%cdata, gfaf%int_size)
        end if
    else
        call gfaf%error%raise_deverror('case', &
            'Unsupported array specification in label')
    end if

    gfaf%rec = gfaf%rec + nrecs

end function gfaf_read_data

! ======================================================================

function gfaf_read_item(gfaf, label) result(dbase)
    !! Read 1 item in the FAF file.
    !!
    !! Reads 1 item corresponding to label in Gaussian FAF.
    class(gfaf_parser), intent(inout) :: gfaf
        !! Gaussian FAF parser instance.
    character(len=*), intent(in) :: label
        !! Label of the quantity of interest.
    type(gfaf_data) :: dbase
        !! Extracted data.

    integer :: IType, LenBuf, N1, N2, N3, N4, N5, NI, NR, NRI, NTot
    logical :: EOF
    character(len=:), allocatable :: key
    character(len=MAXLENLAB) :: CBuf

    key = upcase(trim(label))

    ! Ensure that we are at the beginning of the file.
    if (gfaf%rec /= gfaf%n_base_recs) call gfaf%skip(-gfaf%n_base_recs)

    ! Now let us search for the key
    do
        call Rd_Labl(gfaf%unit, gfaf%version, CBuf, NI, NR, NTot, LenBuf, N1, &
            N2, N3, N4, N5, IType, NRI, EOF, LabVer=gfaf%int_size)
        if (EOF) then
            dbase%dtype = '0'
            return
            ! write(msg, '("Label """,a,""" not found in file.")') trim(label)
            ! call gfaf%error%raise_error('keyword', 'not found', trim(msg))
            ! return
        end if
        gfaf%rec = gfaf%rec + 1
        if (CBuf == key) then
            dbase = gfaf%read_data(CBuf, NI, NR, NTot, LenBuf, N1, N2, N3, &
                                   N4, N5, IType, NRI)
            exit
        else
            call gfaf%skip((NTot+LenBuf-1)/LenBuf)
        end if
    end do

end function gfaf_read_item

! ======================================================================

function gfaf_read_items(gfaf, labels) result(dbase)
    ! Read multiple items in the FAF file.
    !!
    !! Reads several items corresponding to label in Gaussian FAF.
    class(gfaf_parser), intent(inout) :: gfaf
        !! Gaussian FAF parser instance.
    character(len=*), dimension(:), intent(in) :: labels
        !! Labels of the quantities of interest.
    type(gfaf_data), dimension(:), allocatable :: dbase
        !! Extracted data.

    integer :: idx, ikey, IType, LenBuf, N1, N2, N3, N4, N5, NI, nkeys, NR, &
        NRI, NTot
    logical :: EOF
    integer, dimension(:), allocatable :: indexes ! store indexes still to find
    character(len=MAXLENLAB) :: CBuf

    ! Build a list of key indexes that will be updated to check the keys
    ! to search.
    nkeys = size(labels)
    allocate(dbase(nkeys), indexes(nkeys))
    indexes = [(i, i = 1, nkeys)]

    ! Ensure that we are at the beginning of the file.
    if (gfaf%rec /= gfaf%n_base_recs) call gfaf%skip(-gfaf%n_base_recs)

    ! Now let us search for the keys
    main: do
        call Rd_Labl(gfaf%unit, gfaf%version, CBuf, NI, NR, NTot, LenBuf, N1, &
            N2, N3, N4, N5, IType, NRI, EOF, LabVer=gfaf%int_size)
        if (EOF) exit
        gfaf%rec = gfaf%rec + 1
        do ikey = 1, nkeys
            if (CBuf == upcase(labels(indexes(ikey)))) then
                idx = indexes(ikey)
                dbase(idx) = gfaf%read_data(CBuf, NI, NR, NTot, LenBuf, N1, &
                                           N2, N3, N4, N5, IType, NRI)
                indexes(ikey:) = eoshift(indexes(ikey:), 1)
                nkeys = nkeys - 1
                if (nkeys == 0) exit main
                exit
            end if
        end do
    end do main

    if (nkeys > 0) then
        do ikey = 1, nkeys
            idx = indexes(ikey)
            dbase(idx)%dtype = '0'
            dbase(idx)%key = trim(labels(idx))
        end do
    end if

end function gfaf_read_items

! ======================================================================

subroutine gfaf_read_keys(gfaf, count, keys, recs, load)
    !! Read keys stored in FAF and basic information.
    !!
    !! The routine can be used to load information in the parser DB for faster
    !! search (default with no argument), and/or to retrieve information.
    !! If data are requested, the routine avoids reading unnecessarily the
    !! file.
    class(gfaf_parser), intent(inout) :: gfaf
        !! Gaussian FAF parser instance.
    integer, intent(out), optional :: count
        !! Number of labels stored in FAF.
    character(len=MAXLENLAB), dimension(:), allocatable, intent(out), &
        optional :: keys
        !! List of label keys stored in FAF.
    integer, dimension(:), allocatable, intent(out), optional :: recs
        !! List of record positions
    logical, intent(out), optional :: load
        !! Load information in FAF parser, default if no other arguments given.

    integer :: IType, LenBuf, N1, N2, N3, N4, N5, NI, nkeys, NR, NRI, NTot
    logical :: do_load, EOF, get_count, get_keys, get_recs
    character(len=MAXLENLAB) :: CBuf
    
    get_count = present(count)
    get_keys = present(keys)
    get_recs = present(recs)
    if (present(load)) then
        do_load = load
    else
        do_load = .not.(get_count .or. get_keys .or. get_recs)
    end if

    ! Check if nothing to do
    if (.not.(get_count .or. get_keys .or. get_recs .or. do_load)) &
        return

    ! Check if loading has already been done
    if (gfaf%preloaded) then
        if (get_count) count = size(gfaf%label_keys)
        if (get_keys) keys = gfaf%label_keys
        if (get_recs) recs = gfaf%label_recs
    else
        ! Not present, we need to analyse the file.
        ! First run to find number of keys.
        nkeys = 0
        call gfaf%skip(-gfaf%n_base_recs)
        do
            call Rd_Labl(gfaf%unit, gfaf%version, CBuf, NI, NR, NTot, LenBuf, &
                         N1, N2, N3, N4, N5, IType, NRI, EOF, gfaf%int_size)
            if (EOF) exit
            nkeys = nkeys + 1
            gfaf%rec = gfaf%rec + 1
            call gfaf%skip((NTot+LenBuf-1)/LenBuf)
        end do
        ! Number of keys retrieved, allocate
        if (do_load) allocate(gfaf%label_keys(nkeys), gfaf%label_recs(nkeys))
        if (get_count) count = nkeys
        if (get_keys) allocate(keys(nkeys))
        if (get_recs) allocate(recs(nkeys))
        ! Early exit if we only needed number of keys
        if (.not.(do_load .or. get_keys .or. get_recs)) return
        ! True run, extract information
        call gfaf%skip(-gfaf%n_base_recs)
        nkeys = 0
        do
            call Rd_Labl(gfaf%unit, gfaf%version, CBuf, NI, NR, NTot, LenBuf, &
                         N1, N2, N3, N4, N5, IType, NRI, EOF, gfaf%int_size)
            if (EOF) exit
            nkeys = nkeys + 1
            gfaf%rec = gfaf%rec + 1
            if (do_load) then
                gfaf%label_keys(nkeys) = CBuf
                gfaf%label_recs(nkeys) = gfaf%rec
            end if
            if (get_keys) keys(nkeys) = CBuf
            if (get_recs) recs(nkeys) = gfaf%rec
            call gfaf%skip((NTot+LenBuf-1)/LenBuf)
        end do
    end if

end subroutine gfaf_read_keys

! ======================================================================

subroutine gfaf_skip_records(gfaf, skip)
    !! Skip a number of records in a FAF file.
    !!
    !! Skips a number of records, either starting from the current record
    !! (positive values) or starting from the beginning of the file (negative
    !! values).
    class(gfaf_parser), intent(inout) :: gfaf
        !! Gaussian FAF parser instance.
    integer, intent(in) :: skip
        !! Number of records to skip.

    integer :: i, offset

    if (skip < 0) then
        offset = abs(skip) - gfaf%rec
    else
        offset = skip
    end if

    if (offset < 0) then
        rewind gfaf%unit
        do i = 1, abs(skip)
            read(gfaf%unit)
        end do
        gfaf%rec = abs(skip)
    else if (offset > 0) then
        do i = 1, offset
            read(gfaf%unit)
        end do
        gfaf%rec = gfaf%rec + offset
    end if

end subroutine gfaf_skip_records

! ======================================================================

end module gfaf_io