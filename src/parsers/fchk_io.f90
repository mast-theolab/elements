module fchk_io
    !! Gaussian Formatted checkpoint file Input/Ouput module.
    !!
    !! The module provides objects (derived types) to operate on Gaussian's
    !! Formatted checkpoint file (fchk).
    !!
    !! @todo
    !! * add int_kind, real_kind for more versatile parser
    !! * add possibility to store keys for faster search
    !! @endtodo
    use iso_fortran_env, only: int32, int64, real64
    use run_env, only: CoreExecObject, run
    use numeric, only: intwp

    private
    public :: fchk_data, fchk_parser

    integer, parameter :: IPOS_NEQ = 48, IPOS_TYPE = 44, LFMT_C = 12, &
        LFMT_H = 8, &
        NCOLS_C = 5, & ! Number of elements (columns) per line for character
        NCOLS_R = 5, & ! Number of elements (columns) per line for real
        NCOLS_I = 6, & ! Number of elements (columns) per line for integer
        NCOLS_L = 72, & ! Number of elements (columns) per line for logical
        NCOLS_H = 9  ! Number of elements (columns) per line for Hollerith
    integer, parameter, public :: LHEAD = 42

    type :: fchk_data
        character(len=:), allocatable :: key
        character(len=1) :: dtype
            !! Data type, as string.  Possible values are:
            !! * '0': unset (data was not extracted)
            !! * 'I': integer
            !! * 'R': real
            !! * 'C': complex
            !! * 'L': logical
            !! * 'S': character string
        integer :: size = 0
            !! Total number of elements.
        integer, dimension(:), allocatable :: idata
        real(real64), dimension(:), allocatable :: rdata
        complex(real64), dimension(:), allocatable :: cdata
        logical, dimension(:), allocatable :: ldata
        character(len=:), dimension(:), allocatable :: sdata
    end type fchk_data

    type, extends (CoreExecObject) :: fchk_obj
        !! Basic class to handle Gaussian FChk files.
        private
        character(len=:), allocatable :: name
            !! Filename.
        integer :: unit = 0
            !! Fortran unit.
        integer :: rec = 0
            !! Current record read/written.
        character(len=:), allocatable, public :: gaussian
            !! Gaussian version
    contains
        procedure :: filename => fchk_get_name
        procedure :: close => fchk_close

    end type fchk_obj

    type, extends(fchk_obj) :: fchk_parser
        private
        logical :: preloaded_keys = .false.
            !! If true, keys have been preloaded with their record positions.
        character(len=LHEAD), dimension(:), allocatable :: label_keys
            !! Store the label keys contained in FAF.
        integer, dimension(:), allocatable :: label_recs
            !! Store the label positions in FAF for faster retrieval.
    contains
        procedure, private :: fchk_read_item, fchk_read_items
        procedure, private :: read_data => fchk_read_data
        procedure :: skip => fchk_skip_records
        procedure :: keys => fchk_read_keys
        generic :: get => fchk_read_item, fchk_read_items
    end type fchk_parser

    interface fchk_parser
        module procedure init_fchk_parser
    end interface fchk_parser

contains

! ======================================================================

function init_fchk_parser(fname, skip_version, preload, exit_on_error, &
                          silent) result(fchk)
    !! Constructor-like function to create a fchk parser instance.
    character(len=*), intent(in) :: fname
        !! Name of the Fortran Array File.
    logical, intent(in), optional :: skip_version
        !! Skip version parsing if not present.
        !! Otherwise, raise an error if version is missing.
    logical, intent(in), optional :: preload
        !! Preload keys and their positions for faster searches.
    logical, intent(in), optional :: exit_on_error
        !! Exit if an error is encountered.  By default, `.true.`
    logical, intent(in), optional :: silent
        !! Do not print messages.  By default, `.false.`.
    type(fchk_parser) :: fchk
        !! Gaussian FChk parser instance.

    integer, parameter :: len_lab = 64, max_reclen = 10, ikind = kind(1)
    integer :: ios
    logical :: exists, exit_ok, get_version, no_print, preload_keys
    type(fchk_data) :: gvers

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

    call fchk%error%init(exit_on_error=exit_ok, no_printing=no_print)

    ! File operations
    if (present(preload)) then
        preload_keys = preload
    else
        preload_keys = .false.
    end if

    inquire(file=fname, exist=exists)
    if (.not.exists) then
        call fchk%error%raise_error('file', 'missing', 'file not found')
        return
    else
        allocate(character(len=len_trim(fname)) :: fchk%name)
        fchk%name = trim(fname)
        open(newunit=fchk%unit, file=fchk%name, iostat=ios, &
             form='formatted', action='read', status='old')
        if (ios /= 0) then
            call fchk%error%raise_error('file', 'open', 'operation failed')
            return
        end if
    end if

    if (preload_keys) then
        call fchk%keys()
    end if

    if (present(skip_version)) then
        get_version = .not.skip_version
    else
        get_version = .true.
    endif

    if (get_version) then
        gvers = fchk%get('Gaussian Version')
        if (gvers%dtype /= '0') then
            fchk%gaussian = gvers%sdata(1)
        else
            fchk%gaussian = 'Unknown'
        end if
    end if

end function init_fchk_parser

! ======================================================================

function fchk_close(fchk) result(status)
    !! Close connection to Gaussian FChk file.
    class(fchk_obj), intent(inout), target :: fchk
        !! Gaussian generic FChk file handler instance.
    logical :: status
        !! Return if file closed properly.

    integer :: ios
    character(len=256) :: msg

    close(fchk%unit, iostat=ios)
    status = ios == 0
    if (.not.status) then
        msg = ' '
        write(msg, '("Failed to close fchk file: ",a)') trim(fchk%name)
        call fchk%error%raise_error('file', 'close', msg)
    end if

end function fchk_close

! ======================================================================

function fchk_get_name(fchk) result(fname)
    !! Query name of the Gaussian FChk file.
    class(fchk_obj), intent(in), target :: fchk
        !! Gaussian generic FChk file handler instance.
    character(len=:), pointer :: fname
        !! Filename.
    fname => fchk%name
end function fchk_get_name

! ======================================================================

function fchk_read_data(fchk, line, label) result(dbase)
    !! Read 1 dataset from a Gaussian Formatted checkpoint file.
    !!
    !! Reads data associated to a label in the formatted checkpoint file.
    class(fchk_parser), intent(inout) :: fchk
        !! Gaussian FChk parser instance.
    character(len=*), intent(in) :: line
        !! Header line, to process.
    character(len=*), intent(in) :: label
        !! Label associated to data.
    type(fchk_data) :: dbase
        !! Extracted data base.

    integer :: i, nrecs, nvals
    character(len=64) :: sval, fmt

    dbase%key = trim(label)
    if (index(line(LHEAD+1:), 'N=') > 0) then
        read(line(LHEAD+1:), *) dbase%dtype, sval, nvals
        select case (dbase%dtype)
            case ('R')
                dbase%size = nvals
                allocate(dbase%rdata(dbase%size))
                read(fchk%unit, *) (dbase%rdata(i), i=1, dbase%size)
                nrecs = ceiling(dbase%size/real(NCOLS_R))
            case ('I')
                dbase%size = nvals
                allocate(dbase%idata(dbase%size))
                read(fchk%unit, *) (dbase%idata(i), i=1, dbase%size)
                nrecs = ceiling(dbase%size/real(NCOLS_I))
            case ('C')
                dbase%dtype = 'S'
                dbase%size = 1
                allocate(character(len=LFMT_C*nvals):: dbase%sdata(1))
                write(fmt, '("(a",i0,")")') NCOLS_C*LFMT_C
                nrecs = ceiling(nvals/real(NCOLS_C))
                do i = 1, nrecs
                    read(fchk%unit, fmt) dbase%sdata((i-1)*NCOLS_C*LFMT_C+1:)
                end do
            case ('H')
                dbase%dtype = 'S'
                dbase%size = nvals
                allocate(character(len=LFMT_H):: dbase%sdata(nvals))
                ! NOTE: we assume that each element is separated by at least
                ! one blank character to facilitate parsing
                ! If needed, the code below could be adapted to parse touching
                ! strings.
                ! write(fmt, '("(a",i0,")")') LFMT_H
                ! nrecs = ceiling(dbase%size/real(NCOLS_H))
                ! do i = 1, nrecs
                !     read(fchk%unit, *) (dbase%sdata(i), i=1, dbase%size)
                ! end do
                read(fchk%unit, *) (dbase%sdata(i), i=1, dbase%size)
                nrecs = ceiling(dbase%size/real(NCOLS_H))
            case ('L')
                dbase%size = nvals
                allocate(dbase%ldata(nvals))
                read(fchk%unit, *) (dbase%ldata(i), i=1, dbase%size)
                nrecs = ceiling(dbase%size/real(NCOLS_L))
            case default
                call fchk%error%raise_deverror( &
                    'case', 'Unsupported data type in FChk', &
                    source='fchk_read_data')
                return
        end select
        fchk%rec = fchk%rec + nrecs
    else
        read(line(LHEAD+1:), *) dbase%dtype, sval
        dbase%size = 1
        select case (dbase%dtype)
            case ('R')
                allocate(dbase%rdata(dbase%size))
                read(sval, *) dbase%rdata(1)
            case ('I')
                allocate(dbase%idata(dbase%size))
                read(sval, *) dbase%idata(1)
            case default
                call fchk%error%raise_deverror( &
                    'case', 'Unsupported data type in FChk', &
                    source='fchk_read_data')
                return
        end select
    end if

end function fchk_read_data

! ======================================================================

function fchk_read_item(fchk, label) result(dbase)
    !! Read 1 item in the FChk file.
    !!
    !! Reads 1 item corresponding to label in Gaussian FChk file.
    class(fchk_parser), intent(inout) :: fchk
        !! Gaussian Fchk parser instance.
    character(len=*), intent(in) :: label
        !! Label of the quantity of interest.
    type(fchk_data) :: dbase
        !! Extracted data.

    integer :: i, ios, irec, n
    character(len=256) :: line

    if (fchk%preloaded_keys) then
        irec = 0
        do i = 1, size(fchk%label_keys)
            if (label == fchk%label_keys(i)) then
                irec = fchk%label_recs(i)
                exit
            end if
        end do
        if (irec == 0) then
            dbase%dtype = '0'
            dbase%key = trim(label)
            return
        else
            call fchk%skip(-irec)
            read(fchk%unit, '(a)', iostat=ios) line
            fchk%rec = fchk%rec + 1
            dbase = fchk%read_data(line, label)
        end if
    else
        rewind(fchk%unit)
        fchk%rec = 0
        do
            read(fchk%unit, '(a)', iostat=ios) line
            fchk%rec = fchk%rec + 1
            if (ios < 0) then
                dbase%dtype = '0'
                dbase%key = trim(label)
                return
            end if
            if (line(1:LHEAD) == label) then
                dbase = fchk%read_data(line, label)
                exit
            else
                n = nrecs_fchk_block(line)
                if (n > 0) call fchk%skip(n)
            end if
        end do
    end if

end function fchk_read_item

! ======================================================================

function fchk_read_items(fchk, labels) result(dbase)
    ! Read multiple items in the FChk file.
    !!
    !! Reads several items corresponding to label in Gaussian FChk file.
    class(fchk_parser), intent(inout) :: fchk
        !! Gaussian FChk parser instance.
    character(len=*), dimension(:), intent(in) :: labels
        !! Labels of the quantities of interest.
    type(fchk_data), dimension(:), allocatable :: dbase
        !! Extracted data.

    integer :: i, ikey, ilab, imap, ios, nkeys, nlabs, nmap
    integer, dimension(:), allocatable :: irecs, map
    logical :: found
    character(len=256) :: line

    nlabs = size(labels)
    allocate(dbase(nlabs))

    if (fchk%preloaded_keys) then
        nkeys = size(fchk%label_keys)
        allocate(irecs(nlabs))
        irecs = 0

        ! Collect records of available labels/keys
        nmap = 0
        do ilab = 1, nlabs
            do ikey = 1, nkeys
                if (labels(ilab) == fchk%label_keys(ikey)) then
                    irecs(ilab) = fchk%label_recs(ikey)
                    nmap = nmap + 1
                    exit
                end if
            end do
        end do

        ! Now build a mapping to sorted keys so we go linearly within the file
        allocate(map(nmap))
        nmap = 0
        do ilab = 1, nlabs
            irec = irecs(ilab)
            if (irec > 0) then
                if (nmap == 0) then
                    map(1) = ilab
                    nmap = nmap + 1
                else
                    imap = 1
                    do while (irec > irecs(map(imap)))
                        imap = imap + 1
                        if (imap > nmap) exit
                    end do
                    map(imap:) = eoshift(map(imap:), -1)
                    map(imap) = ilab
                    nmap = nmap + 1
                end if
            else
                dbase(ilab)%dtype = '0'
                dbase(ilab)%key = trim(labels(ilab))
            end if
        end do

        ! Now that the mapping is done, read file.
        do imap = 1, nmap
            ilab = map(imap)
            irec = irecs(ilab)
            call fchk%skip(-irec)
            read(fchk%unit, '(a)', iostat=ios) line
            fchk%rec = fchk%rec + 1
            dbase(ilab) = fchk%read_data(line, labels(ilab))
        end do

    else ! no preloading
        allocate(map(nlabs))
        map = [(i, i = 1, nlabs)]
        rewind(fchk%unit)
        fchk%rec = 0
        do
            read(fchk%unit, '(a)', iostat=ios) line
            fchk%rec = fchk%rec + 1
            if (ios < 0) exit
            found = .false.
            do i = 1, nlabs
                ilab = map(i)
                if (line(1:LHEAD) == labels(ilab)) then
                    found = .true.
                    dbase(ilab) = fchk%read_data(line, labels(ilab))
                    map(i:) = eoshift(map(i:), 1)
                    nlabs = nlabs - 1
                    if (nlabs == 0) return
                    exit
                end if
            end do
            if (.not.found) then
                n = nrecs_fchk_block(line)
                if (n > 0) call fchk%skip(n)
            end if
        end do

        if (nlabs > 0) then
            do i = 1, nlabs
                ilab = map(i)
                dbase(ilab)%dtype = '0'
                dbase(ilab)%key = trim(labels(ilab))
            end do
        end if
    end if

end function fchk_read_items

! ======================================================================

subroutine fchk_read_keys(fchk, count, keys, recs, load)
    !! Read keys stored in FChk and basic information.
    !!
    !! The routine can be used to load information in the parser DB for faster
    !! search (default with no argument), and/or to retrieve information.
    !! If data are requested, the routine avoids reading unnecessarily the
    !! file.
    class(fchk_parser), intent(inout) :: fchk
        !! Gaussian FChk parser instance.
    integer, intent(out), optional :: count
        !! Number of labels stored in FAF.
    character(len=LHEAD), dimension(:), allocatable, intent(out), &
        optional :: keys
        !! List of label keys stored in FAF.
    integer, dimension(:), allocatable, intent(out), optional :: recs
        !! List of record positions
    logical, intent(out), optional :: load
        !! Load information in FChk parser, default if no other arguments given.

    integer :: ios, nkeys
    logical :: do_load, get_count, get_keys, get_recs
    character(len=256) :: line

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
    if (fchk%preloaded_keys) then
        if (get_count) count = size(fchk%label_keys)
        if (get_keys) keys = fchk%label_keys
        if (get_recs) recs = fchk%label_recs
    else
        ! Not present, we need to analyse the file.
        ! First run to find number of keys.
        nkeys = 0
        rewind(fchk%unit)
        fchk%rec = 0
        do
            read(fchk%unit, '(a)', iostat=ios) line
            fchk%rec = fchk%rec + 1
            if (ios < 0) exit
            nkeys = nkeys + 1
            call fchk%skip(nrecs_fchk_block(line))
        end do
        ! Number of keys retrieved, allocate
        if (do_load) allocate(fchk%label_keys(nkeys), fchk%label_recs(nkeys))
        if (get_count) count = nkeys
        if (get_keys) allocate(keys(nkeys))
        if (get_recs) allocate(recs(nkeys))
        ! Early exit if we only needed number of keys
        if (.not.(do_load .or. get_keys .or. get_recs)) return
        ! True run, extract information
        rewind(fchk%unit)
        fchk%rec = 0
        nkeys = 0
        do
            read(fchk%unit, '(a)', iostat=ios) line
            if (ios < 0) exit
            nkeys = nkeys + 1
            fchk%rec = fchk%rec + 1
            if (do_load) then
                fchk%label_keys(nkeys) = line(:LHEAD)
                fchk%label_recs(nkeys) = fchk%rec
            end if
            if (get_keys) keys(nkeys) = line(:LHEAD)
            if (get_recs) recs(nkeys) = fchk%rec
            call fchk%skip(nrecs_fchk_block(line))
        end do
        ! Once finished, confirm the keys have been loaded if it has been
        ! requested
        if (do_load) fchk%preloaded_keys = .true.
    end if

end subroutine fchk_read_keys

! ======================================================================

subroutine fchk_skip_records(fchk, skip)
    !! Skip a number of records in a FChk file.
    !!
    !! Skips a number of records, either starting from the current record
    !! (positive values) or starting from the beginning of the file (negative
    !! values).
    !!
    !! @note
    !! For negative values, the routine stops at the preceding record, assuming
    !! the target record is to be read after.
    !! @endnote
    class(fchk_parser), intent(inout) :: fchk
        !! Gaussian FChk parser instance.
    integer, intent(in) :: skip
        !! Number of records to skip.

    integer :: i, offset, ioff

    if (skip < 0) then
        offset = abs(skip) - fchk%rec
        ioff = -1
    else
        offset = skip
        ioff = 0
    end if

    if (offset < 0) then
        rewind fchk%unit
        do i = 1, abs(skip)+ioff
            read(fchk%unit, *)
        end do
        fchk%rec = abs(skip)+ioff
    else if (offset > 0) then
        do i = 1, offset+ioff
            read(fchk%unit, *)
        end do
        fchk%rec = fchk%rec + offset + ioff
    end if

end subroutine fchk_skip_records

! ======================================================================

function nrecs_fchk_block(line) result(nrecs)
    !! Returns the number of records occupied by a block.
    !!
    !! Returns the number of records occupied by a data block based on its
    !! header specifications.
    !! The returned value ignores the header line, so a scalar will occupy 0
    !! records.
    character(len=*), intent(in) :: line
        !! Header specification line.
    integer :: nrecs
        !! Number of records occupied by block

    integer :: ios, ncols, nvals

    if (line(IPOS_NEQ:IPOS_NEQ+1) == 'N=') then
        select case (line(IPOS_TYPE:IPOS_TYPE))
        case ('I')
            ncols = NCOLS_I
        case ('R')
            ncols = NCOLS_R
        case ('C')
            ncols = NCOLS_C
        case ('L')
            ncols = NCOLS_L
        case ('H')
            ncols = NCOLS_H
        case default
            print *, trim(line)
            call run%error%raise_deverror( &
                'case', 'Unsupported data type in FChk', &
                source='fchk_read_data')
                return
        end select
        read(line(IPOS_NEQ+2:), *, iostat=ios) nvals
        nrecs = (nvals+ncols-1)/ncols
    else
        nrecs = 0
    end if

end function nrecs_fchk_block

! ======================================================================

end module fchk_io