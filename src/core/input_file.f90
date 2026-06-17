submodule (input) input_file
    !! Submodule containing the definition of procedures related to the
    !! File/Program instances.
    use string, only: locase
    use fchk_io, only: fchk_data, fchk_parser

    implicit none

contains

! ======================================================================
! PSEUDO-CONSTRUCTORS
! ======================================================================

module procedure init_data_file
    !! Initialize an instance of the DataFile class.
    !!
    !! Initialize a DataFile instance based on the input file and,
    !! optionally the ftype.
    !!
    !! The function tries to open the file and identifies the file type
    !! if not specified, and the program that has generated it, as well
    !! as the version.

    integer :: iu, ios
    logical :: exit_ok, exists, no_print
    type(ErrorHandle) :: err

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

    call dfile%error%init(exit_on_error=exit_ok, no_printing=no_print)

    ! check if file exists
    inquire(file=fname, exist=exists)
    if (.not.exists) then
        call dfile%error%raise_error('file', 'missing', 'file not found')
        return
    end if

    dfile%name = fname
    open(file=dfile%name, newunit=iu, action='read', iostat=ios)
    if (ios /= 0) then
        call dfile%error%raise_error('file', 'open', 'operation failed')
        return
    end if
    close(iu)

    if (present(ftype)) then
        dfile%type = alias_file_type(ftype)
        if (dfile%type == ' ') then
            call dfile%error%raise_error('file', 'type', &
                                         'unsupported file format')
            return
        end if
    else
        call err%init(.false., .true.)
        dfile%type = get_file_type(dfile%name, err=err)
        if (err%raised()) then
            call dfile%error%from(err)
            return
        end if
    end if

    dfile%prog = init_program_version(dfile%name, dfile%type, prog_name, &
        prog_major, prog_minor, no_prog_version_ok, exit_on_error, silent)

end procedure init_data_file

! ======================================================================

module procedure init_program_version
    !! Extracts information on program version from file.
    !!
    !! Extracts information about the program version from file `fname`.
    !! The filetype (`ftype`) must be set.
    !!
    !! It is assumed that the validity of the file has been checked, the
    !! function directly opens the file to parse it.
    !!
    !! @note
    !! If major version is provided, there is no automatic search at all.
    !! If `minor` is not set, it is assumed to be `NA`.
    !!
    !! If `no_version_ok` is set and no version is found, the version is set
    !! to `NA`.
    !!
    !! Since this may not make sense for some programs/formats, the routine
    !! accepts empty major version.
    !! @endnote

    integer :: pos
    logical :: exit_ok, has_progname, has_version, no_print
    character(len=:), allocatable :: ft, gvers
    character(len=100) :: msg
    type(fchk_parser) :: fchk
    type(fchk_data) :: fchk_db

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

    call prog%error%init(exit_on_error=exit_ok, no_printing=no_print)

    has_progname = present(prog_name)
    has_version = present(major)

    if (has_progname) prog%name = trim(prog_name)

    if (has_version) then
        if (len_trim(major) == 0) then
            prog%major = 'NA'
        else
            prog%major = trim(major)
            if (present(minor)) then
                if (len_trim(minor) == 0) then
                    prog%minor = 'NA'
                    prog%version = prog%major
                else
                    prog%minor = trim(minor)
                    prog%version = prog%major // ' ' // prog%minor
                end if
            else
                prog%minor = 'NA'
                prog%version = prog%major
            end if
        end if
    end if

    if (.not.has_progname .or. .not.has_version) then
        ft = alias_file_type(ftype)
        if (ft == ' ') then
            call prog%error%raise_error('file', 'type', 'unrecognized value')
            return
        end if
        ! Note
        select case(ft)
            case ('GFChk')
                if (.not.has_progname) prog%name = 'Gaussian'
                if (.not.has_version) then
                    fchk = fchk_parser(fname)
                    if (allocated(fchk%gaussian)) then
                        gvers = fchk%gaussian
                    else
                        fchk_db = fchk%get('Gaussian Version')
                        if (fchk_db%dtype == '0') then
                            if (no_version_ok) then
                                prog%major = 'NA'
                                prog%minor = 'NA'
                                prog%version = 'NA'
                            else
                                call prog%error%raise_error( &
                                    'data', 'missing', &
                                    'Gaussian version is missing')
                            end if
                        end if
                        gvers = fchk_db%sdata(1)
                    end if
                    pos = index(gvers, 'Rev')
                    if (pos == 0) then
                        call prog%error%raise_error( &
                            'data', 'structure', &
                            'unable to parse the Gaussian version string')
                        return
                    else
                        if (gvers(pos+3:pos+3) == '-') then
                            ! Structure of the type CDVRev-XXX
                            prog%minor = gvers(pos+4:)
                            prog%major = gvers(:pos-1)
                            if (prog%major /= 'CDV') then
                                call prog%error%raise_error( &
                                    'data', 'structure', &
                                    'unable to parse the Gaussian version &
                                    &string')
                                return
                            end if
                            if (prog%major(1:1) == 'C') prog%major(1:1) = 'G'
                        else
                            ! Structure of the type <arch>-GxxRevXXX
                            if (index(gvers(:pos-1), '-') == 0) then
                                call prog%error%raise_error( &
                                    'data', 'structure', &
                                    'unable to parse the Gaussian version &
                                    &string')
                                return
                            end if
                            prog%minor = gvers(pos+3:)
                            prog%major = gvers(pos-3:pos-1)
                        end if
                        prog%version = prog%major // ' ' // prog%minor
                    end if
                end if
            case default
                write(msg, &
                    '("Support of file type """,a,""" not yet implemented")') &
                    ft
                call prog%error%raise_error('file', 'type', &
                    'unsupported file type')
                return
        end select
    end if

end procedure init_program_version

! ======================================================================
! MODULE PROCEDURES
! ======================================================================

module procedure check_prog_version
    !! Check if version in ProgramInfo instance matches query.
    !!
    !! Returns a True/False if the version in a ProgramInfo instance,
    !! matches the major and optionally minor revision.
    !! The format of the version is software dependent.
    !! The procedure can extract it from:
    !!
    !! 1. a ProgramInfo instance directly, provided as `prog_info`.
    !! 2. a DataFile instance, provided as `file_data`.
    !!
    !! The procedure takes the first available in this order.
    class(ProgramInfo), pointer :: prog => null()

    if (present(prog_info)) then
        prog => prog_info
    else if (present(file_data)) then
        prog => file_data%prog
    else
        res = .false.
        return
    end if

    res = locase(trim(major)) == locase(prog%major)
    if (res .and. present(minor)) &
        res = locase(trim(minor)) == locase(prog%minor)

end procedure check_prog_version

! ======================================================================

module procedure get_datafile_name
    !! Gets name of the file in DataFile instance.
    !!
    !! Returns the filename stored in the DataFile instance.

    name = dfile%name

end procedure get_datafile_name

! ======================================================================

module procedure get_datafile_type
    !! Gets type of the ifle in DataFile instance.
    !!
    !! Returns the filetype stored in the DataFile instance.

    dtype = dfile%type

end procedure get_datafile_type

! ======================================================================

module procedure get_file_type
    !! Gets file type.
    !!
    !! Gets file type, based on the filename or the file content.

    integer :: ios, iu, pos
    logical :: always_read, do_read, force_check, req_read
    character(len=10) :: ftype_ext, ftype_file
    character(len=1024) :: line

    if (present(soft_check)) then
        force_check = .not.soft_check
    else
        force_check = .false.
    end if

    req_read = present(read_file)
    if (req_read) then
        always_read = read_file
    else
        always_read = .false.
    end if

    pos = index(fname, '.', back=.true.)
    ftype_ext = alias_file_type(fname(pos+1:), is_ext=.true.)

    if (ftype_ext /= ' ') then
        do_read = force_check .or. always_read
    else
        do_read = .not.(req_read .and. .not.always_read)
    end if

    ftype_file = ' '
    if (do_read) then
        open(file=fname, newunit=iu, action='read', iostat=ios)
        if (ios /= 0) then
            if (present(err)) then
                call err%raise_error('file', 'open', &
                                     'could not open file to find type.')
            else
                call run%error%raise_error('file', 'open', &
                                           'could not open file to find type.')
            end if
            return
        end if
        read(iu, '(a)', iostat=ios) line
        if (ios /= 0) goto 10  ! Empty file, bypass
        ! For log files, 'Entering Gaussian...' on first line, so check first
        if (line(:33) == ' Entering Gaussian System, Link 0') then
            ftype_file = 'GLog'
        else
            do
                read(iu, '(a)', iostat=ios) line
                if (ios == 0) then
                    if (line(:49) == &
                        'Route                                      C   N=' &
                        ) then
                        ftype_file = 'GFChk'
                        exit
                    end if
                end if
            end do
        close(iu)
        end if
    end if
 10 continue

    if (ftype_file /= ' ') then
        ftype = trim(ftype_file)
    else if (ftype_ext /= ' ') then
        ftype = trim(ftype_ext)
    else
        if (present(err)) then
            call err%raise_error('file', 'type', &
                                 'could not define the file type.')
        else
            call run%error%raise_error('file', 'type', &
                                       'could not define the file type.')
        end if
    end if

end procedure get_file_type

! ======================================================================

module procedure get_prog_name
    !! Returns the name of the program stored in a ProgramInfo instance.
    !!
    !! Returns the name of the program stored in a ProgramInfo instance.
    !! The procedure can extract it from:
    !!
    !! 1. a ProgramInfo instance directly, provided as `prog_info`.
    !! 2. a DataFile instance, provided as `file_data`.
    !!
    !! The procedure takes the first available in this order.

    if (present(prog_info)) then
        name = prog_info%name
    else if (present(file_data)) then
        name = file_data%prog%name
    else
        name = ' '
    end if

end procedure get_prog_name

! ======================================================================

module procedure get_prog_version
    !! Returns the version of the program stored in a ProgramInfo instance.
    !!
    !! Returns the version of the program stored in a ProgramInfo
    !! instance.  The procedure can extract it from:
    !!
    !! 1. a ProgramInfo instance directly, provided as `prog_info`.
    !! 2. a DataFile instance, provided as `file_data`.
    !!
    !! The procedure takes the first available in this order.
    !!
    !! NOTE: The format of the version is vendor dependent and may not
    !! be a good basis for a generic parsing.  `check_prog_version` is
    !! a bit more robust for this.

    if (present(prog_info)) then
        version = prog_info%version
    else if (present(file_data)) then
        version = file_data%prog%version
    else
        version = ' '
    end if

end procedure get_prog_version

! ======================================================================

module procedure get_prog_major
    !! Returns the major revision stored in a ProgramInfo instance.
    !!
    !! Returns the major revision of the program stored in a ProgramInfo
    !! instance.  The procedure can extract it from:
    !!
    !! 1. a ProgramInfo instance directly, provided as `prog_info`.
    !! 2. a DataFile instance, provided as `file_data`.
    !!
    !! The procedure takes the first available in this order.

    if (present(prog_info)) then
        version = prog_info%major
    else if (present(file_data)) then
        version = file_data%prog%major
    else
        version = ' '
    end if

end procedure get_prog_major

! ======================================================================

module procedure get_prog_minor
    !! Returns the minor revision stored in a ProgramInfo instance.
    !!
    !! Returns the minor revision of the program stored in a ProgramInfo
    !! instance.  The procedure can extract it from:
    !!
    !! 1. a ProgramInfo instance directly, provided as `prog_info`.
    !! 2. a DataFile instance, provided as `file_data`.
    !!
    !! The procedure takes the first available in this order.
    !!

    if (present(prog_info)) then
        version = prog_info%minor
    else if (present(file_data)) then
        version = file_data%prog%minor
    else
        version = ' '
    end if

end procedure get_prog_minor

! ======================================================================
! SUB-MODULE COMPONENTS
! ======================================================================

function alias_file_type(name, is_ext) result(ftype)
    !! Gets the conventional filetype name based on a generic name.
    !!
    !! Given a generic name of a file extension (without the leading .),
    !! the function returns the conventional filetype name to be used by
    !! other functions

    character(len=*), intent(in) :: name
    !! Name of the file type or file extension to check.
    logical, intent(in), optional :: is_ext
    !! If True, `name` refers to an extension, not a file type.
    character(len=:), allocatable :: ftype
    !! Conventional filetype name.

    logical :: from_ext

    if (present(is_ext)) then
        from_ext = is_ext
    else
        from_ext = .false.
    end if

    if (from_ext) then
        select case (locase(name))
            case ('fchk', 'fck', 'fch')
                ftype = 'GFChk'
            case ('log', 'out')
                ftype = 'GLog'
            case ('baf')
                ftype = 'GBAF'
            case ('faf')
                ftype = 'GFAF'
            case default
                ftype = ' '
        end select
    else
        select case (locase(name))
            case ('fch', 'fchk', 'gfchk')
                ftype = 'GFChk'
            case ('log', 'glog')
                ftype = 'GLog'
            case ('baf', 'gbaf')
                ftype = 'GBAF'
            case ('faf', 'gfaf')
                ftype = 'GFAF'
            case default
                ftype = ' '
        end select
    end if

end function alias_file_type

! ======================================================================

end submodule input_file
