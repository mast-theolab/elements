submodule (input) input_file
    !! Submodule containing the definition of procedures related to the
    !! File/Program instances.
    use string, only: locase
    use parsefchk, only: fchkdata, fchkparser
    use exception, only: BaseException, Error, InitError, RaiseArgError, &
        RaiseFileError, RaiseKeyError

    implicit none

contains

! ======================================================================
! MODULE INTERFACE
! ======================================================================

module procedure init_file
    !! Initialize an instance of the DataFile class.
    !!
    !! Initialize a DataFile instance based on the input file and,
    !! optionally the ftype.
    !!
    !! The function tries to open the file and identifies the file type
    !! if not specified, and the program that has generated it, as well
    !! as the version.

    integer :: iu, ios
    logical :: exists
    class(BaseException), allocatable :: suberr

    ! initialize error handling
    file%error = InitError()

    ! check if file exists
    inquire(file=fname, exist=exists)
    if (.not.exists) then
        call RaiseFileError(file%error, fname, 'searching', 'File not found.')
        return
    end if

    file%name = fname
    open(file=file%name, newunit=iu, action='read', iostat=ios)
    if (ios /= 0) then
        call RaiseFileError(file%error, fname, 'opening', 'Operation failed')
        return
    end if
    close(iu)

    if (present(ftype)) then
        file%type = alias_file_type(ftype)
        if (file%type == ' ') then
            call RaiseArgError(file%error, 'ftype', 'Unsupported file format')
            return
        end if
    else
        file%type = get_file_type(file%name, err=suberr)
        if (suberr%raised()) then
            file%error = suberr
            return
        end if
    end if

    file%prog = get_program_version(file%name, file%type)

end procedure init_file

! ======================================================================

module procedure get_file_type
    !! Gets file type.
    !!
    !! Gets file type, based on the filename or the file content.

    integer :: ios, iu, pos
    logical :: do_read
    character(len=10) :: ftype_ext, ftype_file
    character(len=1024) :: line

    err = InitError()

    if (present(read_file)) then
        do_read = read_file
    else
        do_read = .true.
    end if

    pos = index(fname, '.', back=.true.)
    ftype_ext = alias_file_type(fname(pos+1:), is_ext=.true.)

    ftype_file = ' '
    if (do_read) then
        open(file=fname, newunit=iu, action='read', iostat=ios)
        if (ios /= 0) then
            call RaiseFileError(err, fname, 'opening', &
                                'Could not open file to find type.')
            return
        end if
        do
            read(iu, '(a)', iostat=ios) line
            if (ios == 0) then
                if (line(:49) == &
                    'Gaussian Version                           C   N=') then
                    ftype_file = 'GFChk'
                    exit
                else if (line(:33) == ' Entering Gaussian System, Link 0') then
                    ftype_file = 'GLog'
                end if
            end if
        end do
        close(iu)
    end if

    if (ftype_file /= ' ') then
        ftype = trim(ftype_file)
    else if (ftype_ext /= ' ') then
        ftype = trim(ftype_ext)
    else
        call RaiseKeyError(err, 'file type', 'analysing', &
                           'Could not define the file type.')
    end if

end procedure get_file_type

! ======================================================================

module procedure get_program_version
    !! Extracts information on program version from file.
    !!
    !! Extracts information about the program version from file `fname`.
    !! The filetype (`ftype`) must be set.
    !!
    !! It is assumed that the validity of the file has been checked, the
    !! function directly opens the file to parse it.

    integer :: iu, pos
    character(len=:), allocatable :: ft, key
    type(fchkparser) :: fchk
    type(fchkdata) :: fchk_db

    if (present(err)) err = InitError()

    ft = alias_file_type(ftype)
    if (ft == ' ') then
        if (present(err)) then
            call RaiseArgError(err, 'ftype', 'Unrecognized value.')
            return
        else
            print *, 'Unsupported ftype, stopping.'
            stop 1
        end if
    end if
    ! Note
    open(file=fname, newunit=iu, action='read', status='old')
    select case(ft)
        case ('GFChk')
            prog%name = 'Gaussian'
            fchk = fchkparser(fname)
            fchk_db = fchk%read('Gaussian Version')
            if (fchk_db%dtype == '0') then
                if (present(err)) then
                    call RaiseKeyError(err, 'Gaussian version', &
                        'looking for', &
                        'Could not determine Gaussian version, missing key.')
                    return
                else
                    print *, &
                        'Could not find the Gaussian version key. Stopping.'
                    stop 1
                end if
            end if
            pos = index(fchk_db%cdata, '-')
            if (pos == 0) then
                if (present(err)) then
                    call RaiseKeyError(err, 'Gaussian version', 'parsing', &
                                       'Unknown structure.')
                    return
                else
                    print *, 'Failed to parse the Gaussian version'
                    stop 1
                end if
            else
                key = trim(fchk_db%cdata(pos+1:))
                pos = index(key, 'Rev')
                if (pos == 0) then
                    if (present(err)) then
                        call RaiseKeyError(err, 'Gaussian version', &
                            'parsing', 'Unknown structure.')
                        return
                    else
                        print *, 'Failed to parse the Gaussian version'
                        stop 1
                    end if
                else
                    prog%major = key(:pos-1)
                    ! Take care of the fact of FCHK generated by C86DV
                    ! Assume it is GDV for simplicity
                    if (prog%major(1:1) == 'C') prog%major(1:1) = 'G'
                    prog%minor = key(pos+4:)
                    prog%version = prog%major // ' ' // prog%minor
                end if
            end if
        case default
            print '("Support of file type """,a,""" not yet implemented")', &
            ft
            stop 1
    end select
    close(iu)

end procedure get_program_version

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

    err = dfile%name

end procedure get_datafile_name

! ======================================================================

module procedure get_datafile_type
    !! Gets type of the ifle in DataFile instance.
    !!
    !! Returns the filetype stored in the DataFile instance.

    err = dfile%type

end procedure get_datafile_type

! ======================================================================

module procedure get_error_instance
    !! Gets error instance.
    !!
    !! Returns the error instance stored in the DataFile instance.

    err = dfile%error

end procedure get_error_instance

! ======================================================================

module procedure set_error_instance
    !! Sets error instance.
    !!
    !! Sets the error instance stored in the DataFile instance.

    dfile%error = err

end procedure set_error_instance

! ======================================================================

module procedure get_error_msg
    !! Gets error message.
    !!
    !! Returns the content of the error instance stored in the DataFile
    !! instance.

    err = dfile%error%msg()

end procedure get_error_msg

! ======================================================================

module procedure check_error_status
    !! Checks error status.
    !!
    !! Returns the error status of the internal component of a DataFile
    !! instance.

    raised = dfile%error%raised()

end procedure check_error_status

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