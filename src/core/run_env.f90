module run_env
    !! Provide environment variables and parameters for program executions.
    !!
    !! The module maintains only 1 dependence, to module string, to facilitate
    !! basic operations.

    use iso_fortran_env, only: output_unit
    use string, only: locase

    implicit none

    integer, protected :: iu_out = output_unit
        !! Output unit.
    integer, protected :: verbosity = 0
        !! Default verbosity level.
    integer, private :: def_print_level = 1
        !! Default minimum printing levels for warning/errors.
    logical, private :: def_exit_error = .true.
        !! Default behavior of error instances to force quit on error.

    integer, parameter, private :: WARN_LVL = 1, ERROR_LVL = 3, &
        ERROR_LOW_LVL = 2, ERROR_DEV_LVL = -1

    type, public :: ErrorHandle
        !! Error type
        !!
        !! The type can manage printing and exiting.
        private
        integer :: level = 0
            !! Internal classifier of the error level.
            !! Negative values are used for unrecoverable errors, currently
            !! handled as non-recoverable errors.
            !!
            !! | code | meaning              |
            !! |-----:|----------------------|
            !! |  -1  | development          |
            !! |   0  | no error             |
            !! |   1  | warning              |
            !! |   2  | error (non-blocking) |
            !! |   3  | error                |
        integer :: category = 0
            !! Category of error.
            !! 10^1^ for the main category, 10^0^ for the sub-category
            !!
            !! | code | meaning                                  |
            !! |-----:|------------------------------------------|
            !! |    0 | generic                                  |
            !! |   10 | memory                                   |
            !! |   11 | memory - allocation                      |
            !! |   12 | memory - exceeded                        |
            !! |   20 | file - generic                           |
            !! |   21 | file - not found                         |
            !! |   22 | file - wrong type                        |
            !! |   23 | file - cannot open                       |
            !! |   24 | file - cannot close                      |
            !! |   25 | file - cannot read                       |
            !! |   26 | file - cannot write                      |
            !! |   27 | file - EOF reached                       |
            !! |   30 | keyword - generic                        |
            !! |   31 | keyword - not found                      |
            !! |   32 | keyword - input arguments                |
            !! |   33 | keyword - duplicate/equivalent keys      |
            !! |   40 | data - generic                           |
            !! |   41 | data - missing                           |
            !! |   42 | data - inconsistency                     |
            !! |   43 | data - incorrect number                  |
            !! |   44 | data - unknown structure                 |
            !! |   45 | data - unexpected/unsupported value      |
            !! |   46 | data - unexpected/unsupported type       |
            !! |   50 | value - generic                          |
            !! |   51 | value - failed conversion                |
            !! |   52 | value - incompatibility                  |
            !! |   53 | value - unset                            |
            !! |   54 | value - wrong type                       |
            !! |   55 | value - wrong/unexpected                 |
            !! |   60 | calc - generic                           |
            !! |   61 | calc - NaN/invalid operation             |
            !! |   62 | calc - singularity (or risk)             |
            !! |   63 | calc - inconsistency in results          |
            !! |   70 | option - generic                         |
            !! |   71 | option - unsupported option              |
            !! |   72 | option - unsupported value               |
            !! |   73 | option - conflicting options             |
            !! |   74 | option - missing option/specification    |
            !! |  -10 | dev - generic error                      |
            !! |  -11 | dev - wrong value in call                |
            !! |  -12 | dev - unsupported case                   |
            !! |  -13 | dev - feature NYI                        |
            !! |  -14 | dev - internal limit reached             |
            !! |  -20 | arg - generic problem                    |
            !! |  -21 | arg - wrong/unexpected value             |
            !! |  -22 | arg - wrong type                         |
            !! |  -23 | arg - empty content                      |
            !! |  -24 | arg - missing                            |
            !! |  -25 | arg - array size/number                  |
            !! |  -26 | arg - conflicts between given arguments  |
        ! integer :: label
        !     !! Internal classifier of the error type.
        !     !! -1: undefined
        !     !!  0: generic
        !     !!  1: runtime (standard execution)
        !     !!  2: development (error in call to routine)
        character(len=:), allocatable :: source
            !! Source of the error, typically the calling procedure.
        character(len=:), allocatable :: cause
            !! Cause of the error.
        character(len=:), allocatable :: details
            !! Details on the origin of the error.
        character(len=:), allocatable :: extra
            !! Extra message.
        logical, dimension(2) :: has_been_set = [.false., .false.]
            !! Keep track of error parameters set.
            !! 1. exit_on_error
            !! 2. minimum level of printing
        logical :: error_set = .false.
            !! An error status has been set.
        logical :: exit_on_error = .false.
            !! Exit if an error is found.
        logical :: has_printed = .false.
            !! Set to true if the warning/error message has been printed.
        integer :: print_min_level = 1
            !! Minimum level for printing.
        contains
            procedure :: init => error_init
            procedure, private :: reset => error_reset
            procedure, private :: set => error_set
            procedure, private :: finish => error_finish
            procedure :: from => error_copy_from_error
            procedure :: get_level => error_query_level
            procedure :: get_id => error_query_category
            procedure :: get_type => error_query_type
            procedure :: has_error => error_query_error
            procedure :: has_warning => error_query_warn
            procedure :: has_type => error_check_type
            procedure :: info => error_query_info
            procedure :: is_ok => error_query_ok
            procedure :: print => error_print
            procedure :: raise_error => error_raise_error
            procedure :: raise_warning => error_raise_warn
            procedure :: raise_argerror => error_raise_argerror
            procedure :: raise_deverror => error_raise_deverror
            procedure :: raise_generror => error_raise_generror
            procedure :: raised => error_query_raised
    end type ErrorHandle

    type, public :: CoreExecObject
        !! Basic object providing core functions
        !!
        !! The object can be inherited to provide core features to objects
        !! used in the ELEMENTS library
        type(ErrorHandle) :: error
    end type CoreExecObject

    type, private :: ExecHandle
        private
        integer :: is_set = 0
            !! Non-default parameters (cf. [run_set_params] for details).
        type(ErrorHandle), public :: &
            error = ErrorHandle( &
                exit_on_error=.true., &
                print_min_level=0, &
                has_been_set=[.true., .true.])
            !! Error handling at execution level.
    contains
        procedure :: set => run_set_params
        procedure :: check => run_check_error
    end type ExecHandle

    type(ExecHandle) :: run

contains

! ======================================================================
! TYPE-BOUND PROCEDURES
! ======================================================================

function error_check_type(err, to_check) result(query)
    !! Check if error type matches `to_check`.
    class(ErrorHandle), intent(in) :: err
        !! Instance of `error_base`.
    character(len=*), intent(in) :: to_check
        !! Error type to check.
    logical :: query
        !! Result of the check query.

    if (err%is_ok()) then
        query = .false.
    else
        select case (locase(trim(to_check)))
            case ('gen', 'generic')
                query = err%category/10 == 0
            case ('mem', 'memory')
                query = err%category/10 == 1
            case ('file')
                query = err%category/10 == 2
            case ('key', 'keyword')
                query = err%category/10 == 3
            case ('dat', 'data', 'qty', 'quantity')
                query = err%category/10 == 4
            case ('val', 'value')
                query = err%category/10 == 5
            case ('calc', 'math')
                query = err%category/10 == 6
            case ('opt', 'option')
                query = err%category/10 == 7
            case ('dev', 'devel', 'developer')
                query = err%category/10 == -1
            case ('arg')
                query = err%category/10 == -2
            case default
                query = .false.
        end select
    end if

end function error_check_type

! ======================================================================

subroutine error_copy_from_error(dest, src, check_and_raise)
    !! Copy error information from an existing error.
    !!
    !! Copies error information from an existing error.
    !! The core parameters (error level for printing and raising errors...)
    !! are not affected.
    !! By default, the system test if an error needs to be raised from
    !! these new data.
    class(ErrorHandle), intent(inout) :: dest
        !! Error instance to update.
    class(ErrorHandle), intent(in) :: src
        !! Source error instance.
    logical, intent(in), optional :: check_and_raise
        !! Check error level and raise if necessary.

    logical :: check

    if (present(check_and_raise)) then
        check = check_and_raise
    else
        check = .true.
    end if

    call dest%reset()
    dest%error_set = src%error_set
    dest%level = src%level
    dest%category = src%category
    if (allocated(src%source)) dest%source = src%source
    if (allocated(src%cause)) dest%cause = src%cause
    if (allocated(src%details)) dest%details = src%details
    if (allocated(src%extra)) dest%extra = src%extra

    if (check) call dest%finish()

end subroutine error_copy_from_error

! ======================================================================

subroutine error_finish(err, no_exit)
    !! Finalize error condition and interrupt if needed.
    class(ErrorHandle), intent(inout) :: err
        !! Instance of `error_base`.
    logical, intent(in), optional :: no_exit
        !! Do not exit on error, the default is set based on error level and
        !! the initial setup of the error.

    integer :: min_print
    logical :: do_exit, do_print, exit_def

    if (present(no_exit)) then
        do_exit = .not.no_exit
    else
        do_exit = .false.
    end if

    if (err%has_been_set(2)) then
        min_print = err%print_min_level
    else
        min_print = def_print_level
    end if
    do_print = err%level >= min_print
    if (do_print) call err%print()

    if (err%has_been_set(1)) then
        exit_def = err%exit_on_error
    else
        exit_def = def_exit_error
    end if
    if (exit_def .and. (err%level < 0 .or. err%level >= 3) &
            .or. do_exit) then
        if (.not.do_print) call err%print()
        stop err%level
    end if

end subroutine error_finish

! ======================================================================

subroutine error_init(err, exit_on_error, no_printing, print_level, &
                      force)
    !! Initialize the `error_base` instance.
    !!
    !! Initializes an `error_base` instance, setting basic parameters.
    class(ErrorHandle), intent(inout) :: err
        !! Instance of `error_base`.
    logical, intent(in), optional :: exit_on_error
        !! A raised error will cause the termination of the program.
    logical, intent(in), optional :: no_printing
        !! The error instance should not automatically print messages.
    character, intent(in), optional :: print_level
        !! Any error at or above chosen level will be printed.
    logical, intent(in), optional :: force
        !! Force reinitialization, even if set before.

    logical :: force_init

    if (.not.present(exit_on_error) .and. .not.present(no_printing) &
        .and. .not.present(print_level)) &
        return

    if (present(force)) then
        force_init = force
    else
        force_init = .false.
    end if

    if (present(exit_on_error)) then
        if (err%has_been_set(1) .and. .not.force_init) then
            call err%set(ERROR_DEV_LVL, &
                         'error instance already parametrized', &
                         details='`exit_on_error` already set.', &
                         source='internal')
            return
        end if
        err%exit_on_error = exit_on_error
        err%has_been_set(1) = .true.
    end if
    if (present(no_printing) .or. present(print_level)) then
        if (err%has_been_set(2) .and. .not.force_init) then
            call err%set(ERROR_DEV_LVL, &
                         'error instance already parametrized', &
                         details='`print_min_level` already set.', &
                         source='internal')
            return
        end if
        err%has_been_set(2) = .true.
    end if
    if (present(no_printing)) then
        err%print_min_level = 100
    else if (present(print_level)) then
        select case(locase(trim(print_level)))
        case ('warn', 'warning')
            err%print_min_level = WARN_LVL
        case ('error')
            err%print_min_level = ERROR_LVL
        case ('any error', 'anyerror')
            err%print_min_level = ERROR_LOW_LVL
        case default
            call err%set(ERROR_DEV_LVL, 'unrecognized error level')
            err%has_been_set(2) = .false.
            return
        end select
    end if

end subroutine error_init

! ======================================================================

subroutine error_print(err)
    !! Print error message.
    !!
    !! Print an error message based on the stored parameters
    !!
    !! @note
    !! This function bypasses the test on the minimum level for printing
    !! and will always print.
    !! @endnote
    class(ErrorHandle), intent(inout) :: err
    !! Instance of `error_base`.

    character(len=3) :: err_code
    character(len=:), allocatable :: err_cat
    logical :: has_details, has_extra, has_source

    if (err%level == ERROR_DEV_LVL) then
        err_code = 'dev'
        err_cat = 'n/a'
    else if (err%category/10 > 0) then
        err_code = 'cat'
        select case (err%category/10)
            case(1)
                err_cat = 'Memory'
            case(2)
                err_cat = 'File'
            case(3)
                err_cat = 'Keyword-related'
            case(4)
                err_cat = 'Data processing'
            case(5)
                err_cat = 'Value operations'
            case(6)
                err_cat = 'Calculations'
            case(7)
                err_cat = 'User options'
            case default
                err_cat = 'Unspecified'
        end select
    else
        err_code = 'gen'
        err_cat = 'n/a'
    end if

    has_details = allocated(err%details)
    has_extra = allocated(err%extra)
    has_source = allocated(err%source)
    if (has_details .and. has_extra .and. has_source) then
        call write_err(err_code, err%cause, err%details, err%extra, &
                       err%source, label=err_cat)
    else if (has_details .and. has_extra) then
        call write_err(err_code, err%cause, details=err%details, &
                       extra=err%extra, label=err_cat)
    else if (has_details .and. has_source) then
        call write_err(err_code, err%cause, details=err%details, &
                       source=err%source, label=err_cat)
    else if (has_extra .and. has_source) then
        call write_err(err_code, err%cause, extra=err%extra, &
                       source=err%source, label=err_cat)
    else if (has_details) then
        call write_err(err_code, err%cause, details=err%details, &
                       label=err_cat)
    else if (has_extra) then
        call write_err(err_code, err%cause, extra=err%extra, &
                       label=err_cat)
    else if (has_source) then
        call write_err(err_code, err%cause, source=err%source, &
                       label=err_cat)
    else
        call write_err(err_code, err%cause, label=err_cat)
    end if

    err%has_printed = .true.

end subroutine error_print

! ======================================================================

function error_query_category(err) result(cat)
    !! Return error category.
    class(ErrorHandle), intent(in) :: err
        !! Instance of `error_base`.
    integer :: cat
        !! Category of error.

    cat = err%category

end function error_query_category

! ======================================================================

function error_query_error(err) result(query)
    !! Check if the error level is: error
    class(ErrorHandle), intent(in) :: err
        !! Instance of `error_base`.
    logical :: query
        !! Result of the query

    query = err%level >= 2 .or. err%level < 0

end function error_query_error

! ======================================================================

subroutine error_query_info(err, cause, details, extra, source, &
                            msg_as_format, msg_multiline)
    !! Return error information.
    !!
    !! Returns information on error.
    !! `msg_as_format`, if present, contains the full message as Fortran
    !! format
    !! `msg_multiline` contains the message with C-style new line.
    class(ErrorHandle), intent(in) :: err
        !! Instance of `error_base`.
    character(len=:), allocatable, intent(out), optional :: cause
        !! Cause of the error.
    character(len=:), allocatable, intent(out), optional :: details
        !! Details on the error.
    character(len=:), allocatable, intent(out), optional :: extra
        !! Extra information on the error.
    character(len=:), allocatable, intent(out), optional :: source
        !! Source of the error.
    character(len=:), allocatable, intent(out), optional :: msg_as_format
        !! Full message as a Fortran compatible format.
    character(len=:), allocatable, intent(out), optional :: msg_multiline
        !! Full message as a string, using C-like newline characters.

    integer :: lstr
    character(len=10000) :: line

    if (present(cause)) cause = err%cause
    if (present(details) .and. allocated(err%details)) &
        details = err%details
    if (present(extra) .and. allocated(err%extra)) &
        extra = err%extra
    if (present(source) .and. allocated(err%source)) &
        source = err%source
    if (present(msg_as_format)) then
        if (allocated(err%source)) then
            line = '("Error encountered from ' // err%source // ': ' &
                // err%cause // '"'
            lstr = 25 + len(err%source) + 2 + len(err%cause) + 1
        else
            line = '("Error encountered: ' // err%cause // '"'
            lstr = 21 + len(err%cause) + 1
        end if
        if (allocated(err%details)) then
            line(lstr+1:) = ',/,"Reason: ' // err%details // '"'
            lstr = lstr + 12 + len(err%details) + 1
        end if
        if (allocated(err%extra)) then
            line(lstr+1:) = ',/,"Other information: ' // err%extra // '"'
            lstr = lstr + 23 + len(err%extra) + 1
        end if
        line(lstr+1:) = ')'
        lstr = lstr + 1
        allocate(character(len=lstr) :: msg_as_format)
        msg_as_format = line(:lstr)
    end if
    if (present(msg_multiline)) then
        if (allocated(err%source)) then
            line = 'Error encountered from ' // err%source // ': ' &
                // err%cause
            lstr = 23 + len(err%source) + 2 + len(err%cause)
        else
            line = 'Error encountered: ' // err%cause
            lstr = 19 + len(err%cause)
        end if
        if (allocated(err%details)) then
            line(lstr+1:) = new_line(line) // 'Reason: ' // err%details
            lstr = lstr + 9 + len(err%details)
        end if
        if (allocated(err%extra)) then
            line(lstr+1:) = new_line(line) // 'Other information: ' &
                // err%extra
            lstr = lstr + 20 + len(err%extra)
        end if
        allocate(character(len=lstr) :: msg_multiline)
        msg_multiline = line(:lstr)
    end if

end subroutine error_query_info

! ======================================================================

function error_query_level(err) result(level)
    !! Return error level.
    class(ErrorHandle), intent(in) :: err
        !! Instance of `error_base`.
    integer :: level
        !! Level of error

    level = err%level

end function error_query_level

! ======================================================================

function error_query_ok(err) result(query)
    !! Check if no error has been set.
    class(ErrorHandle), intent(in) :: err
        !! Instance of `error_base`.
    logical :: query
        !! Result of the query
    query = .not.err%error_set .or. err%level == 0

end function error_query_ok

! ======================================================================

function error_query_raised(err) result(query)
    !! Check if any level of error has been raised
    class(ErrorHandle), intent(in) :: err
        !! Instance of `error_base`.
    logical :: query
        !! Result of the query

    query = err%level > 0

end function error_query_raised

! ======================================================================

function error_query_type(err) result(query)
    !! Return error category.
    class(ErrorHandle), intent(in) :: err
        !! Instance of `error_base`.
    character(len=:), allocatable :: query
        !! Type of error.

    select case (err%category/10)
        case (1)
            query = 'mem'
        case (2)
            query = 'file'
        case (3)
            query = 'key'
        case (4)
            query = 'data'
        case (5)
            query = 'value'
        case (6)
            query = 'calc'
        case (7)
            query = 'opt'
        case (-1)
            query = 'dev'
        case (-2)
            query = 'arg'
        case default
            query = 'gen'
    end select

end function error_query_type

! ======================================================================

function error_query_warn(err) result(query)
    !! Check if the error level is: warning
    class(ErrorHandle), intent(in) :: err
        !! Instance of `error_base`.
    logical :: query
        !! Result of the query

    query = err%level == 1

end function error_query_warn

! ======================================================================

subroutine error_raise_argerror(err, op, cause, details, extra, source)
    !! Raise development-related error on input arguments.
    !!
    !! Sets parameters, messages and behavior for an error of level
    !! "DEVERROR", specific to call arguments.
    class(ErrorHandle), intent(inout) :: err
        !! Instance of `error_base`.
    character(len=*), intent(in) :: op
        !! Type of operation attempted: opening, reading, conversion...
        !! The error type may be kept empty (generic)
    character(len=*), intent(in) :: cause
        !! Cause of the error.
    character(len=*), intent(in), optional :: details
        !! Details on the error.
    character(len=*), intent(in), optional :: extra
        !! Extra information.
    character(len=*), intent(in), optional :: source
        !! Source of the error: name of the procedure, unit, method...

    integer :: cat_code

    cat_code = get_error_code('arg', op)

    call err%set(-2, cause, details, extra, cat_code, source)

    call err%finish()

end subroutine error_raise_argerror

! ======================================================================

subroutine error_raise_deverror(err, op, cause, details, extra, source)
    !! Raise development-related error.
    !!
    !! Sets parameters, messages and behavior for an error of level
    !! "DEVERROR".
    class(ErrorHandle), intent(inout) :: err
        !! Instance of `error_base`.
    character(len=*), intent(in) :: op
        !! Type of operation attempted: opening, reading, conversion...
        !! The error type may be kept empty (generic)
    character(len=*), intent(in) :: cause
        !! Cause of the error.
    character(len=*), intent(in), optional :: details
        !! Details on the error.
    character(len=*), intent(in), optional :: extra
        !! Extra information.
    character(len=*), intent(in), optional :: source
        !! Source of the error: name of the procedure, unit, method...

    integer :: cat_code

    cat_code = get_error_code('dev', op)

    call err%set(ERROR_DEV_LVL, cause, details, extra, cat_code, source)

    call err%finish()

end subroutine error_raise_deverror

! ======================================================================

subroutine error_raise_generror(err, cause, details, extra, source)
    !! Raise a generic, uncategorized error.
    !!
    !! Sets parameters, messages and behavior for a generic error with no
    !! specific category.
    class(ErrorHandle), intent(inout) :: err
        !! Instance of `error_base`.
    character(len=*), intent(in) :: cause
        !! Cause of the error.
    character(len=*), intent(in), optional :: details
        !! Details on the error.
    character(len=*), intent(in), optional :: extra
        !! Extra information.
    character(len=*), intent(in), optional :: source
        !! Source of the error: name of the procedure, unit, method...

    integer :: cat_code

    cat_code = get_error_code('gen', 'gen')

    call err%set(ERROR_LVL, cause, details, extra, cat_code, source)

    call err%finish()

end subroutine error_raise_generror

! ======================================================================

subroutine error_raise_error(err, cat, op, cause, details, extra, source, &
                             low_err, no_exit)
    !! Raise error of level "error".
    !!
    !! Sets parameters, messages and behavior for an error of level
    !! "ERROR".
    !! If `low_risk` is true, the error is likely recoverable and the
    !! job may still proceed.
    class(ErrorHandle), intent(inout) :: err
        !! Instance of `error_base`.
    character(len=*), intent(in) :: cat
        !! (Main) category of the error, in singular (except data):
        !! allocate file, keyword, data/quantity, value
        !! Unknown category are assumed to "generic".
    character(len=*), intent(in) :: op
        !! Type of operation attempted: opening, reading, conversion...
        !! The error type may be kept empty (generic)
    character(len=*), intent(in) :: cause
        !! Cause of the error.
    character(len=*), intent(in), optional :: details
        !! Details on the error.
    character(len=*), intent(in), optional :: extra
        !! Extra information.
    character(len=*), intent(in), optional :: source
        !! Source of the error: name of the procedure, unit, method...
    logical, intent(in), optional :: low_err
        !! The error appears recoverable, and is not critical.
    logical, intent(in), optional :: no_exit
        !! Do not exit on error, the default is set based on error level and
        !! the initial setup of the error.

    integer :: cat_code, level

    cat_code = get_error_code(cat, op)

    if (present(low_err)) then
        if (low_err) then
            level = 2
        else
            level = 3
        end if
    else
        level = 3
    end if

    call err%set(level, cause, details, extra, cat_code, source)

    call err%finish(no_exit)

end subroutine error_raise_error

! ======================================================================

subroutine error_raise_warn(err, cat, op, cause, details, extra, source)
    !! Raise error of a level "warning".
    !!
    !! Sets parameters, messages and behavior for an error of level
    !! "WARNING".
    class(ErrorHandle), intent(inout) :: err
        !! Instance of `error_base`.
    character(len=*), intent(in) :: cat
        !! (Main) category of the error, in singular (except data):
        !! allocate file, keyword, data/quantity, value
        !! Unknown category are assumed to "generic".
    character(len=*), intent(in) :: op
        !! Type of operation attempted: opening, reading, conversion...
        !! The error type may be kept empty (generic)
    character(len=*), intent(in) :: cause
        !! Cause of the error.
    character(len=*), intent(in), optional :: details
        !! Details on the error.
    character(len=*), intent(in), optional :: extra
        !! Extra information.
    character(len=*), intent(in), optional :: source
    !! Source of the error: procedure, unit, method...

    integer :: cat_code

    cat_code = get_error_code(cat, op)

    call err%set(1, cause, details, extra, cat_code, source)

    call err%finish()

end subroutine error_raise_warn

! ======================================================================

subroutine error_reset(err)
    !! Reset the error status.
    !!
    !! Resets the attributes of the error instance.
    class(ErrorHandle), intent(inout) :: err
        !! Instance of `error_base`.

    if (err%error_set) then
        err%error_set = .false.
        err%level = 0
        err%category = 0
        if (allocated(err%source)) deallocate(err%source)
        if (allocated(err%cause)) deallocate(err%cause)
        if (allocated(err%details)) deallocate(err%details)
        if (allocated(err%extra)) deallocate(err%extra)
    end if

end subroutine error_reset

! ======================================================================

subroutine error_set(err, lvl, cause, details, extra, cat, source)
    !! Set error parameters.
    !!
    !! General routine to set error parameters.
    !! This is a low-level routine, not intended to be called directly.
    class(ErrorHandle), intent(inout) :: err
        !! Instance of `error_base`.
    integer, intent(in) :: lvl
        !! Level of error.
    character(len=*), intent(in) :: cause
        !! Cause of the error.
    character(len=*), intent(in), optional :: details
        !! Details on the error.
    character(len=*), intent(in), optional :: extra
        !! Extra information.
    integer, intent(in), optional :: cat
        !! Category of the error.
    character(len=*), intent(in), optional :: source
        !! Source of the error: procedure, unit, method...

    call err%reset()
    err%level = lvl
    err%cause = trim(cause)
    if (present(details)) err%details = trim(details)
    if (present(extra)) err%extra = trim(extra)
    if (present(source)) err%source = trim(source)
    if (present(cat)) then
        err%category = cat
    else
        err%category = 0
    end if
    err%error_set = .true.

end subroutine error_set

! ======================================================================

subroutine run_check_error(this_run, err, msg)
    !! Check error status and exit if error met.
    !!
    !! Checks an error instance given in input and exits the program if
    !! an error has been met.
    class(ExecHandle), intent(inout) :: this_run
        !! ExecHandler instance.
    class(ErrorHandle), intent(inout), optional :: err
        !! Error instance to check.
    character(len=*), intent(in), optional :: msg
        !! Lead message to print before the error.

    character(len=:), allocatable :: fixed_msg

    1000 format(/,A,/,'# Details:')
    1001 format(/,'=> ',A)

    if (present(msg)) then
        fixed_msg = trim(msg)
        if (fixed_msg(len(fixed_msg):len(fixed_msg)) == '.') &
            fixed_msg = fixed_msg(:len(fixed_msg)-1)
    end if

    if (present(err)) then
        if (err%has_error()) then
            if (.not.err%has_printed) then
                if (present(msg)) then
                    write(iu_out, 1000) fixed_msg
                    call err%print()
                end if
            else if (present(msg)) then
                write(iu_out, 1001) fixed_msg
            end if
            stop err%level
        end if
    else
        if (this_run%error%has_error()) then
            if (.not.this_run%error%has_printed) then
                if (present(msg)) then
                    write(iu_out, 1000) fixed_msg
                    call this_run%error%print()
                end if
            else if (present(msg)) then
                write(iu_out, 1001) fixed_msg
            end if
            stop this_run%error%level
        end if
    end if

end subroutine run_check_error

! ======================================================================

subroutine run_set_params(this_run, unit_output, verbosity_level, &
                          print_level, exit_on_error)
    !! Set parameters to run a program.
    !!
    !! Sets core runtime parameters for a program execution.
    class(ExecHandle), intent(inout) :: this_run
        !! ExecHandle instance.
    integer, intent(in), optional :: unit_output
        !! Fortran unit for the general output.
    integer, intent(in), optional :: verbosity_level
        !! Verbosity level.
    integer, intent(in), optional :: print_level
        !! Level of printing for notes (0), warning (1), error (>= 2).
    logical, intent(in), optional :: exit_on_error
        !! Errors instance should force exit on full errors.

    if (present(unit_output)) then
        iu_out = unit_output
        this_run%is_set = ibset(this_run%is_set, 0)
    end if
    if (present(verbosity_level)) then
        verbosity = verbosity_level
        this_run%is_set = ibset(this_run%is_set, 1)
    end if
    if (present(print_level)) then
        def_print_level = print_level
        this_run%is_set = ibset(this_run%is_set, 2)
    end if
    if (present(exit_on_error)) then
        def_exit_error = exit_on_error
        this_run%is_set = ibset(this_run%is_set, 3)
    end if

end subroutine run_set_params

! ======================================================================
! MODULE PROCEDURES
! ======================================================================

function get_error_code(cat, op) result (code)
    !! Parse error category and operation and return the right code.
    character(len=*), intent(in) :: cat
        !! (Main) category of the error, in singular (except data):
        !! allocate file, keyword, data/quantity, value
        !! Unknown category are assumed to "generic".
    character(len=*), intent(in) :: op
        !! Type of operation attempted: opening, reading, conversion...
        !! The error type may be kept empty (generic)
    integer :: code
        !! Code for the category.


    integer :: cat_main, cat_sub

    select case (locase(trim(cat)))
        case ('dev')
            cat_main = -1
        case ('arg')
            cat_main = -2
        case ('mem', 'alloc')
            cat_main = 1
        case ('file')
            cat_main = 2
        case ('key', 'keyword')
            cat_main = 3
        case ('dat', 'data', 'qty', 'quantity')
            cat_main = 4
        case ('val', 'value')
            cat_main = 5
        case ('calc', 'math')
            cat_main = 6
        case ('opt', 'option')
            cat_main = 7
        case default
            cat_main = 0
    end select

    cat_sub = 0
    select case(cat_main)
        case (-1)  ! development
            select case(locase(trim(op)))
                case ('argval')
                    cat_sub = -1
                case ('case', 'cond', 'condition')
                    cat_sub = -2
                case ('version')
                    cat_sub = -3
                case ('nyi', 'notimplemented', 'notyet')
                    cat_sub = -4
                case ('lim', 'limit')
                    cat_sub = -5
            end select
        case (-2)  ! call arguments
            select case(locase(trim(op)))
                case ('val', 'value')
                    cat_sub = -1
                case ('type')
                    cat_sub = -2
                case ('empty')
                    cat_sub = -3
                case ('missing')
                    cat_sub = -4
                case ('number', 'size')
                    cat_sub = -5
                case ('conflict', 'conflicting')
                    cat_sub = -6
            end select
        case (1)  ! memory
            select case(locase(trim(op)))
                case ('alloc', 'allocate', 'allocation')
                    cat_sub = 1
                case ('out', 'exceed')
                    cat_sub = 2
            end select
        case (2)  ! file
            select case(locase(trim(op)))
                case ('not found', 'missing')
                    cat_sub = 1
                case ('type')
                    cat_sub = 2
                case ('open')
                    cat_sub = 3
                case ('close')
                    cat_sub = 4
                case ('read')
                    cat_sub = 5
                case ('write')
                    cat_sub = 6
                case ('eof', 'end')
                    cat_sub = 7
            end select
        case (3)  ! keyword
            select case(locase(trim(op)))
                case ('not found', 'missing')
                    cat_sub = 1
                case ('args', 'arguments', 'input')
                    cat_sub = 2
                case ('dupl', 'duplicate', 'equiv', 'equivalent')
                    cat_sub = 3
            end select
        case (4)  ! data
            select case(locase(trim(op)))
                case ('not found', 'missing')
                    cat_sub = 1
                case ('consistent', 'consistency', 'inconsistent', &
                      'inconsistency')
                    cat_sub = 2
                case ('excess', 'too many', 'num', 'number')
                    cat_sub = 3
                case ('struct', 'structure', 'unknown')
                    cat_sub = 4
                case ('val', 'value')
                    cat_sub = 5
                case ('type')
                    cat_sub = 6
            end select
        case (5)  ! value
            select case(locase(trim(op)))
                case ('conv', 'convert', 'conversion')
                    cat_sub = 1
                case ('incompatible', 'struct', 'structure')
                    cat_sub = 2
                case ('unset', 'not set')
                    cat_sub = 3
                case ('type')
                    cat_sub = 4
                case ('wrong')
                    cat_sub = 5
            end select
        case (6)  ! calculations/computations
            select case(locase(trim(op)))
                case ('nan', 'invalid')
                    cat_sub = 1
                case ('singularity', 'divergence')
                    cat_sub = 2
                case ('inconsistency')
                    cat_sub = 3
            end select
        case (7)  ! user-options
            select case(locase(trim(op)))
                case ('unknown', 'incorrect')
                    cat_sub = 1
                case ('value')
                    cat_sub = 2
                case ('conflict', 'conflicting')
                    cat_sub = 3
                case ('missing', 'not found')
                    cat_sub = 4
            end select
    end select

    code = cat_main*10 + cat_sub

end function get_error_code

! ======================================================================

subroutine write_err(nature, cause, details, extra, source, label)
    !! Write error message on default unit.
    !!
    !! Writes an error message.  The formatting depends on the nature
    !!   of the error:
    !! * Generic/gen: generic/basic error
    !! * Specific/cat: specific error (for a given category)
    !! * Developer/dev: coding/development error
    character(len=*), intent(in) :: nature
    !! Nature of the error.
    character(len=*), intent(in) :: cause
    !! Cause of the error.
    character(len=*), intent(in), optional :: details
    !! Details on the error.
    character(len=*), intent(in), optional :: extra
    !! Extra information.
    character(len=*), intent(in), optional :: source
    !! Source of the error: procedure, unit, method...
    character(len=*), intent(in), optional :: label
    !! Label for specific errors.

    character(len=:), allocatable :: err_label

    1000 format(/,'Error encountered: ',a)
    1001 format(/,'Error encountered in [',a,']: ',a)
    1002 format(/,'Internal error encountered: ',a)
    1003 format(/,'Unclassified error: ',a)
    1004 format(/,a,1x,'error encountered: ',a)
    1010 format('-- Reason: ',a)
    1020 format('-- Note: ',a)

    select case (locase(trim(nature)))
        case ('generic', 'gen', 'basic', 'std')
            write(iu_out, 1000) trim(cause)
            if (present(details)) write(iu_out, 1010) trim(details)
            if (present(extra)) write(iu_out, 1020) trim(extra)
        case ('specific', 'cat')
            if (present(label)) then
                err_label = trim(label)
            else
                err_label = ' '
            end if
            if (len_trim(err_label) > 0) then
                write(iu_out, 1004) err_label, trim(cause)
            else
                write(iu_out, 1000) trim(cause)
            end if
            if (present(details)) write(iu_out, 1010) trim(details)
            if (present(extra)) write(iu_out, 1020) trim(extra)
        case ('deverr', 'dev', 'arg')
            if (present(source)) then
                write(iu_out, 1001) trim(source), trim(cause)
            else
                write(iu_out, 1002) trim(cause)
            end if
            if (present(details)) write(iu_out, 1010) trim(details)
            if (present(extra)) write(iu_out, 1020) trim(extra)
        case default
            write(iu_out, 1003) trim(cause)
    end select

end subroutine write_err

! ======================================================================

end module run_env