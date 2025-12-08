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

    type, public :: error_handle
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
            !! | code | meaning                   |
            !! |-----:|---------------------------|
            !! |    0 | generic                   |
            !! |   10 | memory                    |
            !! |   11 | memory - allocation       |
            !! |   12 | memory - exceeded         |
            !! |   20 | file - generic            |
            !! |   21 | file - not found          |
            !! |   22 | file - wrong type         |
            !! |   23 | file - cannot open        |
            !! |   24 | file - cannot close       |
            !! |   25 | file - cannot read        |
            !! |   26 | file - cannot write       |
            !! |   27 | file - EOF reached        |
            !! |   30 | keyword - generic         |
            !! |   31 | keyword - not found       |
            !! |   40 | data - generic            |
            !! |   41 | data - missing            |
            !! |   42 | data - inconsistency      |
            !! |   43 | data - excess             |
            !! |   50 | value - generic           |
            !! |   51 | value - failed conversion |
            !! |  -10 | dev - generic error       |
            !! |  -11 | dev - wrong value in call |
            !! |  -12 | dev - unsupported case    |
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
            procedure :: raise_error => error_raise_error
            procedure :: raise_warning => error_raise_warn
            procedure :: raise_deverror => error_raise_deverror
            procedure :: has_error => error_query_error
            procedure :: has_warning => error_query_warn
            procedure :: is_ok => error_query_ok
            procedure :: info => error_get
            procedure :: print => error_print
    end type error_handle

    type, public :: base_obj
        !! Basic object providing core functions
        !!
        !! The object can be inherited to provide core features to objects
        !! used in the ELEMENTS library
        type(error_handle) :: error
    end type base_obj

    type, private :: run_handle
        private
        integer :: is_set = 0
            !! Non-default parameters (cf. [run_set_params] for details).
        type(error_handle), public :: &
            error = error_handle( &
                exit_on_error=.true., &
                print_min_level=0, &
                has_been_set=[.true., .true.])
            !! Error handling at execution level.
    contains
        procedure :: set => run_set_params
        procedure :: check => run_check_error 
    end type run_handle
    
    type(run_handle) :: run

contains

! ======================================================================

subroutine error_init(err, exit_on_error, no_printing, print_level, &
                      force)
    !! Initialize the `error_base` instance.
    !!
    !! Initializes an `error_base` instance, setting basic parameters.
    class(error_handle), intent(inout) :: err
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
            call err%set(-1, 'Error instance already parametrized', &
                         details='`exit_on_error` already set.', &
                         source='internal')
            return
        end if
        err%exit_on_error = exit_on_error
        err%has_been_set(1) = .true.
    end if
    if (present(no_printing) .or. present(print_level)) then
        if (err%has_been_set(2) .and. .not.force_init) then
            call err%set(-1, 'Error instance already parametrized', &
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
            err%print_min_level = 1
        case ('error')
            err%print_min_level = 3
        case ('any error', 'anyerror')
            err%print_min_level = 2
        case default
            call err%set(-1, 'Unrecognized error level')
            err%has_been_set(2) = .false.
            return
        end select
    end if

end subroutine error_init

! ======================================================================

subroutine error_reset(err)
    !! Reset the error status.
    !!
    !! Resets the attributes of the error instance.
    class(error_handle), intent(inout) :: err
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
    class(error_handle), intent(inout) :: err
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

subroutine error_finish(err, no_exit)
    !! Finalize error condition and interrupt if needed.
    class(error_handle), intent(inout) :: err
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

subroutine error_raise_deverror(err, op, cause, details, extra, source)
    !! Raise development-related error..
    !!
    !! Sets parameters, messages and behavior for an error of level
    !! "DEVERROR".
    class(error_handle), intent(inout) :: err
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

    call err%set(-1, cause, details, extra, cat_code, source)

    call err%finish()

end subroutine error_raise_deverror

! ======================================================================

subroutine error_raise_error(err, cat, op, cause, details, extra, source, &
                             low_err, no_exit)
    !! Raise error of level "error".
    !!
    !! Sets parameters, messages and behavior for an error of level
    !! "ERROR".
    !! If `low_risk` is true, the error is likely recoverable and the
    !! job may still proceed.
    class(error_handle), intent(inout) :: err
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
    class(error_handle), intent(inout) :: err
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

function error_query_ok(err) result(query)
    !! Check if no error has been set.
    class(error_handle), intent(in) :: err
        !! Instance of `error_base`.
    logical :: query
        !! Result of the query
    query = .not.err%error_set .or. err%level == 0

end function error_query_ok

! ======================================================================

function error_query_error(err) result(query)
    !! Check if the error level is: error
    class(error_handle), intent(in) :: err
        !! Instance of `error_base`.
    logical :: query
        !! Result of the query

    query = err%level >= 2 .or. err%level < 0

end function error_query_error

! ======================================================================

function error_query_warn(err) result(query)
    !! Check if the error level is: warning
    class(error_handle), intent(in) :: err
        !! Instance of `error_base`.
    logical :: query
        !! Result of the query

    query = err%level == 1

end function error_query_warn

! ======================================================================

subroutine error_get(err, cause, details, extra, source, msg_as_format, &
                     msg_multiline)
    !! Return error information.
    !!
    !! Returns information on error.
    !! `msg_as_format`, if present, contains the full message as Fortran
    !! format
    !! `msg_multiline` contains the message with C-style new line.
    class(error_handle), intent(in) :: err
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
    !! Full message as a Fortran compatible fortran.
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

end subroutine error_get

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
    class(error_handle), intent(inout) :: err
    !! Instance of `error_base`.

    character(len=3) :: err_code
    character(len=:), allocatable :: err_cat
    logical :: has_details, has_extra, has_source

    if (err%level == -1) then
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
            end select
        case (4)  ! data
            select case(locase(trim(op)))
                case ('not found', 'missing')
                    cat_sub = 1
                case ('consistent', 'consistency', 'inconsistent', &
                      'inconsistency')
                    cat_sub = 2
                case ('excess', 'too many')
                    cat_sub = 3
            end select
        case (5)  ! value
            select case(locase(trim(op)))
                case ('conv', 'convert', 'conversion')
                    cat_sub = 1
            end select
    end select

    code = cat_main*10 + cat_sub

end function get_error_code

! ======================================================================

subroutine run_check_error(this_run, err)
    !! Check error status and exit if error met.
    !!
    !! Checks an error instance given in input and exits the program if
    !! an error has been met.
    class(run_handle), intent(inout) :: this_run
        !! run_handler instance.
    class(error_handle), intent(inout), optional :: err
        !! Error instance to check.

    if (present(err)) then
        if (err%has_error()) then
            if (.not.err%has_printed) call err%print()
            stop err%level
        end if
    else
        if (this_run%error%has_error()) then
            if (.not.this_run%error%has_printed) call this_run%error%print()
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
    class(run_handle), intent(inout) :: this_run
        !! this_run_handler instance.
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
        case ('deverr', 'dev')
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