module exception
    !! Module for exception handling
    !!
    !! Module providing types for exception handling in other modules.
    use output, only: write_err
    use string, only: locase

    implicit none

    private
    public :: InitError, RaiseError, RaiseAllocateError, RaiseArgError, &
        RaiseFileError, RaiseKeyError, RaiseQuantityError, RaiseTermination, &
        RaiseValueError

    type, public, abstract :: BaseException
        private
        logical :: status = .False.
        character(len=1024) :: errmsg = ' '
    contains
        procedure(set_message), private, deferred :: set_msg
        procedure, private :: check_status, get_base_message
        generic :: raise => set_msg
        generic :: raised => check_status
        generic :: msg => get_base_message
    end type BaseException

    abstract interface
        subroutine set_message(this, msg)
            import BaseException
            class(BaseException) :: this
            character(len=*), intent(in), optional :: msg
        end subroutine set_message
    end interface

    ! Derived types
    type, public, extends(BaseException) :: NoError
    contains
        procedure, private :: set_msg => set_noerror
    end type NoError

    type, public, extends(BaseException) :: Termination
    contains
        procedure, private :: set_msg => set_term_message
    end type Termination

    type, public, extends(BaseException) :: Error
    contains
        procedure, private :: set_msg => set_error_message
    end type Error

    type, public, extends(BaseException) :: ValueError
    contains
        procedure, private :: set_msg => set_value_message
    end type ValueError

    type, public, extends(BaseException) :: AllocateError
    contains
        procedure, private :: set_msg => set_allocate_message
    end type AllocateError

    type, public, extends(BaseException) :: FileError
        character(len=256), public :: file = ' '
        character(len=64), public :: action = ' '
    contains
        procedure, private :: set_msg => set_file_message
    end type FileError

    type, public, extends(BaseException) :: KeyError
        character(len=256), public :: key = ' '
        character(len=64), public :: action = ' '
    contains
        procedure, private :: set_msg => set_key_message
    end type KeyError

    type, public, extends(BaseException) :: ArgumentError
        character(len=256), public :: arg = ' '
        character(len=512), public :: reason = ' '
    contains
        procedure, private :: set_msg => set_argerr_message
    end type ArgumentError

    type, public, extends(BaseException) :: QuantityError
        character(len=256), public :: qty = ' '
        character(len=512), public :: issue = ' '
    contains
        procedure, private :: set_msg => set_qtyerr_message
    end type QuantityError

    !! Add nature of error: argument/value...
    type, public :: ErrorHandler
        !! A type to manage error raised by procedures
        !!
        !! The type can manage printing and exiting.
        private
        integer :: level
        !! Internal classifier of the error level.
        !! 0: no error
        !! 1: warning
        !! 2: error (non-blocking)
        !! 3: error
        integer :: label
        !! Internal classifier of the error type.
        !! -1: undefined
        !!  0: generic
        !!  1: runtime (standard execution)
        !!  2: development (error in call to routine)
        character(len=:), allocatable :: source
        !! Source of the error, typically the calling procedure.
        character(len=:), allocatable :: cause
        !! Cause of the error.
        character(len=:), allocatable :: details
        !! Details on the origin of the error.
        character(len=:), allocatable :: extra
        !! Extra message.
        logical :: has_been_set = .false.
        !! The error parameters have been set.
        logical :: error_set = .false.
        !! An error status has been set.
        logical :: exit_on_error = .true.
        !! Exit if an error is found
        integer :: print_min_level = 1
        !! Minimum level for printing
        contains
            procedure :: setup => handler_setup
            procedure, private :: reset => handler_reset
            procedure, private :: set => handler_set_error
            procedure :: raise_error => handler_raise_error
            procedure :: raise_warning => handler_raise_warn
            procedure :: has_error => handler_query_error
            procedure :: has_warning => handler_query_warn
            procedure :: is_ok => handler_query_ok
            procedure :: info => handler_get_error
            procedure :: print => handler_print_error
    end type ErrorHandler

    type(ErrorHandler), public :: runstat

contains

! ======================================================================

function InitError() result(err)
    !! Initializer for the Error system
    !!
    !! The function creates an instance derived from `BaseException` set
    !!   to `NoError`.
    !! This should be typically used at the beginning of a procedure
    !!   that uses the error system provided by this module to initialize
    !!   the associated variable with an "OK" status for later testing.
    type(NoError), allocatable :: err

    allocate(err)

    return
end function InitError

! ======================================================================

subroutine RaiseError(err, msg)
    !! Sets a basic error
    !!
    !! Takes an instance of BaseException and sets it to the default
    !!   Error with the message chosen in input.
    class(BaseException), allocatable, intent(inout) :: err
    !! Original error, updated on return
    character(len=*), intent(in) :: msg
    !! Message to include in the error

    type(Error), allocatable :: newerr

    allocate(newerr)
    call newerr%raise(msg)

    if (allocated(err)) deallocate(err)
    err = newerr

    return
end subroutine RaiseError

! ======================================================================

subroutine RaiseArgError(err, arg, reason)
    !! Sets an argument error
    !!
    !! Takes an instance of BaseException and sets it to the
    !!   argument-related error, with the message built from the input
    !!   argument and optional motive.
    class(BaseException), allocatable, intent(inout) :: err
    !! Original error, updated on return
    character(len=*), intent(in) :: arg
    !! Argument that raised the error
    character(len=*), intent(in), optional :: reason
    !! Reason for the error related to `arg`

    type(ArgumentError), allocatable :: newerr

    allocate(newerr)
    newerr%arg = trim(arg)
    if (present(reason)) newerr%reason = trim(reason)
    call newerr%raise()

    if (allocated(err)) deallocate(err)
    err = newerr

    return
end subroutine RaiseArgError

! ======================================================================

subroutine RaiseTermination(err, msg)
    !! Set a termination signal.
    !!
    !! Takes an instance of BaseException and sets it to Termination
    !!   with an optional message set in input.
    class(BaseException), allocatable, intent(inout) :: err
    !! Original error, updated on return
    character(len=*), intent(in), optional :: msg
    !! Message to include in the error

    type(Termination), allocatable :: newerr

    allocate(newerr)
    if (present(msg)) then
        call newerr%raise(msg)
    else
        call newerr%raise(' ')
    end if

    if (allocated(err)) deallocate(err)
    err = newerr

    return
end subroutine RaiseTermination

! ======================================================================

subroutine RaiseValueError(err, msg, motive)
    !! Sets a value error.
    !!
    !! Takes an instance of BaseException and sets it to the
    !!   value-related error.
    class(BaseException), allocatable, intent(inout) :: err
    !! Original error, updated on return
    character(len=*), intent(in), optional :: msg
    !! General error message.
    character(len=*), intent(in), optional :: motive
    !! Complement or alternative message: reason why value is invalid.

    type(ValueError), allocatable :: newerr
    character(len=1024) :: new_msg

    allocate(newerr)
    if (present(msg) .and. present(motive)) then
        write(new_msg, '(a,a," Reason: ",a)') trim(msg), new_line(' '), &
            trim(motive)
    else if (present(msg)) then
        new_msg = trim(msg)
    else if (present(motive)) then
        write(new_msg, '("Reason: ",a)') trim(motive)
    else
        new_msg = ' '
    end if
    call newerr%raise(new_msg)

    if (allocated(err)) deallocate(err)
    err = newerr

end subroutine RaiseValueError

! ======================================================================

subroutine RaiseAllocateError(err, what, msg)
    !! Sets a memory allocation error.
    !!
    !! Takes an instance of BaseException and sets it to the
    !!   allocate-related error.
    class(BaseException), allocatable, intent(inout) :: err
    !! Original error, updated on return
    character(len=*), intent(in), optional :: what
    !! What was allocated.
    character(len=*), intent(in), optional :: msg
    !! Error message provided by the system.

    type(AllocateError), allocatable :: newerr
    character(len=1024) :: new_msg

    allocate(newerr)
    if (present(what) .and. present(msg)) then
        write(new_msg, '("Allocation of ",a," failed.",a,"Reason: ",a)') &
            trim(what), new_line(' '), trim(msg)
    else if (present(what)) then
        write(new_msg, '("Allocation of ",a," failed.")') trim(what)
    else if (present(msg)) then
        write(new_msg, '("Memory allocation failed. Reason: ",a)') trim(msg)
    else
        new_msg = ' '
    end if
    call newerr%raise(new_msg)

    if (allocated(err)) deallocate(err)
    err = newerr

end subroutine RaiseAllocateError

! ======================================================================

subroutine RaiseFileError(err, file, action, osmsg)
    !! Sets a file-related error.
    !!
    !! Takes an instance of BaseException and sets it to the
    !!   file-related error.
    !!
    !! Note: action should be as gerund (opening, closing...)
    class(BaseException), allocatable, intent(inout) :: err
    !! Original error, updated on return
    character(len=*), intent(in) :: file
    !! File name.
    character(len=*), intent(in) :: action
    !! Action applied to the file.
    character(len=*), intent(in), optional :: osmsg
    !! Optional OS message encountered.

    type(FileError), allocatable :: newerr

    allocate(newerr)
    if (file /= ' ') then
        newerr%file = '"' // trim(file) // '"'
    else
        newerr%file = '<unknown>'
    end if
    if (action /= ' ') then
        newerr%action = trim(action)
    else
        newerr%action = 'operating on'
    end if
    call newerr%raise(osmsg)

    if (allocated(err)) deallocate(err)
    err = newerr

end subroutine RaiseFileError

! ======================================================================

subroutine RaiseKeyError(err, key, action, reason)
    !! Sets a keyword/quantity-related error.
    !!
    !! Takes an instance of BaseException and sets it to the
    !!   key-related error.
    !!
    !! Note: action should be as gerund (parsing, searching...)
    class(BaseException), allocatable, intent(inout) :: err
    !! Original error, updated on return
    character(len=*), intent(in) :: key
    !! Keyword / quantity name.
    character(len=*), intent(in) :: action
    !! Action on the key (search/parse...).
    character(len=*), intent(in), optional :: reason
    !! Reason for the error.

    type(KeyError), allocatable :: newerr

    allocate(newerr)
    if (key /= ' ') then
        newerr%key = '"' // trim(key) // '"'
    else
        newerr%key = '<unknown>'
    end if
    if (action /= ' ') then
        newerr%action = trim(action)
    else
        newerr%action = 'looking for'
    end if
    call newerr%raise(reason)

    if (allocated(err)) deallocate(err)
    err = newerr

end subroutine RaiseKeyError

! ======================================================================

subroutine RaiseQuantityError(err, qty, msg)
    !! Sets an quantity-related error
    !!
    !! Takes an instance of BaseException and sets it to the
    !!   quantity-related error, with the message built from the input
    !!   quantity and optional motive.
    class(BaseException), allocatable, intent(inout) :: err
    !! Original error, updated on return
    character(len=*), intent(in), optional :: qty
    !! Argument that raised the error
    character(len=*), intent(in), optional :: msg
    !! msg for the error related to `qty`

    type(QuantityError), allocatable :: newerr

    allocate(newerr)
    if (present(qty)) newerr%qty = trim(qty)
    if (present(msg)) newerr%issue = trim(msg)
    call newerr%raise()

    if (allocated(err)) deallocate(err)
    err = newerr

    return
end subroutine RaiseQuantityError

! ======================================================================

function check_status(this) result(stat)
    !! Check status of an exception
    !!
    !! Checks the status of an exception and returns True if raised
    class(BaseException) :: this
    logical :: stat

    stat = this%status

    return
end function check_status

! ======================================================================

function get_base_message(this) result(msg)
    !! Check status of an exception
    !!
    !! Checks the status of an exception and returns True if raised
    class(BaseException) :: this
    character(len=len_trim(this%errmsg)) :: msg

    msg = trim(this%errmsg)

    return
end function get_base_message

! ======================================================================

subroutine set_status(this, msg)
    !! Sets status and error message of an exception
    !!
    !! Sets a message and changes the status of an exception
    class(BaseException) :: this
    character(len=*), intent(in) :: msg

    this%status = .True.
    this%errmsg = trim(msg)

    return
end subroutine set_status

! ======================================================================

subroutine set_noerror(this, msg)
    !! No-error case
    !!
    !! Sets a message and makes sure the status is not raised.
    class(NoError) :: this
    character(len=*), intent(in), optional :: msg

    if (present(msg)) then
        this%errmsg = trim(msg)
    else
        this%errmsg = ' '
    end if
    this%status = .False.

    return
end subroutine set_noerror

! ======================================================================

subroutine set_error_message(this, msg)
    !! Sets message of a standard error.
    !!
    !! Sets a message and changes the status of a basic error
    class(Error) :: this
    character(len=*), intent(in), optional :: msg

    character(len=1024) :: new_msg

    if (present(msg)) then
        call set_status(this, msg)
    else
        new_msg = 'Error found.'
        call set_status(this, new_msg)
    end if

    return
end subroutine set_error_message

! ======================================================================

subroutine set_allocate_message(this, msg)
    !! Sets error message for an AllocateError exception.
    !!
    !! Sets the error message and updates the status of an exception.
    !! The error message can be set in two ways:
    !! - by providing directly the msg (no different from BaseException)
    !! - through the attribute `msg`
    class(AllocateError) :: this
    character(len=*), intent(in), optional :: msg

    character(len=1024) :: new_msg

    if (present(msg)) then
        call set_status(this, msg)
    else
        new_msg = 'Memory allocation failed.'
        call set_status(this, new_msg)
    end if

    return
end subroutine set_allocate_message

! ======================================================================

subroutine set_value_message(this, msg)
    !! Sets error message for a ValueError exception.
    !!
    !! Sets the error message and updates the status of an exception.
    !! The error message can be set in two ways:
    !! - by providing directly the msg (no different from BaseException)
    !! - through the attribute `msg
    class(ValueError) :: this
    character(len=*), intent(in), optional :: msg

    character(len=1024) :: new_msg

    if (present(msg)) then
        call set_status(this, msg)
    else
        new_msg = 'Incorrect value.'
        call set_status(this, new_msg)
    end if

    return
end subroutine set_value_message

! ======================================================================

subroutine set_file_message(this, msg)
    !! Sets error message for a FileError exception.
    !!
    !! Sets the error message and updates the status of an exception.
    !! The first part of the message is built from internal attributes
    !! `file` and `action`.
    !! An optional second part (the reason) is included if `msg` is
    !! provided.
    class(FileError) :: this
    character(len=*), intent(in), optional :: msg

    character(len=1024) :: new_msg

    if (present(msg)) then
        write(new_msg, '("Error while ",a," file ",a,".",a,"Reason: ",a)') &
            trim(this%action), trim(this%file), new_line(' '), trim(msg)
    else
        write(new_msg, '("Error encountered while ",a," file ",a,".")') &
            trim(this%action), trim(this%file)
    end if

    call set_status(this, new_msg)

    return
end subroutine set_file_message

! ======================================================================

subroutine set_key_message(this, msg)
    !! Sets error message for a KeyError exception.
    !!
    !! Sets the error message and updates the status of an exception.
    !! The first part of the message is built from internal attributes
    !! `key` and `action`.
    !! An optional second part (the reason) is included if `msg` is
    !! provided.
    class(KeyError) :: this
    character(len=*), intent(in), optional :: msg

    character(len=1024) :: new_msg

    if (present(msg)) then
        write(new_msg, '("Error while ",a,1x,a,".",a,"Reason: ",a)') &
            trim(this%action), trim(this%key), new_line(' '), trim(msg)
    else
        write(new_msg, '("Error encountered while ",a,1x,a,".")') &
            trim(this%action), trim(this%key)
    end if

    call set_status(this, new_msg)

    return
end subroutine set_key_message

! ======================================================================

subroutine set_argerr_message(this, msg)
    !! Set error message for an ArgumentError exception.
    !!
    !! Sets the error message and updates the status of an exception.
    !! The error message can be set in two ways:
    !! - by providing directly the msg (no different from BaseException)
    !! - through the attributes `arg` and `reason`
    class(ArgumentError) :: this
    character(len=*), intent(in), optional :: msg

    character(len=1024) :: new_msg

    1000 format('Error found with argument "',a,'"')
    1001 format('Error in argument "',a,'": ',a)
    if (present(msg)) then
        call set_status(this, msg)
    else
        if (this%arg /= ' ') then
            if (this%reason /= ' ') then
                write(new_msg, 1001) trim(this%arg), trim(this%reason)
            else
                write(new_msg, 1000) trim(this%arg)
            end if
        else
            new_msg = 'Error within the input arguments'
        end if
        call set_status(this, new_msg)
    end if

    return
end subroutine set_argerr_message

! ======================================================================

subroutine set_term_message(this, msg)
    !! Set message of a normal termination.
    !!
    !! Sets a message and changes the status to a special-type error
    !!   for normal termination.
    class(Termination) :: this
    character(len=*), intent(in), optional :: msg

    if (present(msg)) then
        call set_status(this, msg)
    else
        call set_status(this, ' ')
    end if

    return
end subroutine set_term_message

! ======================================================================

subroutine set_qtyerr_message(this, msg)
    !! Set error message for an QuantityError exception.
    !!
    !! Sets the error message and updates the status of an exception.
    !! The error message can be set in two ways:
    !! - by providing directly the msg (no different from BaseException)
    !! - through the attributes `arg` and `reason`
    class(QuantityError) :: this
    character(len=*), intent(in), optional :: msg

    character(len=1024) :: new_msg

    1000 format('Error found while processing quantity "',a,'"')
    1001 format('Error with quantity "',a,'": ',a)
    if (present(msg)) then
        call set_status(this, msg)
    else
        if (this%qty /= ' ') then
            if (this%issue /= ' ') then
                write(new_msg, 1001) trim(this%qty), trim(this%issue)
            else
                write(new_msg, 1000) trim(this%qty)
            end if
        else if (this%issue /= ' ') then
            new_msg = trim(this%issue)
        else
            new_msg = 'Error while processing some quantities'
        end if
        call set_status(this, new_msg)
    end if

    return
end subroutine set_qtyerr_message

! ======================================================================

subroutine handler_setup(this, exit_on_error, no_printing, print_level, &
                         force)
    !! Initialize the ErrorHandler instance
    !!
    !! Initializes an ErrorHandler instance, setting basic parameters.
    class(ErrorHandler), intent(inout) :: this
    !! Instance of ErrorHandler.
    logical, intent(in), optional :: exit_on_error
    !! A raised error will cause the termination of the program.
    logical, intent(in), optional :: no_printing
    !! The error instance should not automatically print messages.
    character, intent(in), optional :: print_level
    !! Any error at or above chosen level will be printed.
    logical, intent(in), optional :: force
    !! Force reinitialization, even if set before.

    logical :: force_setup

    if (.not.present(exit_on_error) .and. .not.present(no_printing) &
        .and. .not.present(print_level)) &
        return

    if (present(force)) then
        force_setup = force
    else
        force_setup = .false.
    end if

    if (this%has_been_set .and. force_setup) then
        call this%raise_error('Error instance already parametrized', &
                              cat='dev')
        return
    end if

    if (present(exit_on_error)) this%exit_on_error = exit_on_error
    if (present(no_printing)) then
        this%print_min_level = 100
    else if (present(print_level)) then
        select case(locase(trim(print_level)))
        case ('warn', 'warning')
            this%print_min_level = 1
        case ('error')
            this%print_min_level = 3
        case ('any error', 'anyerror')
            this%print_min_level = 2
        case default
            call this%raise_error( &
                'Unrecognized error level', cat='dev')
            return
        end select
    end if
    this%has_been_set = .true.

end subroutine handler_setup

! ======================================================================

subroutine handler_reset(this)
    !! Reset the error status.
    !!
    !! Resets the attributes of the error instance.
    class(ErrorHandler), intent(inout) :: this
    !! Instance of ErrorHandler.

    if (this%error_set) then
        this%error_set = .false.
        if (allocated(this%source)) deallocate(this%source)
        if (allocated(this%cause)) deallocate(this%cause)
        if (allocated(this%details)) deallocate(this%details)
        if (allocated(this%extra)) deallocate(this%extra)
        this%level = 0
        this%label = 0
    end if

end subroutine handler_reset

! ======================================================================

subroutine handler_set_error(this, level, cause, details, extra, cat, source)
    !! Set error parameters.
    !!
    !! General routine to set error parameters.
    class(ErrorHandler), intent(inout) :: this
    !! Instance of ErrorHandler.
    integer, intent(in) :: level
    !! Level of error.
    character(len=*), intent(in) :: cause
    !! Cause of the error.
    character(len=*), intent(in), optional :: details
    !! Details on the error.
    character(len=*), intent(in), optional :: extra
    !! Extra information.
    character(len=*), intent(in), optional :: cat
    !! Category of the error.
    character(len=*), intent(in), optional :: source
    !! Source of the error: procedure, unit, method...

    call this%reset()
    this%level = level
    this%cause = trim(cause)
    if (present(details)) this%details = trim(details)
    if (present(extra)) this%extra = trim(extra)
    if (present(source)) this%source = trim(source)
    if (present(cat)) then
        select case (locase(cat(:3)))
        case ('run', 'std')
            this%label = 1
        case ('dev')
            this%label = 2
        case ('n/a', 'und')
            this%label = -1
        case default
            this%label = 0
        end select
    else
        this%label = 1
    end if
    this%error_set = .true.

end subroutine handler_set_error

! ======================================================================

subroutine handler_raise_error(this, cause, details, extra, cat, source, &
                               low_risk, no_exit)
    !! Raise error of level "error".
    !!
    !! Sets parameters, messages and behavior for an error of level
    !! "ERROR".
    !! If `low_risk` is true, the error is likely recoverable and the
    !! job may still proceed.
    class(ErrorHandler), intent(inout) :: this
    !! Instance of ErrorHandler.
    character(len=*), intent(in) :: cause
    !! Cause of the error.
    character(len=*), intent(in), optional :: details
    !! Details on the error.
    character(len=*), intent(in), optional :: extra
    !! Extra information.
    character(len=*), intent(in), optional :: cat
    !! Category of the error.
    character(len=*), intent(in), optional :: source
    !! Source of the error: procedure, unit, method...
    logical, intent(in), optional :: low_risk
    !! The risk posed by the error is low.
    logical, intent(in), optional :: no_exit
    !! Do not exit on error.

    integer :: level
    logical :: do_exit

    if (present(low_risk)) then
        if (low_risk) then
            level = 2
        else
            level = 3
        end if
    else
        level = 3
    end if
    call this%set(level, cause, details, extra, cat, source)

    if (this%level >= this%print_min_level) call this%print()

    do_exit = this%exit_on_error .and. this%level >= 3
    if (present(no_exit)) do_exit = do_exit .and. .not.no_exit

    if (do_exit) stop this%level

end subroutine handler_raise_error

! ======================================================================

subroutine handler_raise_warn(this, cause, details, extra, cat, source)
    !! Raise error of a level "warning".
    !!
    !! Sets parameters, messages and behavior for an error of level
    !! "WARNING".
    class(ErrorHandler), intent(inout) :: this
    !! Instance of ErrorHandler.
    character(len=*), intent(in) :: cause
    !! Cause of the error.
    character(len=*), intent(in), optional :: details
    !! Details on the error.
    character(len=*), intent(in), optional :: extra
    !! Extra information.
    character(len=*), intent(in), optional :: cat
    !! Category of the error.
    character(len=*), intent(in), optional :: source
    !! Source of the error: procedure, unit, method...

    call this%set(1, cause, details, extra, cat, source)

    if (this%level >= this%print_min_level) call this%print()

end subroutine handler_raise_warn

! ======================================================================

function handler_query_ok(this) result(query)
    !! Check if no error has been set.
    class(ErrorHandler), intent(in) :: this
    !! Instance of ErrorHandler.
    logical :: query
    !! Result of the query
    query = .not.this%error_set .or. this%level == 0
end function handler_query_ok

! ======================================================================

function handler_query_error(this) result(query)
    !! Check if the error level is: error
    class(ErrorHandler), intent(in) :: this
    !! Instance of ErrorHandler.
    logical :: query
    !! Result of the query
    query = this%level >= 2
end function handler_query_error

! ======================================================================

function handler_query_warn(this) result(query)
    !! Check if the error level is: warning
    class(ErrorHandler), intent(in) :: this
    !! Instance of ErrorHandler.
    logical :: query
    !! Result of the query
    query = this%level == 1
end function handler_query_warn

! ======================================================================

subroutine handler_get_error(this, cause, details, extra, source, &
                             msg_as_format, msg_multiline)
    !! Return error information.
    !!
    !! Returns information on error.
    !! `msg_as_format`, if present, contains the full message as Fortran
    !! format
    !! `msg_multiline` contains the message with C new line.
    class(ErrorHandler), intent(in) :: this
    !! Instance of ErrorHandler.
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

    if (present(cause)) cause = this%cause
    if (present(details) .and. allocated(this%details)) &
        details = this%details
    if (present(extra) .and. allocated(this%extra)) &
        extra = this%extra
    if (present(source) .and. allocated(this%source)) &
        source = this%source
    if (present(msg_as_format)) then
        if (allocated(this%source)) then
            line = '("Error encountered from ' // this%source // ': ' &
                // this%cause // '"'
            lstr = 25 + len(this%source) + 2 + len(this%cause) + 1
        else
            line = '("Error encountered: ' // this%cause // '"'
            lstr = 21 + len(this%cause) + 1
        end if
        if (allocated(this%details)) then
            line(lstr+1:) = ',/,"Reason: ' // this%details // '"'
            lstr = lstr + 12 + len(this%details) + 1
        end if
        if (allocated(this%extra)) then
            line(lstr+1:) = ',/,"Other information: ' // this%extra // '"'
            lstr = lstr + 23 + len(this%extra) + 1
        end if
        line(lstr+1:) = ')'
        lstr = lstr + 1
        allocate(character(len=lstr) :: msg_as_format)
        msg_as_format = line(:lstr)
    end if
    if (present(msg_multiline)) then
        if (allocated(this%source)) then
            line = 'Error encountered from ' // this%source // ': ' &
                // this%cause
            lstr = 23 + len(this%source) + 2 + len(this%cause)
        else
            line = 'Error encountered: ' // this%cause
            lstr = 19 + len(this%cause)
        end if
        if (allocated(this%details)) then
            line(lstr+1:) = new_line(line) // 'Reason: ' // this%details
            lstr = lstr + 9 + len(this%details)
        end if
        if (allocated(this%extra)) then
            line(lstr+1:) = new_line(line) // 'Other information: ' &
                // this%extra
            lstr = lstr + 20 + len(this%extra)
        end if
        allocate(character(len=lstr) :: msg_multiline)
        msg_multiline = line(:lstr)
    end if

end subroutine handler_get_error

! ======================================================================

subroutine handler_print_error(this)
    !! Print error message.
    !!
    !! Print an error message based on the stored parameters
    !!
    !! @note
    !! This function bypasses the test on the minimum level for printing
    !! and will always print.
    !! @endnote
    class(ErrorHandler), intent(in) :: this
    !! Instance of ErrorHandler.

    character(len=3) :: label
    logical :: has_details, has_extra, has_source

    select case(this%label)
    case (0)
        label = 'gen'
    case (1)
        label = 'std'
    case (2)
        label = 'dev'
    case default
        label = 'n/a'
    end select

    has_details = allocated(this%details)
    has_extra = allocated(this%extra)
    has_source = allocated(this%source)
    if (has_details .and. has_extra .and. has_source) then
        call write_err(label, this%cause, this%details, this%extra, &
                       this%source)
    else if (has_details .and. has_extra) then
        call write_err(label, this%cause, details=this%details, &
                       extra=this%extra)
    else if (has_details .and. has_source) then
        call write_err(label, this%cause, details=this%details, &
                       source=this%source)
    else if (has_extra .and. has_source) then
        call write_err(label, this%cause, extra=this%extra, source=this%source)
    else if (has_details) then
        call write_err(label, this%cause, details=this%details)
    else if (has_extra) then
        call write_err(label, this%cause, extra=this%extra)
    else if (has_source) then
        call write_err(label, this%cause, source=this%source)
    else
        call write_err(label, this%cause)
    end if

end subroutine handler_print_error

! ======================================================================

end module exception
