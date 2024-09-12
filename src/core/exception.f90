module exception
    !! Module for exception handling
    !!
    !! Module providing types for exception handling in other modules.

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

end module exception
