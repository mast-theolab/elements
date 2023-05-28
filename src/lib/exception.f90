module exception
    !! Module for exception handling
    !!
    !! Module providing types for exception handling in other modules.
    private
    public :: InitError, RaiseError, RaiseArgError, RaiseTermination, &
        RaiseValueError
    ! public :: InitError!, RaiseArgError, RaiseError

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

    type, public, extends(BaseException) :: ArgumentError
        character(len=256), public :: arg = ' '
        character(len=512), public :: reason = ' '
    contains
        procedure, private :: set_msg => set_argerr_message
    end type ArgumentError

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

subroutine RaiseTermination(err, msg_)
    !! Set a termination signal.
    !!
    !! Takes an instance of BaseException and sets it to Termination
    !!   with an optional message set in input.
    class(BaseException), allocatable, intent(inout) :: err
    !! Original error, updated on return
    character(len=*), intent(in), optional :: msg_
    !! Message to include in the error

    type(Termination), allocatable :: newerr
    
    allocate(newerr)
    if (present(msg_)) then
        call newerr%raise(msg_)
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
    !! Reason why value is invalid.
    character(len=*), intent(in), optional :: motive
    !! Reason why value is invalid.

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
    !! Set message of a standard error
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

subroutine set_value_message(this, msg)
    !! Set message of a standard error
    !!
    !! Sets a message and changes the status of a basic error
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

subroutine set_argerr_message(this, msg)
    !! Set error message for an ArgumentError exception
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

end module exception
