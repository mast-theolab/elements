submodule (parse_cmdline) parse_cmdline_setval

    implicit none

contains

! ======================================================================

function set_real_bound_msg(val_min, val_max) result(msg)
    !! Build a suitable message to state the accepted values.
    !!
    !! Builds and returns a string containing a message on accepted
    !!   values.
    !! The message is adapted based on provided bounds.

    real(real64), intent(in), optional :: val_min
    !! Minimum accepted value.
    real(real64), intent(in), optional :: val_max
    !! Maximum accepted value.
    character(len=:), allocatable :: msg
    !! Statement.

    integer :: N
    character(len=:), allocatable :: txt, txt2

    if (present(val_min) .and. present(val_max)) then
        txt = val_txt(val_min)
        txt2 = val_txt(val_max)
        N = 28 + len_trim(txt) + len_trim(txt2)
        allocate(character(len=N) :: msg)
        write(msg, '("value must be between ",a," and ",a,".")') trim(txt), &
            trim(txt2)
    else if (present(val_min)) then
        txt = val_txt(val_min)
        N = 27 + len_trim(txt)
        allocate(character(len=N) :: msg)
        write(msg, '("minimum accepted value is ",a,".")') trim(txt)
    else if (present(val_max)) then
        txt = val_txt(val_max)
        N = 27 + len_trim(txt)
        allocate(character(len=N) :: msg)
        write(msg, '("maximum accepted valuminimum accepted value is e is ",a,".")') trim(txt)
    else
        error stop "val_min and val_max are both missing"
    end if

contains
    function val_txt(value) result(text)
        !! Write a text suitable for the format.

        real(real64), intent(in) :: value
        character(len=:), allocatable :: text

        integer :: i
        character(len=64) :: fmt, test

        if (abs(value) < tiny(value)) then
            text = '0.0'
        else if (abs(value) >= 0.01_real64 .and. abs(value) <= 5000.0_real64) &
            then
            write(test, '(f0.16)') value
            i = verify(test, ' 0', back=.True.)
            text = test(:i)
        else
            if (value < 0.0_real64) then
                i = 11
            else
                i = 10
            end if
            allocate(character(len=i) :: text)
            write(fmt, '("(es",i0,".4)")') i
            write(test, fmt) value
            text = trim(test)
        end if
    end function val_txt

end function set_real_bound_msg

! ======================================================================

module procedure upd_argint_value
    !! Update the value of an integer-type argument.
    !!
    !! Updates the value of argument `this` based, either based on
    !!   internal constants or the content of string.
    !! The update can be an increment or an assignment based on
    !!   the internal logic of the argument.

    integer :: ios
    integer(int64) :: value

    character(len=256) :: msg

    err = InitError()

    if (this%max_num == 0) then
        if (this%add_value .and. this%is_set.eq.2) then
            this%value = this%value + this%const
        else
            this%value = this%const
            this%is_set = 2
        end if
    else
        if (.not.present(string)) then
            call RaiseError(err, 'Missing value to update ' // this%label)
            return
        end if
        if (len_trim(string) > 0) then
            read(string, *, iostat=ios) value
            if (ios /= 0) then
                call RaiseValueError(err, &
                    'Expected integer value for ' // this%label)
                return
            end if
            if (value < this%min_ok .or. value > this%max_ok) then
                if (value < this%min_ok) then
                    write(msg, '("minimum accepted value is, ",i0)') &
                        this%min_ok
                else
                    write(msg, '("maximum accepted value is, ",i0)') &
                        this%max_ok
                end if
                call RaiseValueError(err, &
                    'Value out of range for ' // this%label, &
                    msg)
                return
            end if
            if (this%add_value .and. this%is_set.eq.2) then
                this%value = this%value + value
            else
                this%value = value
                this%is_set = 2
            end if
        end if
    end if

end procedure upd_argint_value

! ======================================================================

module procedure set_argint_list_str
    !! Build a list of integers from the content of `string`.
    !!
    !! Builds a list of integers contained in `string` and separated
    !!   by `sep`.
    !! The list needs to be built in one step to allow the proper
    !!   allocation of the space.

    integer :: i0, i1, ios, N
    integer(int64) :: value
    character(len=256) :: msg, msg2
    character(len=:), allocatable :: seps

    err = InitError()

    if (.not.present(string)) then
        call RaiseArgError(err, 'string', 'Missing string containing values.')
        return
    end if

    if (allocated(this%values)) deallocate(this%values)

    if (present(sep)) then
        ! allocate(character(len=len(sep)) :: seps)
        seps = sep
    else
        seps = ','
    end if

    ! First, let us find how many blocks are present
    N = 1
    do i0 = 1, len(string)
        if (index(seps, string(i0:i0)) > 0) N = N + 1
    end do
    if (N < this%min_num .or. N > this%max_num) then
        write(msg, '("Incorrect number of arguments for ",a)') this%label
        call RaiseError(err, trim(msg))
        return
    end if

    allocate(this%values(N))

    ! Now let us parse the content of string.
    N = 1
    i0 = 1
    do
        i1 = i0 + 1
        do while (index(seps, string(i1:i1)) == 0)
            i1 = i1 + 1
            if (i1 == len(string)) exit
        end do
        read(string(i0:i1-1), *, iostat=ios) value
        if (ios /= 0) then
            write(msg, '("Incorrect value in position ",i0," for ",a)') N, &
                this%label
            call RaiseValueError(err, trim(msg))
            return
        end if
        if (value < this%min_ok .or. value > this%max_ok) then
            write(msg, '("Value in position ",i0," out of range for ",a)') N, &
                this%label
            if (value < this%min_ok) then
                write(msg, '("minimum accepted value is, ",i0)') this%min_ok
            else
                write(msg, '("maximum accepted value is, ",i0)') this%max_ok
            end if
            call RaiseValueError(err, msg, msg2)
            return
        end if
        this%values(N) = value
        if (i1 < len(string)) then
            i0 = i1 + 1
            N = N + 1
        else
            exit
        end if
    end do

    ! Sanitary check, this part should not happen
    if (N < size(this%values)) then
        write(msg, '("Parser has failed to identify all components in ",a)') &
            this%label
        call RaiseError(err, trim(msg))
        return
    end if

    this%is_set = 2

end procedure set_argint_list_str

! ======================================================================

module procedure set_argint_list_arr
    !! Build a list of integers from a list of `strings`.
    !!
    !! Builds a list of integers stored as character strings `strings`.
    !! The list needs to be built in one step to allow the proper
    !!   allocation of the space.
    !!
    !! @note
    !!   For simplicity, the function does not support list of strings
    !!   that would be themselves character-separated lists.
    !!   Each string must be a separate number.
    !! @endnote

    integer :: i, ios, N
    integer(int64) :: value
    character(len=256) :: msg, msg2

    err = InitError()

    if (allocated(this%values)) deallocate(this%values)

    N = size(strings)
    if (N < this%min_num .or. N > this%max_num) then
        write(msg, '("Incorrect number of arguments for ",a)') this%label
        call RaiseError(err, msg)
        return
    end if

    allocate(this%values(N))

    ! Now let us parse the strings.
    do i = 1, N
        read(strings(i), *, iostat=ios) value
        if (ios /= 0) then
            write(msg, '("Incorrect value in position ",i0," for ",a)') i, &
                this%label
            call RaiseValueError(err, trim(msg))
            return
        end if
        if (value < this%min_ok .or. value > this%max_ok) then
            write(msg, '("Value in position ",i0," out of range for ",a)') i, &
                this%label
            if (value < this%min_ok) then
                write(msg, '("minimum accepted value is, ",i0)') this%min_ok
            else
                write(msg, '("maximum accepted value is, ",i0)') this%max_ok
            end if
            call RaiseValueError(err, msg, msg2)
            return
        end if
        this%values(i) = value
    end do

    this%is_set = 2

end procedure set_argint_list_arr

! ======================================================================

module procedure upd_argreal_value
    !! Update the value of a real-type argument.
    !!
    !! Updates the value of argument `this` based, either based on
    !!   internal constants or the content of string.
    !! The update can be an increment or an assignment based on
    !!   the internal logic of the argument.

    integer :: ios
    real(real64) :: value

    character(len=:), allocatable :: msg

    err = InitError()

    if (this%max_num == 0) then
        if (this%add_value .and. this%is_set.eq.2) then
            this%value = this%value + this%const
        else
            this%value = this%const
            this%is_set = 2
        end if
    else
        if (.not.present(string)) then
            call RaiseError(err, 'Missing value to update ' // this%label)
            return
        end if
        if (len_trim(string) > 0) then
            read(string, *, iostat=ios) value
            if (ios /= 0) then
                call RaiseValueError(err, &
                    'Expected real value for ' // this%label)
                return
            end if
            if (value < this%min_ok .or. value > this%max_ok) then
                if (value < this%min_ok) then
                    msg = set_real_bound_msg(val_min=this%min_ok)
                else
                    msg = set_real_bound_msg(val_max=this%max_ok)
                end if
                call RaiseValueError(err, &
                    'Value out of range for ' // this%label, &
                    msg)
                return
            end if
            if (this%add_value .and. this%is_set.eq.2) then
                this%value = this%value + value
            else
                this%value = value
                this%is_set = 2
            end if
        end if
    end if

end procedure upd_argreal_value

! ======================================================================

module procedure set_argreal_list_str
    !! Build a list of reals from the content of `string`.
    !!
    !! Builds a list of reals contained in `string` and separated
    !!   by `sep`.
    !! The list needs to be built in one step to allow the proper
    !!   allocation of the space.

    integer :: i0, i1, ios, N
    real(real64) :: value
    character(len=256) :: msg
    character(len=:), allocatable :: msg2
    character(len=:), allocatable :: seps

    err = InitError()

    if (.not.present(string)) then
        call RaiseArgError(err, 'string', 'Missing string containing values.')
        return
    end if

    if (allocated(this%values)) deallocate(this%values)

    if (present(sep)) then
        ! allocate(character(len=len(sep)) :: seps)
        seps = sep
    else
        seps = ','
    end if

    ! First, let us find how many blocks are present
    N = 1
    do i0 = 1, len(string)
        if (index(seps, string(i0:i0)) > 0) N = N + 1
    end do
    if (N < this%min_num .or. N > this%max_num) then
        write(msg, '("Incorrect number of arguments for ",a)') this%label
        call RaiseError(err, trim(msg))
        return
    end if

    allocate(this%values(N))

    ! Now let us parse the content of string.
    N = 1
    i0 = 1
    do
        i1 = i0 + 1
        do while (index(seps, string(i1:i1)) == 0)
            i1 = i1 + 1
            if (i1 == len(string)) exit
        end do
        read(string(i0:i1-1), *, iostat=ios) value
        if (ios /= 0) then
            write(msg, '("Incorrect value in position ",i0," for ",a)') N, &
                this%label
            call RaiseValueError(err, trim(msg))
            return
        end if
        if (value < this%min_ok .or. value > this%max_ok) then
            write(msg, '("Value in position ",i0," out of range for ",a)') N, &
                this%label
            if (value < this%min_ok) then
                msg2 = set_real_bound_msg(val_min=this%min_ok)
            else
                msg2 = set_real_bound_msg(val_max=this%max_ok)
            end if
            call RaiseValueError(err, msg, msg2)
            return
        end if
        this%values(N) = value
        if (i1 < len(string)) then
            i0 = i1 + 1
            N = N + 1
        else
            exit
        end if
    end do

    ! Sanitary check, this part should not happen
    if (N < size(this%values)) then
        write(msg, '("Parser has failed to identify all components in ",a)') &
            this%label
        call RaiseError(err, trim(msg))
        return
    end if

    this%is_set = 2

end procedure set_argreal_list_str

! ======================================================================

module procedure set_argreal_list_arr
    !! Build a list of reals from a list of `strings`.
    !!
    !! Builds a list of reals stored as character strings `strings`.
    !! The list needs to be built in one step to allow the proper
    !!   allocation of the space.
    !!
    !! @note
    !!   For simplicity, the function does not support list of strings
    !!   that would be themselves character-separated lists.
    !!   Each string must be a separate number.
    !! @endnote

    integer :: i, ios, N
    real(real64) :: value
    character(len=256) :: msg
    character(len=:), allocatable :: msg2

    err = InitError()

    if (allocated(this%values)) deallocate(this%values)

    N = size(strings)
    if (N < this%min_num .or. N > this%max_num) then
        write(msg, '("Incorrect number of arguments for ",a)') this%label
        call RaiseError(err, msg)
        return
    end if

    allocate(this%values(N))

    ! Now let us parse the strings.
    do i = 1, N
        read(strings(i), *, iostat=ios) value
        if (ios /= 0) then
            write(msg, '("Incorrect value in position ",i0," for ",a)') i, &
                this%label
            call RaiseValueError(err, msg)
            return
        end if
        if (value < this%min_ok .or. value > this%max_ok) then
            write(msg, '("Value in position ",i0," out of range for ",a)') i, &
                this%label
            if (value < this%min_ok) then
                msg2 = set_real_bound_msg(val_min=this%min_ok)
            else
                msg2 = set_real_bound_msg(val_max=this%max_ok)
            end if
            call RaiseValueError(err, msg, msg2)
            return
        end if
        this%values(i) = value
    end do

    this%is_set = 2

end procedure set_argreal_list_arr

! ======================================================================

module procedure upd_argbool_value
    !! Update the value of an logical-type argument.
    !!
    !! Updates the value of argument `this` based, either based on
    !!   internal constants or the content of string.
    !! The update can be an increment or an assignment based on
    !!   the internal logic of the argument.

    err = InitError()

    this%value = this%const
    this%is_set = 2

end procedure upd_argbool_value

! ======================================================================

module procedure upd_argchar_value
    !! Update the value of a character-type argument.
    !!
    !! Updates the value of argument `this` based, either based on
    !!   internal constants or the content of string.
    !! The update can be an increment or an assignment based on
    !!   the internal logic of the argument.

    err = InitError()

    if (.not.present(string)) then
        call RaiseError(err, 'Missing value to update ' // this%label)
        return
    end if
    if (len_trim(string) > 0) then
        if (allocated(this%value)) deallocate(this%value)
        allocate(character(len=len_trim(string)) :: this%value)
        this%value = trim(string)
    end if

    this%is_set = 2

end procedure upd_argchar_value

! ======================================================================

module procedure set_argchar_list_str
    !! Store `string` in list of character arguments.
    !!
    !! Constructs a 1-item list for character-style arguments from
    !!   `string`.
    !! The list needs to be built in one step to allow the proper
    !!   allocation of the space.

    integer :: N
    character(len=256) :: msg

    err = InitError()

    if (.not.present(string)) then
        call RaiseArgError(err, 'string', 'Missing string containing values.')
        return
    end if

    if (allocated(this%values)) deallocate(this%values)

    ! First, let us find how many blocks are present
    N = 1
    if (N < this%min_num .or. N > this%max_num) then
        write(msg, '("Incorrect number of arguments for ",a)') this%label
        call RaiseError(err, trim(msg))
        return
    end if

    allocate(character(len=len_trim(string)) :: this%values(N))

    this%values(1) = trim(string)
    this%is_set = 2

end procedure set_argchar_list_str

! ======================================================================

module procedure set_argchar_list_arr
    !! Build a list of strings from a list of `strings`.
    !!
    !! Builds a list of strings stored as character strings `strings`.
    !! The list needs to be built in one step to allow the proper
    !!   allocation of the space.

    integer :: i, lmax, lstr, N
    character(len=256) :: msg

    err = InitError()

    if (allocated(this%values)) deallocate(this%values)

    N = size(strings)
    if (N < this%min_num .or. N > this%max_num) then
        write(msg, '("Incorrect number of arguments for ",a)') this%label
        call RaiseError(err, msg)
        return
    end if

    lmax = 0
    do i = 1, N
        lstr = len_trim(strings(i))
        if (lstr > lmax) lmax = lstr
    end do
    allocate(character(len=lmax) :: this%values(N))

    ! Now let us parse the strings.
    do i = 1, N
        this%values(i) = trim(strings(i))
    end do

    this%is_set = 2

end procedure set_argchar_list_arr

! ======================================================================

module procedure upd_noarrayI
    !! Update value from array - dummy version
    !!
    !! Dummy version of a procedure to build value from array of strings.

    err = InitError()

    call RaiseError(err, 'Values cannot be set from array of strings.')

end procedure upd_noarrayI

! ======================================================================

module procedure upd_noarrayR
    !! Update value from array - dummy version
    !!
    !! Dummy version of a procedure to build value from array of strings.

    err = InitError()

    call RaiseError(err, 'Values cannot be set from array of strings.')

end procedure upd_noarrayR

! ======================================================================

module procedure upd_noarrayB
    !! Update value from array - dummy version
    !!
    !! Dummy version of a procedure to build value from array of strings.

    err = InitError()

    call RaiseError(err, 'Values cannot be set from array of strings.')

end procedure upd_noarrayB

! ======================================================================

module procedure upd_noarrayC
    !! Update value from array - dummy version
    !!
    !! Dummy version of a procedure to build value from array of strings.

    err = InitError()

    call RaiseError(err, 'Values cannot be set from array of strings.')

end procedure upd_noarrayC

! ======================================================================

module procedure upd_noscalar
    !! Update value from scalar - dummy version
    !!
    !! Dummy version of a procedure to build value from string.

    err = InitError()

    call RaiseError(err, 'Values cannot be set from a string.')

end procedure upd_noscalar

! ======================================================================

module procedure upd_noarray
    !! Update value from array - dummy version
    !!
    !! Dummy version of a procedure to build value from string.

    err = InitError()

    call RaiseError(err, 'Values cannot be set from a string.')

end procedure upd_noarray

! ======================================================================


end submodule parse_cmdline_setval
