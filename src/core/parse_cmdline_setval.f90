submodule (parse_cmdline) parse_cmdline_setval

    implicit none

contains

! ======================================================================
! TYPE-BOUND PROCEDURES (DEFINITIONS)
! ======================================================================

module procedure argbool_from_list_no
    !! Update value from array - dummy version
    !!
    !! Dummy version of a procedure to build value from array of strings.

    done = .false.
    call arg%error%raise_argerror('type', &
        'scalar logical argument cannot be set from array of strings')

end procedure argbool_from_list_no

! ======================================================================

module procedure argbool_from_str
    !! Update the value of an logical-type argument.
    !!
    !! Updates the value of argument `this` based, either based on
    !!   internal constants or the content of string.
    !! The update can be an increment or an assignment based on
    !!   the internal logic of the argument.

    arg%value = arg%const
    arg%is_set = 2
    done = .true.

end procedure argbool_from_str

! ======================================================================

module procedure argchar_from_list_no
    !! Update value from array - dummy version
    !!
    !! Dummy version of a procedure to build value from array of strings.

    done = .false.
    call arg%error%raise_argerror('type', &
        'scalar character argument cannot be set from array of strings')

end procedure argchar_from_list_no

! ======================================================================

module procedure argchar_from_str
    !! Update the value of a character-type argument.
    !!
    !! Updates the value of argument `this` based, either based on
    !!   internal constants or the content of string.
    !! The update can be an increment or an assignment based on
    !!   the internal logic of the argument.

    done = .false.
    if (.not.present(argval)) then
        call arg%error%raise_argerror('missing', &
            'missing value to update ' // arg%label)
        return
    end if
    if (len_trim(argval) > 0) then
        if (allocated(arg%value)) deallocate(arg%value)
        allocate(character(len=len_trim(argval)) :: arg%value)
        arg%value = trim(argval)
    end if

    arg%is_set = 2
    done = .true.

end procedure argchar_from_str

! ======================================================================

module procedure argchars_from_list
    !! Build a list of strings from a list of `strings`.
    !!
    !! Builds a list of strings stored as character strings `strings`.
    !! The list needs to be built in one step to allow the proper
    !!   allocation of the space.

    integer :: i, lmax, lstr, N
    character(len=256) :: msg

    done = .false.

    if (allocated(arg%values)) deallocate(arg%values)

    N = size(argvals)
    if (N < arg%min_num .or. N > arg%max_num) then
        write(msg, '("incorrect number of arguments for ",a)') arg%label
        call arg%error%raise_argerror('number', msg)
        return
    end if

    lmax = 0
    do i = 1, N
        lstr = len_trim(argvals(i))
        if (lstr > lmax) lmax = lstr
    end do
    allocate(character(len=lmax) :: arg%values(N))

    ! Now let us parse the strings.
    do i = 1, N
        arg%values(i) = trim(argvals(i))
    end do

    arg%is_set = 2
    done = .true.

end procedure argchars_from_list

! ======================================================================

module procedure argchars_from_str
    !! Store `string` in list of character arguments.
    !!
    !! Constructs a 1-item list for character-style arguments from
    !!   `string`.
    !! The list needs to be built in one step to allow the proper
    !!   allocation of the space.

    integer :: N
    character(len=256) :: msg

    done = .false.

    if (.not.present(argval)) then
        call arg%error%raise_argerror('missing', &
            'missing string containing values')
        return
    end if

    if (allocated(arg%values)) deallocate(arg%values)

    ! First, let us find how many blocks are present
    N = 1
    if (N < arg%min_num .or. N > arg%max_num) then
        write(msg, '("incorrect number of arguments for ",a)') arg%label
        call arg%error%raise_argerror('number', msg)
        return
    end if

    allocate(character(len=len_trim(argval)) :: arg%values(N))

    arg%values(1) = trim(argval)
    arg%is_set = 2
    done = .true.

end procedure argchars_from_str

! ======================================================================

module procedure arggen_from_list_no
    !! Update value from array - dummy version
    !!
    !! Dummy version of a procedure to build value from string.

    done = .false.
    call arg%error%raise_argerror('type', &
        'scalar generic argument cannot be set from array of strings')

end procedure arggen_from_list_no

! ======================================================================

module procedure arggen_from_str_no
    !! Update value from scalar - dummy version
    !!
    !! Dummy version of a procedure to build value from string.

    done = .false.
    call arg%error%raise_argerror('type', &
        'scalar generic argument cannot be set from string')

end procedure arggen_from_str_no

! ======================================================================

module procedure argint_from_list_no
    !! Update value from array - dummy version
    !!
    !! Dummy version of a procedure to build value from array of strings.

    done = .false.
    call arg%error%raise_argerror('type', &
        'scalar integer argument cannot be set from array of string')

end procedure argint_from_list_no

! ======================================================================

module procedure argint_from_str
    !! Update the value of an integer-type argument.
    !!
    !! Updates the value of argument `this` based, either based on
    !!   internal constants or the content of string.
    !! The update can be an increment or an assignment based on
    !!   the internal logic of the argument.

    integer :: ios
    integer(int64) :: ival
    character(len=256) :: msg

    done = .false.

    if (arg%max_num == 0) then
        if (arg%add_value .and. arg%is_set.eq.2) then
            arg%value = arg%value + arg%const
        else
            arg%value = arg%const
            arg%is_set = 2
        end if
    else
        if (.not.present(argval)) then
            call arg%error%raise_argerror('missing', &
                'missing value to update ' // arg%label)
            return
        end if
        if (len_trim(argval) > 0) then
            read(argval, *, iostat=ios) ival
            if (ios /= 0) then
                call arg%error%raise_argerror('type', &
                    'expected integer value for ' // arg%label)
                return
            end if
            if (ival < arg%min_ok .or. ival > arg%max_ok) then
                if (ival < arg%min_ok) then
                    write(msg, '("minimum accepted value is, ",i0)') arg%min_ok
                else
                    write(msg, '("maximum accepted value is, ",i0)') arg%max_ok
                end if
                call arg%error%raise_argerror('value', &
                    'value out of range for ' // arg%label, msg)
                return
            end if
            if (arg%add_value .and. arg%is_set.eq.2) then
                arg%value = arg%value + ival
            else
                arg%value = ival
                arg%is_set = 2
            end if
        end if
    end if

    done = .true.

end procedure argint_from_str

! ======================================================================

module procedure argints_from_list
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
    integer(int64) :: ival
    character(len=256) :: msg, msg2

    done = .false.

    if (allocated(arg%values)) deallocate(arg%values)

    N = size(argvals)
    if (N < arg%min_num .or. N > arg%max_num) then
        write(msg, '("incorrect number of arguments for ",a)') arg%label
        call arg%error%raise_argerror('number', msg)
        return
    end if

    allocate(arg%values(N))

    ! Now let us parse the strings.
    do i = 1, N
        read(argvals(i), *, iostat=ios) ival
        if (ios /= 0) then
            write(msg, '("incorrect value in position ",i0," for ",a)') N, &
                arg%label
            call arg%error%raise_argerror('value', msg)
            return
        end if
        if (ival < arg%min_ok .or. ival > arg%max_ok) then
            write(msg, '("value in position ",i0," out of range for ",a)') i, &
                arg%label
            if (ival < arg%min_ok) then
                write(msg, '("minimum accepted value is, ",i0)') arg%min_ok
            else
                write(msg, '("maximum accepted value is, ",i0)') arg%max_ok
            end if
            call arg%error%raise_argerror('value', msg, msg2)
            return
        end if
        arg%values(i) = ival
    end do

    arg%is_set = 2
    done = .true.

end procedure argints_from_list

! ======================================================================

module procedure argints_from_str
    !! Build a list of integers from the content of `string`.
    !!
    !! Builds a list of integers contained in `string` and separated
    !!   by `sep`.
    !! The list needs to be built in one step to allow the proper
    !!   allocation of the space.

    integer :: i0, i1, ios, N
    integer(int64) :: ival
    character(len=256) :: msg, msg2
    character(len=:), allocatable :: seps

    done = .false.

    if (.not.present(argval)) then
        call arg%error%raise_argerror('missing', &
            'missing string containing values')
        return
    end if

    if (allocated(arg%values)) deallocate(arg%values)

    if (present(sep)) then
        ! allocate(character(len=len(sep)) :: seps)
        seps = sep
    else
        seps = ','
    end if

    ! First, let us find how many blocks are present
    N = 1
    do i0 = 1, len(argval)
        if (index(seps, argval(i0:i0)) > 0) N = N + 1
    end do
    if (N < arg%min_num .or. N > arg%max_num) then
        write(msg, '("incorrect number of arguments for ",a)') arg%label
        call arg%error%raise_argerror('number', msg)
        return
    end if

    allocate(arg%values(N))

    ! Now let us parse the content of string.
    N = 1
    i0 = 1
    do
        i1 = i0 + 1
        do while (index(seps, argval(i1:i1)) == 0)
            i1 = i1 + 1
            if (i1 == len(argval)) exit
        end do
        read(argval(i0:i1-1), *, iostat=ios) ival
        if (ios /= 0) then
            write(msg, '("incorrect value in position ",i0," for ",a)') N, &
                arg%label
            call arg%error%raise_argerror('value', msg)
            return
        end if
        if (ival < arg%min_ok .or. ival > arg%max_ok) then
            write(msg, '("value in position ",i0," out of range for ",a)') N, &
                arg%label
            if (ival < arg%min_ok) then
                write(msg, '("minimum accepted value is, ",i0)') arg%min_ok
            else
                write(msg, '("maximum accepted value is, ",i0)') arg%max_ok
            end if
            call arg%error%raise_argerror('value', msg, msg2)
            return
        end if
        arg%values(N) = ival
        if (i1 < len(argval)) then
            i0 = i1 + 1
            N = N + 1
        else
            exit
        end if
    end do

    ! Sanitary check, this part should not happen
    if (N < size(arg%values)) then
        write(msg, '("parser has failed to identify all components in ",a)') &
            arg%label
        call arg%error%raise_deverror('argval', msg)
        return
    end if

    arg%is_set = 2
    done = .true.

end procedure argints_from_str

! ======================================================================

module procedure argreal_from_list_no
    !! Update value from array - dummy version
    !!
    !! Dummy version of a procedure to build value from array of strings.

    done = .false.
    call arg%error%raise_argerror('type', &
        'scalar real argument cannot be set from array of strings.')

end procedure argreal_from_list_no

! ======================================================================

module procedure argreal_from_str
    !! Update the value of a real-type argument.
    !!
    !! Updates the value of argument `this` based, either based on
    !!   internal constants or the content of string.
    !! The update can be an increment or an assignment based on
    !!   the internal logic of the argument.

    integer :: ios
    real(real64) :: rval
    character(len=:), allocatable :: msg

    done = .false.

    if (arg%max_num == 0) then
        if (arg%add_value .and. arg%is_set.eq.2) then
            arg%value = arg%value + arg%const
        else
            arg%value = arg%const
            arg%is_set = 2
        end if
    else
        if (.not.present(argval)) then
            call arg%error%raise_argerror('missing', &
                'missing value to update ' // arg%label)
            return
        end if
        if (len_trim(argval) > 0) then
            read(argval, *, iostat=ios) rval
            if (ios /= 0) then
                call arg%error%raise_argerror('type', &
                    'expected real value for ' // arg%label)
                return
            end if
            if (rval < arg%min_ok .or. rval > arg%max_ok) then
                if (rval < arg%min_ok) then
                    msg = set_real_bound_msg(val_min=arg%min_ok)
                else
                    msg = set_real_bound_msg(val_max=arg%max_ok)
                end if
                call arg%error%raise_argerror('value', &
                    'value out of range for ' // arg%label, msg)
                return
            end if
            if (arg%add_value .and. arg%is_set.eq.2) then
                arg%value = arg%value + rval
            else
                arg%value = rval
                arg%is_set = 2
            end if
        end if
    end if

    done = .true.

end procedure argreal_from_str

! ======================================================================

module procedure argreals_from_list
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
    real(real64) :: rval
    character(len=256) :: msg
    character(len=:), allocatable :: msg2

    done = .false.

    if (allocated(arg%values)) deallocate(arg%values)

    N = size(argvals)
    if (N < arg%min_num .or. N > arg%max_num) then
        write(msg, '("incorrect number of arguments for ",a)') arg%label
        call arg%error%raise_argerror('number', msg)
        return
    end if

    allocate(arg%values(N))

    ! Now let us parse the strings.
    do i = 1, N
        read(argvals(i), *, iostat=ios) rval
        if (ios /= 0) then
            write(msg, '("incorrect value in position ",i0," for ",a)') N, &
                arg%label
            call arg%error%raise_argerror('value', msg)
            return
        end if
        if (rval < arg%min_ok .or. rval > arg%max_ok) then
            write(msg, '("value in position ",i0," out of range for ",a)') i, &
                arg%label
            if (rval < arg%min_ok) then
                msg2 = set_real_bound_msg(val_min=arg%min_ok)
            else
                msg2 = set_real_bound_msg(val_max=arg%max_ok)
            end if
            call arg%error%raise_argerror('value', msg, msg2)
            return
        end if
        arg%values(i) = rval
    end do

    arg%is_set = 2

end procedure argreals_from_list

! ======================================================================

module procedure argreals_from_str
    !! Build a list of reals from the content of `string`.
    !!
    !! Builds a list of reals contained in `string` and separated
    !!   by `sep`.
    !! The list needs to be built in one step to allow the proper
    !!   allocation of the space.

    integer :: i0, i1, ios, N
    real(real64) :: rval
    character(len=256) :: msg
    character(len=:), allocatable :: msg2
    character(len=:), allocatable :: seps

    done = .true.

    if (.not.present(argval)) then
        call arg%error%raise_argerror('missing', &
            'missing string containing values')
        return
    end if

    if (allocated(arg%values)) deallocate(arg%values)

    if (present(sep)) then
        ! allocate(character(len=len(sep)) :: seps)
        seps = sep
    else
        seps = ','
    end if

    ! First, let us find how many blocks are present
    N = 1
    do i0 = 1, len(argval)
        if (index(seps, argval(i0:i0)) > 0) N = N + 1
    end do
    if (N < arg%min_num .or. N > arg%max_num) then
        write(msg, '("incorrect number of arguments for ",a)') arg%label
        call arg%error%raise_argerror('number', msg)
        return
    end if

    allocate(arg%values(N))

    ! Now let us parse the content of string.
    N = 1
    i0 = 1
    do
        i1 = i0 + 1
        do while (index(seps, argval(i1:i1)) == 0)
            i1 = i1 + 1
            if (i1 == len(argval)) exit
        end do
        read(argval(i0:i1-1), *, iostat=ios) rval
        if (ios /= 0) then
            write(msg, '("incorrect value in position ",i0," for ",a)') N, &
                arg%label
            call arg%error%raise_argerror('value', msg)
            return
        end if
        if (rval < arg%min_ok .or. rval > arg%max_ok) then
            write(msg, '("value in position ",i0," out of range for ",a)') N, &
                arg%label
            if (rval < arg%min_ok) then
                msg2 = set_real_bound_msg(val_min=arg%min_ok)
            else
                msg2 = set_real_bound_msg(val_max=arg%max_ok)
            end if
            call arg%error%raise_argerror('value', msg, msg2)
            return
        end if
        arg%values(N) = rval
        if (i1 < len(argval)) then
            i0 = i1 + 1
            N = N + 1
        else
            exit
        end if
    end do

    ! Sanitary check, this part should not happen
    if (N < size(arg%values)) then
        write(msg, '("parser has failed to identify all components in ",a)') &
            arg%label
        call arg%error%raise_deverror('argval', msg)
        return
    end if

    arg%is_set = 2
    done = .true.

end procedure argreals_from_str

! ======================================================================
! INTERNAL PROCEDURES
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
        call run%error%raise_deverror('argval', &
            'val_min and val_max are both missing')
        return
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

end submodule parse_cmdline_setval
