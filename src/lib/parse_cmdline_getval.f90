submodule (parse_cmdline) parse_cmdline_getval

    implicit none

contains

! ======================================================================

module procedure argval_set_by_user
    !! Return True is value set by user.
    !!
    !! Checks if value set by user and returns True in that case.

    integer :: iarg

    err = InitError()

    iarg = this%get_argname_id(argname)
    if (iarg == 0) then
        call RaiseError(err, 'Unknown argument name.')
        return
    end if

    res = this%args(iarg)%arg%is_set == 2

end procedure argval_set_by_user

! ======================================================================

module procedure get_value_int32val
    !! Get value for a scalar integer-type option.
    !!
    !! Gets value related to a scalar integer-type option.

    integer :: iarg

    err = InitError()

    iarg = this%get_argname_id(argname)
    if (iarg == 0) then
        call RaiseError(err, 'Unknown argument name.')
        return
    end if

    select type (opt => this%args(iarg)%arg)
        class is (ArgInt)
            if (opt%is_set == 0) then
                call RaiseError(err, 'Unset argument')
                return
            end if
            if (abs(opt%value) > huge(1_int32)) then
                call RaiseError(err, &
                    'Cannot represent value with current precision')
                return
            else
                result = opt%value
            end if
        class default
            call RaiseError(err, 'Argument is not a scalar integer.')
            return
    end select

end procedure get_value_int32val

! ======================================================================

module procedure get_value_int64val
    !! Get value for a scalar integer-type option.
    !!
    !! Gets value related to a scalar integer-type option.

    integer :: iarg

    err = InitError()

    iarg = this%get_argname_id(argname)
    if (iarg == 0) then
        call RaiseError(err, 'Unknown argument name.')
        return
    end if

    select type (opt => this%args(iarg)%arg)
        class is (ArgInt)
            if (opt%is_set == 0) then
                call RaiseError(err, 'Unset argument')
                return
            end if
            result = opt%value
        class default
            call RaiseError(err, 'Argument is not a scalar integer.')
            return
    end select

end procedure get_value_int64val

! ======================================================================

module procedure get_value_int32arr
    !! Get list of integer-type values from argname.
    !!
    !! Gets values related to a scalar integer-type option.

    integer :: i, iarg

    err = InitError()

    iarg = this%get_argname_id(argname)
    if (iarg == 0) then
        call RaiseError(err, 'Unknown argument name.')
        return
    end if

    select type (opt => this%args(iarg)%arg)
        class is (ArgIntList)
            if (opt%is_set == 0) then
                call RaiseError(err, 'Unset argument')
                return
            end if
            do i = 1, size(opt%values)
                if (abs(opt%values(i)) > huge(1_int32)) then
                    call RaiseError(err, &
                        'Cannot represent value with current precision')
                    return
                end if
            end do
            result = opt%values
        class default
            call RaiseError(err, 'Argument is not a list of integers.')
            return
    end select

end procedure get_value_int32arr

! ======================================================================

module procedure get_value_int64arr
    !! Get list of integer-type values from argname.
    !!
    !! Gets values related to a scalar integer-type option.

    integer :: iarg

    err = InitError()

    iarg = this%get_argname_id(argname)
    if (iarg == 0) then
        call RaiseError(err, 'Unknown argument name.')
        return
    end if

    select type (opt => this%args(iarg)%arg)
        class is (ArgIntList)
            if (opt%is_set == 0) then
                call RaiseError(err, 'Unset argument')
                return
            end if
            result = opt%values
        class default
            call RaiseError(err, 'Argument is not a list of integers.')
            return
    end select

end procedure get_value_int64arr

! ======================================================================

module procedure get_value_real32val
    !! Get value for a scalar real-type option.
    !!
    !! Gets value related to a scalar real-type option.

    integer :: iarg

    err = InitError()

    iarg = this%get_argname_id(argname)
    if (iarg == 0) then
        call RaiseError(err, 'Unknown argument name.')
        return
    end if

    select type (opt => this%args(iarg)%arg)
        class is (ArgReal)
            if (opt%is_set == 0) then
                call RaiseError(err, 'Unset argument')
                return
            end if
            if (abs(opt%value) > huge(1_real32)) then
                call RaiseError(err, &
                    'Cannot represent value with current precision')
                return
            else
                result = opt%value
            end if
        class default
            call RaiseError(err, 'Argument is not a scalar integer.')
            return
    end select

end procedure get_value_real32val

! ======================================================================

module procedure get_value_real64val
    !! Get value for a scalar real-type option.
    !!
    !! Gets value related to a scalar real-type option.

    integer :: iarg

    err = InitError()

    iarg = this%get_argname_id(argname)
    if (iarg == 0) then
        call RaiseError(err, 'Unknown argument name.')
        return
    end if

    select type (opt => this%args(iarg)%arg)
        class is (ArgReal)
            if (opt%is_set == 0) then
                call RaiseError(err, 'Unset argument')
                return
            end if
            result = opt%value
        class default
            call RaiseError(err, 'Argument is not a scalar real.')
            return
    end select

end procedure get_value_real64val

! ======================================================================

module procedure get_value_real32arr
    !! Get list of real-type values from argname.
    !!
    !! Gets values related to a scalar real-type option.

    integer :: i, iarg

    err = InitError()

    iarg = this%get_argname_id(argname)
    if (iarg == 0) then
        call RaiseError(err, 'Unknown argument name.')
        return
    end if

    select type (opt => this%args(iarg)%arg)
        class is (ArgRealList)
            if (opt%is_set == 0) then
                call RaiseError(err, 'Unset argument')
                return
            end if
            do i = 1, size(opt%values)
                if (abs(opt%values(i)) > huge(1_real32)) then
                    call RaiseError(err, &
                        'Cannot represent value with current precision')
                    return
                end if
            end do
            result = opt%values
        class default
            call RaiseError(err, 'Argument is not a list of integers.')
            return
    end select

end procedure get_value_real32arr

! ======================================================================

module procedure get_value_real64arr
    !! Get list of real-type values from argname.
    !!
    !! Gets values related to a scalar real-type option.

    integer :: iarg

    err = InitError()

    iarg = this%get_argname_id(argname)
    if (iarg == 0) then
        call RaiseError(err, 'Unknown argument name.')
        return
    end if

    select type (opt => this%args(iarg)%arg)
        class is (ArgRealList)
            if (opt%is_set == 0) then
                call RaiseError(err, 'Unset argument')
                return
            end if
            result = opt%values
        class default
            call RaiseError(err, 'Argument is not a list of integers.')
            return
    end select

end procedure get_value_real64arr

! ======================================================================

module procedure get_value_boolval
    !! Get value for a scalar logical-type option.
    !!
    !! Gets value related to a scalar logical-type option.

    integer :: iarg

    err = InitError()

    iarg = this%get_argname_id(argname)
    if (iarg == 0) then
        call RaiseError(err, 'Unknown argument name.')
        return
    end if

    select type (opt => this%args(iarg)%arg)
        class is (ArgBool)
            if (opt%is_set == 0) then
                call RaiseError(err, 'Unset argument')
                return
            end if
            result = opt%value
        class default
            call RaiseError(err, 'Argument is not a scalar logical.')
            return
    end select

end procedure get_value_boolval

! ======================================================================

module procedure get_value_charval
    !! Get value for a scalar character-type option.
    !!
    !! Gets value related to a scalar character-type option.

    integer :: iarg

    err = InitError()

    iarg = this%get_argname_id(argname)
    if (iarg == 0) then
        call RaiseError(err, 'Unknown argument name.')
        return
    end if

    select type (opt => this%args(iarg)%arg)
        class is (ArgChar)
            if (opt%is_set == 0) then
                call RaiseError(err, 'Unset argument')
                return
            end if
            result = opt%value
        class default
            call RaiseError(err, 'Argument is not a scalar character.')
            return
    end select

end procedure get_value_charval

! ======================================================================

module procedure get_value_chararr
    !! Get list of character-type values from argname.
    !!
    !! Gets values related to a list character-type option.

    integer :: iarg

    err = InitError()

    iarg = this%get_argname_id(argname)
    if (iarg == 0) then
        call RaiseError(err, 'Unknown argument name.')
        return
    end if

    select type (opt => this%args(iarg)%arg)
        class is (ArgCharList)
            if (opt%is_set == 0) then
                call RaiseError(err, 'Unset argument')
                return
            end if
            result = opt%values
        class default
            call RaiseError(err, 'Argument is not a list of strings.')
            return
    end select

end procedure get_value_chararr

! ======================================================================

end submodule parse_cmdline_getval
