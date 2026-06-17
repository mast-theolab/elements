submodule (parse_cmdline) parse_cmdline_getval

    implicit none

contains

! ======================================================================

module procedure argsdb_getval_bool
    !! Get value for a scalar logical-type option.
    !!
    !! Gets value related to a scalar logical-type option.

    integer :: iarg

    iarg = argsDB%get_argname_id(argname)
    if (iarg == 0) then
        call argsDB%error%raise_deverror('argval', 'unknown argument name')
        return
    end if

    select type (opt => argsDB%args(iarg)%arg)
        class is (arg_bool)
            if (opt%is_set == 0) then
                call argsDB%error%raise_error('value', 'unset', &
                    'value of argument is not set')
                return
            end if
            result = opt%value
        class default
            call argsDB%error%raise_error('value', 'type', &
                'argument is not a logical scalar')
            return
    end select

end procedure argsdb_getval_bool

! ======================================================================

module procedure argsdb_getval_char
    !! Get value for a scalar character-type option.
    !!
    !! Gets value related to a scalar character-type option.

    integer :: iarg

    iarg = argsDB%get_argname_id(argname)
    if (iarg == 0) then
        call argsDB%error%raise_deverror('argval', 'unknown argument name')
        return
    end if

    select type (opt => argsDB%args(iarg)%arg)
        class is (arg_char)
            if (opt%is_set == 0) then
                call argsDB%error%raise_error('value', 'unset', &
                    'value of argument is not set')
                return
            end if
            result = opt%value
        class default
            call argsDB%error%raise_error('value', 'type', &
                'argument is not a character scalar')
            return
    end select

end procedure argsdb_getval_char

! ======================================================================

module procedure argsdb_getval_int32
    !! Get value for a scalar integer-type option.
    !!
    !! Gets value related to a scalar integer-type option.

    integer :: iarg

    iarg = argsDB%get_argname_id(argname)
    if (iarg == 0) then
        call argsDB%error%raise_deverror('argval', 'unknown argument name')
        return
    end if

    select type (opt => argsDB%args(iarg)%arg)
        class is (arg_int)
            if (opt%is_set == 0) then
                call argsDB%error%raise_error('value', 'unset', &
                    'value of argument is not set')
                return
            end if
            if (abs(opt%value) > huge(1_int32)) then
                call argsDB%error%raise_error('value', 'convert', &
                    'range insufficient to represent value of argument')
                return
            else
                result = int(opt%value, int32)
            end if
        class default
            call argsDB%error%raise_error('value', 'type', &
                'argument is not an integer scalar')
            return
    end select

end procedure argsdb_getval_int32

! ======================================================================

module procedure argsdb_getval_int64
    !! Get value for a scalar integer-type option.
    !!
    !! Gets value related to a scalar integer-type option.

    integer :: iarg

    iarg = argsDB%get_argname_id(argname)
    if (iarg == 0) then
        call argsDB%error%raise_deverror('argval', 'unknown argument name')
        return
    end if

    select type (opt => argsDB%args(iarg)%arg)
        class is (arg_int)
            if (opt%is_set == 0) then
                call argsDB%error%raise_error('value', 'unset', &
                    'value of argument is not set')
                return
            end if
            result = opt%value
        class default
            call argsDB%error%raise_error('value', 'type', &
                'argument is not an integer scalar')
            return
    end select

end procedure argsdb_getval_int64

! ======================================================================

module procedure argsdb_getval_real32
    !! Get value for a scalar real-type option.
    !!
    !! Gets value related to a scalar real-type option.

    integer :: iarg

    iarg = argsDB%get_argname_id(argname)
    if (iarg == 0) then
        call argsDB%error%raise_deverror('argval', 'unknown argument name')
        return
    end if

    select type (opt => argsDB%args(iarg)%arg)
        class is (arg_real)
            if (opt%is_set == 0) then
                call argsDB%error%raise_error('value', 'unset', &
                    'value of argument is not set')
                return
            end if
            if (abs(opt%value) > huge(1_real32)) then
                call argsDB%error%raise_error('value', 'convert', &
                    'precision insufficient to represent value of argument')
                return
            else
                result = real(opt%value, real32)
            end if
        class default
            call argsDB%error%raise_error('value', 'type', &
                'argument is not a real scalar')
            return
    end select

end procedure argsdb_getval_real32

! ======================================================================

module procedure argsdb_getval_real64
    !! Get value for a scalar real-type option.
    !!
    !! Gets value related to a scalar real-type option.

    integer :: iarg

    iarg = argsDB%get_argname_id(argname)
    if (iarg == 0) then
        call argsDB%error%raise_deverror('argval', 'unknown argument name')
        return
    end if

    select type (opt => argsDB%args(iarg)%arg)
        class is (arg_real)
            if (opt%is_set == 0) then
                call argsDB%error%raise_error('value', 'unset', &
                    'value of argument is not set')
                return
            end if
            result = opt%value
        class default
            call argsDB%error%raise_error('value', 'type', &
                'argument is not a real scalar')
            return
    end select

end procedure argsdb_getval_real64

! ======================================================================

module procedure argsdb_getvals_char
    !! Get list of character-type values from argname.
    !!
    !! Gets values related to a list character-type option.

    integer :: iarg

    iarg = argsDB%get_argname_id(argname)
    if (iarg == 0) then
        call argsDB%error%raise_deverror('argval', 'unknown argument name')
        return
    end if

    select type (opt => argsDB%args(iarg)%arg)
        class is (arg_chars)
            if (opt%is_set == 0) then
                call argsDB%error%raise_error('value', 'unset', &
                    'value of argument is not set')
                return
            end if
            result = opt%values
        class default
            call argsDB%error%raise_error('value', 'type', &
                'argument is not a list of character strings')
            return
    end select

end procedure argsdb_getvals_char

! ======================================================================

module procedure argsdb_getvals_int32
    !! Get list of integer-type values from argname.
    !!
    !! Gets values related to a scalar integer-type option.

    integer :: i, iarg

    iarg = argsDB%get_argname_id(argname)
    if (iarg == 0) then
        call argsDB%error%raise_deverror('argval', 'unknown argument name')
        return
    end if

    select type (opt => argsDB%args(iarg)%arg)
        class is (arg_ints)
            if (opt%is_set == 0) then
                call argsDB%error%raise_error('value', 'unset', &
                    'value of argument is not set')
                return
            end if
            do i = 1, size(opt%values)
                if (abs(opt%values(i)) > huge(1_int32)) then
                    call argsDB%error%raise_error('value', 'convert', &
                        'range insufficient to represent values of argument')
                    return
                end if
            end do
            result = int(opt%values, int32)
        class default
            call argsDB%error%raise_error('value', 'type', &
                'argument is not a list of integers')
            return
    end select

end procedure argsdb_getvals_int32

! ======================================================================

module procedure argsdb_getvals_int64
    !! Get list of integer-type values from argname.
    !!
    !! Gets values related to a scalar integer-type option.

    integer :: iarg

    iarg = argsDB%get_argname_id(argname)
    if (iarg == 0) then
        call argsDB%error%raise_deverror('argval', 'unknown argument name')
        return
    end if

    select type (opt => argsDB%args(iarg)%arg)
        class is (arg_ints)
            if (opt%is_set == 0) then
                call argsDB%error%raise_error('value', 'unset', &
                    'value of argument is not set')
                return
            end if
            result = opt%values
        class default
            call argsDB%error%raise_error('value', 'type', &
                'argument is not a list of integers')
            return
    end select

end procedure argsdb_getvals_int64

! ======================================================================

module procedure argsdb_getvals_real32
    !! Get list of real-type values from argname.
    !!
    !! Gets values related to a scalar real-type option.

    integer :: i, iarg

    iarg = argsDB%get_argname_id(argname)
    if (iarg == 0) then
        call argsDB%error%raise_deverror('argval', 'unknown argument name')
        return
    end if

    select type (opt => argsDB%args(iarg)%arg)
        class is (arg_reals)
            if (opt%is_set == 0) then
                call argsDB%error%raise_error('value', 'unset', &
                    'value of argument is not set')
                return
            end if
            do i = 1, size(opt%values)
                if (abs(opt%values(i)) > huge(1_real32)) then
                    call argsDB%error%raise_error('value', 'convert', &
                        'precision insufficient to represent value of &
                        &argument')
                    return
                end if
            end do
            result = real(opt%values, real32)
        class default
            call argsDB%error%raise_error('value', 'type', &
                'argument is not a list of real numbers')
            return
    end select

end procedure argsdb_getvals_real32

! ======================================================================

module procedure argsdb_getvals_real64
    !! Get list of real-type values from argname.
    !!
    !! Gets values related to a scalar real-type option.

    integer :: iarg

    iarg = argsDB%get_argname_id(argname)
    if (iarg == 0) then
        call argsDB%error%raise_deverror('argval', 'unknown argument name')
        return
    end if

    select type (opt => argsDB%args(iarg)%arg)
        class is (arg_reals)
            if (opt%is_set == 0) then
                call argsDB%error%raise_error('value', 'unset', &
                    'value of argument is not set')
                return
            end if
            result = opt%values
        class default
            call argsDB%error%raise_error('value', 'type', &
                'argument is not a list of real numbers')
            return
    end select

end procedure argsdb_getvals_real64

! ======================================================================

module procedure argsdb_val_is_userset
    !! Return True is value set by user.
    !!
    !! Checks if value set by user and returns True in that case.

    integer :: iarg

    res = .false.

    iarg = argsDB%get_argname_id(argname)
    if (iarg == 0) then
        call argsDB%error%raise_deverror('argval', 'unknown argument name')
        return
    end if

    res = argsDB%args(iarg)%arg%is_set == 2

end procedure argsdb_val_is_userset

! ======================================================================

end submodule parse_cmdline_getval
