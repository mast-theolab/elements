submodule (parse_cmdline) parse_cmdline_addarg

    implicit none

contains

! ======================================================================

module procedure set_arg_names
    !! Set argument forms and labels for the parser.
    !!
    !! Sets the forms and names of the arguments.
    !! To facilitate parsing, the argument names are expected to follow
    !!   UNIX/GNU-like conventions:
    !! - short argument names must be a single letter, they can be
    !!   preceded by the argument symbol (typically: -/+).
    !!   Otherwise, '-' is added.
    !!   It is expected that short arguments can then be combined in the
    !!     commandline.
    !! - long argument names must be preceded by two "prefixes"
    !! - short and long arguments can be omitted, in which case it is
    !!   expected that the argument is assumed positional and any value
    !!   not related to any optional argument is assigned to it.
    !! - `label` is used internally for reference.  It is mandatory for
    !!   positional arguments.  For optional arguments, the label is based
    !!   on the long name, short if no long name provided, and overridden
    !!   if `label` is provided.
    !!
    !! Hence, either `short`/`long` must be provided or `label`.
    !! Optional arguments (starting with "prefixes") can be set as
    !!   required with `required` set to True (it is automatic by
    !!   default).
    integer :: i, larg
    character(len=1) :: prefix
    character(len=:), allocatable :: label_
    character(len=:), allocatable :: arg

    res = InitError()

    if (present(short) .or. present(long)) then
        if (present(long)) then
            if (len_trim(long) == 0) then
                call RaiseError(res, 'Empty version of the long form')
                return
            end if
            arg = trim(long)
            larg = len(arg)
            i = 1
            if (findstr(prefixes, arg(1:1)) > 0) then
                prefix = arg(1:1)
                do
                    i = i + 1
                    if (i > larg) then
                        call RaiseError(res, 'Empty long form')
                        return
                    else if (i > 3) then
                        call RaiseError(res, 'Too many prefixes')
                        return
                    end if
                    if (findstr(prefixes, arg(i:i)) > 0) then
                        if (arg(i:i) /= prefix) then
                            call RaiseError(res, &
                                'Inconsistency in the prefixes')
                            return
                        end if
                    else if (i == 2) then
                        call RaiseError(res, &
                            'Two prefixes needed for long forms')
                        return
                    else
                        exit
                    end if
                end do
                this%long_name = arg
                label_ = this%long_name(i:)
            else
                allocate(character(len=larg+2) :: this%long_name)
                label_ = arg
                this%long_name = '--' // arg
            end if
        end if
        if (present(short)) then
            arg = trim(short)
            larg = len(arg)
            select case (larg)
                case(0)
                    call RaiseError(res, 'Empty version of the short form')
                    return
                case(1)
                    if (findstr(prefixes, arg(1:1)) > 0) then
                        call RaiseError(res, 'Empty short form')
                        return
                    end if
                    allocate(character(len=2) :: this%short_name)
                    this%short_name = '-' // arg
                    if (.not.allocated(label_)) label_ = arg
                case(2)
                    if (findstr(prefixes, arg(1:1)) == 0) then
                        call RaiseError(res, 'Unsupported prefix')
                        return
                    end if
                    this%short_name = arg
                    if (.not.allocated(label_)) label_ = arg(2:2)
                case default
                    call RaiseError(res, &
                        'Short version should be 1 character long.')
                    return
            end select
        end if
        if (present(label)) label_ = trim(label)
        if (.not.allocated(label_)) then
            call RaiseError(res, 'DEVERR: `label` ended up uninitialized.')
            return
        else
            this%label = label_
        end if
        if (present(required)) then
            this%is_req = required
        else
            this%is_req = .False.
        end if
        this%is_pos = .False.
    else
        if (.not.present(label)) then
            call RaiseError(res, 'Missing label for positional argument.')
            return
        end if
        this%label = trim(label)
        if (present(required)) then
            if (.not.required) then
                call RaiseError(res, 'Positional arguments cannot be optional')
                return
            end if
        end if
        this%is_req = .True.
        this%is_pos = .True.
    end if
end procedure set_arg_names

! ======================================================================

module procedure set_arg_helpmsg
    !! Set the help message for the argument.
    !!
    !! Sets the help message to be displayed for a given argument.

    err = InitError()

    if (len_trim(msg) == 0) then
        call RaiseError(err, 'Missing help message')
        return
    end if

    this%help_msg = trim(msg)

end procedure set_arg_helpmsg

! ======================================================================

module procedure add_argument_int
    !! Add an integer-type argument to the arguments DB.
    !!
    !! Adds an integer-type argument to the database of arguments.
    !! The argument can expect one or more values, the routine will
    !!   properly assign the right internal type based on the
    !!   parameters.
    !! Supported type of arguments:
    !! - 'scalar'/'int': scalar integer
    !! - 'counter': scalar integer, incremented each time
    !! - 'list': list of integers
    !! - 'range': special list, which expect 1 to 3 elements

    integer(int64) :: ival
    logical :: is_scalar
    character(len=1024) :: errmsg
    character(len=MAX_ARGLEN) :: argname
    type(ArgInt), target :: arg_val
    type(ArgIntList), target :: arg_arr
    class(ArgObj), pointer :: new_arg
    class(BaseException), allocatable :: suberr

    ! Initialize error
    this%error = InitError()

    ! Check argument type and basic pre-processing
    select case (locase(argtype))
        case('int', 'scalar')
            is_scalar = .True.
            new_arg => arg_val
        case('counter')
            arg_val%add_value = .True.
            arg_val%min_num = 0
            arg_val%max_num = 0
            if (present(const_value)) then
                if (is_number(const_value)) then
                    arg_val%const = to_int64(const_value)
                else
                    call RaiseArgError(this%error, 'const_value', &
                                       'Not a number')
                    return
                end if
            else
                arg_val%const = 1
            end if
            is_scalar = .True.
            new_arg => arg_val
        case('list')
            arg_arr%max_num = huge(1)
            is_scalar = .False.
            new_arg => arg_arr
        case('range')
            arg_arr%max_num = 3
            is_scalar = .False.
            new_arg => arg_arr
        case default
            call RaiseArgError(this%error, 'argtype', &
                               'Unrecognized argument type')
            return
    end select

    ! Definition of the argument name/identifier
    suberr = new_arg%set_arg(this%prefixes, shortname, longname, label, &
                             required)
    if (suberr%raised()) then
        select type(suberr)
            class is (Error)
                write(errmsg, '(a,a," Source error: ",a)') &
                    'Error encountered in definition of argument', &
                    new_line(' '), trim(suberr%msg())
                call RaiseError(this%error, errmsg)
                return
            class default
                call RaiseError(this%error, &
                    'Unknown error met in argument definition')
                return
        end select
    end if

    if (this%chk_name_overlap(new_arg, argname)) then
        write(errmsg, '("''",a,"'' has already been defined.")') trim(argname)
        call RaiseArgError(this%error, 'names', errmsg)
        return
    end if

    ! Check there is space in DB
    this%nargs = this%nargs + 1
    if (this%nargs > MAX_ARGS) then
        errmsg = 'Number of arguments exceeds available storage in DB.'
        call RaiseError(this%error, errmsg)
        return
    end if

    ! Definition of the help message
    if (present(help)) then
        suberr = new_arg%set_help(help)
        if (suberr%raised()) then
            select type(suberr)
                class is (Error)
                    write(errmsg, '(a,a,"Source error:",a)') &
                        'Error encountered in setup of help message', &
                        new_line(' '), trim(suberr%msg())
                    call RaiseArgError(this%error, errmsg)
                    return
                class default
                    call RaiseError(this%error, &
                        'Unknown error met while setting up the help message')
                    return
            end select
        end if
    else
        suberr = new_arg%set_help('Information not available.')
    end if

    ! Setup of the parameters
    if (is_scalar) then
        if (present(def_value)) then
            if (is_number(def_value)) then
                arg_val%value = to_int64(def_value)
                arg_val%is_set = 1
            else
                call RaiseArgError(this%error, 'def_value', 'Not a number')
                return
            end if
        end if
        if (present(min_value)) then
            if (is_number(min_value)) then
                arg_val%min_ok = to_int64(min_value)
            else
                call RaiseArgError(this%error, 'min_value', 'Not a number')
                return
            end if
        end if
        if (present(max_value)) then
            if (is_number(max_value)) then
                arg_val%max_ok = to_int64(max_value)
            else
                call RaiseArgError(this%error, 'max_value', 'Not a number')
                return
            end if
        end if
        if (present(const_value)) then
            if (.not.is_number(const_value)) then
                call RaiseArgError(this%error, 'const_value', 'Not a number')
                return
            end if
            ival = to_int64(const_value)
            if (ival < arg_val%min_ok &
                .or. ival > arg_val%max_ok) then
                errmsg = 'Constant value inconsistent with values permitted &
                    &for argument'
                call RaiseArgError(this%error, 'const_value', errmsg)
                return
            end if
            arg_val%const = ival
            arg_val%min_num = 0
            arg_val%max_num = 0
        end if
        if (present(add_value)) arg_val%add_value = add_value
        this%args(this%nargs)%arg = arg_val
        if (new_arg%is_pos) call updDB_pos_arg(1, 1)
    else
        if (present(min_value)) then
            if (is_number(min_value)) then
                arg_arr%min_ok = to_int64(min_value)
            else
                call RaiseArgError(this%error, 'min_value', 'Not a number')
                return
            end if
        end if
        if (present(max_value)) then
            if (is_number(max_value)) then
                arg_arr%max_ok = to_int64(max_value)
            else
                call RaiseArgError(this%error, 'max_value', 'Not a number')
                return
            end if
        end if
        if (present(min_nvals)) then
            if (.not.is_number(min_nvals)) then
                call RaiseArgError(this%error, 'min_nvals', 'Not a number')
                return
            end if
            ival = to_int64(min_nvals)
            if (ival <= 0) then
                errmsg = 'Minimum number of values cannot be negative.'
                call RaiseArgError(this%error, 'min_nvals', errmsg)
                return
            end if
            arg_arr%min_num = ival
        end if
        if (present(max_nvals)) then
            if (is_number(max_nvals)) then
                arg_arr%max_num = to_int64(max_nvals)
            else
                call RaiseArgError(this%error, 'max_nvals', 'Not a number')
                return
            end if
        end if
        if (arg_arr%min_num > arg_arr%max_num) then
            errmsg = 'Minimum number of values expected larger than the &
                &maximum'
            call RaiseArgError(this%error, 'nvals', errmsg)
            return
        end if
        this%args(this%nargs)%arg = arg_arr
        if (new_arg%is_pos) &
            call updDB_pos_arg(arg_arr%min_num, arg_arr%max_num)
    end if

contains
    subroutine updDB_pos_arg(min_nvals, max_nvals)
        !! Update DB with info regarding positional arguments
        integer, intent(in) :: min_nvals
        !! Minimum number of values.
        integer, intent(in) :: max_nvals
        !! Maximum number of values.

        integer :: iarg

        this%nargs_pos = this%nargs_pos + 1
        if (this%nargs_pos > MAX_ARGS) then
            call RaiseError(this%error, 'Too many positional arguments')
            return
        end if

        this%iargs_pos(1,this%nargs_pos) = this%nargs
        if (this%nargs_pos == 1) then
            this%iargs_pos(2,this%nargs_pos) = 1
        else
            this%iargs_pos(2,this%nargs_pos) = &
                abs(this%iargs_pos(3,this%nargs_pos-1)) + 1
        end if

        if (min_nvals < max_nvals) then
            do iarg = 1, this%nargs_pos-1
                if (this%iargs_pos(3,iarg) < 0) then
                    errmsg = 'Only one positional argument with a variable &
                        &number of values allowed in the DB.'
                    call RaiseError(this%error, errmsg)
                    return
                end if
            end do
            this%iargs_pos(3,this%nargs_pos) = &
                -(this%iargs_pos(2,this%nargs_pos) + min_nvals - 1)
        else
            this%iargs_pos(3,this%nargs_pos) = &
                this%iargs_pos(2,this%nargs_pos) + min_nvals - 1
        end if
    end subroutine updDB_pos_arg

end procedure add_argument_int

! ======================================================================

module procedure add_argument_real
    !! Add an real-type argument to the arguments DB.
    !!
    !! Adds an real-type argument to the database of arguments.
    !! The argument can expect one or more values, the routine will
    !!   properly assign the right internal type based on the
    !!   parameters.
    !! Supported type of arguments:
    !! - 'scalar'/'real': scalar real
    !! - 'counter': scalar real, incremented each time
    !! - 'list': list of reals
    !! - 'range': special list, which expect 1 to 3 elements

    integer :: ival
    real(real64) :: rval
    logical :: is_scalar
    character(len=1024) :: errmsg
    character(len=MAX_ARGLEN) :: argname
    type(ArgReal), target :: arg_val
    type(ArgRealList), target :: arg_arr
    class(ArgObj), pointer :: new_arg
    class(BaseException), allocatable :: suberr

    ! Initialize error
    this%error = InitError()

    ! Check argument type and basic pre-processing
    select case (locase(argtype))
        case('real', 'scalar')
            is_scalar = .True.
            new_arg => arg_val
        case('counter')
            arg_val%add_value = .True.
            arg_val%min_num = 0
            arg_val%max_num = 0
            if (present(const_value)) then
                if (is_number(const_value)) then
                    arg_val%const = to_real64(const_value)
                else
                    call RaiseArgError(this%error, 'const_value', &
                                       'Not a number')
                    return
                end if
            else
                arg_val%const = 1
            end if
            is_scalar = .True.
            new_arg => arg_val
        case('list')
            arg_arr%max_num = huge(1)
            is_scalar = .False.
            new_arg => arg_arr
        case('range')
            arg_arr%max_num = 3
            is_scalar = .False.
            new_arg => arg_arr
        case default
            call RaiseArgError(this%error, 'argtype', &
                               'Unrecognized argument type')
            return
    end select

    ! Definition of the argument name/identifier
    suberr = new_arg%set_arg(this%prefixes, shortname, longname, label, &
                             required)
    if (suberr%raised()) then
        select type(suberr)
            class is (Error)
                write(errmsg, '(a,a," Source error: ",a)') &
                    'Error encountered in definition of argument', &
                    new_line(' '), trim(suberr%msg())
                call RaiseError(this%error, errmsg)
                return
            class default
                call RaiseError(this%error, &
                    'Unknown error met in argument definition')
                return
        end select
    end if

    if (this%chk_name_overlap(new_arg, argname)) then
        write(errmsg, '("''",a,"'' has already been defined.")') trim(argname)
        call RaiseArgError(this%error, 'names', errmsg)
        return
    end if

    ! Check there is space in DB
    this%nargs = this%nargs + 1
    if (this%nargs > MAX_ARGS) then
        errmsg = 'Number of arguments exceeds available storage in DB.'
        call RaiseError(this%error, errmsg)
        return
    end if

    ! Definition of the help message
    if (present(help)) then
        suberr = new_arg%set_help(help)
        if (suberr%raised()) then
            select type(suberr)
                class is (Error)
                    write(errmsg, '(a,a,"Source error:",a)') &
                        'Error encountered in setup of help message', &
                        new_line(' '), trim(suberr%msg())
                    call RaiseArgError(this%error, errmsg)
                    return
                class default
                    call RaiseError(this%error, &
                        'Unknown error met while setting up the help message')
                    return
            end select
        end if
    else
        suberr = new_arg%set_help('Information not available.')
    end if

    ! Setup of the parameters
    if (is_scalar) then
        if (present(def_value)) then
            if (is_number(def_value)) then
                arg_val%value = to_real64(def_value)
                arg_val%is_set = 1
            else
                call RaiseArgError(this%error, 'def_value', 'Not a number')
                return
            end if
        end if
        if (present(min_value)) then
            if (is_number(min_value)) then
                arg_val%min_ok = to_real64(min_value)
            else
                call RaiseArgError(this%error, 'min_value', 'Not a number')
                return
            end if
        end if
        if (present(max_value)) then
            if (is_number(max_value)) then
                arg_val%max_ok = to_real64(max_value)
            else
                call RaiseArgError(this%error, 'max_value', 'Not a number')
                return
            end if
        end if
        if (present(const_value)) then
            if (.not.is_number(const_value)) then
                call RaiseArgError(this%error, 'const_value', 'Not a number')
                return
            end if
            rval = to_real64(const_value)
            if (rval < arg_val%min_ok &
                .or. rval > arg_val%max_ok) then
                errmsg = 'Constant value inconsistent with values permitted &
                    &for argument'
                call RaiseArgError(this%error, 'const_value', errmsg)
                return
            end if
            arg_val%const = rval
            arg_val%min_num = 0
            arg_val%max_num = 0
        end if
        if (present(add_value)) arg_val%add_value = add_value
        this%args(this%nargs)%arg = arg_val
        if (new_arg%is_pos) call updDB_pos_arg(1, 1)
    else
        if (present(min_value)) then
            if (is_number(min_value)) then
                arg_arr%min_ok = to_real64(min_value)
            else
                call RaiseArgError(this%error, 'min_value', 'Not a number')
                return
            end if
        end if
        if (present(max_value)) then
            if (is_number(max_value)) then
                arg_arr%max_ok = to_real64(max_value)
            else
                call RaiseArgError(this%error, 'max_value', 'Not a number')
                return
            end if
        end if
        if (present(min_nvals)) then
            if (.not.is_number(min_nvals)) then
                call RaiseArgError(this%error, 'min_nvals', 'Not a number')
                return
            end if
            ival = to_int64(min_nvals)
            if (ival <= 0) then
                errmsg = 'Minimum number of values cannot be negative.'
                call RaiseArgError(this%error, 'min_nvals', errmsg)
                return
            end if
            arg_arr%min_num = ival
        end if
        if (present(max_nvals)) then
            if (is_number(max_nvals)) then
                arg_arr%max_num = to_int64(max_nvals)
            else
                call RaiseArgError(this%error, 'max_nvals', 'Not a number')
                return
            end if
        end if
        if (arg_arr%min_num > arg_arr%max_num) then
            errmsg = 'Minimum number of values expected larger than the &
                &maximum'
            call RaiseArgError(this%error, 'nvals', errmsg)
            return
        end if
        this%args(this%nargs)%arg = arg_arr
        if (new_arg%is_pos) &
            call updDB_pos_arg(arg_arr%min_num, arg_arr%max_num)
    end if

contains
    subroutine updDB_pos_arg(min_nvals, max_nvals)
        !! Update DB with info regarding positional arguments
        integer, intent(in) :: min_nvals
        !! Minimum number of values.
        integer, intent(in) :: max_nvals
        !! Maximum number of values.

        integer :: iarg

        this%nargs_pos = this%nargs_pos + 1
        if (this%nargs_pos > MAX_ARGS) then
            call RaiseError(this%error, 'Too many positional arguments')
            return
        end if

        this%iargs_pos(1,this%nargs_pos) = this%nargs
        if (this%nargs_pos == 1) then
            this%iargs_pos(2,this%nargs_pos) = 1
        else
            this%iargs_pos(2,this%nargs_pos) = &
                abs(this%iargs_pos(3,this%nargs_pos-1)) + 1
        end if

        if (min_nvals < max_nvals) then
            do iarg = 1, this%nargs_pos-1
                if (this%iargs_pos(3,iarg) < 0) then
                    errmsg = 'Only one positional argument with a variable &
                        &number of values allowed in the DB.'
                    call RaiseError(this%error, errmsg)
                    return
                end if
            end do
            this%iargs_pos(3,this%nargs_pos) = &
                -(this%iargs_pos(2,this%nargs_pos) + min_nvals - 1)
        else
            this%iargs_pos(3,this%nargs_pos) = &
                this%iargs_pos(2,this%nargs_pos) + min_nvals - 1
        end if
    end subroutine updDB_pos_arg

end procedure add_argument_real

! ======================================================================

module procedure add_argument_bool
    !! Add an logical-type argument to the arguments DB.
    !!
    !! Adds an logical-type argument to the database of arguments.
    !! The argument can expect one or more values, the routine will
    !!   properly assign the right internal type based on the
    !!   parameters.
    !! Supported type of arguments:
    !! - 'store_true': store True when found.
    !! - 'store_false': store False when found.

    character(len=1024) :: errmsg
    character(len=MAX_ARGLEN) :: argname
    type(ArgBool), target :: arg_val
    class(ArgObj), pointer :: new_arg
    class(BaseException), allocatable :: suberr

    ! Initialize error
    this%error = InitError()

    ! Check argument type and basic pre-processing
    select case (locase(argtype))
        case('store_true')
            arg_val%const = .True.
        case('store_false')
            arg_val%const = .False.
        case default
            call RaiseArgError(this%error, 'argtype', &
                               'Unrecognized argument type')
            return
    end select
    new_arg => arg_val
    arg_val%min_num = 0
    arg_val%max_num = 0
    arg_val%is_pos = .False.

    if (.not.present(shortname) .and. .not.present(longname)) then
        call RaiseArgError(this%error, 'names', &
                           'Boolean arguments cannot be positional')
        return
    end if

    ! Definition of the argument name/identifier
    suberr = new_arg%set_arg(this%prefixes, shortname, longname, label, &
                             required)
    if (suberr%raised()) then
        select type(suberr)
            class is (Error)
                write(errmsg, '(a,a," Source error: ",a)') &
                    'Error encountered in definition of argument', &
                    new_line(' '), trim(suberr%msg())
                call RaiseError(this%error, errmsg)
                return
            class default
                call RaiseError(this%error, &
                    'Unknown error met in argument definition')
                return
        end select
    end if

    if (this%chk_name_overlap(new_arg, argname)) then
        write(errmsg, '("''",a,"'' has already been defined.")') trim(argname)
        call RaiseArgError(this%error, 'names', errmsg)
        return
    end if

    ! Check there is space in DB
    ! Finalization: insert in DB
    this%nargs = this%nargs + 1
    if (this%nargs > MAX_ARGS) then
        errmsg = 'Number of arguments exceeds available storage in DB.'
        call RaiseError(this%error, errmsg)
        return
    end if

    if (new_arg%label == 'help') this%iarg_help = this%nargs

    ! Definition of the help message
    if (present(help)) then
        suberr = new_arg%set_help(help)
        if (suberr%raised()) then
            select type(suberr)
                class is (Error)
                    write(errmsg, '(a,a,"Source error:",a)') &
                        'Error encountered in setup of help message', &
                        new_line(' '), trim(suberr%msg())
                    call RaiseArgError(this%error, errmsg)
                    return
                class default
                    call RaiseError(this%error, &
                        'Unknown error met while setting up the help message')
                    return
            end select
        end if
    else
        suberr = new_arg%set_help('Information not available.')
    end if

    ! Setup of the parameters
    if (present(def_value)) then
        arg_val%value = def_value
        arg_val%is_set = 1
    end if

    this%args(this%nargs)%arg = arg_val

end procedure add_argument_bool

! ======================================================================

module procedure add_argument_char
    !! Add an character-type argument to the arguments DB.
    !!
    !! Adds an character-type argument to the database of arguments.
    !! The argument can expect one or more values, the routine will
    !!   properly assign the right internal type based on the
    !!   parameters.
    !! Supported type of arguments:
    !! - 'string': scalar character string
    !! - 'list': list of strings

    integer :: ival
    logical :: is_scalar
    character(len=1024) :: errmsg
    character(len=MAX_ARGLEN) :: argname
    type(ArgChar), target :: arg_val
    type(ArgCharList), target :: arg_arr
    class(ArgObj), pointer :: new_arg
    class(BaseException), allocatable :: suberr

    ! Initialize error
    this%error = InitError()

    ! Check argument type and basic pre-processing
    select case (locase(argtype))
        case('string', 'scalar')
            is_scalar = .True.
            new_arg => arg_val
        case('list')
            arg_arr%max_num = huge(1)
            is_scalar = .False.
            new_arg => arg_arr
        case default
            call RaiseArgError(this%error, 'argtype', &
                               'Unrecognized argument type')
            return
    end select

    ! Definition of the argument name/identifier
    suberr = new_arg%set_arg(this%prefixes, shortname, longname, label, &
                             required)
    if (suberr%raised()) then
        select type(suberr)
            class is (Error)
                write(errmsg, '(a,a," Source error: ",a)') &
                    'Error encountered in definition of argument', &
                    new_line(' '), trim(suberr%msg())
                call RaiseError(this%error, errmsg)
                return
            class default
                call RaiseError(this%error, &
                    'Unknown error met in argument definition')
                return
        end select
    end if

    if (this%chk_name_overlap(new_arg, argname)) then
        write(errmsg, '("''",a,"'' has already been defined.")') trim(argname)
        call RaiseArgError(this%error, 'names', errmsg)
        return
    end if

    ! Check there is space in DB
    ! Finalization: insert in DB
    this%nargs = this%nargs + 1
    if (this%nargs > MAX_ARGS) then
        errmsg = 'Number of arguments exceeds available storage in DB.'
        call RaiseError(this%error, errmsg)
        return
    end if

    ! Definition of the help message
    if (present(help)) then
        suberr = new_arg%set_help(help)
        if (suberr%raised()) then
            select type(suberr)
                class is (Error)
                    write(errmsg, '(a,a,"Source error:",a)') &
                        'Error encountered in setup of help message', &
                        new_line(' '), trim(suberr%msg())
                    call RaiseArgError(this%error, errmsg)
                    return
                class default
                    call RaiseError(this%error, &
                        'Unknown error met while setting up the help message')
                    return
            end select
        end if
    else
        suberr = new_arg%set_help('Information not available.')
    end if

    ! Setup of the parameters
    if (is_scalar) then
        if (present(def_value)) arg_val%value = trim(def_value)
        this%args(this%nargs)%arg = arg_val
        if (new_arg%is_pos) call updDB_pos_arg(1, 1)
    else
        if (present(min_nvals)) then
            if (.not.is_number(min_nvals)) then
                call RaiseArgError(this%error, 'min_nvals', 'Not a number')
                return
            end if
            ival = to_int64(min_nvals)
            if (ival <= 0) then
                errmsg = 'Minimum number of values cannot be negative.'
                call RaiseArgError(this%error, 'min_nvals', errmsg)
                return
            end if
            arg_arr%min_num = ival
        end if
        if (present(max_nvals)) then
            if (is_number(max_nvals)) then
                arg_arr%max_num = to_int64(max_nvals)
            else
                call RaiseArgError(this%error, 'max_nvals', 'Not a number')
                return
            end if
        end if
        if (arg_arr%min_num > arg_arr%max_num) then
            errmsg = 'Minimum number of values expected larger than the &
                &maximum'
            call RaiseArgError(this%error, 'nvals', errmsg)
            return
        end if
        this%args(this%nargs)%arg = arg_arr
        if (new_arg%is_pos) &
            call updDB_pos_arg(arg_arr%min_num, arg_arr%max_num)
    end if

contains
    subroutine updDB_pos_arg(min_nvals, max_nvals)
        !! Update DB with info regarding positional arguments
        integer, intent(in) :: min_nvals
        !! Minimum number of values.
        integer, intent(in) :: max_nvals
        !! Maximum number of values.

        integer :: iarg

        this%nargs_pos = this%nargs_pos + 1
        if (this%nargs_pos > MAX_ARGS) then
            call RaiseError(this%error, 'Too many positional arguments')
            return
        end if

        this%iargs_pos(1,this%nargs_pos) = this%nargs
        if (this%nargs_pos == 1) then
            this%iargs_pos(2,this%nargs_pos) = 1
        else
            this%iargs_pos(2,this%nargs_pos) = &
                abs(this%iargs_pos(3,this%nargs_pos-1)) + 1
        end if

        if (min_nvals < max_nvals) then
            do iarg = 1, this%nargs_pos-1
                if (this%iargs_pos(3,iarg) < 0) then
                    errmsg = 'Only one positional argument with a variable &
                        &number of values allowed in the DB.'
                    call RaiseError(this%error, errmsg)
                    return
                end if
            end do
            this%iargs_pos(3,this%nargs_pos) = &
                -(this%iargs_pos(2,this%nargs_pos) + min_nvals - 1)
        else
            this%iargs_pos(3,this%nargs_pos) = &
                this%iargs_pos(2,this%nargs_pos) + min_nvals - 1
        end if
    end subroutine updDB_pos_arg

end procedure add_argument_char

! ======================================================================

end submodule parse_cmdline_addarg