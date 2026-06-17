submodule (parse_cmdline) parse_cmdline_addarg

    use numeric, only: is_number, to_int, to_int64, to_real64

    implicit none

contains

! ======================================================================
! TYPE-BOUND PROCEDURES (DEFINITIONS) - CMDLINEARGSDB
! ======================================================================

module procedure argsdb_add_bool
    !! Add a logical-type argument to the arguments DB.
    !!
    !! Adds a logical-type argument to the database of arguments.
    !! The argument can expect one or more values, the routine will
    !!   properly assign the right internal type based on the
    !!   parameters.
    !! Supported type of arguments:
    !! - 'store_true': store True when found.
    !! - 'store_false': store False when found.

    logical :: add_ok
    character(len=1024) :: errmsg
    character(len=MAX_ARGLEN) :: argname
    character(len=:), allocatable :: err_cause
    type(arg_bool), target :: arg_val
    class(arg_obj), pointer :: new_arg

    ! Check argument type and basic pre-processing
    select case (locase(argtype))
        case('store_true')
            arg_val%const = .true.
        case('store_false')
            arg_val%const = .false.
        case default
            call argsDB%error%raise_deverror('argval', &
                'unrecognized type of logical argument')
            return
    end select
    new_arg => arg_val
    arg_val%min_num = 0
    arg_val%max_num = 0
    arg_val%is_pos = .false.

    if (.not.present(shortname) .and. .not.present(longname)) then
        call argsDB%error%raise_deverror('argval', &
            'logical arguments cannot be positional')
        return
    end if

    ! Definition of the argument name/identifier
    add_ok = new_arg%set_arg(argsDB%prefixes, shortname, longname, label, &
                             required)
    if (.not.add_ok) then
        call new_arg%error%info(err_cause)
        write(errmsg, '("unable to add argument num. ",i0)') argsDB%nargs+1
        call argsDB%error%raise_deverror('argval', errmsg, err_cause)
        return
    end if

    if (argsDB%chk_name_overlap(new_arg, argname)) then
        write(errmsg, '("''",a,"'' has already been defined.")') trim(argname)
        call argsDB%error%raise_deverror('argval', errmsg)
        return
    end if

    ! Check there is space in DB
    ! Finalization: insert in DB
    argsDB%nargs = argsDB%nargs + 1
    if (argsDB%nargs > MAX_ARGS) then
        call argsDB%error%raise_deverror('limit', &
            'number of arguments exceeds available storage in DB.')
        return
    end if

    if (new_arg%label == 'help') argsDB%iarg_help = argsDB%nargs

    ! Definition of the help message
    if (present(help)) then
        add_ok = new_arg%set_help(help)
        if (.not.add_ok) then
            call new_arg%error%info(err_cause)
            write(errmsg, '("unable to set up help of argument num. ",i0)') &
                argsDB%nargs
            call argsDB%error%raise_deverror('argval', errmsg, err_cause)
            return
        end if
    else
        add_ok = new_arg%set_help('Information not available.')
    end if

    ! Setup of the parameters
    if (present(def_value)) then
        arg_val%value = def_value
        arg_val%is_set = 1
    else
        arg_val%is_set = 0
    end if

    argsDB%args(argsDB%nargs)%arg = arg_val

end procedure argsdb_add_bool

! ======================================================================

module procedure argsdb_add_char
    !! Add a character-type argument to the arguments DB.
    !!
    !! Adds a character-type argument to the database of arguments.
    !! The argument can expect one or more values, the routine will
    !!   properly assign the right internal type based on the
    !!   parameters.
    !! Supported type of arguments:
    !! - 'string': scalar character string
    !! - 'list': list of strings

    integer :: ival
    logical :: add_ok, is_scalar
    character(len=1024) :: errmsg
    character(len=MAX_ARGLEN) :: argname
    character(len=:), allocatable :: err_cause
    type(arg_char), target :: arg_val
    type(arg_chars), target :: arg_arr
    class(arg_obj), pointer :: new_arg

    ! Check argument type and basic pre-processing
    select case (locase(argtype))
        case('string', 'scalar')
            is_scalar = .true.
            new_arg => arg_val
        case('list')
            arg_arr%max_num = huge(1)
            is_scalar = .false.
            new_arg => arg_arr
        case default
            call argsDB%error%raise_deverror('argval', &
                'unrecognized type of character argument')
            return
    end select

    ! Definition of the argument name/identifier
    add_ok = new_arg%set_arg(argsDB%prefixes, shortname, longname, label, &
                             required)
    if (.not.add_ok) then
        call new_arg%error%info(err_cause)
        write(errmsg, '("unable to add argument num. ",i0)') argsDB%nargs+1
        call argsDB%error%raise_deverror('argval', errmsg, err_cause)
        return
    end if

    if (argsDB%chk_name_overlap(new_arg, argname)) then
        write(errmsg, '("''",a,"'' has already been defined.")') trim(argname)
        call argsDB%error%raise_deverror('argval', errmsg)
        return
    end if

    ! Check there is space in DB
    argsDB%nargs = argsDB%nargs + 1
    if (argsDB%nargs > MAX_ARGS) then
        call argsDB%error%raise_deverror('limit', &
            'number of arguments exceeds available storage in DB.')
        return
    end if

    ! Definition of the help message
    if (present(help)) then
        add_ok = new_arg%set_help(help)
        if (.not.add_ok) then
            call new_arg%error%info(err_cause)
            write(errmsg, '("unable to set up help of argument num. ",i0)') &
                argsDB%nargs
            call argsDB%error%raise_deverror('argval', errmsg, err_cause)
            return
        end if
    else
        add_ok = new_arg%set_help('Information not available.')
    end if

    ! Setup of the parameters
    if (is_scalar) then
        if (present(def_value)) arg_val%value = trim(def_value)
        argsDB%args(argsDB%nargs)%arg = arg_val
        if (new_arg%is_pos) call argsdb_upd_pos_dbase(argsDB, 1, 1)
    else
        if (present(min_nvals)) then
            if (.not.is_number(min_nvals)) then
                call argsDB%error%raise_deverror('argval', &
                    '"min_nvals" is not a number')
                return
            end if
            ival = to_int(min_nvals)
            if (ival <= 0) then
                errmsg = 'minimum number of values cannot be negative.'
                call argsDB%error%raise_deverror('argval', errmsg)
                return
            end if
            arg_arr%min_num = ival
        end if
        if (present(max_nvals)) then
            if (is_number(max_nvals)) then
                arg_arr%max_num = to_int(max_nvals)
            else
                call argsDB%error%raise_deverror('argval', &
                    '"max_nvals" is not a number')
                return
            end if
        end if
        if (arg_arr%min_num > arg_arr%max_num) then
            errmsg = 'minimum number of values expected larger than the &
                &maximum'
            call argsDB%error%raise_deverror('argval', errmsg)
            return
        end if
        argsDB%args(argsDB%nargs)%arg = arg_arr
        if (new_arg%is_pos) &
            call argsdb_upd_pos_dbase(argsDB, arg_arr%min_num, arg_arr%max_num)
    end if

end procedure argsdb_add_char

! ======================================================================

module procedure argsdb_add_int
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

    integer :: ival
    integer(int64) :: ival64
    logical :: add_ok, is_scalar
    character(len=1024) :: errmsg
    character(len=MAX_ARGLEN) :: argname
    character(len=:), allocatable :: err_cause
    type(arg_int), target :: arg_val
    type(arg_ints), target :: arg_arr
    class(arg_obj), pointer :: new_arg

    ! Check argument type and basic pre-processing
    select case (locase(argtype))
        case('int', 'scalar')
            is_scalar = .true.
            new_arg => arg_val
        case('counter')
            arg_val%add_value = .true.
            arg_val%min_num = 0
            arg_val%max_num = 0
            if (present(const_value)) then
                if (is_number(const_value)) then
                    arg_val%const = to_int(const_value)
                else
                    call argsDB%error%raise_deverror('argval', &
                        '"const_value" is not a number.')
                    return
                end if
            else
                arg_val%const = 1
            end if
            is_scalar = .true.
            new_arg => arg_val
        case('list')
            arg_arr%max_num = huge(1)
            is_scalar = .false.
            new_arg => arg_arr
        case('range')
            arg_arr%max_num = 3
            is_scalar = .false.
            new_arg => arg_arr
        case default
            call argsDB%error%raise_deverror('argval', &
                'unrecognized type of integer argument')
            return
    end select

    ! Definition of the argument name/identifier
    add_ok = new_arg%set_arg(argsDB%prefixes, shortname, longname, label, &
                             required)
    if (.not.add_ok) then
        call new_arg%error%info(err_cause)
        write(errmsg, '("unable to add argument num. ",i0)') argsDB%nargs+1
        call argsDB%error%raise_deverror('argval', errmsg, err_cause)
        return
    end if

    if (argsDB%chk_name_overlap(new_arg, argname)) then
        write(errmsg, '("''",a,"'' has already been defined.")') trim(argname)
        call argsDB%error%raise_deverror('argval', errmsg)
        return
    end if

    ! Check there is space in DB
    argsDB%nargs = argsDB%nargs + 1
    if (argsDB%nargs > MAX_ARGS) then
        call argsDB%error%raise_deverror('limit', &
            'number of arguments exceeds available storage in DB.')
        return
    end if

    ! Definition of the help message
    if (present(help)) then
        add_ok = new_arg%set_help(help)
        if (.not.add_ok) then
            call new_arg%error%info(err_cause)
            write(errmsg, '("unable to set up help of argument num. ",i0)') &
                argsDB%nargs
            call argsDB%error%raise_deverror('argval', errmsg, err_cause)
            return
        end if
    else
        add_ok = new_arg%set_help('Information not available.')
    end if

    ! Setup of the parameters
    if (is_scalar) then
        if (present(def_value)) then
            if (is_number(def_value)) then
                arg_val%value = to_int(def_value)
                arg_val%is_set = 1
            else
                call argsDB%error%raise_deverror('argval', &
                    'default value is not a number')
                return
            end if
        else
            arg_val%is_set = 0
        end if
        if (present(min_value)) then
            if (is_number(min_value)) then
                arg_val%min_ok = to_int(min_value)
            else
                call argsDB%error%raise_deverror('argval', &
                    'minimum value is not a number')
                return
            end if
        end if
        if (present(max_value)) then
            if (is_number(max_value)) then
                arg_val%max_ok = to_int(max_value)
            else
                call argsDB%error%raise_deverror('argval', &
                    'maximum value is not a number')
                return
            end if
        end if
        if (present(const_value)) then
            if (.not.is_number(const_value)) then
                call argsDB%error%raise_deverror('argval', &
                    'constant value is not a number')
                return
            end if
            ival64 = to_int64(const_value)
            if (ival64 < arg_val%min_ok &
                .or. ival64 > arg_val%max_ok) then
                errmsg = 'constant value inconsistent with values permitted &
                    &for argument'
                call argsDB%error%raise_deverror('argval', errmsg)
                return
            end if
            arg_val%const = ival64
            arg_val%min_num = 0
            arg_val%max_num = 0
        end if
        if (present(add_value)) arg_val%add_value = add_value
        argsDB%args(argsDB%nargs)%arg = arg_val
        if (new_arg%is_pos) call argsdb_upd_pos_dbase(argsDB, 1, 1)
    else
        if (present(min_value)) then
            if (is_number(min_value)) then
                arg_arr%min_ok = to_int(min_value)
            else
                call argsDB%error%raise_deverror('argval', &
                    'minimum value is not a number')
                return
            end if
        end if
        if (present(max_value)) then
            if (is_number(max_value)) then
                arg_arr%max_ok = to_int(max_value)
            else
                call argsDB%error%raise_deverror('argval', &
                    'maximum value is not a number')
                return
            end if
        end if
        if (present(min_nvals)) then
            if (.not.is_number(min_nvals)) then
                call argsDB%error%raise_deverror('argval', &
                    '"min_nvals" is not a number')
                return
            end if
            ival = to_int(min_nvals)
            if (ival <= 0) then
                errmsg = 'minimum number of values cannot be negative.'
                call argsDB%error%raise_deverror('argval', errmsg)
                return
            end if
            arg_arr%min_num = ival
        end if
        if (present(max_nvals)) then
            if (is_number(max_nvals)) then
                arg_arr%max_num = to_int(max_nvals)
            else
                call argsDB%error%raise_deverror('argval', &
                    '"max_nvals" is not a number')
                return
            end if
        end if
        if (arg_arr%min_num > arg_arr%max_num) then
            errmsg = 'minimum number of values expected larger than the &
                &maximum'
            call argsDB%error%raise_deverror('argval', errmsg)
            return
        end if
        argsDB%args(argsDB%nargs)%arg = arg_arr
        if (new_arg%is_pos) &
            call argsdb_upd_pos_dbase(argsDB, arg_arr%min_num, arg_arr%max_num)
    end if

end procedure argsdb_add_int

! ======================================================================

module procedure argsdb_add_real
    !! Add a real-type argument to the arguments DB.
    !!
    !! Adds a real-type argument to the database of arguments.
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
    logical :: add_ok, is_scalar
    character(len=1024) :: errmsg
    character(len=MAX_ARGLEN) :: argname
    character(len=:), allocatable :: err_cause
    type(arg_real), target :: arg_val
    type(arg_reals), target :: arg_arr
    class(arg_obj), pointer :: new_arg

    ! Check argument type and basic pre-processing
    select case (locase(argtype))
        case('real', 'scalar')
            is_scalar = .true.
            new_arg => arg_val
        case('counter')
            arg_val%add_value = .true.
            arg_val%min_num = 0
            arg_val%max_num = 0
            if (present(const_value)) then
                if (is_number(const_value)) then
                    arg_val%const = to_real64(const_value)
                else
                    call argsDB%error%raise_deverror('argval', &
                        '"const_value" is not a number.')
                    return
                end if
            else
                arg_val%const = 1
            end if
            is_scalar = .true.
            new_arg => arg_val
        case('list')
            arg_arr%max_num = huge(1)
            is_scalar = .false.
            new_arg => arg_arr
        case('range')
            arg_arr%max_num = 3
            is_scalar = .false.
            new_arg => arg_arr
        case default
            call argsDB%error%raise_deverror('argval', &
                'unrecognized type of real argument')
            return
    end select

    ! Definition of the argument name/identifier
    add_ok = new_arg%set_arg(argsDB%prefixes, shortname, longname, label, &
                             required)
    if (.not.add_ok) then
        call new_arg%error%info(err_cause)
        write(errmsg, '("unable to add argument num. ",i0)') argsDB%nargs+1
        call argsDB%error%raise_deverror('argval', errmsg, err_cause)
        return
    end if

    if (argsDB%chk_name_overlap(new_arg, argname)) then
        write(errmsg, '("''",a,"'' has already been defined.")') trim(argname)
        call argsDB%error%raise_deverror('argval', errmsg)
        return
    end if

    ! Check there is space in DB
    argsDB%nargs = argsDB%nargs + 1
    if (argsDB%nargs > MAX_ARGS) then
        call argsDB%error%raise_deverror('limit', &
            'number of arguments exceeds available storage in DB.')
        return
    end if

    ! Definition of the help message
    if (present(help)) then
        add_ok = new_arg%set_help(help)
        if (.not.add_ok) then
            call new_arg%error%info(err_cause)
            write(errmsg, '("unable to set up help of argument num. ",i0)') &
                argsDB%nargs
            call argsDB%error%raise_deverror('argval', errmsg, err_cause)
            return
        end if
    else
        add_ok = new_arg%set_help('Information not available.')
    end if

    ! Setup of the parameters
    if (is_scalar) then
        if (present(def_value)) then
            if (is_number(def_value)) then
                arg_val%value = to_real64(def_value)
                arg_val%is_set = 1
            else
                call argsDB%error%raise_deverror('argval', &
                    'default value is not a number')
                return
            end if
        else
            arg_val%is_set = 0
        end if
        if (present(min_value)) then
            if (is_number(min_value)) then
                arg_val%min_ok = to_real64(min_value)
            else
                call argsDB%error%raise_deverror('argval', &
                    'minimum value is not a number')
                return
            end if
        end if
        if (present(max_value)) then
            if (is_number(max_value)) then
                arg_val%max_ok = to_real64(max_value)
            else
                call argsDB%error%raise_deverror('argval', &
                    'maximum value is not a number')
                return
            end if
        end if
        if (present(const_value)) then
            if (.not.is_number(const_value)) then
                call argsDB%error%raise_deverror('argval', &
                    'constant value is not a number')
                return
            end if
            rval = to_real64(const_value)
            if (rval < arg_val%min_ok &
                .or. rval > arg_val%max_ok) then
                errmsg = 'constant value inconsistent with values permitted &
                    &for argument'
                call argsDB%error%raise_deverror('argval', errmsg)
                return
            end if
            arg_val%const = rval
            arg_val%min_num = 0
            arg_val%max_num = 0
        end if
        if (present(add_value)) arg_val%add_value = add_value
        argsDB%args(argsDB%nargs)%arg = arg_val
        if (new_arg%is_pos) call argsdb_upd_pos_dbase(argsDB, 1, 1)
    else
        if (present(min_value)) then
            if (is_number(min_value)) then
                arg_arr%min_ok = to_real64(min_value)
            else
                call argsDB%error%raise_deverror('argval', &
                    'minimum value is not a number')
                return
            end if
        end if
        if (present(max_value)) then
            if (is_number(max_value)) then
                arg_arr%max_ok = to_real64(max_value)
            else
                call argsDB%error%raise_deverror('argval', &
                    'maximum value is not a number')
                return
            end if
        end if
        if (present(min_nvals)) then
            if (.not.is_number(min_nvals)) then
                call argsDB%error%raise_deverror('argval', &
                    '"min_nvals" is not a number')
                return
            end if
            ival = to_int(min_nvals)
            if (ival <= 0) then
                errmsg = 'minimum number of values cannot be negative.'
                call argsDB%error%raise_deverror('argval', errmsg)
                return
            end if
            arg_arr%min_num = ival
        end if
        if (present(max_nvals)) then
            if (is_number(max_nvals)) then
                arg_arr%max_num = to_int(max_nvals)
            else
                call argsDB%error%raise_deverror('argval', &
                    '"max_nvals" is not a number')
                return
            end if
        end if
        if (arg_arr%min_num > arg_arr%max_num) then
            errmsg = 'minimum number of values expected larger than the &
                &maximum'
            call argsDB%error%raise_deverror('argval', errmsg)
            return
        end if
        argsDB%args(argsDB%nargs)%arg = arg_arr
        if (new_arg%is_pos) &
            call argsdb_upd_pos_dbase(argsDB, arg_arr%min_num, arg_arr%max_num)
    end if

end procedure argsdb_add_real

! ======================================================================
! TYPE-BOUND PROCEDURES (DEFINITIONS) - ARG_OBJ
! ======================================================================

module procedure argobj_set_helpmsg
    !! Set the help message for the argument.
    !!
    !! Sets the help message to be displayed for a given argument.
    done = .false.

    if (len_trim(msg) == 0) then
        call arg%error%raise_deverror('argval', 'missing help message')
        return
    end if

    arg%help_msg = trim(msg)

    done = .true.

end procedure argobj_set_helpmsg

! ======================================================================

module procedure argobj_set_names
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
    integer :: i, lname
    character(len=1) :: prefix
    character(len=:), allocatable :: arglabel, argname

    done = .false.

    if (present(short) .or. present(long)) then
        if (present(long)) then
            if (len_trim(long) == 0) then
                call arg%error%raise_deverror('argval', &
                    'long form of argument is empty')
                return
            end if
            argname = trim(long)
            lname = len(argname)
            i = 1
            if (findstr(prefixes, argname(1:1)) > 0) then
                prefix = argname(1:1)
                do
                    i = i + 1
                    if (i > lname) then
                        call arg%error%raise_deverror('argval', &
                            'long form of argument is empty')
                        return
                    else if (i > 3) then
                        call arg%error%raise_deverror('argval', &
                            'too many prefixes in long form of argument')
                        return
                    end if
                    if (findstr(prefixes, argname(i:i)) > 0) then
                        if (argname(i:i) /= prefix) then
                            call arg%error%raise_deverror('argval', &
                                'inconsistency in prefixes to argument')
                            return
                        end if
                    else if (i == 2) then
                        call arg%error%raise_deverror('argval', &
                            'two prefixes are needed for long forms')
                        return
                    else
                        exit
                    end if
                end do
                arg%long_name = argname
                arglabel = arg%long_name(i:)
            else
                allocate(character(len=lname+2) :: arg%long_name)
                arglabel = argname
                arg%long_name = '--' // argname
            end if
        end if
        if (present(short)) then
            argname = trim(short)
            lname = len(argname)
            select case (lname)
                case(0)
                    call arg%error%raise_deverror('argval', &
                        'short form of argument is empty')
                    return
                case(1)
                    if (findstr(prefixes, argname(1:1)) > 0) then
                        call arg%error%raise_deverror('argval', &
                            'short form of argument is empty')
                        return
                    end if
                    allocate(character(len=2) :: arg%short_name)
                    arg%short_name = '-' // argname
                    if (.not.allocated(arglabel)) arglabel = argname
                case(2)
                    if (findstr(prefixes, argname(1:1)) == 0) then
                        call arg%error%raise_deverror('argval', &
                            'unsupported prefix in short form of argument')
                        return
                    end if
                    arg%short_name = argname
                    if (.not.allocated(arglabel)) arglabel = argname(2:2)
                case default
                    call arg%error%raise_deverror('argval', &
                            'short version should be 1 character long')
                    return
            end select
        end if
        if (present(label)) arglabel = trim(label)
        if (.not.allocated(arglabel)) then
            call arg%error%raise_deverror('argval', &
                '"label" ended up uninitialized')
            return
        else
            arg%label = arglabel
        end if
        if (present(required)) then
            arg%is_req = required
        else
            arg%is_req = .false.
        end if
        arg%is_pos = .false.
    else
        if (.not.present(label)) then
            call arg%error%raise_deverror('argval', &
                'missing label for positional argument')
            return
        end if
        arg%label = trim(label)
        if (present(required)) then
            if (.not.required) then
                call arg%error%raise_deverror('argval', &
                    'positional arguments cannot be optional')
                return
            end if
        end if
        arg%is_req = .true.
        arg%is_pos = .true.
    end if
    done = .true.

end procedure argobj_set_names

! ======================================================================
! INTERNAL PROCEDURES
! ======================================================================

subroutine argsdb_upd_pos_dbase(argsDB, min_nvals, max_nvals)
    !! Update information on positional arguments in database.
    !!
    !! Updates database of positional arguments in `CmdLineArgsDB` instance.
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! `CmdLineArgsDB` instance.
    integer, intent(in) :: min_nvals
        !! Minimum number of values.
    integer, intent(in) :: max_nvals
        !! Maximum number of values.

    integer :: iarg

    argsDB%nargs_pos = argsDB%nargs_pos + 1
    if (argsDB%nargs_pos > MAX_ARGS) then
        call argsDB%error%raise_deverror('limit', &
            'too many positional arguments')
        return
    end if

    argsDB%iargs_pos(1,argsDB%nargs_pos) = argsDB%nargs
    if (argsDB%nargs_pos == 1) then
        argsDB%iargs_pos(2,argsDB%nargs_pos) = 1
    else
        argsDB%iargs_pos(2,argsDB%nargs_pos) = &
            abs(argsDB%iargs_pos(3,argsDB%nargs_pos-1)) + 1
    end if

    if (min_nvals < max_nvals) then
        do iarg = 1, argsDB%nargs_pos-1
            if (argsDB%iargs_pos(3,iarg) < 0) then
                call argsDB%error%raise_deverror('argval', &
                    'only one positional argument with a variable number of &
                    &values allowed in the DB')
                return
            end if
        end do
        argsDB%iargs_pos(3,argsDB%nargs_pos) = &
            -(argsDB%iargs_pos(2,argsDB%nargs_pos) + min_nvals - 1)
    else
        argsDB%iargs_pos(3,argsDB%nargs_pos) = &
            argsDB%iargs_pos(2,argsDB%nargs_pos) + min_nvals - 1
    end if

end subroutine argsdb_upd_pos_dbase

! ======================================================================

end submodule parse_cmdline_addarg
