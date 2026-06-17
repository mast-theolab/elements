submodule (parse_cmdline) parse_cmdline_parser

    implicit none

contains

! ======================================================================

module procedure argsdb_parse_args
    !! Parse a list of arguments in input or from commandline.
    !!
    !! Parses a list of arguments from the commandline or alternatively
    !!   from a list given in argument to the subroutine.

    integer :: i, iarg, iarg_pos, iargDB, ierr, ioff, ios, isep, j, nargs, &
        nargs_opt, nargs_pos, nvals
    real :: rval
    logical :: do_append, found, ok
    logical, dimension(:), allocatable :: todo
    character(len=MAX_ARGLEN) :: arg, key, val
    character(len=2) :: shortarg
    character(len=256) :: errmsg
    character(len=MAX_ARGLEN), dimension(:), allocatable :: arg_vals
    character(len=:), allocatable :: err_cause

    ! Find number of arguments and initialize arrays
    if (present(arglist)) then
        nargs = size(arglist)
    else
        nargs = command_argument_count()
    end if
    allocate(todo(nargs), arg_vals(nargs))
    todo = .true.
    arg_vals = ' '

    ! Preprocessing: check some information in the argument list:
    ! - presence of '--' to stop optional arguments
    ! - presence of short options, which may be fused. For this one, we
    !   check the first option, check if it expects value. If yes, not fused.
    !   Otherwise, we check all options are ok and set them.
    nargs_opt = 0
    do iarg = 1, nargs
        arg = get_arg(iarg, ierr)
        if (ierr /= 0) then
            write(errmsg, '("Error encountered while reading command-line &
                &argument num. ",i0)') iarg
            call argsDB%error%raise_error('data', 'unknown', errmsg)
            return
        end if
        if (argsDB%iarg_help > 0) then
            if (allocated(argsDB%args(argsDB%iarg_help)%arg%short_name)) then
                if (arg == argsDB%args(argsDB%iarg_help)%arg%short_name) &
                    call argsDB%print_help()
            end if
            if (allocated(argsDB%args(argsDB%iarg_help)%arg%long_name)) then
                if (arg == argsDB%args(argsDB%iarg_help)%arg%long_name) &
                    call argsDB%print_help()
            end if
        end if
        if (len_trim(arg) > 2 &
            .and. findstr(argsDB%prefixes, arg(1:1)) > 0 &
            .and. arg(2:2) /= arg(1:1) &
            .and. nargs_opt == 0) then
            ! first check if not a negative number.
            ! if this is the case, proceed to next argument.
            read(arg, *, iostat=ios) rval
            if (ios == 0) cycle
            ! not a number, continue
            shortarg = arg(:2)
            iargDB = chk_shortarg(shortarg)
            if (iargDB == 0) return
            if (argsDB%args(iargDB)%arg%max_num == 0) then
                ok = argsDB%args(iargDB)%arg%set_value()
                if (.not.ok) then
                    call argsDB%args(iargDB)%arg%error%info(err_cause)
                    if (argsDB%args(iargDB)%arg%error%has_type('arg')) then
                        write(errmsg, '("error while parsing ''",a,"''")') &
                            shortarg
                        call argsDB%error%raise_error('data', 'unknown', &
                            errmsg, err_cause)
                    else
                        write(errmsg, &
                              '("unexpected error met with ''",a,"''")') &
                            shortarg
                        call argsDB%error%raise_error('data', ' ', errmsg, &
                            err_cause)
                    end if
                    return
                end if
                do i = 3, len_trim(arg)
                    shortarg = arg(1:1) // arg(i:i)
                    iargDB = chk_shortarg(shortarg)
                    if (iargDB == 0) return
                    associate(opt => argsDB%args(iargDB)%arg)
                    if (opt%max_num == 0) then
                        ok = opt%set_value()
                        if (.not.ok) then
                            call opt%error%info(err_cause)
                            if (opt%error%has_type('arg')) then
                                write(errmsg, &
                                      '("error while parsing ''",a,"''")') &
                                    shortarg
                                call argsDB%error%raise_error('data', &
                                    'unknown', errmsg, err_cause)
                            else
                                write(errmsg, &
                                      '("unexpected error met with &
                                      &''",a,"''")') &
                                    shortarg
                                call argsDB%error%raise_error('data', ' ', &
                                    errmsg, err_cause)
                            end if
                            return
                        end if
                    else
                        write(errmsg, &
                              '("argument num. ",i0," mixes incompatible &
                              &options with and without values.")') iarg
                        call argsDB%error%raise_error('data', 'structure', &
                            errmsg)
                        return
                    end if
                    end associate
                end do
                todo(iarg) = .false.
            end if
        else if (arg == '--' .and. nargs_opt == 0) then
            nargs_opt = iarg - 1
        end if
    end do
    if (nargs_opt == 0) nargs_opt = nargs

    ! Now main processing
    ! We proceed in a matrix way: for each option, check the argument list
    ! This makes it easier to store once basic information on the parameters
    !   of an option and to collect multiple values associated to the same
    !   option.
    do iargDB = 1, argsDB%nargs
        associate(opt => argsDB%args(iargDB)%arg)
        if(opt%is_pos) cycle
        select type (opt)
            class is (arg_ints)
                do_append = .true.
            class is (arg_reals)
                do_append = .true.
            class is (arg_chars)
                do_append = .true.
            class default
                do_append = .false.
        end select
        nvals = 0
        do iarg = 1, nargs_opt
            if (todo(iarg)) then
                arg = get_arg(iarg, ierr)
                key = ' '
                val = ' '
                found = .false.
                if (findstr(argsDB%prefix_long, arg(:2)) > 0) then
                    isep = index(arg, '=')
                    if (isep > 0) then
                        key = arg(:isep-1)
                    else
                        key = arg
                    end if
                    found = opt%long_name == key
                else if (findstr(argsDB%prefix_short, arg(:1)) > 0) then
                    if (len_trim(arg) == 1) then
                        write(errmsg, &
                              '("incorrect argument in position ",i0)') &
                            iarg
                        call argsDB%error%raise_error('data', 'structure', &
                            errmsg)
                        return
                    end if
                    isep = index(arg, '=')
                    if (isep > 0) then
                        key = arg(:isep-1)
                    else
                        if (len_trim(arg) > 2) then
                            isep = 2
                            key = arg(:isep)
                        else
                            isep = 0
                            key = arg
                        end if
                    end if
                    found = opt%short_name == key
                end if
                if (found) then
                    todo(iarg) = .false.
                    if (isep > 0) then
                        if (opt%max_num == 0) then
                            write(errmsg, '("argument ''",a,"'' does not &
                                &accept values")') trim(key)
                            call argsDB%error%raise_error('value', &
                                'incompatible', errmsg)
                            return
                        end if
                        val = arg(isep+1:)
                    else if(opt%max_num > 0) then
                        val = get_arg(iarg+1, ierr)
                        todo(iarg+1) = .false.
                    else
                        val = ' '
                    end if
                    if (do_append) then
                        nvals = nvals + 1
                        arg_vals(nvals) = val
                    else
                        ok = opt%set_value(val)
                        if (.not.ok) then
                            call opt%error%info(err_cause)
                            if (opt%error%has_type('arg')) then
                                write(errmsg, &
                                      '("error while parsing ''",a,"''")') &
                                    trim(key)
                                call argsDB%error%raise_error('data', &
                                    'unknown', errmsg, err_cause)
                            else
                                write(errmsg, &
                                      '("unexpected error met with &
                                      &''",a,"''")') &
                                    trim(key)
                                call argsDB%error%raise_error('data', ' ', &
                                    errmsg, err_cause)
                            end if
                            return
                        end if
                    end if
                end if
            end if
        end do
        if(nvals > 0 .and. do_append) then
            if (nvals == 1) then
                ok = opt%set_value(arg_vals(1))
            else
                ok = opt%set_value(arg_vals(:nvals))
            end if
            if (.not.ok) then
                call opt%error%info(err_cause)
                if (opt%error%has_type('arg')) then
                    write(errmsg,  '("error while parsing ''",a,"''")') &
                        trim(opt%label)
                    call argsDB%error%raise_error('data', &
                        'unknown', errmsg, err_cause)
                else
                    write(errmsg, &
                          '("unexpected error met with ''",a,"''")') &
                        trim(opt%label)
                    call argsDB%error%raise_error('data', ' ', errmsg, &
                        err_cause)
                end if
                return
            end if
        end if
        end associate
    end do

    ! Finalization: check unprocessed arguments
    nargs_pos = 0
    ! - Positional arguments or unknown keywords among optional arguments
    do iarg = 1, nargs_opt
        if (todo(iarg)) then
            arg = get_arg(iarg, ierr)
            if (findstr(argsDB%prefixes, arg(:1)) > 0) then
                if (findstr(argsDB%prefix_long, arg(:2)) > 0) then
                    key = arg
                else
                    key = arg(:2)
                end if
                write(errmsg, '("unknown argument: ",a)') trim(key)
                call argsDB%error%raise_error('data', 'unknown', errmsg)
                return
            else
                nargs_pos = nargs_pos + 1
                arg_vals(nargs_pos) = arg
            end if
        end if
    end do
    ! - Positional arguments after delimiter ('--')
    do iarg = nargs_opt+2, nargs
        nargs_pos = nargs_pos + 1
        arg_vals(nargs_pos) = arg
    end do
    ! - Check positional arguments
    if (nargs_pos > 0) then
        if (argsDB%nargs_pos == 0) then
            errmsg = ' '
            if (nargs_pos == 1) then
                write(errmsg, '("unexpected positional argument: ",a)') &
                    trim(arg_vals(1))
            else
                write(errmsg, '("unexpected positional arguments: ",a)') &
                    trim(arg_vals(1))
                do iarg = 2, nargs_pos
                    i = len_trim(errmsg) + 1
                    write(errmsg(i:), '(", ",a)') trim(arg_vals(iarg))
                end do
            end if
            call argsDB%error%raise_error('data', 'unknown', errmsg)
            return
        else if (nargs_pos < abs(argsDB%iargs_pos(3,argsDB%nargs_pos))) then
            write(errmsg, '("not enough positional arguments: at least ",i0,&
                &" expected, ",i0," given")') &
                abs(argsDB%iargs_pos(3,argsDB%nargs_pos)), nargs_pos
            call argsDB%error%raise_error('data', 'missing', errmsg)
            return
        else
            found = .false. ! Found variable length of arguments
            do i = 1, argsDB%nargs_pos
                if (argsDB%iargs_pos(3,i) < 0) then
                    found = .true.
                    exit
                end if
            end do
            if (.not.found) then
                if (nargs_pos /= argsDB%iargs_pos(3,argsDB%nargs_pos)) then
                    write(errmsg, '("mismatch in number of positional &
                        &arguments: at least ",i0," expected, ",i0," &
                        &given")') argsDB%iargs_pos(3,argsDB%nargs_pos), &
                        nargs_pos
                    call argsDB%error%raise_error('data', 'inconsistency', &
                        errmsg)
                    return
                end if
            end if
        end if
        ioff = 0
        do iarg_pos = 1, argsDB%nargs_pos
            iargDB = argsDB%iargs_pos(1,iarg_pos)
            i = ioff + argsDB%iargs_pos(2,iarg_pos)
            j = ioff + argsDB%iargs_pos(3,iarg_pos)
            if (j < 0) then
                ioff = nargs_pos - abs(argsDB%iargs_pos(3,argsDB%nargs_pos))
                ! Number of variable elements: ioff+1
                if (ioff > argsDB%args(iargDB)%arg%max_num - &
                    argsDB%args(iargDB)%arg%min_num) then
                    write(errmsg, '("too many arguments for option: ",a,"; &
                        &up to ",i0," can be given")') &
                        trim(argsDB%args(iargDB)%arg%label), &
                        argsDB%args(iargDB)%arg%max_num
                    call argsDB%error%raise_error('data', 'excess', errmsg)
                    return
                end if
                j = ioff + abs(j)
            end if
            if (j > i) then
                ok = argsDB%args(iargDB)%arg%set_value(arg_vals(i:j))
            else
                ok = argsDB%args(iargDB)%arg%set_value(arg_vals(i))
            end if
            if (.not.ok) then
                call argsDB%args(iargDB)%arg%error%info(err_cause)
                if (argsDB%args(iargDB)%arg%error%has_type('arg')) then
                    write(errmsg,  '("error while parsing ''",a,"''")') &
                        trim(argsDB%args(iargDB)%arg%label)
                    call argsDB%error%raise_error('data', &
                        'unknown', errmsg, err_cause)
                else
                    write(errmsg, &
                          '("unexpected error met with ''",a,"''")') &
                        trim(argsDB%args(iargDB)%arg%label)
                    call argsDB%error%raise_error('data', ' ', errmsg, &
                        err_cause)
                end if
                return
            end if
        end do
    end if

    ! Check if we are missing arguments
    do iargDB = 1, argsDB%nargs
        associate(opt => argsDB%args(iargDB)%arg)
        if (opt%is_req .and. opt%is_set /=2) then
            write(errmsg, &
                  '("option ''",a,"'' is required and missing.")') &
                trim(opt%label)
            call argsDB%error%raise_error('data', 'missing', errmsg)
        end if
        end associate
    end do

contains
    function get_arg(index_arg, index_err) result(res)
        !! Get argument from argument list.
        integer, intent(in) :: index_arg
            !! Argument index.
        integer, intent(out) :: index_err
            !! Error status on getting argument.
            !! 0: no error, -1: not found.
        character(len=MAX_ARGLEN) :: res

        index_err = 0
        if (present(arglist)) then
            if (index_arg <= size(arglist)) then
                res = arglist(index_arg)
            else
                index_err = -1
            end if
        else
            call get_command_argument(index_arg, res, status=index_err)
            if (index_err > 0) index_err = -1
        end if
    end function get_arg

    function chk_shortarg(name) result(index_argDB)
        !! Check if short arg. `name` corresponds to a defined option.
        character(len=2), intent(in) :: name
            !! Short argument name.
        integer :: index_argDB
            !! Index of argument corresponding to `shortarg`, 0 otherwise.

        logical :: is_found

        is_found = .false.
        do index_argDB = 1, argsDB%nargs
            if (argsDB%args(index_argDB)%arg%short_name == name) then
                is_found = .true.
                if (index_argDB == argsDB%iarg_help) call argsDB%print_help()
                exit
            end if
        end do
        if (.not.is_found) then
            write(errmsg, '("unknown option as argument num. ",i0,": ",a)') &
                iarg
            call argsDB%error%raise_error('data', 'unknown', errmsg)
            index_argDB = 0
        end if
    end function chk_shortarg

end procedure argsdb_parse_args

! ======================================================================

end submodule parse_cmdline_parser
