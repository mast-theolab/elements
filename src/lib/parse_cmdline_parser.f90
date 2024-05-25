submodule (parse_cmdline) parse_cmdline_parser

contains

! ======================================================================

module procedure parse_args_list
    !! Parse a list of arguments in input or from commandline
    !!
    !! Parses a list of arguments from the commandline or alternatively
    !!   from a list given in argument to the subroutine.
    implicit none

    integer :: i, iarg, iarg_pos, iargDB, ierr, ioff, ios, isep, j, nargs, &
        nargs_opt, nargs_pos, nvals
    real :: rval
    logical :: do_append, found, ignore_next
    logical, dimension(:), allocatable :: todo
    character(len=MAX_ARGLEN) :: arg, key, val
    character(len=2) :: shortarg
    character(len=256) :: errmsg
    character(len=MAX_ARGLEN), dimension(:), allocatable :: arg_vals
    class(BaseException), allocatable :: suberr

    err = InitError()

    ! Find number of arguments and initialize arrays
    if (present(arglist)) then
        nargs = size(arglist)
    else
        nargs = command_argument_count()
    end if
    allocate(todo(nargs), arg_vals(nargs))
    todo = .True. ; arg_vals = ' '

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
            call RaiseError(err, errmsg)
            return
        end if
        if (this%iarg_help > 0) then
            if (allocated(this%args(this%iarg_help)%arg%short_name)) then
                if (arg == this%args(this%iarg_help)%arg%short_name) &
                    call this%print_help()
            end if
            if (allocated(this%args(this%iarg_help)%arg%long_name)) then
                if (arg == this%args(this%iarg_help)%arg%long_name) &
                    call this%print_help()
            end if
        end if
        if (len_trim(arg) > 2 &
            .and. findstr(this%prefixes, arg(1:1)) > 0 &
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
            if (this%args(iargDB)%arg%max_num == 0) then
                suberr = this%args(iargDB)%arg%set_value()
                if (suberr%raised()) then
                    select type(suberr)
                        class is (Error)
                            write(errmsg, &
                                  '("Error while parsing ''",a,"'': ",a)') &
                                shortarg, trim(suberr%msg())
                            call RaiseError(err, errmsg)
                            return
                        class is (ValueError)
                            call RaiseValueError(err, suberr%msg())
                            return
                        class default
                            write(errmsg, '("Unknown error while parsing &
                                &'',a,'': ",a)') shortarg, trim(suberr%msg())
                            call RaiseError(err, errmsg)
                            return
                    end select
                end if
                do i = 3, len_trim(arg)
                    shortarg = arg(1:1) // arg(i:i)
                    iargDB = chk_shortarg(shortarg)
                    if (iargDB == 0) return
                    if (this%args(iargDB)%arg%max_num == 0) then
                        suberr = this%args(iargDB)%arg%set_value()
                        if (suberr%raised()) then
                            select type(suberr)
                                class is (Error)
                                    write(errmsg, &
                                        '("Error while parsing ''",a,"'':&
                                        & ",a)') &
                                        shortarg, trim(suberr%msg())
                                    call RaiseError(err, errmsg)
                                    return
                                class is (ValueError)
                                    call RaiseValueError(err, suberr%msg())
                                    return
                                class default
                                    write(errmsg, &
                                        '("Unknown error while parsing &
                                        &'',a,'': ", a)') shortarg, &
                                        trim(suberr%msg())
                                    call RaiseError(err, errmsg)
                                    return
                            end select
                        end if
                    else
                        write(errmsg, '("Argument num. ",i0," mixes options &
                            &that require values and others not.",a,&
                            &"They cannot be merged.")') iarg, new_line(' ')
                        call RaiseError(err, errmsg)
                        return
                    end if
                end do
                todo(iarg) = .False.
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
    do iargDB = 1, this%nargs
        associate(opt => this%args(iargDB)%arg)
        if(opt%is_pos) cycle
        select type (opt)
            class is (ArgIntList)
                do_append = .True.
            class is (ArgRealList)
                do_append = .True.
            class is (ArgCharList)
                do_append = .True.
            class default
                do_append = .False.
        end select
        nvals = 0
        do iarg = 1, nargs_opt
            if (todo(iarg)) then
                arg = get_arg(iarg, ierr)
                key = ' '
                val = ' '
                found = .False.
                if (findstr(this%prefix_long, arg(:2)) > 0) then
                    isep = index(arg, '=')
                    if (isep > 0) then
                        key = arg(:isep-1)
                    else
                        key = arg
                    end if
                    found = opt%long_name == key
                else if (findstr(this%prefix_short, arg(:1)) > 0) then
                    if (len_trim(arg) == 1) then
                        write(errmsg, '("Wrong argument in position ",i0)') &
                            iarg
                        call RaiseError(err, errmsg)
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
                    todo(iarg) = .False.
                    if (isep > 0) then
                        if (opt%max_num == 0) then
                            write(errmsg, '("Argument ''",a,"'' does not &
                                &accept values")') trim(key)
                            call RaiseError(err, errmsg)
                            return
                        end if
                        val = arg(isep+1:)
                    else if(opt%max_num > 0) then
                        val = get_arg(iarg+1, ierr)
                        todo(iarg+1) = .False.
                    else
                        val = ' '
                    end if
                    if (do_append) then
                        nvals = nvals + 1
                        arg_vals(nvals) = val
                    else
                        suberr = opt%set_value(val)
                        if (suberr%raised()) then
                            select type(suberr)
                                class is (Error)
                                    write(errmsg, &
                                        '("Error while parsing ''",a,"'':&
                                        & ",a)') trim(key), trim(suberr%msg())
                                    call RaiseError(err, errmsg)
                                    return
                                class is (ValueError)
                                    call RaiseValueError(err, suberr%msg())
                                    return
                                class default
                                    write(errmsg, &
                                        '("Unknown error while parsing &
                                        &''",a,"'': ",a)') trim(key), &
                                        trim(suberr%msg())
                                    call RaiseError(err, errmsg)
                                    return
                            end select
                        end if
                    end if
                end if
            end if
        end do
        if(nvals > 0 .and. do_append) then
            if (nvals == 1) then
                suberr = opt%set_value(arg_vals(1))
            else
                suberr = opt%set_value(arg_vals(:nvals))
            end if
            if (suberr%raised()) then
                select type(suberr)
                    class is (Error)
                        write(errmsg, &
                            '("Error while parsing values for ''",a,"''")') &
                                trim(opt%label)
                        call RaiseError(err, errmsg)
                        return
                    class is (ValueError)
                        call RaiseValueError(err, suberr%msg())
                        return
                    class default
                        write(errmsg, &
                            '("Unknown error while parsing values for &
                            &''",a,"''")') trim(opt%label)
                        call RaiseError(err, errmsg)
                        return
                end select
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
            if (findstr(this%prefixes, arg(:1)) > 0) then
                if (findstr(this%prefix_long, arg(:2)) > 0) then
                    key = arg
                else
                    key = arg(:2)
                end if
                write(errmsg, '("Unknown argument: ",a)') trim(key)
                call RaiseError(err, errmsg)
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
        if (this%nargs_pos == 0) then
            errmsg = ' '
            if (nargs_pos == 1) then
                write(errmsg, '("Unexpected positional argument: ",a)') &
                    trim(arg_vals(1))
            else
                write(errmsg, '("Unexpected positional arguments: ",a)') &
                    trim(arg_vals(1))
                do iarg = 2, nargs_pos
                    i = len_trim(errmsg) + 1
                    write(errmsg(i:), '(", ",a)') trim(arg_vals(iarg))
                end do
            end if
            call RaiseError(err, errmsg)
            return
        else if (nargs_pos < abs(this%iargs_pos(3,this%nargs_pos))) then
            write(errmsg, '("Not enough positional arguments: at least ",i0,&
                &" expected, ",i0," given.")') &
                abs(this%iargs_pos(3,this%nargs_pos)), nargs_pos
            call RaiseError(err, errmsg)
            return
        else
            found = .False. ! Found variable length of arguments
            do i = 1, this%nargs_pos
                if (this%iargs_pos(3,i) < 0) then
                    found = .True.
                    exit
                end if
            end do
            if (.not.found) then
                if (nargs_pos /= this%iargs_pos(3,this%nargs_pos)) then
                    write(errmsg, '("Mismatch in number of positional &
                        &arguments: at least ",i0," expected, ",i0," &
                        &given.")') this%iargs_pos(3,this%nargs_pos), nargs_pos
                    call RaiseError(err, errmsg)
                    return
                end if
            end if
        end if
        ioff = 0
        do iarg_pos = 1, this%nargs_pos
            iargDB = this%iargs_pos(1,iarg_pos)
            i = ioff + this%iargs_pos(2,iarg_pos)
            j = ioff + this%iargs_pos(3,iarg_pos)
            if (j < 0) then
                ioff = nargs_pos - abs(this%iargs_pos(3,this%nargs_pos))
                ! Number of variable elements: ioff+1
                if (ioff > this%args(iargDB)%arg%max_num - &
                    this%args(iargDB)%arg%min_num) then
                    write(errmsg, '("Too many arguments for option: ",a,". &
                        &Up to ",i0," can be given.")') &
                        trim(this%args(iargDB)%arg%label), &
                        this%args(iargDB)%arg%max_num
                    Call RaiseError(err, errmsg)
                    return
                end if
                j = ioff + abs(j)
            end if
            if (j > i) then
                suberr = this%args(iargDB)%arg%set_value(arg_vals(i:j))
            else
                suberr = this%args(iargDB)%arg%set_value(arg_vals(i))
            end if
            if (suberr%raised()) then
                select type(suberr)
                    class is (Error)
                        write(errmsg, &
                            '("Error while parsing values for ''",a,"''")') &
                                this%args(iargDB)%arg%label
                        call RaiseError(err, errmsg)
                        return
                    class is (ValueError)
                        call RaiseValueError(err, suberr%msg())
                        return
                    class default
                        write(errmsg, &
                            '("Unknown error while parsing &
                            &''",a,"''")') this%args(iargDB)%arg%label
                        call RaiseError(err, errmsg)
                        return
                end select
            end if
        end do
    end if

    ! Check if we are missing arguments
    do iargDB = 1, this%nargs
        associate(opt => this%args(iargDB)%arg)
        if (opt%is_req .and. opt%is_set /=2) then
            write(errmsg, &
                  '("Error: Option ''",a,"'' is required and missing.")') &
                trim(opt%label)
            call RaiseError(err, errmsg)
        end if
        end associate
    end do

contains
    function get_arg(iarg, ierr) result(res)
        !! Get argument from argument list.
        integer, intent(in) :: iarg
        !! Argument index.
        integer, intent(out) :: ierr
        !! Error instance.
        character(len=MAX_ARGLEN) :: res

        ierr = 0
        if (present(arglist)) then
            if (iarg <= size(arglist)) then
                res = arglist(iarg)
            else
                ierr = -1
            end if
        else
            call get_command_argument(iarg, res, status=ierr)
            if (ierr > 0) ierr = -1
        end if
    end function get_arg

    function chk_shortarg(shortarg) result(iargDB)
        !! Check if shortarg corresponds to a defined option.
        character(len=2), intent(in) :: shortarg
        !! Short argument name.
        integer :: iargDB
        !! Index of argument corresponding to `shortarg`, 0 otherwise.

        logical :: found

        found = .False.
        do iargDB = 1, this%nargs
            if (this%args(iargDB)%arg%short_name == shortarg) then
                found = .True.
                if (iargDB == this%iarg_help) call this%print_help()
                exit
            end if
        end do
        if (.not.found) then
            write(errmsg, '("Unknown option as argument num. ",i0,": ",a)') &
                iarg
            call RaiseError(err, errmsg)
            iargDB = 0
        end if
    end function chk_shortarg

end procedure parse_args_list

! ======================================================================

end submodule parse_cmdline_parser