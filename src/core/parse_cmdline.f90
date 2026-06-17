module parse_cmdline

    use iso_fortran_env, only: int32, int64, real32, real64
    use run_env, only: CoreExecObject, ErrorHandle, run
    use string, only: findstr, locase, upcase

    implicit none

    integer, parameter, private :: &
        ! High number of arguments (e.g. 400) risks raising a warning message
        ! from GFortran about a potential problem of stack allocation, 200
        ! should be enough already.
        MAX_ARGS = 200, &
        MAX_ARGLEN = 256

    ! character(len=*), dimension(*), parameter :: &
    !     prefixes = ['-'], &
    !     prefixes_short = [(prefixes(i), i=1,size(prefixes))], &
    !     prefixes_long = [(prefixes(i)//prefixes(i), i=1,size(prefixes))]

    type, extends(CoreExecObject), private, abstract :: arg_obj
        !! Basic argument object
        character(len=:), allocatable :: short_name
        character(len=:), allocatable :: long_name
        character(len=:), allocatable :: help_msg
        character(len=:), private, allocatable :: label
        integer :: min_num = 1
        integer :: max_num = 1
        integer :: is_set = 0
        !! is_set: value has been set: 0: no, 1: default, 2: user
        logical :: is_req  ! Argument is not optional
        logical :: is_pos  ! Argument is positional
    contains
        procedure(set_arg_value1), deferred, private :: set_arg_scalar
        procedure(set_arg_valueN), deferred, private :: set_arg_array
        generic :: set_value => set_arg_scalar, set_arg_array
        procedure :: set_arg => argobj_set_names
        procedure :: set_help => argobj_set_helpmsg
    end type arg_obj

    abstract interface
        function set_arg_value1(arg, argval, sep) result(done)
            import arg_obj, ErrorHandle
            class(arg_obj), intent(inout) :: arg
            character(len=*), intent(in), optional :: argval
            character(len=*), intent(in), optional :: sep
            logical :: done
        end function set_arg_value1
    end interface

    abstract interface
        function set_arg_valueN(arg, argvals) result(done)
            import arg_obj, ErrorHandle
            class(arg_obj), intent(inout) :: arg
            character(len=*), dimension(:), intent(in) :: argvals
            logical :: done
        end function set_arg_valueN
    end interface

    type, private, extends(arg_obj) :: arg_int
        integer(int64) :: value
        ! For some reason, GFortran assumes HUGE(value) is real.
        ! To bypass the problem, we pass a constant value.
        integer(int64) :: min_ok = -huge(1_int64)
        integer(int64) :: max_ok = +huge(1_int64)
        integer(int64) :: const
        logical :: add_value = .false.
    contains
        procedure :: set_arg_scalar => argint_from_str
        procedure :: set_arg_array => argint_from_list_no
    end type arg_int

    type, private, extends(arg_obj) :: arg_ints
        integer(int64), dimension(:), allocatable :: values
        integer(int64) :: min_ok = -huge(1_int64)
        integer(int64) :: max_ok = +huge(1_int64)
        logical :: append = .false.
        !! if append is True, the argument can be called multiple time.
        !!   each time, the new value is appended.  In this case,
        !!   character-separated values cannot be used.
    contains
        procedure :: set_arg_scalar => argints_from_str
        procedure :: set_arg_array => argints_from_list
    end type arg_ints

    type, private, extends(arg_obj) :: arg_real
        real(real64) :: value
        real(real64) :: min_ok = -huge(1.0_real64)
        real(real64) :: max_ok = +huge(1.0_real64)
        real(real64) :: const
        logical :: add_value = .false.
    contains
        procedure :: set_arg_scalar => argreal_from_str
        procedure :: set_arg_array => argreal_from_list_no
    end type arg_real

    type, private, extends(arg_obj) :: arg_reals
        real(real64), dimension(:), allocatable :: values
        real(real64) :: min_ok = -huge(1.0_real64)
        real(real64) :: max_ok = +huge(1.0_real64)
        logical :: append = .false.
        !! if append is True, the argument can be called multiple time.
        !!   each time, the new value is appended.  In this case,
        !!   character-separated values cannot be used.
    contains
        procedure :: set_arg_scalar => argreals_from_str
        procedure :: set_arg_array => argreals_from_list
    end type arg_reals

    type, private, extends(arg_obj) :: arg_bool
        logical :: value
        logical :: const
    contains
        procedure :: set_arg_scalar => argbool_from_str
        procedure :: set_arg_array => argbool_from_list_no
    end type arg_bool

    type, private, extends(arg_obj) :: arg_char
        character(len=:), allocatable :: value
    contains
        procedure :: set_arg_scalar => argchar_from_str
        procedure :: set_arg_array => argchar_from_list_no
    end type arg_char

    type, private, extends(arg_obj) :: arg_chars
        character(len=:), dimension(:), allocatable :: values
        logical :: append = .false.
        !! if append is True, the argument can be called multiple time.
        !!   each time, the new value is appended.  In this case,
        !!   character-separated values cannot be used.
    contains
        procedure :: set_arg_scalar => argchars_from_str
        procedure :: set_arg_array => argchars_from_list
    end type arg_chars

    type, extends(arg_obj), private :: arg_gen
        !! A dummy container to build list of different arguments
        class(arg_obj), allocatable :: arg
    contains
        procedure :: set_arg_scalar => arggen_from_str_no
        procedure :: set_arg_array => arggen_from_list_no
    end type arg_gen

    type, public, extends(CoreExecObject) :: CmdLineArgsDB
        private
        type(arg_gen), dimension(MAX_ARGS) :: args
        integer, dimension(3,MAX_ARGS) :: iargs_pos
        character(len=1), dimension(:), allocatable :: prefixes
        character(len=1), dimension(:), allocatable :: prefix_short
        character(len=2), dimension(:), allocatable :: prefix_long
        character(len=:), allocatable :: progname
        integer :: iarg_help = 0
        !! Stores indexes of positional arguments. <0 for arbitrary number
        integer :: nargs = 0, nargs_pos = 0
    contains
        procedure, private :: argsdb_getval_int32, argsdb_getval_int64, &
            argsdb_getvals_int32, argsdb_getvals_int64, &
            argsdb_getval_real32, argsdb_getval_real64, &
            argsdb_getvals_real32, argsdb_getvals_real64, &
            argsdb_getval_bool, &
            argsdb_getval_char, argsdb_getvals_char
        procedure, private :: get_argname_id
        procedure :: chk_name_overlap => chk_argname_overlap
        procedure :: add_arg_int => argsdb_add_int
        procedure :: add_arg_real => argsdb_add_real
        procedure :: add_arg_bool => argsdb_add_bool
        procedure :: add_arg_char => argsdb_add_char
        procedure :: is_user_set => argsdb_val_is_userset
        generic :: get_value => argsdb_getval_int32, argsdb_getval_int64, &
            argsdb_getvals_int32, argsdb_getvals_int64, &
            argsdb_getval_real32, argsdb_getval_real64, &
            argsdb_getvals_real32, argsdb_getvals_real64, &
            argsdb_getval_bool, &
            argsdb_getval_char, argsdb_getvals_char
        procedure :: parse_args => argsdb_parse_args
        procedure :: print_help => print_help
    end type CmdLineArgsDB

    interface CmdLineArgsDB
        module procedure init_args_db
    end interface CmdLineArgsDB

interface

! ----------------------------------------------------------------------
! TYPE-BOUND PROCEDURES (INTERFACES) - CMDLINEARGSDB
! ----------------------------------------------------------------------

module subroutine argsdb_add_bool(argsDB, argtype, shortname, longname, &
                                  label, help, required, def_value)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argtype
        !! Type of argument.
    character(len=*), intent(in), optional :: shortname
        !! Short argument name (as `-o`).
    character(len=*), intent(in), optional :: longname
        !! Long argument name (as `--option`).
    character(len=*), intent(in), optional :: label
        !! Internal argument name, which may be used for help.
    character(len=*), intent(in), optional :: help
        !! Help message to be displayed in the help
    logical, intent(in), optional :: required
        !! Argument must be given. By default, based on presence of prefix(es).
    logical, intent(in), optional :: def_value
        !! Default value, to use if not set by user.
end subroutine argsdb_add_bool

! ----------------------------------------------------------------------

module subroutine argsdb_add_char(argsDB, argtype, shortname, longname, &
                                  label, help, required, def_value, &
                                  min_nvals, max_nvals)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argtype
        !! Type of argument.
    character(len=*), intent(in), optional :: shortname
        !! Short argument name (as `-o`).
    character(len=*), intent(in), optional :: longname
        !! Long argument name (as `--option`).
    character(len=*), intent(in), optional :: label
        !! Internal argument name, which may be used for help.
    character(len=*), intent(in), optional :: help
        !! Help message to be displayed in the help
    logical, intent(in), optional :: required
        !! Argument must be given; by default, based on presence of prefix(es).
    character(len=*), intent(in), optional :: def_value
        !! Default value, to use if not set by user.
    class(*), intent(in), optional :: min_nvals
        !! For lists of values, minimum size of the list
    class(*), intent(in), optional :: max_nvals
        !! For lists of values, maixmum size of the list.
end subroutine argsdb_add_char

! ----------------------------------------------------------------------

module subroutine argsdb_add_int(argsDB, argtype, shortname, longname, label, &
                                 help, required, def_value, min_value, &
                                 max_value, const_value, add_value, &
                                 min_nvals, max_nvals)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argtype
        !! Type of argument.
    character(len=*), intent(in), optional :: shortname
        !! Short argument name (as `-o`).
    character(len=*), intent(in), optional :: longname
        !! Long argument name (as `--option`).
    character(len=*), intent(in), optional :: label
        !! Internal argument name, which may be used for help.
    character(len=*), intent(in), optional :: help
        !! Help message to be displayed in the help
    logical, intent(in), optional :: required
        !! Argument must be given; by default, based on presence of prefix(es).
    class(*), intent(in), optional :: def_value
        !! Default value, to use if not set by user.
    class(*), intent(in), optional :: min_value
        !! Minimum accepted value
    class(*), intent(in), optional :: max_value
        !! Maximum accepted value
    class(*), intent(in), optional :: const_value
        !! Constant value, to use when argument found (no value expected).
    logical, intent(in), optional :: add_value
        !! Increment value associated to argument each time it is encountered.
    class(*), intent(in), optional :: min_nvals
        !! For lists of values, minimum size of the list.
    class(*), intent(in), optional :: max_nvals
        !! For lists of values, maixmum size of the list.
end subroutine argsdb_add_int

! ----------------------------------------------------------------------

module subroutine argsdb_add_real(argsDB, argtype, shortname, longname, &
                                  label, help, required, def_value, &
                                  min_value, max_value, const_value, &
                                  add_value, min_nvals, max_nvals)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object
    character(len=*), intent(in) :: argtype
        !! Type of argument.
    character(len=*), intent(in), optional :: shortname
        !! Short argument name (as `-o`).
    character(len=*), intent(in), optional :: longname
        !! Long argument name (as `--option`).
    character(len=*), intent(in), optional :: label
        !! Internal argument name, which may be used for help.
    character(len=*), intent(in), optional :: help
        !! Help message to be displayed in the help
    logical, intent(in), optional :: required
        !! Argument must be given; by default, based on presence of prefix(es).
    class(*), intent(in), optional :: def_value
        !! Default value, to use if not set by user.
    class(*), intent(in), optional :: min_value
        !! Minimum accepted value
    class(*), intent(in), optional :: max_value
        !! Maximum accepted value
    class(*), intent(in), optional :: const_value
        !! Constant value, to use when argument found (no value expected).
    logical, intent(in), optional :: add_value
        !! Increment value associated to argument each time it is encountered.
    class(*), intent(in), optional :: min_nvals
        !! For lists of values, minimum size of the list
    class(*), intent(in), optional :: max_nvals
        !! For lists of values, maixmum size of the list.
end subroutine argsdb_add_real

! ----------------------------------------------------------------------

module subroutine argsdb_getval_bool(argsDB, argname, result)
    class(CmdLineArgsDB), intent(inout) :: argsDB
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    logical, intent(out) :: result
    !! Associated value.
end subroutine argsdb_getval_bool

! ----------------------------------------------------------------------

module subroutine argsdb_getval_char(argsDB, argname, result)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argname
        !! Argumemt name/label.
    character(len=*), intent(out) :: result
        !! Associated value.
end subroutine argsdb_getval_char

! ----------------------------------------------------------------------

module subroutine argsdb_getval_int32(argsDB, argname, result)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argname
        !! Argumemt name/label.
    integer(int32), intent(out) :: result
        !! Associated value.
end subroutine argsdb_getval_int32

! ----------------------------------------------------------------------

module subroutine argsdb_getval_int64(argsDB, argname, result)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argname
        !! Argumemt name/label.
    integer(int64), intent(out) :: result
        !! Associated value.
end subroutine argsdb_getval_int64

! ----------------------------------------------------------------------

module subroutine argsdb_getval_real32(argsDB, argname, result)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argname
        !! Argumemt name/label.
    real(real32), intent(out) :: result
        !! Associated value.
end subroutine argsdb_getval_real32

! ----------------------------------------------------------------------

module subroutine argsdb_getval_real64(argsDB, argname, result)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argname
        !! Argumemt name/label.
    real(real64), intent(out) :: result
        !! Associated value.
end subroutine argsdb_getval_real64

! ----------------------------------------------------------------------

module subroutine argsdb_getvals_char(argsDB, argname, result)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argname
        !! Argumemt name/label.
    character(len=:), dimension(:), allocatable, intent(out) :: result
        !! Associated value.
end subroutine argsdb_getvals_char

! ----------------------------------------------------------------------

module subroutine argsdb_getvals_int32(argsDB, argname, result)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argname
        !! Argumemt name/label.
    integer(int32), dimension(:), allocatable, intent(out) :: result
        !! Associated value.
end subroutine argsdb_getvals_int32

! ----------------------------------------------------------------------

module subroutine argsdb_getvals_int64(argsDB, argname, result)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argname
        !! Argumemt name/label.
    integer(int64), dimension(:), allocatable, intent(out) :: result
        !! Associated value.
end subroutine argsdb_getvals_int64

! ----------------------------------------------------------------------

module subroutine argsdb_getvals_real32(argsDB, argname, result)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argname
        !! Argumemt name/label.
    real(real32), dimension(:), allocatable, intent(out) :: result
        !! Associated value.
end subroutine argsdb_getvals_real32

! ----------------------------------------------------------------------

module subroutine argsdb_getvals_real64(argsDB, argname, result)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argname
        !! Argumemt name/label.
    real(real64), dimension(:), allocatable, intent(out) :: result
        !! Associated value.
end subroutine argsdb_getvals_real64

! ----------------------------------------------------------------------

module subroutine argsdb_parse_args(argsDB, arglist)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), dimension(:), intent(in), optional :: arglist
        !! List of arguments to parse.
end subroutine argsdb_parse_args

! ----------------------------------------------------------------------

module function argsdb_val_is_userset(argsDB, argname) result(res)
    class(CmdLineArgsDB), intent(inout) :: argsDB
        !! Arguments database object.
    character(len=*), intent(in) :: argname
        !! Argumemt name/label.
    logical :: res
        !! Boolean stating if value set by user.
end function argsdb_val_is_userset


! ----------------------------------------------------------------------
! TYPE-BOUND PROCEDURES (INTERFACES) - ARGUMENTS OBJECTS
! ----------------------------------------------------------------------

module function argbool_from_list_no(arg, argvals) result(done)
    class(arg_bool), intent(inout) :: arg
        !! Argument object.
    character(len=*), dimension(:), intent(in) :: argvals
        !! List of values stored as strings.
    logical :: done
        !! Success status of the assignment operation.
end function argbool_from_list_no

! ----------------------------------------------------------------------

module function argbool_from_str(arg, argval, sep) result(done)
    class(arg_bool), intent(inout) :: arg
        !! Argument object.
    character(len=*), intent(in), optional :: argval
        !! Value stored as string.
    character(len=*), intent(in), optional :: sep
        !! Dummy argument.
    logical :: done
        !! Success status of the assignment operation.
end function argbool_from_str

! ----------------------------------------------------------------------

module function argchar_from_list_no(arg, argvals) result(done)
    class(arg_char), intent(inout) :: arg
        !! Argument object.
    character(len=*), dimension(:), intent(in) :: argvals
        !! List of values stored as strings.
    logical :: done
        !! Success status of the assignment operation.
end function argchar_from_list_no

! ----------------------------------------------------------------------

module function argchar_from_str(arg, argval, sep) result(done)
    class(arg_char), intent(inout) :: arg
        !! Argument object.
    character(len=*), intent(in), optional :: argval
        !! Value stored as string.
    character(len=*), intent(in), optional :: sep
        !! Dummy argument.
    logical :: done
        !! Success status of the assignment operation.
end function argchar_from_str

! ----------------------------------------------------------------------

module function argchars_from_list(arg, argvals) result(done)
    class(arg_chars), intent(inout) :: arg
        !! Argument object.
    character(len=*), dimension(:), intent(in) :: argvals
        !! List of values stored as strings.
    logical :: done
        !! Success status of the assignment operation.
end function argchars_from_list

! ----------------------------------------------------------------------

module function argchars_from_str(arg, argval, sep) result(done)
    class(arg_chars), intent(inout) :: arg
        !! Argument object.
    character(len=*), intent(in), optional :: argval
        !! String containing the values.
    character(len=*), intent(in), optional :: sep
        !! List of one-character separators, given as a string.
    logical :: done
        !! Success status of the assignment operation.
end function argchars_from_str

! ----------------------------------------------------------------------

module function arggen_from_list_no(arg, argvals) result(done)
    class(arg_gen), intent(inout) :: arg
        !! Argument object.
    character(len=*), dimension(:), intent(in) :: argvals
        !! Dummy argument.
    logical :: done
        !! Success status of the assignment operation.
end function arggen_from_list_no

! ----------------------------------------------------------------------

module function arggen_from_str_no(arg, argval, sep) result(done)
    class(arg_gen), intent(inout) :: arg
        !! Argument object.
    character(len=*), intent(in), optional :: argval
        !! Dummy argument.
    character(len=*), intent(in), optional :: sep
        !! Dummy argument.
    logical :: done
        !! Success status of the assignment operation.
end function arggen_from_str_no

! ----------------------------------------------------------------------

module function argint_from_list_no(arg, argvals) result(done)
    class(arg_int), intent(inout) :: arg
        !! Argument object.
    character(len=*), dimension(:), intent(in) :: argvals
        !! List of values stored as strings.
    logical :: done
        !! Success status of the assignment operation.
end function argint_from_list_no

! ----------------------------------------------------------------------

module function argint_from_str(arg, argval, sep) result(done)
    class(arg_int), intent(inout) :: arg
        !! Argument object.
    character(len=*), intent(in), optional :: argval
        !! Value stored as string.
    character(len=*), intent(in), optional :: sep
        !! Dummy argument.
    logical :: done
        !! Success status of the assignment operation.
end function argint_from_str

! ----------------------------------------------------------------------

module function argints_from_list(arg, argvals) result(done)
    class(arg_ints), intent(inout) :: arg
        !! Argument object.
    character(len=*), dimension(:), intent(in) :: argvals
        !! List of values stored as strings.
    logical :: done
        !! Success status of the assignment operation.
end function argints_from_list

! ----------------------------------------------------------------------

module function argints_from_str(arg, argval, sep) result(done)
    class(arg_ints), intent(inout) :: arg
        !! Argument object.
    character(len=*), intent(in), optional :: argval
        !! String containing the values.
    character(len=*), intent(in), optional :: sep
        !! list of one-character separators, given as a string.
    logical :: done
        !! Success status of the assignment operation.
end function argints_from_str

! ----------------------------------------------------------------------

module function argobj_set_helpmsg(arg, msg) result(done)
    class(arg_obj), intent(inout) :: arg
        !! Argument object.
    character(len=*), intent(in) :: msg
        !! Help message.
    logical :: done
        !! Success status of the operation.
end function argobj_set_helpmsg

! ----------------------------------------------------------------------

module function argobj_set_names(arg, prefixes, short, long, label, &
                                 required) result(done)
    class(arg_obj), intent(inout) :: arg
        !! Argument object.
    character(len=*), dimension(:), intent(in) :: prefixes
        !! List of prefix characters.
    character(len=*), intent(in), optional :: short
        !! Short form of the argument.
    character(len=*), intent(in), optional :: long
        !! Long form of the argument.
    character(len=*), intent(in), optional :: label
        !! Internal and reference name of the argument.
    logical, intent(in), optional :: required
        !! Set the argument as required; only for optional argument.
    logical :: done
        !! Success status of the operation.
end function argobj_set_names

! ----------------------------------------------------------------------

module function argreal_from_list_no(arg, argvals) result(done)
    class(arg_real), intent(inout) :: arg
        !! Argument object.
    character(len=*), dimension(:), intent(in) :: argvals
        !! List of values stored as strings.
    logical :: done
        !! Success status of the assignment operation.
end function argreal_from_list_no

! ----------------------------------------------------------------------

module function argreal_from_str(arg, argval, sep) result(done)
    class(arg_real), intent(inout) :: arg
        !! Argument object.
    character(len=*), intent(in), optional :: argval
        !! Value stored as string.
    character(len=*), intent(in), optional :: sep
        !! Dummy argument.
    logical :: done
        !! Success status of the assignment operation.
end function argreal_from_str

! ----------------------------------------------------------------------

module function argreals_from_list(arg, argvals) result(done)
    class(arg_reals), intent(inout) :: arg
        !! Argument object.
    character(len=*), dimension(:), intent(in) :: argvals
        !! List of values stored as strings.
    logical :: done
        !! Success status of the assignment operation.
end function argreals_from_list

! ----------------------------------------------------------------------

module function argreals_from_str(arg, argval, sep) result(done)
    class(arg_reals), intent(inout) :: arg
        !! Argument object.
    character(len=*), intent(in), optional :: argval
        !! String containing the values.
    character(len=*), intent(in), optional :: sep
        !! List of one-character separators, given as a string.
    logical :: done
        !! Success status of the assignment operation.
end function argreals_from_str

! ----------------------------------------------------------------------

end interface

contains

! ======================================================================
! PSEUDO-CONSTRUCTORS
! ======================================================================

function init_args_db(add_help, prefixes, progname) result(argsDB)
    !! Initialize database to store arguments parameters
    !!
    !! Initializes a database object to store information on supported
    !!   options in the commandline and extract/process data.
    !! The system expects short-name optional argument to be preceded by
    !!   1 prefix character, long names by 2 identical characters.
    logical, intent(in), optional :: add_help
        !! Add help keywords: -h/--help.
    character(len=*), intent(in), optional :: prefixes
        !! List of accepted 1-char prefixes for optional arguments, as a string.
    character(len=*), intent(in), optional :: progname
        !! Name of the program.
    type(CmdLineArgsDB) :: argsDB
        !! Database of command-line arguments.

    integer :: i, istat, n_prefix
    logical :: add_help_
    character(len=1) :: prefix

    if (present(prefixes)) then
        n_prefix = len(prefixes)
        allocate(argsDB%prefixes(n_prefix), argsDB%prefix_short(n_prefix), &
                 argsDB%prefix_long(n_prefix), stat=istat)
        if (istat /= 0) then
            call argsDB%error%raise_error('mem', 'allocate', &
                'failed to set up memory for user-defined argument prefixes')
            return
        end if
        do i = 1, n_prefix
            prefix = prefixes(i:i)
            argsDB%prefixes(i) = prefix
            argsDB%prefix_short(i) = prefix
            argsDB%prefix_long(i) = prefix // prefix
        end do
    else
        allocate(argsDB%prefixes(1), argsDB%prefix_short(1), argsDB%prefix_long(1), &
                 stat=istat)
        if (istat /= 0) then
            call argsDB%error%raise_error('mem', 'allocate', &
                'failed to set up default arguments prefixes')
            return
        end if
        argsDB%prefixes(1) = '-'
        argsDB%prefix_short(1) = '-'
        argsDB%prefix_long(1) = '--'
    end if

    if (present(progname)) then
        argsDB%progname = trim(progname)
    else
        call get_command_argument(0, length=i)
        allocate(character(len=i) :: argsDB%progname)
        call get_command_argument(0, argsDB%progname)
    end if

    if (present(add_help)) then
        add_help_ = add_help
    else
        add_help_ = .true.
    end if

    if (add_help_) then
        call argsDB%add_arg_bool('store_true', '-h', '--help', &
                                 help='Print this help message', &
                                 def_value=.false.)
    end if

end function init_args_db

! ======================================================================
! TYPE-BOUND PROCEDURES
! ======================================================================

function chk_argname_overlap(this, new_arg, argname) result(res)
    !! Check if new argument overlaps with registered names in DB.
    !!
    !! Checks if any of the identifiers of an argument (short name,
    !!   long name, label) overlaps with arguments already inserted into
    !!   the arguments DB.
    !! Returns True if an overlap exist.

    class(CmdLineArgsDB), intent(in) :: this
    !! Arguments database object
    class(arg_obj), intent(in) :: new_arg
    !! New argument to be inserted in database
    character(len=*), intent(out) :: argname
    !! Overlapping argument name.
    logical :: res
    !! Result of the check

    integer :: iarg

    argname = ' '
    res = .false.
    if (allocated(new_arg%short_name)) then
        do iarg = 1, this%nargs
            associate (opt => this%args(iarg)%arg)
            if (allocated(opt%short_name)) then
                if (opt%short_name == new_arg%short_name) then
                    res = .true.
                    argname = new_arg%short_name
                    return
                end if
            end if
            end associate
        end do
    end if

    if (allocated(new_arg%long_name)) then
        do iarg = 1, this%nargs
            associate (opt => this%args(iarg)%arg)
            if (allocated(opt%long_name)) then
                if (opt%long_name == new_arg%long_name) then
                    res = .true.
                    argname = new_arg%long_name
                    return
                end if
            end if
            end associate
        end do
    end if

    if (allocated(new_arg%label)) then
        do iarg = 1, this%nargs
            if (this%args(iarg)%arg%label == new_arg%label) then
                res = .true.
                argname = new_arg%label
                return
            end if
        end do
    end if

end function chk_argname_overlap

! ======================================================================

function get_argname_id(this, argname) result(ind)
    !! Get index of argname in database.
    !!
    !! Returns the index of `argname` in database `this`.
    !! Returns 0 if not found.

    class(CmdLineArgsDB), intent(in) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    integer :: ind
    !! Associated value.

    integer :: iarg

    ind = 0
    do iarg = 1, this%nargs
        if (this%args(iarg)%arg%label == argname) then
            ind = iarg
            exit
        end if
    end do

end function get_argname_id

! ======================================================================

subroutine print_help(this)
    !! Print help message
    !!
    !! Builds and prints the help message.
    !! The subroutine then terminates any execution

    class(CmdLineArgsDB), intent(in) :: this
    !! Arguments database object.

    integer :: iargDB, lopt, lline, lline0, N
    character(len=1024) :: line, tmpline

    write(line, '("Usage: ",a)') trim(this%progname)
    lline0 = len_trim(line)
    lline = lline0
    do iargDB = 1, this%nargs
        if (iargDB == this%iarg_help) cycle
        associate(opt => this%args(iargDB)%arg)
        if (.not.opt%is_req) then
            tmpline = ' ['
            lopt = 2
        else
            tmpline = ' '
            lopt = 1
        end if
        if (allocated(opt%short_name)) then
            write(tmpline(lopt+1:), '(a,1x)') opt%short_name
            lopt = lopt + len(opt%short_name)
        else if (allocated(opt%long_name)) then
            write(tmpline(lopt+1:), '(a,1x)') opt%long_name
            lopt = lopt + len(opt%long_name)
        end if
        call build_val_list(opt, tmpline, lopt)
        if (.not.opt%is_req) then
            write(tmpline(lopt+1:), '("]")')
            lopt = lopt + 1
        end if
        end associate
        if (lline + lopt > 80) then
            write(*, '(a)') line(:lline)
            line = ' '
            lline = lline0
        end if
        write(line(lline+1:), '(a)') tmpline(:lopt)
        lline = lline + lopt
    end do
    write(*, '(a)') line(:lline)

    ! Here add program description

    ! Add full list of options
    ! First consider positional arguments
    write(*, '(/,"Positional arguments:")')
    do iargDB = 1, this%nargs
        associate(opt => this%args(iargDB)%arg)
        if (opt%is_pos) then
            line = ' '
            lline = 2
            call build_val_list(opt, line, lline)
            write(*, '(a)') line(:lline)
            if (allocated(opt%help_msg)) then
                write(*, '(24x,a)') trim(opt%help_msg)
            end if
        end if
        end associate
    end do

    ! Then optional arguments
    write(*, '(/,"Optional arguments:")')
    do iargDB = 1, this%nargs
        associate(opt => this%args(iargDB)%arg)
        if (.not.opt%is_pos) then
            line = ' '
            lline = 2
            if (allocated(opt%short_name)) then
                write(line(lline+1:), '(a)') opt%short_name
                lline = lline + len(opt%short_name)
            end if
            if (allocated(opt%long_name)) then
                if (lline > 2) then
                    write(line(lline+1:), '(", ")')
                    lline = lline + 2
                end if
                write(line(lline+1:), '(a)') opt%long_name
                lline = lline + len(opt%long_name)
            end if
            call build_val_list(opt, line, lline)
            write(*, '(a)') line(:lline)
            if (allocated(opt%help_msg)) then
                write(*, '(24x,a)') trim(opt%help_msg)
            end if
        end if
        end associate
    end do

    ! Add footer

    stop

contains

subroutine build_val_list(arg, argline, l_argline)
    !! Build list of values in vals

    interface
        function tostring(text)
            character(len=*), intent(in) :: text
            character(len=len(text)) :: tostring
        end function tostring
    end interface

    class(arg_obj), intent(in) :: arg
    !! Argument object.
    character(len=*), intent(inout) :: argline
    !! Line with description of argument call, already containing argument.
    integer, intent(inout) :: l_argline
    !! Length actually used in `argline`, updated on output.

    integer :: i
    procedure(tostring), pointer :: convert

    if (arg%is_pos) then
        convert => nochange
    else
        convert => upcase
    end if

    do i = 1, arg%min_num
        write(argline(l_argline+1:), '(1x,a)') convert(arg%label)
        l_argline = l_argline + 1 + len(arg%label)
    end do
    if (arg%max_num > arg%min_num) then
        write(argline(l_argline+1:), '(" [")')
        l_argline = l_argline + 2
        N = arg%max_num - arg%min_num
        if (N > 3) N = -3
        do i = 1, abs(N)
            write(argline(l_argline+1:), '(a,1x)') convert(arg%label)
            l_argline = l_argline + len(arg%label) + 1
        end do
        if (N < 0) then
            write(argline(l_argline:), '("...")')
            l_argline = l_argline + 4
        end if
        write(argline(l_argline:), '("]")')
        l_argline = l_argline + 1
    end if

end subroutine build_val_list

function nochange(text)
    ! do nothing (used as target for conversion)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: nochange

    nochange = text
end function

end subroutine print_help

! ======================================================================

end module parse_cmdline
