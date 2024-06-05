module parse_cmdline

    use iso_fortran_env, only: int32, int64, real32, real64
    use exception, only: ArgumentError, BaseException, Error, InitError, &
        RaiseAllocateError, RaiseArgError, RaiseError, RaiseTermination, &
        RaiseValueError, ValueError
    use numeric, only: is_number, to_int64, to_real64
    use string, only: findstr, locase, upcase

    implicit none

    integer, parameter, private :: &
        MAX_ARGS = 400, &
        MAX_ARGLEN = 256

    ! character(len=*), dimension(*), parameter :: &
    !     prefixes = ['-'], &
    !     prefixes_short = [(prefixes(i), i=1,size(prefixes))], &
    !     prefixes_long = [(prefixes(i)//prefixes(i), i=1,size(prefixes))]

    type, private, abstract :: ArgObj
        !! Basic argument object
        character(len=:), allocatable :: short_name
        character(len=:), allocatable :: long_name
        character(len=:), allocatable :: help_msg
        character(len=:), private, allocatable :: label
        integer :: min_num = 1
        integer :: max_num = 1
        integer :: is_set
        !! is_set: value has been set: 0: no, 1: default, 2: user
        logical :: is_req  ! Argument is not optional
        logical :: is_pos  ! Argument is positional
    contains
        procedure(set_arg_value1), deferred, private :: set_arg_scalar
        procedure(set_arg_valueN), deferred, private :: set_arg_array
        generic :: set_value => set_arg_scalar, set_arg_array
        procedure :: set_arg => set_arg_names
        procedure :: set_help => set_arg_helpmsg
    end type ArgObj

    abstract interface
        function set_arg_value1(this, string, sep) result(err)
            import ArgObj, BaseException
            class(ArgObj), intent(inout) :: this
            character(len=*), intent(in), optional :: string
            character(len=*), intent(in), optional :: sep
            class(BaseException), allocatable :: err
        end function set_arg_value1
    end interface

    abstract interface
        function set_arg_valueN(this, strings) result(err)
            import ArgObj, BaseException
            class(ArgObj), intent(inout) :: this
            character(len=*), dimension(:), intent(in) :: strings
            class(BaseException), allocatable :: err
        end function set_arg_valueN
    end interface

    type, private, extends(ArgObj) :: ArgInt
        integer(int64) :: value
        ! For some reason, GFortran assumes HUGE(value) is real.
        ! To bypass the problem, we pass a constant value.
        integer(int64) :: min_ok = -huge(1_int64)
        integer(int64) :: max_ok = +huge(1_int64)
        integer(int64) :: const
        logical :: add_value = .False.
    contains
        procedure :: set_arg_scalar => upd_argint_value
        procedure :: set_arg_array => upd_noarrayI
    end type ArgInt

    type, private, extends(ArgObj) :: ArgIntList
        integer(int64), dimension(:), allocatable :: values
        integer(int64) :: min_ok = -huge(1_int64)
        integer(int64) :: max_ok = +huge(1_int64)
        logical :: append = .False.
        !! if append is True, the argument can be called multiple time.
        !!   each time, the new value is appended.  In this case,
        !!   character-separated values cannot be used.
    contains
        procedure :: set_arg_scalar => set_argint_list_str
        procedure :: set_arg_array => set_argint_list_arr
    end type ArgIntList

    type, private, extends(ArgObj) :: ArgReal
        real(real64) :: value
        real(real64) :: min_ok = -huge(1.0_real64)
        real(real64) :: max_ok = +huge(1.0_real64)
        real(real64) :: const
        logical :: add_value = .False.
    contains
        procedure :: set_arg_scalar => upd_argreal_value
        procedure :: set_arg_array => upd_noarrayR
    end type ArgReal

    type, private, extends(ArgObj) :: ArgRealList
        real(real64), dimension(:), allocatable :: values
        real(real64) :: min_ok = -huge(1.0_real64)
        real(real64) :: max_ok = +huge(1.0_real64)
        logical :: append = .False.
        !! if append is True, the argument can be called multiple time.
        !!   each time, the new value is appended.  In this case,
        !!   character-separated values cannot be used.
    contains
        procedure :: set_arg_scalar => set_argreal_list_str
        procedure :: set_arg_array => set_argreal_list_arr
    end type ArgRealList

    type, private, extends(ArgObj) :: ArgBool
        logical :: value
        logical :: const
    contains
        procedure :: set_arg_scalar => upd_argbool_value
        procedure :: set_arg_array => upd_noarrayB
    end type ArgBool

    type, private, extends(ArgObj) :: ArgChar
        character(len=:), allocatable :: value
    contains
        procedure :: set_arg_scalar => upd_argchar_value
        procedure :: set_arg_array => upd_noarrayC
    end type ArgChar

    type, private, extends(ArgObj) :: ArgCharList
        character(len=:), dimension(:), allocatable :: values
        logical :: append = .False.
        !! if append is True, the argument can be called multiple time.
        !!   each time, the new value is appended.  In this case,
        !!   character-separated values cannot be used.
    contains
        procedure :: set_arg_scalar => set_argchar_list_str
        procedure :: set_arg_array => set_argchar_list_arr
    end type ArgCharList

    type, extends(ArgObj), private :: GenArg
        !! A dummy container to build list of different arguments
        class(ArgObj), allocatable :: arg
    contains
        procedure :: set_arg_scalar => upd_noscalar
        procedure :: set_arg_array => upd_noarray
    end type GenArg

    type, public :: CmdArgDB
        private
        type(GenArg), dimension(MAX_ARGS) :: args
        integer, dimension(3,MAX_ARGS) :: iargs_pos
        character(len=1), dimension(:), allocatable :: prefixes
        character(len=1), dimension(:), allocatable :: prefix_short
        character(len=2), dimension(:), allocatable :: prefix_long
        character(len=:), allocatable :: progname
        integer :: iarg_help = 0
        !! Stores indexes of positional arguments. <0 for arbitrary number
        integer :: nargs = 0, nargs_pos = 0
        ! Stores error status, that can be updated if needed
        ! The first time, contains the initialization status
        class(BaseException), allocatable :: error
    contains
        procedure, private :: get_value_int32val, get_value_int64val, &
            get_value_int32arr, get_value_int64arr, &
            get_value_real32val, get_value_real64val, &
            get_value_real32arr, get_value_real64arr, &
            get_value_boolval, &
            get_value_charval, get_value_chararr
        procedure, private :: get_argname_id
        procedure :: chk_name_overlap => chk_argname_overlap
        procedure :: add_arg_int => add_argument_int
        procedure :: add_arg_real => add_argument_real
        procedure :: add_arg_bool => add_argument_bool
        procedure :: add_arg_char => add_argument_char
        procedure :: is_user_set => argval_set_by_user
        generic :: get_value => get_value_int32val, get_value_int64val, &
            get_value_int32arr, get_value_int64arr, &
            get_value_real32val, get_value_real64val, &
            get_value_real32arr, get_value_real64arr, &
            get_value_boolval, &
            get_value_charval, get_value_chararr
        procedure :: parse_args => parse_args_list
        procedure :: print_help => print_help
        procedure :: has_error => check_argDB_error_status
        procedure :: get_error => get_argDB_error_msg
        procedure :: exception => get_argDB_error
    end type CmdArgDB

    interface CmdArgDB
        module procedure init_args_db
    end interface CmdArgDB

interface

! ----------------------------------------------------------------------

module subroutine parse_args_list(this, arglist)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), dimension(:), intent(in), optional :: arglist
end subroutine parse_args_list

! ----------------------------------------------------------------------

module subroutine add_argument_int(this, argtype, shortname, longname, label, &
    help, required, def_value, min_value, max_value, const_value, add_value, &
    min_nvals, max_nvals)
    class(CmdArgDB), intent(inout) :: this
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
    !! Argument must be given. By default, based on presence of prefix(es).
    class(*), intent(in), optional :: def_value
    !! Default value, to use if not set by user.
    class(*), intent(in), optional :: min_value
    !! Minimum accepted value
    class(*), intent(in), optional :: max_value
    !! Maximum accepted value
    class(*), intent(in), optional :: const_value
    !! Constant value, to use whenever the option is given (no value expected)
    logical, intent(in), optional :: add_value
    !! Increment value associated to argument each time it is encountered.
    class(*), intent(in), optional :: min_nvals
    !! For lists of values, minimum size of the list
    class(*), intent(in), optional :: max_nvals
    !! For lists of values, maixmum size of the list.
end subroutine add_argument_int

! ----------------------------------------------------------------------

module subroutine add_argument_real(this, argtype, shortname, longname, &
    label, help, required, def_value, min_value, max_value, const_value, &
    add_value, min_nvals, max_nvals)
    class(CmdArgDB), intent(inout) :: this
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
    !! Argument must be given. By default, based on presence of prefix(es).
    class(*), intent(in), optional :: def_value
    !! Default value, to use if not set by user.
    class(*), intent(in), optional :: min_value
    !! Minimum accepted value
    class(*), intent(in), optional :: max_value
    !! Maximum accepted value
    class(*), intent(in), optional :: const_value
    !! Constant value, to use whenever the option is given (no value expected)
    logical, intent(in), optional :: add_value
    !! Increment value associated to argument each time it is encountered.
    class(*), intent(in), optional :: min_nvals
    !! For lists of values, minimum size of the list
    class(*), intent(in), optional :: max_nvals
    !! For lists of values, maixmum size of the list.
end subroutine add_argument_real

! ----------------------------------------------------------------------

module subroutine add_argument_bool(this, argtype, shortname, longname, &
    label, help, required, def_value)
    class(CmdArgDB), intent(inout) :: this
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
    !! Argument must be given. By default, based on presence of prefix(es).
    logical, intent(in), optional :: def_value
    !! Default value, to use if not set by user.
end subroutine add_argument_bool

! ----------------------------------------------------------------------

module subroutine add_argument_char(this, argtype, shortname, longname, &
    label, help, required, def_value, min_nvals, max_nvals)
    class(CmdArgDB), intent(inout) :: this
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
    !! Argument must be given. By default, based on presence of prefix(es).
    character(len=*), intent(in), optional :: def_value
    !! Default value, to use if not set by user.
    class(*), intent(in), optional :: min_nvals
    !! For lists of values, minimum size of the list
    class(*), intent(in), optional :: max_nvals
    !! For lists of values, maixmum size of the list.
end subroutine add_argument_char

! ----------------------------------------------------------------------

module function set_arg_names(this, prefixes, short, long, label, &
                              required) result(res)
    class(ArgObj), intent(inout) :: this
    !! Argument
    character(len=*), dimension(:), intent(in) :: prefixes
    !! List of prefix characters.
    character(len=*), intent(in), optional :: short
    !! Short form of the argument.
    character(len=*), intent(in), optional :: long
    !! Long fprm of the argument.
    character(len=*), intent(in), optional :: label
    !! Internal and reference name of the argument.
    logical, intent(in), optional :: required
    !! If True, the argument is required, only for optional argument.
    class(BaseException), allocatable :: res
    !! Return error status.
end function set_arg_names

! ----------------------------------------------------------------------

module function set_arg_helpmsg(this, msg) result(err)
    class(ArgObj), intent(inout) :: this
    !! Argument
    character(len=*), intent(in) :: msg
    !! Help message
    class(BaseException), allocatable :: err
    !! Return error status.
end function set_arg_helpmsg

! ----------------------------------------------------------------------

module function upd_argint_value(this, string, sep) result(err)
    class(ArgInt), intent(inout) :: this
    !! Argument
    character(len=*), intent(in), optional :: string
    !! Value stored as string
    character(len=*), intent(in), optional :: sep
    !! Dummy argument.
    class(BaseException), allocatable :: err
    !! Return error status.
end function upd_argint_value

! ----------------------------------------------------------------------

module function set_argint_list_str(this, string, sep) result(err)
    class(ArgIntList), intent(inout) :: this
    !! Argument
    character(len=*), intent(in), optional :: string
    !! String containing the values
    character(len=*), intent(in), optional :: sep
    !! list of one-character separators, given as a string.
    class(BaseException), allocatable :: err
    !! Return error status.
end function set_argint_list_str

! ----------------------------------------------------------------------

module function set_argint_list_arr(this, strings) result(err)
    class(ArgIntList), intent(inout) :: this
    !! Argument
    character(len=*), dimension(:), intent(in) :: strings
    !! List of values stored as strings
    class(BaseException), allocatable :: err
    !! Return error status.
end function set_argint_list_arr

! ----------------------------------------------------------------------

module function upd_argreal_value(this, string, sep) result(err)
    class(ArgReal), intent(inout) :: this
    !! Argument
    character(len=*), intent(in), optional :: string
    !! Value stored as string
    character(len=*), intent(in), optional :: sep
    !! Dummy argument.
    class(BaseException), allocatable :: err
    !! Return error status.
end function upd_argreal_value

! ----------------------------------------------------------------------

module function set_argreal_list_str(this, string, sep) result(err)
    class(ArgRealList), intent(inout) :: this
    !! Argument
    character(len=*), intent(in), optional :: string
    !! String containing the values
    character(len=*), intent(in), optional :: sep
    !! list of one-character separators, given as a string.
    class(BaseException), allocatable :: err
    !! Return error status.
end function set_argreal_list_str

! ----------------------------------------------------------------------

module function set_argreal_list_arr(this, strings) result(err)
    class(ArgRealList), intent(inout) :: this
    !! Argument
    character(len=*), dimension(:), intent(in) :: strings
    !! List of values stored as strings
    class(BaseException), allocatable :: err
    !! Return error status.
end function set_argreal_list_arr

! ----------------------------------------------------------------------

module function upd_argbool_value(this, string, sep) result(err)
    class(ArgBool), intent(inout) :: this
    !! Argument
    character(len=*), intent(in), optional :: string
    !! Value stored as string
    character(len=*), intent(in), optional :: sep
    !! Dummy argument.
    class(BaseException), allocatable :: err
    !! Return error status.
end function upd_argbool_value

! ----------------------------------------------------------------------

module function upd_argchar_value(this, string, sep) result(err)
    class(ArgChar), intent(inout) :: this
    !! Argument
    character(len=*), intent(in), optional :: string
    !! Value stored as string
    character(len=*), intent(in), optional :: sep
    !! Dummy argument.
    class(BaseException), allocatable :: err
    !! Return error status.
end function upd_argchar_value

! ----------------------------------------------------------------------

module function set_argchar_list_str(this, string, sep) result(err)
    class(ArgCharList), intent(inout) :: this
    !! Argument
    character(len=*), intent(in), optional :: string
    !! String containing the values
    character(len=*), intent(in), optional :: sep
    !! list of one-character separators, given as a string.
    class(BaseException), allocatable :: err
    !! Return error status.
end function set_argchar_list_str

! ----------------------------------------------------------------------

module function set_argchar_list_arr(this, strings) result(err)
    class(ArgCharList), intent(inout) :: this
    !! Argument
    character(len=*), dimension(:), intent(in) :: strings
    !! List of values stored as strings
    class(BaseException), allocatable :: err
    !! Return error status.
end function set_argchar_list_arr

! ----------------------------------------------------------------------

module function argval_set_by_user(this, argname) result(res)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    logical :: res
    !! Boolean stating if value set by user.
end function argval_set_by_user

! ----------------------------------------------------------------------

module subroutine get_value_int32val(this, argname, result)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    integer(int32), intent(out) :: result
    !! Associated value.
end subroutine get_value_int32val

! ----------------------------------------------------------------------

module subroutine get_value_int64val(this, argname, result)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    integer(int64), intent(out) :: result
    !! Associated value.
end subroutine get_value_int64val

! ----------------------------------------------------------------------

module subroutine get_value_int32arr(this, argname, result)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    integer(int32), dimension(:), allocatable, intent(out) :: result
    !! Associated value.
end subroutine get_value_int32arr

! ----------------------------------------------------------------------

module subroutine get_value_int64arr(this, argname, result)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    integer(int64), dimension(:), allocatable, intent(out) :: result
    !! Associated value.
end subroutine get_value_int64arr

! ----------------------------------------------------------------------

module subroutine get_value_real32val(this, argname, result)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    real(real32), intent(out) :: result
    !! Associated value.
end subroutine get_value_real32val

! ----------------------------------------------------------------------

module subroutine get_value_real64val(this, argname, result)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    real(real64), intent(out) :: result
    !! Associated value.
end subroutine get_value_real64val

! ----------------------------------------------------------------------

module subroutine get_value_real32arr(this, argname, result)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    real(real32), dimension(:), allocatable, intent(out) :: result
    !! Associated value.
end subroutine get_value_real32arr

! ----------------------------------------------------------------------

module subroutine get_value_real64arr(this, argname, result)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    real(real64), dimension(:), allocatable, intent(out) :: result
    !! Associated value.
end subroutine get_value_real64arr

! ----------------------------------------------------------------------

module subroutine get_value_boolval(this, argname, result)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    logical, intent(out) :: result
    !! Associated value.
end subroutine get_value_boolval

! ----------------------------------------------------------------------

module subroutine get_value_charval(this, argname, result)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    character(len=*), intent(out) :: result
    !! Associated value.
end subroutine get_value_charval

! ----------------------------------------------------------------------

module subroutine get_value_chararr(this, argname, result)
    class(CmdArgDB), intent(inout) :: this
    !! Arguments database object.
    character(len=*), intent(in) :: argname
    !! Argumemt name/label.
    character(len=:), dimension(:), allocatable, intent(out) :: result
    !! Associated value.
end subroutine get_value_chararr

! ----------------------------------------------------------------------

module function upd_noarrayI(this, strings) result(err)
    class(ArgInt), intent(inout) :: this
    !! Argument
    character(len=*), dimension(:), intent(in) :: strings
    !! List of values stored as strings
    class(BaseException), allocatable :: err
    !! Return error status.
end function upd_noarrayI

! ----------------------------------------------------------------------

module function upd_noarrayR(this, strings) result(err)
    class(ArgReal), intent(inout) :: this
    !! Argument
    character(len=*), dimension(:), intent(in) :: strings
    !! List of values stored as strings
    class(BaseException), allocatable :: err
    !! Return error status.
end function upd_noarrayR

! ----------------------------------------------------------------------

module function upd_noarrayB(this, strings) result(err)
    class(ArgBool), intent(inout) :: this
    !! Argument
    character(len=*), dimension(:), intent(in) :: strings
    !! List of values stored as strings
    class(BaseException), allocatable :: err
    !! Return error status.
end function upd_noarrayB

! ----------------------------------------------------------------------

module function upd_noarrayC(this, strings) result(err)
    class(ArgChar), intent(inout) :: this
    !! Argument
    character(len=*), dimension(:), intent(in) :: strings
    !! List of values stored as strings
    class(BaseException), allocatable :: err
    !! Return error status.
end function upd_noarrayC

! ----------------------------------------------------------------------

module function upd_noscalar(this, string, sep) result(err)
    class(GenArg), intent(inout) :: this
    !! Argument
    character(len=*), intent(in), optional :: string
    !! Dummy argument.
    character(len=*), intent(in), optional :: sep
    !! Dummy argument.
    class(BaseException), allocatable :: err
    !! Return error status.
end function upd_noscalar

! ----------------------------------------------------------------------

module function upd_noarray(this, strings) result(err)
    class(GenArg), intent(inout) :: this
    !! Argument
    character(len=*), dimension(:), intent(in) :: strings
    !! Dummy argument.
    class(BaseException), allocatable :: err
    !! Return error status.
end function upd_noarray

! ----------------------------------------------------------------------

end interface

contains

! ======================================================================

function init_args_db(add_help, prefixes, progname) result(db)
    !! Initialize database to store arguments parameters
    !!
    !! Initializes a database object to store information on supported
    !!   options in the commandline and extract/process data.
    !! The system expects short-name optional argument to be preceded by
    !!   1 prefix character, long names by 2 identical characters.

    logical, intent(in), optional :: add_help
    !! Add help keywords: -h/--help
    character(len=*), intent(in), optional :: prefixes
    !! List of accepted 1-char prefixes for optional arguments, as a string.
    character(len=*), intent(in), optional :: progname
    !! Name of the program.
    type(CmdArgDB) :: db
    !! Database of command-line arguments.

    integer :: i, istat, n_prefix
    logical :: add_help_
    character(len=1) :: prefix

    db%error = InitError()

    if (present(prefixes)) then
        n_prefix = len(prefixes)
        allocate(db%prefixes(n_prefix), db%prefix_short(n_prefix), &
                 db%prefix_long(n_prefix), stat=istat)
        if (istat /= 0) then
            call RaiseAllocateError(db%error, &
                                    'user-defined arguments prefixes')
            return
        end if
        do i = 1, n_prefix
            prefix = prefixes(i:i)
            db%prefixes(i) = prefix
            db%prefix_short(i) = prefix
            db%prefix_long(i) = prefix // prefix
        end do
    else
        allocate(db%prefixes(1), db%prefix_short(1), db%prefix_long(1), &
                 stat=istat)
        if (istat /= 0) then
            Call RaiseAllocateError(db%error, 'default arguments prefixes')
            return
        end if
        db%prefixes(1) = '-'
        db%prefix_short(1) = '-'
        db%prefix_long(1) = '--'
    end if

    if (present(progname)) then
        db%progname = trim(progname)
    else
        call get_command_argument(0, length=i)
        allocate(character(len=i) :: db%progname)
        call get_command_argument(0, db%progname)
    end if

    if (present(add_help)) then
        add_help_ = add_help
    else
        add_help_ = .True.
    end if

    if (add_help_) then
        call db%add_arg_bool('store_true', '-h', '--help', &
                             help='Print this help message', &
                             def_value=.false.)
    end if

end function init_args_db

! ======================================================================

function get_argDB_error(argDB) result(error)
    !! Return error attribute stored in CmdArgDB instance.
    !!
    !! Returns the error attribute currently stored in CmdArgDB.
    !! This is mostly intended for testing, for instance with select
    !! type constructs.
    class(CmdArgDB), intent(in) :: argDB
    !! CmdArgDB instance
    class(BaseException), allocatable :: error
    !! Error instance

    error = argDB%error

end function get_argDB_error

! ======================================================================

function check_argDB_error_status(argDB) result(raised)
    class(CmdArgDB), intent(in) :: argDB
    !! CmdArgDB instance.
    logical :: raised
    !! Status of the error

    raised = argDB%error%raised()

end function check_argDB_error_status

! ======================================================================

function get_argDB_error_msg(argDB) result(msg)
    class(CmdArgDB), intent(in) :: argDB
    !! CmdArgDB instance.
    character(len=:), allocatable :: msg
    !! Error instance.

    msg = argDB%error%msg()

end function get_argDB_error_msg

! ======================================================================

function chk_argname_overlap(this, new_arg, argname) result(res)
    !! Check if new argument overlaps with registered names in DB.
    !!
    !! Checks if any of the identifiers of an argument (short name,
    !!   long name, label) overlaps with arguments already inserted into
    !!   the arguments DB.
    !! Returns True if an overlap exist.

    class(CmdArgDB), intent(in) :: this
    !! Arguments database object
    class(ArgObj), intent(in) :: new_arg
    !! New argument to be inserted in database
    character(len=*), intent(out) :: argname
    !! Overlapping argument name.
    logical :: res
    !! Result of the check

    integer :: iarg

    argname = ' '
    res = .False.
    if (allocated(new_arg%short_name)) then
        do iarg = 1, this%nargs
            associate (opt => this%args(iarg)%arg)
            if (allocated(opt%short_name)) then
                if (opt%short_name == new_arg%short_name) then
                    res = .True.
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
                    res = .True.
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
                res = .True.
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

    class(CmdArgDB), intent(in) :: this
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

    class(CmdArgDB), intent(in) :: this
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
        function tostring(string)
            character(len=*), intent(in) :: string
            character(len=len(string)) :: tostring
        end function tostring
    end interface

    class(ArgObj), intent(in) :: arg
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

function nochange(string)
    ! do nothing (used as target for conversion)
    character(len=*), intent(in) :: string
    character(len=len(string)) :: nochange

    nochange = string
end function


end subroutine print_help

! ======================================================================

end module parse_cmdline
