module numeric
    use iso_fortran_env, only: int16, int32, int64, real32, real64

    implicit none

contains

! ======================================================================

function is_integer(arg) result(res)
    !! Check that argument is an integer.
    !!
    !! Checks that argument arg is a valid integer.
    class(*), intent(in) :: arg
    !! Argument to check.
    logical :: res
    !! Result of check.

    select type (arg)
        type is (integer(int16))
            res = .True.
        type is (integer(int32))
            res = .True.
        type is (integer(int64))
            res = .True.
        class default
            res = .False.
    end select

end function is_integer

! ======================================================================

function is_real(arg) result(res)
    !! Check that argument is a real number.
    !!
    !! Checks that argument arg is a valid real number.
    class(*), intent(in) :: arg
    !! Argument to check.
    logical :: res
    !! Result of check.

    select type (arg)
        type is (real(real32))
            res = .True.
        type is (real(real64))
            res = .True.
        class default
            res = .False.
    end select

end function is_real

! ======================================================================

function is_number(arg) result(res)
    !! Check that argument is a number.
    !!
    !! Checks that argument arg is a valid number.
    class(*), intent(in) :: arg
    !! Argument to check.
    logical :: res
    !! Result of check.

    select type (arg)
        type is (integer(int16))
            res = .True.
        type is (integer(int32))
            res = .True.
        type is (integer(int64))
            res = .True.
        type is (real(real32))
            res = .True.
        type is (real(real64))
            res = .True.
        class default
            res = .False.
    end select

end function is_number

! ======================================================================

function to_int64(arg) result(res)
    !! Convert to integer.
    !!
    !! Converts argument `arg` to integer.
    class(*), intent(in) :: arg
    !! Argument to check.
    integer(int64) :: res
    !! Converted number.

    select type (arg)
        type is (integer(int16))
            res = int(arg)
        type is (integer(int32))
            res = int(arg)
        type is (integer(int64))
            res = int(arg)
        type is (real(real32))
            res = int(arg)
        type is (real(real64))
            res = int(arg)
        class default
            error stop 'Cannot convert to integer'
    end select

end function to_int64

! ======================================================================

function to_real64(arg) result(res)
    !! Convert to real.
    !!
    !! Converts argument `arg` to real.
    class(*), intent(in) :: arg
    !! Argument to check.
    real(real64) :: res
    !! Converted number.

    select type (arg)
        type is (integer(int16))
            res = real(arg)
        type is (integer(int32))
            res = real(arg)
        type is (integer(int64))
            res = real(arg)
        type is (real(real32))
            res = real(arg)
        type is (real(real64))
            res = real(arg)
        class default
            error stop 'Cannot convert to real'
    end select

end function to_real64

! ======================================================================

end module numeric