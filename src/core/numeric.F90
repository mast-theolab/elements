module numeric
    use iso_fortran_env, only: int16, int32, int64, real32, real64

    implicit none

#ifdef USE_R8
    integer, parameter :: realwp = real64
#elif defined USE_R4
    integer, parameter :: realwp = real32
#else
    integer, parameter :: realwp = real64
#endif
    !! Working precision for real, can be changed through 

    real(realwp), parameter :: &
        f0 = 0.0_realwp, f1 = 1.0_realwp, f2 = 2.0_realwp, &
        f3 = 3.0_realwp, f4 = 4.0_realwp, f5 = 5.0_realwp, &
        f6 = 6.0_realwp, f7 = 7.0_realwp, f8 = 8.0_realwp, &
        f9 = 9.0_realwp, f10 = 10.0_realwp, &
        fhalf = 0.5_realwp, fquart = 0.25_realwp, &
        f8th = 0.125_realwp, f16th = 0.0625_realwp

    real(realwp), parameter :: &
        small = 1.0e-6_realwp, &
        !! Values can be considered with respect to typical precision.
        near0 = epsilon(0.0_realwp)*f10
        !! Lower values are negligible for a suitable unit of the quantity.

    real(realwp), parameter :: &
        pi = f4*atan(f1)

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
            res = .true.
        type is (real(real64))
            res = .true.
        class default
            res = .false.
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
            res = .true.
        type is (integer(int32))
            res = .true.
        type is (integer(int64))
            res = .true.
        type is (real(real32))
            res = .true.
        type is (real(real64))
            res = .true.
        class default
            res = .false.
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
            res = int(arg, kind=real32)
        type is (integer(int32))
            res = int(arg, kind=real32)
        type is (integer(int64))
            res = int(arg, kind=real32)
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
            res = real(arg, kind=real64)
        type is (integer(int32))
            res = real(arg, kind=real64)
        type is (integer(int64))
            res = real(arg, kind=real64)
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
