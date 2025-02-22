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

    ! "Sensible" thresholds for scientific calculations
    ! These thresholds are used internally for closeness tests.
    ! near0: to test numbers in terms of numeric precision (epsilon)
    ! small: for final quantities with values typically of magnitude around 1
    ! null0: values below this thresholds can be safely considered null.
    real(real32), parameter, private :: &
        near0_r32 = 1.0e-6_real32, &
        small_r32 = 1.0e-4_real32, &
        null0_r32 = 1.0e-32_real32
    real(real64), parameter, private :: &
        near0_r64 = 1.0e-10_real64, &
        small_r64 = 1.0e-6_real64, &
        null0_r64 = 1.0e-32_real64
        
    interface operator(.iscloseto.)
        module procedure :: is_close_to_r32_generic, is_close_to_r64_generic
    end interface

    interface is_close_to
        module procedure :: is_close_to_r32, is_close_to_r64
    end interface

contains

! ======================================================================

elemental function is_close_to_r32(value, target, rel_tol, abs_tol) result(res)
    !! Check if value is close to target within chosen tolerance(s).
    !!
    !! Checks if a given value is close to a target value within numeric
    !! tolerances.  The tolerances can be defined as absolute or relative.
    !! The test is the following:
    !!
    !! |value-target| <= max(rel_tol * max(|value|, |target|), abs_tol)
    !!
    !! @note
    !! Negative numbers in tolerances values are ignored.
    !! The default values are used in this case.
    !! @endnote
    real(real32), intent(in) :: value
    !! Value to check.
    real(real32), intent(in) :: target
    !! Target for the closeness check.
    real(real32), intent(in), optional :: rel_tol
    !! Relative tolerance threshold, by default as epsilon(value)*10
    real(real32), intent(in), optional :: abs_tol
    !! Absolute tolerance threshold, by default as 10^-32.
    !! The default value is chosen sensible in scientific applications.
    logical :: res
    !! result of the closeness test.

    real(real32) :: tol_abs, tol_rel

    if (present(rel_tol)) then
        if (rel_tol < 0.0_real32) then
            tol_rel = near0_r32
        else
            tol_rel = rel_tol
        end if
    else
        tol_rel = near0_r32
    end if
    if (present(abs_tol)) then
        if (abs_tol < 0.0_real32) then
            tol_abs = null0_r32
        end if
        tol_abs = abs_tol
    else
        tol_abs = null0_r32
    end if

    res = abs(value-target) <= max(tol_rel * max(abs(value), abs(target)), &
                                   tol_abs)
end function is_close_to_r32

! ======================================================================

elemental function is_close_to_r32_generic(value, target) result(res)
    !! Check if value is close to target within chosen tolerance(s).
    !!
    !! Checks if a given value is close to a target value within numeric
    !! standard tolerances for scientific applications.
    !!
    !! |value-target| <= max(rel_tol * max(|value|, |target|), abs_tol)
    real(real32), intent(in) :: value
    !! Value to check.
    real(real32), intent(in) :: target
    !! Target for the closeness check.
    logical :: res
    !! result of the closeness test.

    res = abs(value-target) <= max(near0_r32 * max(abs(value), abs(target)), &
                                   null0_r32)
end function is_close_to_r32_generic

! ======================================================================

elemental function is_close_to_r64(value, target, rel_tol, abs_tol) result(res)
    !! Check if value is close to target within chosen tolerance(s).
    !!
    !! Checks if a given value is close to a target value within numeric
    !! tolerances.  The tolerances can be defined as absolute or relative.
    !! The test is the following:
    !!
    !! |value-target| <= max(rel_tol * max(|value|, |target|), abs_tol)
    !!
    !! @note
    !! Negative numbers in tolerances values are ignored.
    !! The default values are used in this case.
    !! @endnote
    real(real64), intent(in) :: value
    !! Value to check.
    real(real64), intent(in) :: target
    !! Target for the closeness check.
    real(real64), intent(in), optional :: rel_tol
    !! Relative tolerance threshold, by default as epsilon(value)*10
    real(real64), intent(in), optional :: abs_tol
    !! Absolute tolerance threshold, by default as 10^-32.
    !! The default value is chosen sensible in scientific applications.
    logical :: res
    !! result of the closeness test.

    real(real64) :: tol_abs, tol_rel

    if (present(rel_tol)) then
        if (rel_tol < 0.0_real64) then
            tol_rel = near0_r64
        else
            tol_rel = rel_tol
        end if
    else
        tol_rel = near0_r64
    end if
    if (present(abs_tol)) then
        if (abs_tol < 0.0_real64) then
            tol_abs = null0_r64
        else
            tol_abs = abs_tol
        end if
    else
        tol_abs = null0_r64
    end if

    res = abs(value-target) <= max(tol_rel * max(abs(value), abs(target)), &
                                   tol_abs)
end function is_close_to_r64

! ======================================================================

elemental function is_close_to_r64_generic(value, target) result(res)
    !! Check if value is close to target within chosen tolerance(s).
    !!
    !! Checks if a given value is close to a target value within numeric
    !! standard tolerances for scientific applications.
    !!
    !! |value-target| <= max(rel_tol * max(|value|, |target|), abs_tol)
    real(real64), intent(in) :: value
    !! Value to check.
    real(real64), intent(in) :: target
    !! Target for the closeness check.
    logical :: res
    !! result of the closeness test.

    res = abs(value-target) <= max(near0_r64 * max(abs(value), abs(target)), &
                                   null0_r64)
end function is_close_to_r64_generic

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
