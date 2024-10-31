program test_blas_ops
    !! Program to test the BLAS interface with some basic operations.

    use iso_fortran_env, only: real32, real64, real128
    use blas_drv
    use output, only: sec_header, iu_out

    integer, parameter :: N = 10
    real(real32), dimension(N) :: XR32, YR32, ZR32, ZZR32
    real(real32), dimension(N, N) :: AR32, BR32, CR32, CCR32
    real(real64), dimension(N) :: XR64, YR64, ZR64, ZZR64
    real(real64), dimension(N, N) :: AR64, BR64, CR64, CCR64
    complex(real32), dimension(N) :: XC32, YC32, ZC32, ZZC32
    complex(real32), dimension(N, N) :: AC32, BC32, CC32, CCC32
    complex(real64), dimension(N) :: XC64, YC64, ZC64, ZZC64
    complex(real64), dimension(N, N) :: AC64, BC64, CC64, CCC64

    real(real32) :: aval32, rmin32, rmax32
    real(real64) :: aval64, rmin64, rmax64
    complex(real32) :: cval32
    complex(real64) :: cval64

    real(real32) :: f0r32 = 0.0_real32, f1r32 = 1.0_real32
    real(real64) :: f0r64 = 0.0_real64, f1r64 = 1.0_real64
    complex(real32) :: f0c32 = (0.0_real32, 0.0_real32), &
        f1c32 = (1.0_real32, 0.0_real32)
    complex(real64) :: f0c64 = (0.0_real64, 0.0_real64), &
        f1c64 = (1.0_real64, 0.0_real64)

    real(real32), parameter :: zero_r32 = epsilon(f0r32)*N
    real(real64), parameter :: zero_r64 = epsilon(f0r64)*N
    complex(real32), parameter :: zero_c32 = epsilon(f0r32)*N
    complex(real64), parameter :: zero_c64 = epsilon(f0r64)*N

    character(len=:), allocatable :: what

    1000 format(3x,a," -> ",a)

    call sec_header(-1, 'Test Program for BLAS Interface')

    call random_number(XR32)
    call random_number(YR32)
    call random_number(AR32)
    call random_number(BR32)
    call random_number(XR64)
    call random_number(YR64)
    call random_number(AR64)
    call random_number(BR64)
    call random_number(XC32%re)
    call random_number(YC32%re)
    call random_number(AC32%re)
    call random_number(BC32%re)
    call random_number(XC64%re)
    call random_number(YC64%re)
    call random_number(AC64%re)
    call random_number(BC64%re)
    call random_number(XC32%im)
    call random_number(YC32%im)
    call random_number(AC32%im)
    call random_number(BC32%im)
    call random_number(XC64%im)
    call random_number(YC64%im)
    call random_number(AC64%im)
    call random_number(BC64%im)

    call sec_header(1, 'Testing xGEMM')

    call sec_header(2, 'Real 32 bits')

    what = 'C = A * B'

    call xgemm('N', 'N', N, N, N, f1r32, AR32, N, BR32, N, f0r32, CR32, N)
    if (any(abs(CR32-matmul(AR32, BR32)) > zero_r32)) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    what = 'C = A^T * B'

    call xgemm('T', 'N', N, N, N, f1r32, AR32, N, BR32, N, f0r32, CR32, N)
    if (any(abs(CR32-matmul(transpose(AR32), BR32)) > zero_r32)) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    call sec_header(2, 'Real 64 bits')

    what = 'C = A * B'

    call xgemm('N', 'N', N, N, N, f1r64, AR64, N, BR64, N, f0r64, CR64, N)
    if (any(abs(CR64-matmul(AR64, BR64)) > zero_r64)) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    what = 'C = A^T * B'

    call xgemm('T', 'N', N, N, N, f1r64, AR64, N, BR64, N, f0r64, CR64, N)
    if (any(abs(CR64-matmul(transpose(AR64), BR64)) > zero_r64)) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    call sec_header(2, 'complex 32 bits')

    what = 'C = A * B'

    call xgemm('N', 'N', N, N, N, f1c32, AC32, N, BC32, N, f0c32, CC32, N)
    CCC32 = CC32-matmul(AC32, BC32)
    CCC32 = conjg(CCC32) * CCC32
    rmin32 = minval(CCC32%re)
    rmax32 = maxval(CCC32%re)
    if (abs(rmin32) > zero_r32 .or. rmax32 > zero_r32) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    what = 'C = A^T * B'

    call xgemm('T', 'N', N, N, N, f1c32, AC32, N, BC32, N, f0c32, CC32, N)
    CCC32 = CC32-matmul(transpose(AC32), BC32)
    CCC32 = conjg(CCC32) * CCC32
    rmin32 = minval(CCC32%re)
    rmax32 = maxval(CCC32%re)
    if (abs(rmin32) > zero_r32 .or. rmax32 > zero_r32) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    call sec_header(2, 'complex 64 bits')

    what = 'C = A * B'

    call xgemm('N', 'N', N, N, N, f1c64, AC64, N, BC64, N, f0c64, CC64, N)
    CCC64 = CC64-matmul(AC64, BC64)
    CCC64 = conjg(CCC64) * CCC64
    rmin64 = minval(CCC64%re)
    rmax64 = maxval(CCC64%re)
    if (abs(rmin64) > zero_r64 .or. rmax64 > zero_r64) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    what = 'C = A^T * B'

    call xgemm('T', 'N', N, N, N, f1c64, AC64, N, BC64, N, f0c64, CC64, N)
    CCC64 = CC64-matmul(transpose(AC64), BC64)
    CCC64 = conjg(CCC64) * CCC64
    rmin64 = minval(CCC64%re)
    rmax64 = maxval(CCC64%re)
    if (abs(rmin64) > zero_r64 .or. rmax64 > zero_r64) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if


    call sec_header(1, 'Testing xGEMV')

    call sec_header(2, 'Real 32 bits')

    what = 'Z = A * X'

    call xgemv('N', N, N, f1r32, AR32, N, XR32, 1, f0r32, ZR32, 1)
    if (any(abs(ZR32-matmul(AR32, XR32)) > zero_r32)) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    what = 'Z = A^T * X'

    call xgemv('T', N, N, f1r32, AR32, N, XR32, 1, f0r32, ZR32, 1)
    if (any(abs(ZR32-matmul(transpose(AR32), XR32)) > zero_r32)) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    call sec_header(2, 'Real 64 bits')

    what = 'Z = A * X'

    call xgemv('N', N, N, f1r64, AR64, N, XR64, 1, f0r64, ZR64, 1)
    if (any(abs(ZR64-matmul(AR64, XR64)) > zero_r64)) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    what = 'Z = A^T * X'

    call xgemv('T', N, N, f1r64, AR64, N, XR64, 1, f0r64, ZR64, 1)
    if (any(abs(ZR64-matmul(transpose(AR64), XR64)) > zero_r64)) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    call sec_header(2, 'complex 32 bits')

    what = 'Z = A * X'

    call xgemv('N', N, N, f1c32, AC32, N, XC32, 1, f0c32, ZC32, 1)
    ZZC32 = ZC32-matmul(AC32, XC32)
    ZZC32 = conjg(ZZC32) * ZZC32
    rmin32 = minval(ZZC32%re)
    rmax32 = maxval(ZZC32%re)
    if (abs(rmin32) > zero_r32 .or. rmax32 > zero_r32) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    what = 'Z = A^T * X'

    call xgemv('T', N, N, f1c32, AC32, N, XC32, 1, f0c32, ZC32, 1)
    ZZC32 = ZC32-matmul(transpose(AC32), XC32)
    ZZC32 = conjg(ZZC32) * ZZC32
    rmin32 = minval(ZZC32%re)
    rmax32 = maxval(ZZC32%re)
    if (abs(rmin32) > zero_r32 .or. rmax32 > zero_r32) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    call sec_header(2, 'complex 64 bits')

    what = 'Z = A * X'

    call xgemv('N', N, N, f1c64, AC64, N, XC64, 1, f0c64, ZC64, 1)
    ZZC64 = ZC64-matmul(AC64, XC64)
    ZZC64 = conjg(ZZC64) * ZZC64
    rmin64 = minval(ZZC64%re)
    rmax64 = maxval(ZZC64%re)
    if (abs(rmin64) > zero_r64 .or. rmax64 > zero_r64) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    what = 'Z = A^T * X'

    call xgemv('T', N, N, f1c64, AC64, N, XC64, 1, f0c64, ZC64, 1)
    ZZC64 = ZC64-matmul(transpose(AC64), XC64)
    ZZC64 = conjg(ZZC64) * ZZC64
    rmin64 = minval(ZZC64%re)
    rmax64 = maxval(ZZC64%re)
    if (abs(rmin64) > zero_r64 .or. rmax64 > zero_r64) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if


    call sec_header(1, 'Testing xDOT')

    call sec_header(2, 'Real 32 bits')

    what = 'C = A * B'

    aval32 = xdot(N, XR32, 1, YR32, 1)
    if (abs(aval32 - dot_product(XR32, YR32)) > zero_r32) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    call sec_header(2, 'Real 64 bits')

    what = 'C = A * B'

    aval64 = xdot(N, XR64, 1, YR64, 1)
    if (abs(aval64 - dot_product(XR64, YR64)) > zero_r64) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    call sec_header(2, 'complex 32 bits')

    what = 'Z = X^T * Y'

    cval32 = xdot('T', N, XC32, 1, YC32, 1)
    cval32 = cval32 - dot_product(conjg(XC32), YC32)
    cval32 = conjg(cval32) * cval32
    if (cval32%re > zero_r32) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    what = 'Z = X^H * Y'

    cval32 = xdot('H', N, XC32, 1, YC32, 1)
    cval32 = cval32 - dot_product(XC32, YC32)
    cval32 = conjg(cval32) * cval32
    if (cval32%re > zero_r32) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    call sec_header(2, 'complex 64 bits')

    what = 'Z = X^T * Y'

    cval64 = xdot('T', N, XC64, 1, YC64, 1)
    cval64 = cval64 - dot_product(conjg(XC64), YC64)
    cval64 = conjg(cval64) * cval64
    if (cval64%re > zero_r64) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if

    what = 'Z = X^H * Y'

    cval64 = xdot('H', N, XC64, 1, YC64, 1)
    cval64 = cval64 - dot_product(XC64, YC64)
    cval64 = conjg(cval64) * cval64
    if (cval64%re > zero_r64) then
        print 1000, what, 'FAILED'
        stop 1
    else
        print 1000, what, 'PASSED'
    end if



    ! print '(10F10.7)', CR32
    ! print '(10F10.7)', CR32-matmul(AR32, BR32)
    ! call xgemv('N', N, N, 1.0_real32, AR32X N, B132, N, 0.Z_real12, CR32, N)

end program test_blas_ops
