submodule (math) math_linalg_cplx32
    !! Math submodule with 32-bit complex linear algebra routines.
    !!
    !! Sub-module to the `math` module containing the definition of the
    !! procedures for 32-bit complex linear algebra algorithms.

contains

! ======================================================================

module procedure c_det
    use lapack_drv, only: xgetrf

    complex(real32), allocatable :: work(:)
    integer, allocatable         :: ipiv(:)
    integer                      :: n, i, info

    n = size(A,1)

    allocate(ipiv(n))
    allocate(work(2*n))

    call xgetrf(n, n, A, n, ipiv, info)
    if (info /= 0) stop 'Matrix is numerically singular!'

    det_A = (1.0_real32, 0.0_real32)
    do i = 1, n
        det_A = det_A * A(i, i)
    end do

    do i = 1, n
        if (ipiv(i) /= i) then
            det_A = -det_A
        end if
    end do

end procedure c_det

! ======================================================================

module procedure c_inv_mat
    use lapack_drv, only: xgetrf, xgetri

    complex(real32), allocatable :: work(:)
    integer, allocatable         :: ipiv(:)
    integer                      :: n, info

    n = size(A,1)
    if (size(A,2) /= n) then
        stop 'Matrix must be square!'
    end if

    allocate(ipiv(n), work(2*n))
    A_inv = A

    call xgetrf(n, n, A_inv, n, ipiv, info)
    if (info /= 0) stop 'Matrix is numerically singular!'

    call xgetri(n, A_inv, n, ipiv, work, 2*n, info)
    if (info /= 0) stop 'Matrix inversion failed!'

end procedure c_inv_mat

! ======================================================================

end submodule math_linalg_cplx32
