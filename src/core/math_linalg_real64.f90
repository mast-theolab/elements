submodule (math) math_linalg_real64
    !! Math submodule with 64-bit real linear algebra routines.
    !!
    !! Sub-module to the `math` module containing the definition of the
    !! procedures for 64-bit real linear algebra algorithms.

contains

! ======================================================================

module procedure d_cross

    vecC(1) = vecA(2)*vecB(3) - vecA(3)*vecB(2)
    vecC(2) = vecA(3)*vecB(1) - vecA(1)*vecB(3)
    vecC(3) = vecA(1)*vecB(2) - vecA(2)*vecB(1)

end procedure d_cross

! ======================================================================

module procedure d_cross_mat_mat
    integer :: i, j, m, n
    real(real64), dimension(3) :: vec

    m = size(matA, 2) ; n = size(matB, 2)
    allocate(matC(3,m,n))

    do i = 1, n
        vec = matB(:,i)
        do j = 1, m
            matC(1,j,i) = matA(2,j)*vec(3) - matA(3,j)*vec(2)
            matC(2,j,i) = matA(3,j)*vec(1) - matA(1,j)*vec(3)
            matC(3,j,i) = matA(1,j)*vec(2) - matA(2,j)*vec(1)
        end do
    end do

end procedure d_cross_mat_mat

! ======================================================================

module procedure d_cross_mat_vec

    integer :: i, m

    m = size(matA,2)

    allocate(matC(3,m))

    do i = 1, m
        matC(1,i) = matA(2,i)*vecB(3) - matA(3,i)*vecB(2)
        matC(2,i) = matA(3,i)*vecB(1) - matA(1,i)*vecB(3)
        matC(3,i) = matA(1,i)*vecB(2) - matA(2,i)*vecB(1)
    end do

end procedure d_cross_mat_vec

! ======================================================================

module procedure d_cross_vec_mat

    integer :: i, n

    n = size(matB,2)

    allocate(matC(3,n))

    do i = 1, n
        matC(1,i) = vecA(2)*matB(3,i) - VecA(3)*matB(2,i)
        matC(2,i) = vecA(3)*matB(1,i) - VecA(1)*matB(3,i)
        matC(3,i) = vecA(1)*matB(2,i) - VecA(2)*matB(1,i)
    end do

end procedure d_cross_vec_mat

! ======================================================================

module procedure d_det
    real(real64), allocatable :: work(:)
    integer, allocatable      :: ipiv(:)
    integer                   :: n, i, info

    n = size(A,1)

    allocate(ipiv(n))
    allocate(work(2*n))

    call xgetrf(n, n, A, n, ipiv, info)
    if (info /= 0) stop 'Matrix is numerically singular!'

    det_A = 0.0_real64
    do i = 1, n
        det_A = det_A * A(i, i)
    end do

    do i = 1, n
        if (ipiv(i) /= i) then
            det_A = -det_A
        end if
    end do

end procedure d_det

! ======================================================================

module procedure d_inv_mat
    real(real64), allocatable :: work(:)
    integer, allocatable     :: ipiv(:)
    integer                  :: n, info

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

end procedure d_inv_mat

! ======================================================================

end submodule math_linalg_real64