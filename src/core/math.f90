module math
    !! Module containing mathematics-related elements
    !!
    !! Module with additional mathematics-related procedures
    use iso_fortran_env, only: real32, real64

    implicit none

    integer, dimension(:,:), allocatable :: itri_pa
    real(real64), dimension(:,:), allocatable :: rtri_pa
    real(real64), parameter :: pi_r64 = 4.0_real64*atan(1.0_real64)
    real(real32), parameter :: pi_r32 = 4.0_real32*atan(1.0_real32)

    interface cross
        module procedure s_cross, d_cross
    end interface cross

    interface inv_mat
        module procedure s_inv_mat, d_inv_mat, c_inv_mat, z_inv_mat
    end interface inv_mat

    interface det
        module procedure s_det, d_det, c_det, z_det
    end interface det

contains

! ======================================================================

subroutine s_inv_mat(A, A_inv)
    !! Invert a real (32 bits) matrix using LAPACK.
    use lapack_drv, only: xgetrf, xgetri

    real(real32), dimension(:, :), intent(in)  :: A
    !! Input matrix A of size (n, n)
    real(real32), dimension(:, :), intent(out) :: A_inv
    !! Inverse of the input matrix A

    real(real32), allocatable :: work(:)
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

end subroutine s_inv_mat

! ======================================================================

subroutine d_inv_mat(A, A_inv)
    !! Invert a real (64 bits) matrix using LAPACK.
    use lapack_drv, only: xgetrf, xgetri

    real(real64), dimension(:, :), intent(in)  :: A
    !! Input matrix A of size (n, n)
    real(real64), dimension(:, :), intent(out) :: A_inv
    !! Inverse of the input matrix A

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

end subroutine d_inv_mat

! ======================================================================

subroutine c_inv_mat(A, A_inv)
    !! Invert a complex (32 bits) matrix using LAPACK.
    use lapack_drv, only: xgetrf, xgetri

    complex(real32), dimension(:, :), intent(in)  :: A
    !! Input complex matrix A of size (n, n)
    complex(real32), dimension(:, :), intent(out) :: A_inv
    !! Inverse of the input complex matrix A

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

end subroutine c_inv_mat

! ======================================================================

subroutine z_inv_mat(A, A_inv)
    !! Invert a complex (64 bits) matrix using LAPACK.
    use lapack_drv, only: xgetrf, xgetri

    complex(real64), dimension(:, :), intent(in)  :: A
    !! Input complex matrix A of size (n, n)
    complex(real64), dimension(:, :), intent(out) :: A_inv
    !! Inverse of the input complex matrix A

    complex(real64), allocatable :: work(:)
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

end subroutine z_inv_mat

! ======================================================================

function s_det(A) result(det_A)
    !! Compute the determinant of a real (32 bits) matrix.
    use lapack_drv, only: xgetrf

    real(real32), dimension(:, :), intent(in) :: A
    !! Input matrix A of size (n, n)
    real(real32) :: det_A
    !! Output determinant det_A

    real(real32), allocatable :: work(:)
    integer, allocatable      :: ipiv(:)
    integer                   :: n, i, info

    n = size(A,1)

    allocate(ipiv(n))
    allocate(work(2*n))

    call xgetrf(n, n, A, n, ipiv, info)
    if (info /= 0) stop 'Matrix is numerically singular!'

    det_A = 0.0_real32
    do i = 1, n
        det_A = det_A * A(i, i)
    end do

    do i = 1, n
        if (ipiv(i) /= i) then
            det_A = -det_A
        end if
    end do

end function s_det

! ======================================================================

function d_det(A) result(det_A)
    !! Compute the determinant of a real (64 bits) matrix.
    use lapack_drv, only: xgetrf

    real(real64), dimension(:, :), intent(in) :: A
    !! Input matrix A of size (n, n)
    real(real64) :: det_A
    !! Output determinant det_A

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

end function d_det

! ======================================================================

function c_det(A) result(det_A)
    !! Compute the determinant of a complex (32 bits) matrix.
    use lapack_drv, only: xgetrf

    complex(real32), dimension(:, :), intent(in) :: A
    !! Input complex matrix A of size (n, n)
    complex(real32) :: det_A
    !! Output determinant det_A

    complex(real32), allocatable :: work(:)
    integer, allocatable         :: ipiv(:)
    integer                      :: n, i, info

    n = size(A,1)

    allocate(ipiv(n))
    allocate(work(2*n))

    call xgetrf(n, n, A, n, ipiv, info)
    if (info /= 0) stop 'Matrix is numerically singular!'

    det_A = (0.0_real32, 0.0_real32)
    do i = 1, n
        det_A = det_A * A(i, i)
    end do

    do i = 1, n
        if (ipiv(i) /= i) then
            det_A = -det_A
        end if
    end do

end function c_det

! ======================================================================

function z_det(A) result(det_A)
    !! Compute the determinant of a complex (64 bits) matrix.
    use lapack_drv, only: xgetrf

    complex(real64), dimension(:, :), intent(in) :: A
    !! Input complex matrix A of size (n, n)
    complex(real64) :: det_A
    !! Output determinant det_A

    complex(real64), allocatable :: work(:)
    integer, allocatable         :: ipiv(:)
    integer                      :: n, i, info

    n = size(A,1)

    allocate(ipiv(n))
    allocate(work(2*n))

    call xgetrf(n, n, A, n, ipiv, info)
    if (info /= 0) stop 'Matrix is numerically singular!'

    det_A = (0.0_real64, 0.0_real64)
    do i = 1, n
        det_A = det_A * A(i, i)
    end do

    do i = 1, n
        if (ipiv(i) /= i) then
            det_A = -det_A
        end if
    end do

end function z_det

! ======================================================================

recursive function factorial(n) result(n1)
    !! Compute the factorial n!
    !!
    !! Given a value n, computes the corresponding factorial
    integer :: n1
    !! factorial
    integer, intent(in) :: n
    !! positive number number

    if (n < 0) then
        n1 = -1
    else if (n <= 1) then
        n1 = n
    else
        n1 = n*factorial(n-1)
    end if

    return
end function factorial

! ======================================================================

function d_cross(vecA, vecB) result(vecC)
    !! Compute the cross product: C = A x B
    !!
    !! Computes the cross vector between 2 Cartesian vectors.
    !! @note: double precision version
    real(real64), dimension(3) :: vecA
    !! vector A
    real(real64), dimension(3) :: vecB
    !! vector B
    real(real64), dimension(3) :: vecC
    !! vector C

    VecC(1) = vecA(2)*VecB(3) - VecA(3)*VecB(2)
    VecC(2) = vecA(3)*VecB(1) - VecA(1)*VecB(3)
    VecC(3) = vecA(1)*VecB(2) - VecA(2)*VecB(1)

    return
end function d_cross

! ======================================================================

function s_cross(vecA, vecB) result(vecC)
    !! Compute the cross product: C = A x B
    !!
    !! Computes the cross vector between 2 Cartesian vectors.
    !! @note: simple precision version
    real(real32), dimension(3) :: vecA
    !! vector A
    real(real32), dimension(3) :: vecB
    !! vector B
    real(real32), dimension(3) :: vecC
    !! vector C

    VecC(1) = vecA(2)*VecB(3) - VecA(3)*VecB(2)
    VecC(2) = vecA(3)*VecB(1) - VecA(1)*VecB(3)
    VecC(3) = vecA(1)*VecB(2) - VecA(2)*VecB(1)

    return
end function s_cross

! ======================================================================

subroutine build_PascalTriangle(n, do_real)
    !! Build and store Pascal's triangle
    !!
    !! Builds a 2D array with the coefficients of Pascal's triangle
    !!   up to a chosen limit.
    integer, intent(in) :: n
    !! Size of the triangle
    logical, intent(in), optional :: do_real
    !! If present and true, build a version with real coefficients

    integer :: i, j

    allocate(itri_pa(0:n,0:n))
    do i = 0, n
        itri_pa(i,0) = 1
        itri_pa(i,i) = 1
    end do
    do i = 2, n
        do j = 1, i-1
            itri_pa(i,j) = itri_pa(i-1,j-1) + itri_pa(i-1,j)
        end do
    end do

    if (present(do_real)) then
        if (do_real) rtri_pa = real(itri_pa, kind=real64)
    end if

    return
end subroutine build_PascalTriangle

! ======================================================================

pure function int_xn_e2ax2(n, a) result(res)
    !! Compute the integral int(x^n*exp(-2*a*x^2),x=-inf..inf)
    !!
    !! Computes and returns the integral:
    !! \[ \int_{-\infty}^{+\infty} x^n e^{-2 a x^2} dx \]
    integer, intent(in) :: n
    !! Power of x
    real(real64), intent(in) :: a
    real(real64) :: res
    !! Exponential coefficient

    integer :: i
    real(real64) :: b, as

    if (mod(n, 2) /= 0) then
        res = 0.0_real64
        return
    end if
    b = 1.0_real64
    i = 1
    do
        b = b*real(i, kind=real64)
        i = i + 2
        if (i > n-1) exit
    end do
    as = sqrt(a)
    res = b*sqrt(2.0_real64*pi_r64)/(as*2.0_real64)**(n+1)
    return
end function int_xn_e2ax2

! ======================================================================

elemental function phii_xn_phij(only_R, ni, xi, ai, nj, xj, aj, nk) result(res)
    !! Integral between Gaussian-type AOs, int(phi_j(x) x^n phi_j(x) dx)
    !!
    !! Computes the overlap integral between two Gaussian-type atomic
    !!   orbitals, in the form:
    !! \[ I = \int_{-\infty}^{+\infty} dx
    !!      (x-x_i)^{n_i} e^{-a_i (x-x_i)^2}
    !!      x^{n_k}
    !!      (x-x_j)^{n_j} e^{-a_j (x-x_j)^2} \]
    !! which can also be rewritten,
    !! \[ I = \sqrt{\pi} * R * e^{-a_i a_j * {\delta_{ij}}^2/(a_i + a_j)} \]
    !! with \(\delta_{ij} = x_j - x_i\)
    !! Alternatively, only the R part can be returned.
    logical, intent(in) :: only_R
    !! Return only R instead of the whole integral
    integer, intent(in) :: ni
    !! Power of the x coordinate in orbital i
    real(real64), intent(in) :: xi
    !! Center of orbital i
    real(real64), intent(in) :: ai
    !! Exponential coefficient for orbital i
    integer, intent(in) :: nj
    !! Power of the x coordinate in orbital j
    real(real64), intent(in) :: xj
    !! Center of orbital j
    real(real64), intent(in) :: aj
    !! Exponential coefficient for orbital j
    integer, intent(in) :: nk
    !! Power nk in \(x^{n_k}\)
    real(real64) :: res

    integer :: ii, jj, kk
    real(real64) :: ci, cj, ck, fac, faci, facj, ovaiaj, xij

    ovaiaj = 1.0_real64/sqrt(ai + aj)
    xij = xj - xi

    if (only_R) then
        fac = ovaiaj**(ni+nj+nk+1)/sqrt(pi_r64)
    else
        fac = ovaiaj**(ni+nj+nk+1)
    end if
    ci = aj*xij*ovaiaj
    cj = -ai*xij*ovaiaj
    ck = (ai*xi + aj*xj)*ovaiaj
    res = 0.0_real64
    do ii = 0, ni
        faci = itri_pa(ni,ii)*ci**(ni-ii)
        do jj = 0, nj
            facj = faci*itri_pa(nj,jj)*cj**(nj-jj)
            do kk = 0, nk
                res = res + &
                    facj*itri_pa(nk,kk)*ck**(nk-kk)*int_xn_e2ax2(ii+jj+kk, &
                                                                 0.5_real64)
            end do
        end do
    end do
    res = fac*res
    if (.not. only_R) &
        res = res * exp(-ai*aj*xij**2/(ai+aj))

    return
end function phii_xn_phij

! ======================================================================

end module math
