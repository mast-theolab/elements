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

! ----------------------------------------------------------------------

interface
    module function d_cross(vecA, vecB) result(vecC)
        !! Cross vector between 2 64-bits real Cartesian vectors.
        real(real64), dimension(3), intent(in) :: vecA
            !! vector A
        real(real64), dimension(3), intent(in) :: vecB
            !! vector B
        real(real64), dimension(3) :: vecC
            !! vector C
    end function d_cross

    module function d_cross_mat_mat(matA, matB) result(matC)
        !! Computes the cross vector between 2 lists of Cartesian vectors.
        real(real64), dimension(:,:), intent(in) :: matA
            !! list of vectors A.
        real(real64), dimension(:,:), intent(in) :: matB
            !! list of vectors B.
        real(real64), dimension(:,:,:), allocatable :: matC
            !! list of lists of vectors C.
    end function d_cross_mat_mat

    module function d_cross_mat_vec(matA, vecB) result(matC)
        !! Computes the cross product between a list of 64-bits real
        !! Cartesian vectors and 1 64-bits real vector.
        real(real64), dimension(:,:), intent(in) :: matA
            !! list of vectors A (dimension: 3,N).
        real(real64), dimension(3), intent(in) :: vecB
            !! vector B.
        real(real64), dimension(:,:), allocatable :: matC
            !! list of vectors C.
    end function d_cross_mat_vec

    module function d_cross_vec_mat(vecA, matB) result(matC)
        !! Computes the cross product between 1 64-bits real Cartesian
        !! vector and a list of 64-bits real vectors.
        real(real64), dimension(3), intent(in) :: vecA
            !! vector A.
        real(real64), dimension(:,:), intent(in) :: matB
            !! list of vectors B (dimension: 3,N).
        real(real64), dimension(:,:), allocatable :: matC
            !! list of vectors C.
    end function d_cross_vec_mat

    module function s_cross(vecA, vecB) result(vecC)
        !! Cross vector between 2 32-bits real Cartesian vectors.
        real(real32), dimension(3), intent(in) :: vecA
            !! vector A
        real(real32), dimension(3), intent(in) :: vecB
            !! vector B
        real(real32), dimension(3) :: vecC
            !! vector C
    end function s_cross

    module function s_cross_mat_mat(matA, matB) result(matC)
        !! Computes the cross vector between 2 lists of Cartesian vectors.
        real(real32), dimension(:,:), intent(in) :: matA
            !! list of vectors A.
        real(real32), dimension(:,:), intent(in) :: matB
            !! list of vectors B.
        real(real32), dimension(:,:,:), allocatable :: matC
            !! list of lists of vectors C.
    end function s_cross_mat_mat

    module function s_cross_mat_vec(matA, vecB) result(matC)
        !! Computes the cross product between a list of 32-bits real
        !! Cartesian vectors and 1 32-bits real vector.
        real(real32), dimension(:,:), intent(in) :: matA
            !! list of vectors A (dimension: 3,N).
        real(real32), dimension(3), intent(in) :: vecB
            !! vector B.
        real(real32), dimension(:,:), allocatable :: matC
            !! list of vectors C.
    end function s_cross_mat_vec

    module function s_cross_vec_mat(vecA, matB) result(matC)
        !! Computes the cross product between 1 32-bits real Cartesian
        !! vector and a list of 32-bits real vectors.
        real(real32), dimension(3), intent(in) :: vecA
            !! vector A.
        real(real32), dimension(:,:), intent(in) :: matB
            !! list of vectors B (dimension: 3,N).
        real(real32), dimension(:,:), allocatable :: matC
            !! list of vectors C.
    end function s_cross_vec_mat

end interface

! ----------------------------------------------------------------------

    interface operator(.x.)
        !! Operator to compute the cross product: C = A x B
        !! @note
        !! `A` or `B` can be list of vectors, in which case `C` is a 
        !! list of vectors
        !! @endnote
        module procedure s_cross, d_cross, &
            s_cross_vec_mat, d_cross_vec_mat, &
            s_cross_mat_vec, d_cross_mat_vec, &
            s_cross_mat_mat, d_cross_mat_mat
    end interface operator(.x.)

! ----------------------------------------------------------------------

    interface cross
        !! Compute the cross product: C = A x B
        module procedure s_cross, d_cross, &
            s_cross_vec_mat, d_cross_vec_mat, &
            s_cross_mat_vec, d_cross_mat_vec, &
            s_cross_mat_mat, d_cross_mat_mat
    end interface cross

! ----------------------------------------------------------------------

    interface inv_mat
        !! Invert matrices using LAPACK drivers.
        module subroutine s_inv_mat(A, A_inv)
            !! Invert a real (32 bits) matrix using LAPACK.
            use lapack_drv, only: xgetrf, xgetri

            real(real32), dimension(:, :), intent(in)  :: A
                !! Input matrix A of size (n, n)
            real(real32), dimension(:, :), intent(out) :: A_inv
                !! Inverse of the input matrix A
        end subroutine s_inv_mat

        module subroutine d_inv_mat(A, A_inv)
            !! Invert a real (64 bits) matrix using LAPACK.
            use lapack_drv, only: xgetrf, xgetri

            real(real64), dimension(:, :), intent(in)  :: A
                !! Input matrix A of size (n, n)
            real(real64), dimension(:, :), intent(out) :: A_inv
                !! Inverse of the input matrix A
        end subroutine d_inv_mat

        module subroutine c_inv_mat(A, A_inv)
            !! Invert a complex (32 bits) matrix using LAPACK.
            use lapack_drv, only: xgetrf, xgetri

            complex(real32), dimension(:, :), intent(in)  :: A
                !! Input complex matrix A of size (n, n)
            complex(real32), dimension(:, :), intent(out) :: A_inv
                !! Inverse of the input complex matrix A
        end subroutine c_inv_mat

        module subroutine z_inv_mat(A, A_inv)
            !! Invert a complex (64 bits) matrix using LAPACK.
            use lapack_drv, only: xgetrf, xgetri

            complex(real64), dimension(:, :), intent(in)  :: A
                !! Input complex matrix A of size (n, n)
            complex(real64), dimension(:, :), intent(out) :: A_inv
                !! Inverse of the input complex matrix A
        end subroutine z_inv_mat
    end interface inv_mat

! ----------------------------------------------------------------------

    interface det
        !! Compute the determinant of a matrix
        module function s_det(A) result(det_A)
            !! Compute the determinant of a real (32 bits) matrix.
            use lapack_drv, only: xgetrf

            real(real32), dimension(:, :), intent(in) :: A
                !! Input matrix A of size (n, n)
            real(real32) :: det_A
                !! Output determinant det_A
        end function s_det

        module function d_det(A) result(det_A)
            !! Compute the determinant of a real (64 bits) matrix.
            use lapack_drv, only: xgetrf

            real(real64), dimension(:, :), intent(in) :: A
                !! Input matrix A of size (n, n)
            real(real64) :: det_A
                !! Output determinant det_A
        end function d_det

        module function c_det(A) result(det_A)
            !! Compute the determinant of a complex (32 bits) matrix.
            use lapack_drv, only: xgetrf

            complex(real32), dimension(:, :), intent(in) :: A
            !! Input complex matrix A of size (n, n)
            complex(real32) :: det_A
            !! Output determinant det_A
        end function c_det

        module function z_det(A) result(det_A)
            !! Compute the determinant of a complex (64 bits) matrix.
            use lapack_drv, only: xgetrf

            complex(real64), dimension(:, :), intent(in) :: A
            !! Input complex matrix A of size (n, n)
            complex(real64) :: det_A
            !! Output determinant det_A
        end function z_det

    end interface det

contains

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

recursive function factorial(n) result(n1)
    !! Compute the factorial n!
    !!
    !! Given a value n, computes the corresponding factorial
    integer, intent(in) :: n
        !! Positive number.
    integer :: n1
        !! Factorial result.

    if (n < 0) then
        n1 = -1
    else if (n <= 1) then
        n1 = 1
    else
        n1 = n*factorial(n-1)
    end if

    return
end function factorial

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
    do i = 1, n, 2
        b = b*real(i, kind=real64)
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
