module blas_drv
    
    implicit none

! ----------------------------------------------------------------------

    interface xgemv
        !! Generic interface to BLAS xGEMV
        subroutine sgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            character :: trans
            integer :: m, n, lda, incx, incy
            real :: alpha, beta
            real, dimension(lda, *) :: A
            real, dimension(*) :: x, y
        end subroutine sgemv
        subroutine dgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            character :: trans
            integer :: m, n, lda, incx, incy
            double precision :: alpha, beta
            double precision, dimension(lda, *) :: A
            double precision, dimension(*) :: x, y
        end subroutine dgemv
        subroutine cgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            character :: trans
            integer :: m, n, lda, incx, incy
            complex :: alpha, beta
            complex, dimension(lda, *) :: A
            complex, dimension(*) :: x, y
        end subroutine cgemv
        subroutine zgemv(trans, m, n, alpha, A, lda, x, incx, beta, y, incy)
            character :: trans
            integer :: m, n, lda, incx, incy
            complex*16 :: alpha, beta
            complex*16, dimension(lda, *) :: A
            complex*16, dimension(*) :: x, y
        end subroutine zgemv
    end interface xgemv

! ----------------------------------------------------------------------

    interface xgemm
        !! Generic interface to BLAS xGEMM
        subroutine sgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta,&
                         C, ldc)
            character :: transa, transb
            integer :: m, n, k, lda, ldb, ldc
            real :: alpha, beta
            real, dimension(lda, *) :: A
            real, dimension(ldb, *) :: B
            real, dimension(ldc, *) :: C
        end subroutine sgemm
        subroutine dgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta,&
                         C, ldc)
            character :: transa, transb
            integer :: m, n, k, lda, ldb, ldc
            double precision :: alpha, beta
            double precision, dimension(lda, *) :: A
            double precision, dimension(ldb, *) :: B
            double precision, dimension(ldc, *) :: C
        end subroutine dgemm
        subroutine cgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta,&
                         C, ldc)
            character :: transa, transb
            integer :: m, n, k, lda, ldb, ldc
            complex :: alpha, beta
            complex, dimension(lda, *) :: A
            complex, dimension(ldb, *) :: B
            complex, dimension(ldc, *) :: C
        end subroutine cgemm
        subroutine zgemm(transa, transb, m, n, k, alpha, A, lda, B, ldb, beta,&
                         C, ldc)
            character :: transa, transb
            integer :: m, n, k, lda, ldb, ldc
            complex*16 :: alpha, beta
            complex*16, dimension(lda, *) :: A
            complex*16, dimension(ldb, *) :: B
            complex*16, dimension(ldc, *) :: C
        end subroutine zgemm
    end interface xgemm

    interface xdot
        !! Generic interface to BLAS xDOT
        function ddot(n, x, incx, y, incy)
            integer :: n, incx, incy
            double precision, dimension(*) :: x, y
            double precision :: ddot
        end function ddot
        function sdot(n, x, incx, y, incy)
            integer :: n, incx, incy
            real, dimension(*) :: x, y
            real :: sdot
        end function sdot
        function sdsdot(n, b, x, incx, y, incy)
            integer :: n, incx, incy
            real :: b
            real, dimension(*) :: x, y
            real :: sdsdot
        end function sdsdot
        module procedure cdot, zdot
    end interface

contains

! ======================================================================

function cdot(trans, n, cx, incx, cy, incy)
    character :: trans
    integer :: n, incx, incy
    complex, dimension(*) :: cx, cy
    complex :: cdot

    interface
        function cdotc(n_, cx_, incx_, cy_, incy_)
            integer :: n_, incx_, incy_
            complex, dimension(*) :: cx_, cy_
            complex :: cdotc
        end function cdotc
        function cdotu(n_, cx_, incx_, cy_, incy_)
            integer :: n_, incx_, incy_
            complex, dimension(*) :: cx_, cy_
            complex :: cdotu
        end function cdotu
    end interface

    select case(trans)
    case ('H', 'h')
        cdot = cdotc(n, cx, incx, cy, incy)
    case ('T', 't')
        cdot = cdotu(n, cx, incx, cy, incy)
    case default
        stop 1
    end select

end function cdot

! ======================================================================

function zdot(trans, n, cx, incx, cy, incy)
    character :: trans
    integer :: n, incx, incy
    complex*16, dimension(*) :: cx, cy
    complex*16 :: zdot

    interface
        function zdotc(n_, cx_, incx_, cy_, incy_)
            integer :: n_, incx_, incy_
            complex*16, dimension(*) :: cx_, cy_
            complex*16 :: zdotc
        end function zdotc
        function zdotu(n_, cx_, incx_, cy_, incy_)
            integer :: n_, incx_, incy_
            complex*16, dimension(*) :: cx_, cy_
            complex*16 :: zdotu
        end function zdotu
    end interface

    select case(trans)
    case ('H', 'h')
        zdot = zdotc(n, cx, incx, cy, incy)
    case ('T', 't')
        zdot = zdotu(n, cx, incx, cy, incy)
    case default
        stop 1
    end select

end function zdot

end module blas_drv
