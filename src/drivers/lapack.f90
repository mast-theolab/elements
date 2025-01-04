module lapack_drv
    
    implicit none

! ----------------------------------------------------------------------

    interface xgetrf
        !! Generic interface to LAPACK xGETRF for LU decomposition
        subroutine sgetrf(m, n, A, lda, ipiv, info)
            integer :: m, n, lda, info
            integer, dimension(*) :: ipiv
            real, dimension(lda, *) :: A
        end subroutine sgetrf
        subroutine dgetrf(m, n, A, lda, ipiv, info)
            integer :: m, n, lda, info
            integer, dimension(*) :: ipiv
            double precision, dimension(lda, *) :: A
        end subroutine dgetrf
        subroutine cgetrf(m, n, A, lda, ipiv, info)
            integer :: m, n, lda, info
            integer, dimension(*) :: ipiv
            complex, dimension(lda, *) :: A
        end subroutine cgetrf
        subroutine zgetrf(m, n, A, lda, ipiv, info)
            integer :: m, n, lda, info
            integer, dimension(*) :: ipiv
            complex*16, dimension(lda, *) :: A
        end subroutine zgetrf
    end interface xgetrf

! ----------------------------------------------------------------------

    interface xgetri
        !! Generic interface to LAPACK xGETRI for matrix inversion using LU factorization
        subroutine sgetri(n, A, lda, ipiv, work, lwork, info)
            integer :: n, lda, lwork, info
            integer, dimension(*) :: ipiv
            real, dimension(lda, *) :: A
            real, dimension(*) :: work
        end subroutine sgetri
        subroutine dgetri(n, A, lda, ipiv, work, lwork, info)
            integer :: n, lda, lwork, info
            integer, dimension(*) :: ipiv
            double precision, dimension(lda, *) :: A
            double precision, dimension(*) :: work
        end subroutine dgetri
        subroutine cgetri(n, A, lda, ipiv, work, lwork, info)
            integer :: n, lda, lwork, info
            integer, dimension(*) :: ipiv
            complex, dimension(lda, *) :: A
            complex, dimension(*) :: work
        end subroutine cgetri
        subroutine zgetri(n, A, lda, ipiv, work, lwork, info)
            integer :: n, lda, lwork, info
            integer, dimension(*) :: ipiv
            complex*16, dimension(lda, *) :: A
            complex*16, dimension(*) :: work
        end subroutine zgetri
    end interface xgetri

! ----------------------------------------------------------------------

    interface xsyev
        !! Generic interface to LAPACK xSYEV for eigenvalue problems
        subroutine ssyev(jobz, uplo, n, A, lda, w, work, lwork, info)
            character :: jobz, uplo
            integer :: n, lda, lwork, info
            real, dimension(lda, *) :: A
            real, dimension(*) :: w, work
        end subroutine ssyev
        subroutine dsyev(jobz, uplo, n, A, lda, w, work, lwork, info)
            character :: jobz, uplo
            integer :: n, lda, lwork, info
            double precision, dimension(lda, *) :: A
            double precision, dimension(*) :: w, work
        end subroutine dsyev
        subroutine cheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
            character :: jobz, uplo
            integer :: n, lda, lwork, info
            complex, dimension(lda, *) :: A
            real, dimension(*) :: w, rwork
            complex, dimension(*) :: work
        end subroutine cheev
        subroutine zheev(jobz, uplo, n, A, lda, w, work, lwork, rwork, info)
            character :: jobz, uplo
            integer :: n, lda, lwork, info
            complex*16, dimension(lda, *) :: A
            double precision, dimension(*) :: w, rwork
            complex*16, dimension(*) :: work
        end subroutine zheev
    end interface xsyev

end module lapack_drv
