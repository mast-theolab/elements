module arrays
    use iso_fortran_env, only: int32, int64, real32, real64

    implicit none

    interface symm_tri_array
        module procedure symm_tri_arr_r32, symm_tri_arr_r64
    end interface symm_tri_array

contains

! ======================================================================

subroutine symm_tri_arr_r32(N, arr, linear, lower, anti_symm, arr_new)
    !! Symmetrize/antisymmetrize triangular array (real32)
    !!
    !! Symmetrizes/antisymmetrizes array `arr` based on the lower or
    !! upper triangular content, stored in linear form or not.
    !! If arr_new is provided, the new array is returned in arr_new.
    integer, intent(in) :: N
    !! number of rows/columns of the final square matrix
    real(real32), dimension(:), intent(inout) :: arr
    !! array containing the data to symmetrize/antisymmetrize.
    logical, intent(in), optional :: linear
    !! data are stored in linear form (default: True).
    logical, intent(in), optional :: lower
    !! relevant data is the lower triangular, otherwise the upper.
    logical, intent(in), optional :: anti_symm
    !! Antisymmetrize the array (default: False)
    real(real32), dimension(:,:), intent(out), optional :: arr_new

    integer :: i, j, ij, ioff, k, koff
    logical :: do_symm, is_lin, is_LT

    if (present(linear)) then
        is_lin = linear
    else
        is_lin = .True.
    end if
    if (present(lower)) then
        is_LT = lower
    else
        is_LT = .True.
    end if
    if (present(anti_symm)) then
        do_symm = .not.anti_symm
    else
        do_symm = .True.
    end if

    if (present(arr_new)) then
        if (is_lin) then
            if (is_LT) then
                k = 0
                do i = 1, N
                    do j = 1, i
                        k = k + 1
                        arr_new(i,j) = arr(k)
                    end do
                end do
            else
                k = 0
                do i = 1, N
                    do j = 1, i
                        k = k + 1
                        arr_new(j,i) = arr(k)
                    end do
                end do
            end if
        else
            if (is_LT) then
                do j = N, 1, -1
                    ioff = (j-1)*N
                    arr_new(j,j) = arr(ioff+j)
                    do i = N, j+1, -1
                        ij = ioff + i
                        arr_new(i,j) = arr(ij)
                    end do
                end do
            else
                do j = N, 1, -1
                    ioff = (j-1)*N
                    arr_new(j,j) = arr(ioff+j)
                    do i = j-1, 1, -1
                        ij = ioff + i
                        arr_new(i,j) = arr(ij)
                    end do
                end do
            end if
        end if
        if (is_LT) then
            if (do_symm) then
                do i = 1, N
                    do j = i+1, N
                        arr_new(i,j) = arr_new(j,i)
                    end do
                end do
            else
                do i = 1, N
                    do j = i+1, N
                        arr_new(i,j) = -arr_new(j,i)
                    end do
                end do
            end if
        else
            if (do_symm) then
                do i = 1, N
                    do j = 1, i-1
                        arr_new(i,j) = arr_new(j,i)
                    end do
                end do
            else
                do i = 1, N
                    do j = 1, i-1
                        arr_new(i,j) = -arr_new(j,i)
                    end do
                end do
            end if
        end if
    else
        if (is_lin) then
            if (is_LT) then
                ! The most efficient way to compute the indexes in linear form
                ! is to build the upper triangular first.
                ! In the case of lower-triangular, we need to set the sign
                ! correctly, so we need to differentiate symm and antisymm.
                if (do_symm) then
                    do j = N, 1, -1
                        koff = j*(j-1)/2
                        ioff = (j-1)*N
                        do i = j, 1, -1
                            arr(ioff+i) = arr(koff+i)
                        end do
                    end do
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = 1, j-1
                            arr((i-1)*N+j) = arr(ioff+i)
                        end do
                    end do
                else
                    do j = N, 1, -1
                        koff = j*(j-1)/2
                        ioff = (j-1)*N
                        do i = j, 1, -1
                            arr(ioff+i) = -arr(koff+i)
                        end do
                    end do
                    ! we need to correct the diagonal as well
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = 1, j
                            arr((i-1)*N+j) = -arr(ioff+i)
                        end do
                    end do
                end if
            else
                k = 0
                do j = N, 1, -1
                    koff = j*(j-1)/2
                    ioff = (j-1)*N
                    do i = j, 1, -1
                        arr(ioff+i) = arr(koff+i)
                    end do
                end do
                ! Now let us build the missing lower triangular.
                if (do_symm) then
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = 1, j-1
                            arr((i-1)*N+j) = arr(ioff+i)
                        end do
                    end do
                else
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = 1, j-1
                            arr((i-1)*N+j) = -arr(ioff+i)
                        end do
                    end do
                end if
            end if
        else
            ! Array in full shape, we need only build the mising part
            if (is_LT) then
                if (do_symm) then
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = j+1, N
                            arr((i-1)*N+j) = arr(ioff+i)
                        end do
                    end do
                else
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = j+1, N
                            arr((i-1)*N+j) = -arr(ioff+i)
                        end do
                    end do
                end if
            else
                if (do_symm) then
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = 1, j-1
                            arr((i-1)*N+j) = arr(ioff+i)
                        end do
                    end do
                else
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = 1, j-1
                            arr((i-1)*N+j) = -arr(ioff+i)
                        end do
                    end do
                end if
            end if
        end if
    end if
end subroutine symm_tri_arr_r32

! ======================================================================

subroutine symm_tri_arr_r64(N, arr, linear, lower, anti_symm, arr_new)
    !! Symmetrize/antisymmetrize triangular array (real64)
    !!
    !! Symmetrizes/antisymmetrizes array `arr` based on the lower or
    !! upper triangular content, stored in linear form or not.
    !! If arr_new is provided, the new array is returned in arr_new.
    integer, intent(in) :: N
    !! number of rows/columns of the final square matrix
    real(real64), dimension(:), intent(inout) :: arr
    !! array containing the data to symmetrize/antisymmetrize.
    logical, intent(in), optional :: linear
    !! data are stored in linear form (default: True).
    logical, intent(in), optional :: lower
    !! relevant data is the lower triangular, otherwise the upper.
    logical, intent(in), optional :: anti_symm
    !! Antisymmetrize the array (default: False)
    real(real64), dimension(:,:), intent(out), optional :: arr_new

    integer :: i, j, ij, ioff, k, koff
    logical :: do_symm, is_lin, is_LT

    if (present(linear)) then
        is_lin = linear
    else
        is_lin = .True.
    end if
    if (present(lower)) then
        is_LT = lower
    else
        is_LT = .True.
    end if
    if (present(anti_symm)) then
        do_symm = .not.anti_symm
    else
        do_symm = .True.
    end if

    if (present(arr_new)) then
        if (is_lin) then
            if (is_LT) then
                k = 0
                do i = 1, N
                    do j = 1, i
                        k = k + 1
                        arr_new(i,j) = arr(k)
                    end do
                end do
            else
                k = 0
                do i = 1, N
                    do j = 1, i
                        k = k + 1
                        arr_new(j,i) = arr(k)
                    end do
                end do
            end if
        else
            if (is_LT) then
                do j = N, 1, -1
                    ioff = (j-1)*N
                    arr_new(j,j) = arr(ioff+j)
                    do i = N, j+1, -1
                        ij = ioff + i
                        arr_new(i,j) = arr(ij)
                    end do
                end do
            else
                do j = N, 1, -1
                    ioff = (j-1)*N
                    arr_new(j,j) = arr(ioff+j)
                    do i = j-1, 1, -1
                        ij = ioff + i
                        arr_new(i,j) = arr(ij)
                    end do
                end do
            end if
        end if
        if (is_LT) then
            if (do_symm) then
                do i = 1, N
                    do j = i+1, N
                        arr_new(i,j) = arr_new(j,i)
                    end do
                end do
            else
                do i = 1, N
                    do j = i+1, N
                        arr_new(i,j) = -arr_new(j,i)
                    end do
                end do
            end if
        else
            if (do_symm) then
                do i = 1, N
                    do j = 1, i-1
                        arr_new(i,j) = arr_new(j,i)
                    end do
                end do
            else
                do i = 1, N
                    do j = 1, i-1
                        arr_new(i,j) = -arr_new(j,i)
                    end do
                end do
            end if
        end if
    else
        if (is_lin) then
            if (is_LT) then
                ! The most efficient way to compute the indexes in linear form
                ! is to build the upper triangular first.
                ! In the case of lower-triangular, we need to set the sign
                ! correctly, so we need to differentiate symm and antisymm.
                if (do_symm) then
                    do j = N, 1, -1
                        koff = j*(j-1)/2
                        ioff = (j-1)*N
                        do i = j, 1, -1
                            arr(ioff+i) = arr(koff+i)
                        end do
                    end do
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = 1, j-1
                            arr((i-1)*N+j) = arr(ioff+i)
                        end do
                    end do
                else
                    do j = N, 1, -1
                        koff = j*(j-1)/2
                        ioff = (j-1)*N
                        do i = j, 1, -1
                            arr(ioff+i) = -arr(koff+i)
                        end do
                    end do
                    ! we need to correct the diagonal as well
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = 1, j
                            arr((i-1)*N+j) = -arr(ioff+i)
                        end do
                    end do
                end if
            else
                k = 0
                do j = N, 1, -1
                    koff = j*(j-1)/2
                    ioff = (j-1)*N
                    do i = j, 1, -1
                        arr(ioff+i) = arr(koff+i)
                    end do
                end do
                ! Now let us build the missing lower triangular.
                if (do_symm) then
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = 1, j-1
                            arr((i-1)*N+j) = arr(ioff+i)
                        end do
                    end do
                else
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = 1, j-1
                            arr((i-1)*N+j) = -arr(ioff+i)
                        end do
                    end do
                end if
            end if
        else
            ! Array in full shape, we need only build the mising part
            if (is_LT) then
                if (do_symm) then
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = j+1, N
                            arr((i-1)*N+j) = arr(ioff+i)
                        end do
                    end do
                else
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = j+1, N
                            arr((i-1)*N+j) = -arr(ioff+i)
                        end do
                    end do
                end if
            else
                if (do_symm) then
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = 1, j-1
                            arr((i-1)*N+j) = arr(ioff+i)
                        end do
                    end do
                else
                    do j = 1, N
                        ioff = (j-1)*N
                        do i = 1, j-1
                            arr((i-1)*N+j) = -arr(ioff+i)
                        end do
                    end do
                end if
            end if
        end if
    end if
end subroutine symm_tri_arr_r64

! ======================================================================

end module arrays