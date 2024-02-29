module output
    !! Output-related module
    !!
    !! Different output/printing-related procedures
    use iso_fortran_env, only: real32, real64, int32, int64, output_unit
    use string, only: locase

    integer :: iu_out = output_unit

    interface prt_mat
        module procedure :: prt_mat_r32, prt_mat_r64
    end interface prt_mat

    interface len_int
        module procedure :: len_int32, len_int64
    end interface len_int

    interface num_digits_int
        module procedure num_digits_int32, num_digits_int64
    end interface

contains

! ======================================================================

subroutine write_err(nature, msg, error_)
    !! Write error message on default unit
    !!
    !! Writes an error message.  The formatting depends on the nature
    !!   of the error:
    !! * Generic/gen: generic error
    !! * Basic/std: basic/standard error
    implicit none

    character(len=*), intent(in) :: nature
    !! Nature of the error
    character(len=*), intent(in) :: msg
    !! General error message
    character(len=*), intent(in), optional :: error_

    select case (locase(trim(nature)))
        case ('generic', 'gen')
            write(iu_out, '(a)') trim(msg)
            write(iu_out, '(a)') 'Stopping'
        case ('basic', 'std')
            write(iu_out, '(a)') trim(msg)
            write(iu_out, '("Reason:",/,4x,a)') trim(error_)
        case default
            write(iu_out, '("Uncategorized error:",4x,a)') trim(msg)
    end select
end subroutine write_err

! ======================================================================

subroutine prt_mat_r32(mat, n_row, n_col, iunit_, ncol_by_blk_, prec_, thresh_)
    !! Print a real32 matrix
    !!
    !! Prints a real matrix with precision real32.
    implicit none

    real(real32), dimension(:,:), intent(in) :: mat
    !! Matrix to display
    integer, intent(in) :: n_row
    !! Number of rows to display
    integer, intent(in) :: n_col
    !! Number of columns to display
    integer, intent(in), optional :: iunit_
    !! Output unit
    integer, intent(in), optional :: ncol_by_blk_
    !! Number of columns per block
    integer, intent(in), optional :: prec_
    !! Number of digits for precision (<0 for fixed-point notation)
    real(real32), intent(in), optional :: thresh_
    !! If set, elements below thresh_ are not printed.
    integer :: icol, icol0, irow, iu, len_num, n, ncol_by_blk
    integer, parameter :: len_id = 8
    real(real32), dimension(:), allocatable :: vec
    character(len=10) :: num_fmt, id_fmt
    character(len=80) :: dfmt, hfmt

    ! Set output
    if (present(iunit_)) then
        iu = iunit_
    else
        iu = iu_out
    end if

    ! Build formats
    write(id_fmt, '("i",i0)') len_id
    if (.not.present(prec_)) then
        num_fmt = 'es14.6'
        len_num = 14
    else
        if (prec_ > 0) then
            len_num = prec_ + 8
            write(num_fmt, '("es",i0,".",i0)') len_num, prec_
        else
            len_num = -prec_ + 10
            write(num_fmt, '("f",i0,".",i0)') len_num, -prec_
        end if
    end if

    if (.not.present(ncol_by_blk_)) then
        ncol_by_blk = 5
    else
        ncol_by_blk = abs(ncol_by_blk_)
    end if
    allocate(vec(ncol_by_blk))

    write(hfmt, '("(",i0,"x,",i0,"(",i0,"x,",a,"))")') 1+len_id/2, &
        ncol_by_blk, len_num-len_id, id_fmt
    write(dfmt, '("(",a,",1x,",i0,"(",a,"))")') id_fmt, ncol_by_blk, num_fmt

    ! Now print the matrix
    do icol0 = 0, n_col-1, ncol_by_blk
        n = min(n_col-icol0, ncol_by_blk)
        write(iu, hfmt) (icol+icol0, icol=1,n)
        do irow = 1, n_row
            if (present(thresh_)) then
                do icol = 1, n
                    if (abs(mat(irow,icol0+icol)) > thresh_) then
                        vec(icol) = mat(irow,icol0+icol)
                    else
                        vec(icol) = 0.0_real32
                    end if
                end do
            else
                vec = mat(irow,icol0+1:icol0+n)
            end if
            write(iu, dfmt) irow, vec(:n)
        end do
        ! do irow = 1, n_row
        !     write(iu, dfmt) irow, (mat(irow,icol), icol=icol0,maxcol)
        ! end do
    end do

end subroutine prt_mat_r32

! ======================================================================

subroutine prt_mat_r64(mat, n_row, n_col, iunit_, ncol_by_blk_, prec_, thresh_)
    !! Print a real64 matrix
    !!
    !! Prints a real matrix with precision real64.
    implicit none

    real(real64), dimension(:,:), intent(in) :: mat
    !! Matrix to display
    integer, intent(in) :: n_row
    !! Number of rows to display
    integer, intent(in) :: n_col
    !! Number of columns to display
    integer, intent(in), optional :: iunit_
    !! Output unit
    integer, intent(in), optional :: ncol_by_blk_
    !! Number of columns per block
    integer, intent(in), optional :: prec_
    !! Number of digits for precision (<0 for fixed-point notation)
    real(real64), intent(in), optional :: thresh_
    !! If set, elements below thresh_ are not printed.
    integer :: icol, icol0, irow, iu, len_num, n, ncol_by_blk
    integer, parameter :: len_id = 8
    real(real64), dimension(:), allocatable :: vec
    character(len=10) :: num_fmt, id_fmt
    character(len=80) :: dfmt, hfmt

    ! Set output
    if (present(iunit_)) then
        iu = iunit_
    else
        iu = output_unit
    end if

    ! Build formats
    write(id_fmt, '("i",i0)') len_id
    if (.not.present(prec_)) then
        num_fmt = 'es14.6'
        len_num = 14
    else
        if (prec_ > 0) then
            len_num = prec_ + 8
            write(num_fmt, '("es",i0,".",i0)') len_num, prec_
        else
            len_num = -prec_ + 10
            write(num_fmt, '("f",i0,".",i0)') len_num, -prec_
        end if
    end if

    if (.not.present(ncol_by_blk_)) then
        ncol_by_blk = 5
    else
        ncol_by_blk = abs(ncol_by_blk_)
    end if
    allocate(vec(ncol_by_blk))

    write(hfmt, '("(",i0,"x,",i0,"(",i0,"x,",a,"))")') 1+len_id/2, &
        ncol_by_blk, len_num-len_id, id_fmt
    write(dfmt, '("(",a,",1x,",i0,"(",a,"))")') id_fmt, ncol_by_blk, num_fmt

    ! Now print the matrix
    do icol0 = 0, n_col-1, ncol_by_blk
        n = min(n_col-icol0, ncol_by_blk)
        write(iu, hfmt) (icol+icol0, icol=1,n)
        do irow = 1, n_row
            if (present(thresh_)) then
                do icol = 1, n
                    if (abs(mat(irow,icol0+icol)) > thresh_) then
                        vec(icol) = mat(irow,icol0+icol)
                    else
                        vec(icol) = 0.0_real64
                    end if
                end do
            else
                vec = mat(irow,icol0+1:icol0+n)
            end if
            write(iu, dfmt) irow, vec(:n)
        end do
    end do

end subroutine prt_mat_r64

! ======================================================================

integer function len_int32(num) result(lnum)
    !! Length of 32-bit integer num
    !!
    !! Returns the minimum number of characters necessary to store
    !!    `num`.
    implicit none

    integer(int32), intent(in) :: num

    lnum = floor(log10(abs(real(num, kind=real64)))) + 1
    if (num < 0) lnum = lnum + 1

end function len_int32

! ======================================================================

integer function len_int64(num) result(lnum)
    !! Length of 64-bit integer num
    !!
    !! Returns the minimum number of characters necessary to store
    !!    `num`.
    implicit none

    integer(int64), intent(in) :: num

    lnum = floor(log10(abs(real(num, kind=real64)))) + 1
    if (num < 0) lnum = lnum + 1

end function len_int64

! ======================================================================

function num_digits_int32(number)
    !! Finds the number of digits in an integer `number`.
    !!
    !! This functions is useful for instance to set a proper format for
    !! strings.
    implicit none

    integer(int32) :: number
    !! Number of interest.
    integer(int32) :: num_digits_int32
    !! Number of digits in `Number`.

    integer(int32) :: ioff

    if (number < 0) then
        ioff = 2
    else
        ioff = 1
    end if

    num_digits_int32 = floor(log10(real(abs(number), kind=real32))) + ioff

end function num_digits_int32

! ======================================================================

function num_digits_int64(number)
    !! Finds the number of digits in an integer `number`.
    !!
    !! This functions is useful for instance to set a proper format for
    !! strings.
    implicit none

    integer(int64) :: number
    !! Number of interest.
    integer(int64) :: num_digits_int64
    !! Number of digits in `Number`.

    integer(int64) :: ioff

    if (number < 0) then
        ioff = 2
    else
        ioff = 1
    end if

    num_digits_int64 = floor(log10(real(abs(number), kind=real64))) + ioff

end function num_digits_int64

! ======================================================================

end module output
