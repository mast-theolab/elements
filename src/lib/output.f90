module output
    !! Output-related module for GenTensor
    !!
    !! Different output/printing constants and procedures for GenTensor
    use iso_fortran_env, only: real32, real64, int32, int64, output_unit

    interface prt_mat
        module procedure :: prt_mat_r32, prt_mat_r64
    end interface prt_mat

contains

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

end module output
