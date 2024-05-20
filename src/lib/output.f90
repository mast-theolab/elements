module output
    !! Output-related module
    !!
    !! Different output/printing-related procedures to:
    !! - print a real matrix: prt_mat -> prt_mat_r32, prt_mat64
    !! - write an error message -> write_err
    !! - print a section header or title in consistent format -> sec_header
    !! - print atomic coordinates -> prt_coord
    use iso_fortran_env, only: real32, real64, int32, int64, output_unit
    use string, only: locase
    use physic, only: PhysFact

    integer :: iu_out = output_unit
    type(PhysFact), private :: phys

    interface prt_mat
        module procedure :: prt_mat_r32, prt_mat_r64
    end interface prt_mat

    interface prt_coord
        module procedure :: prt_coord_r32, prt_coord_r64
    end interface prt_coord

    interface len_int
        module procedure :: len_int32, len_int64
    end interface len_int

contains

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

subroutine prt_coord_r32(n_at, at_lab, at_crd, at_mass_, iunit_)
    !! Prints the atomic coordinates in a formatted way (32bit version).
    !!
    !! Prints the atomic coordinates as a table, optionally with the
    !! atomic masses included.
    implicit none

    integer, intent(in) :: n_at
    !! Number of atoms
    character(len=*), dimension(n_at), intent(in) :: at_lab
    !! Atomic labels
    real(real32), dimension(3, n_at), intent(in) :: at_crd
    !! Atomic coordinates (in au)
    real(real32), dimension(n_at), intent(in), optional :: at_mass_
    !! Atomic masses (in u)
    integer, intent(in), optional :: iunit_
    !! Output unit

    integer :: ia, iu

    1000 format(" Atom |    Mass    |        X            Y            Z")
    1001 format(" -----+------------+----------------------------------------")
    1010 format(2x,a,2x,"|",2x,f8.4,2x,"|",3(1x,f12.6))
    1100 format(" Atom |        X            Y            Z")
    1101 format(" -----+----------------------------------------")
    1110 format(2x,a,2x,"|",3(1x,f12.6))

    ! Set output
    if (present(iunit_)) then
        iu = iunit_
    else
        iu = iu_out
    end if

    if (present(at_mass_)) then
        write(iu, 1000)
        write(iu, 1001)
        do ia = 1, n_at
            write(iu, 1010) at_lab(ia), at_mass_(ia), &
                phys%bohr2Ang(at_crd(:,ia))
        end do
        write(iu, 1001)
    else
        write(iu, 1100)
        write(iu, 1101)
        do ia = 1, n_at
            write(iu, 1110) at_lab(ia), at_mass_(ia), &
                phys%bohr2Ang(at_crd(:,ia))
        end do
        write(iu, 1101)
    end if

end subroutine prt_coord_r32

! ======================================================================

subroutine prt_coord_r64(n_at, at_lab, at_crd, at_mass_, iunit_)
    !! Prints the atomic coordinates in a formatted way (64bit version).
    !!
    !! Prints the atomic coordinates as a table, optionally with the
    !! atomic masses included.
    implicit none

    integer, intent(in) :: n_at
    !! Number of atoms
    character(len=*), dimension(n_at), intent(in) :: at_lab
    !! Atomic labels
    real(real64), dimension(3, n_at), intent(in) :: at_crd
    !! Atomic coordinates (in au)
    real(real64), dimension(n_at), intent(in), optional :: at_mass_
    !! Atomic masses (in u)
    integer, intent(in), optional :: iunit_
    !! Output unit

    integer :: ia, iu

    1000 format(" Atom      Mass             X            Y            Z")
    1001 format(" ===== ============ ========================================")
    1002 format(" ----- ------------ ----------------------------------------")
    1010 format(2x,a,5x,f8.4,3x,3(1x,f12.6))
    1100 format(" Atom |        X            Y            Z")
    1101 format(" ===== ========================================")
    1102 format(" ----- ----------------------------------------")
    1110 format(2x,a,3x,3(1x,f12.6))

    ! Set output
    if (present(iunit_)) then
        iu = iunit_
    else
        iu = iu_out
    end if

    if (present(at_mass_)) then
        write(iu, 1001)
        write(iu, 1000)
        write(iu, 1002)
        do ia = 1, n_at
            write(iu, 1010) at_lab(ia), at_mass_(ia), &
                phys%bohr2Ang(at_crd(:,ia))
        end do
        write(iu, 1001)
    else
        write(iu, 1101)
        write(iu, 1100)
        write(iu, 1102)
        do ia = 1, n_at
            write(iu, 1110) at_lab(ia), at_mass_(ia), &
                phys%bohr2Ang(at_crd(:,ia))
        end do
        write(iu, 1101)
    end if

end subroutine prt_coord_r64

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

subroutine sec_header(level, title)
    !! Formats and writes a section header in default unit.
    !!
    !! Formats and prints a header with title `title`.
    !! Available header levels are:
    !!
    !! -1. Main title
    !!  0. Chapter
    !!  1. Section / Header1
    !!  2. Subsection / Header2
    !!  3. Subsubsection / Header3
    !!  4. Paragraph / Header4
    implicit none

    integer, intent(in) :: level
    !! Header level
    character(len=*), intent(in) :: title
    !! Header title

    integer :: lblc, lshft, ltitle
    character(len=256) :: fmt

    1000 format('(1x,"/",',i0,'("-"),"\",/,1x,"|",',i0,'x,"|",/,1x,"|",',i0, &
        'x,a,',i0,'x,"|",/,1x,"|",',i0,'x,"|",/,1x,"\",',i0,'("-"),"/")')
    1010 format('(//,1x,',i0,'("*"),/,1x,',i0,'x,a,/,1x,',i0,'("*"))')
    1020 Format('(//,1x,a,/,1x,',i0,'("="))')
    1030 Format('(/,1x,a,/,1x,',i0,'("-"))')
    1040 Format('(/,1x,a,/,1x,',i0,'("^"))')
    1050 Format(/,1x,'### ',a,' ###')

    ltitle = len_trim(title)

    select case(level)
    case(-1)
        if (ltitle > 68) then ! include at least 4 spaces before/after border
            lblc = ltitle + 8
            lshft = 4
        else
            lblc = 76
            lshft = (lblc-ltitle)/2
        endif
        write(fmt, 1000) lblc, lblc, lshft, lblc - ltitle - lshft, lblc, lblc
        write(iu_out, fmt) trim(title)
    case(0)
        if (ltitle > 78) then
            lblc = ltitle + 4
            lshft = 2
        else
            lblc = 78
            lshft = (lblc-ltitle)/2
        endif
        write(fmt, 1010) lblc, lshft, lblc - ltitle - lshft, lblc
        write(iu_out, fmt) trim(title)
    case(1)
        write(fmt, 1020) ltitle
        write(iu_out, fmt) trim(title)
    case(2)
        write(fmt, 1030) ltitle
        write(iu_out, fmt) trim(title)
    case(3)
        write(fmt, 1040) ltitle
        write(iu_out, fmt) trim(title)
    case(4:)
        write(iu_out, 1050) trim(title)
    case default
        print *, 'Unknown header level.  Stopping.'
        stop

    end select

end subroutine sec_header

! ======================================================================

subroutine write_err(nature, msg, error_)
    !! Writes error message on default unit.
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

end module output
