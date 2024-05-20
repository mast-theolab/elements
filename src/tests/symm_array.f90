program test_symm_array
    use arrays, only: symm_tri_array
    use output, only: sec_header, iu_out, prt_mat

    integer :: ilin, ilow, isym, N, NSq, NTT
    real, dimension(:), allocatable, target :: A_LT, A_2D
    real, dimension(:), pointer :: A
    real, dimension(:,:), allocatable :: B
    logical, dimension(2) :: truefalse = [.True., .False.], falsetrue = [.False., .True.]
    character(len=80) :: arg
    character(len=*), parameter :: outfile = 'test_symm_array.txt'
    character(len=*), dimension(2), parameter :: low_yn = ['lower', 'upper'], &
        sym_yn = ['symmetric    ', 'antisymmetric'], &
        lin_yn = ['LT vector ' , 'full array']
    character(len=80) :: sectitle

    open(file=outfile, newunit=iu_out, action='write')
    print '(" Output written in ",a)', outfile
    call sec_header(1, 'Init')
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        read(arg, *) N
        write(iu_out, '(" Reading number of columns/rows: ",i0)') N
    else
        N = 3
        write(iu_out, '(" Using the default value: ",i0)') N
    end if
    NTT = N*(N+1)/2
    NSq = N*N

    write(iu_out, '(1x,a)') "Allocating arrays"
    allocate(A_LT(NTT), A_2D(NSq), B(N,N))

    call sec_header(1, 'Test with suite of integers')

    ! == Data stored in different array: A -> B
    do ilin = 1, 2
        write(sectitle, '("A ",a," - data in new array: B")') &
            trim(lin_yn(ilin))
        call sec_header(2, trim(sectitle))
        if (ilin == 1) then
            call fill_A('int', NTT, A_LT)
            call sec_header(4, 'Initial array A')
            call prt_mat(reshape(A_LT, [NTT, 1]), NTT, 1)
            A => A_LT
        else
            call fill_A('int', NSq, A_2D)
            call sec_header(4, 'Initial array A')
            call prt_mat(reshape(A_2D, [N, N]), N, N)
            A => A_2D
        end if
        do ilow = 1, 2
            do isym = 1, 2
                write(sectitle, '("A ",a," triangular -> B ",a)') &
                    trim(low_yn(ilow)), trim(sym_yn(isym))
                call sec_header(3, trim(sectitle))
                B = 0.0
                call symm_tri_array(N, A, linear=truefalse(ilin), &
                                    lower=truefalse(ilow), &
                                    anti_symm=falsetrue(isym), arr_new=B)
                call prt_mat(B, N, N)
            end do
        end do
    end do

    ! == Data stored in same array: A -> A
    do ilin = 1, 2
        write(sectitle, '("A ",a," - output in A")') &
            trim(lin_yn(ilin))
        call sec_header(2, trim(sectitle))
        if (ilin == 1) then
            call fill_A('int', NTT, A_2D)
            call sec_header(4, 'Initial array A')
            call prt_mat(reshape(A_2D, [NTT, 1]), NTT, 1)
        else
            call fill_A('int', NSq, A_2D)
            call sec_header(4, 'Initial array A')
            call prt_mat(reshape(A_2D, [N, N]), N, N)
        end if
        do ilow = 1, 2
            do isym = 1, 2
                write(sectitle, '("A ",a," triangular -> A ",a)') &
                    trim(low_yn(ilow)), trim(sym_yn(isym))
                call sec_header(3, trim(sectitle))
                B = 0.0
                call symm_tri_array(N, A_2D, linear=truefalse(ilin), &
                                    lower=truefalse(ilow), &
                                    anti_symm=falsetrue(isym))
                call prt_mat(reshape(A_2D, [N, N]), N, N)
                ! Reset A
                if (ilin == 1) then
                    call fill_A('int', NTT, A_2D)
                else
                    call fill_A('int', NSq, A_2D)
                end if
            end do
        end do
    end do



! ======================================================================

contains

subroutine fill_A(opmode, N, A)
    character(len=*), intent(in) :: opmode
    integer, intent(in) :: N
    real, dimension(:), intent(out) :: A

    integer :: i, ij, j

    A = 0.0
    select case(opmode)
    case('int')
        do i = 1, N
            A(i) = real(i)
        end do
    case default
        print *, 'Unsupported opmode'
        stop
    end select
end subroutine fill_A
end program test_symm_array