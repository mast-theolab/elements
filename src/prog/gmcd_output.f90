module gmcd_output
    !! Output-related module for GenTensor
    !!
    !! Different output/printing constants and procedures for GenTensor
    use string, only: locase
    use exception, only: BaseException, ArgumentError, InitError, RaiseError, &
        RaiseArgError

    implicit none

contains

! ======================================================================

function shell_lmxyz(shelltype, spherical, err) result(comps)
    !! Return the types of function for a given shell
    !!
    !! Returns a list of all types of function for a given shell in the
    !!   formats: l-m (pure), xyz (cartesian)
    character(len=*), intent(in) :: shelltype
    !! Type of shell
    logical, intent(in) :: spherical
    !! True for pure (spherical harmonic functions), False for Cart.
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance
    character(len=:), dimension(:), allocatable :: comps

    integer :: i, ix, iy, iz, j, L
    character(len=1) :: lang
    character(len=60) :: fmt
    character(len=2), parameter :: pm_sign = '+-'

    err = InitError()

    select case (shelltype)
        case ('S')
            allocate(character(len=1):: comps(1))
            comps(1) = 's'
        case ('SP')
            allocate(character(len=2):: comps(4))
            comps(1) = 's '
            comps(2) = 'px'
            comps(3) = 'py'
            comps(4) = 'pz'
        case ('P')
            allocate(character(len=2):: comps(3))
            comps(1) = 'px'
            comps(2) = 'py'
            comps(3) = 'pz'
        case ('D')
            if (spherical) then
                allocate(character(len=3):: comps(5))
                comps(1) = 'd0'
                comps(2) = 'd1a'
                comps(3) = 'd1b'
                comps(4) = 'd2b'
                comps(5) = 'd2a'
            else
                allocate(character(len=3):: comps(6))
                comps(1) = 'dxx'
                comps(2) = 'dyy'
                comps(3) = 'dzz'
                comps(4) = 'dxy'
                comps(5) = 'dxz'
                comps(6) = 'dyz'
            end if
        case ('F')
            if (spherical) then
                allocate(character(len=3):: comps(7))
                comps(1) = 'f 0'
                comps(2) = 'f+1'
                comps(3) = 'f-1'
                comps(4) = 'f+2'
                comps(5) = 'f-2'
                comps(6) = 'f+3'
                comps(7) = 'f-3'
            else
                allocate(character(len=4):: comps(10))
                comps( 1) = 'fxxx'
                comps( 2) = 'fyyy'
                comps( 3) = 'fzzz'
                comps( 4) = 'fxxy'
                comps( 5) = 'fxxz'
                comps( 6) = 'fyyx'
                comps( 7) = 'fyyz'
                comps( 8) = 'fzzx'
                comps( 9) = 'fzzy'
                comps(10) = 'fxyz'
            end if
        case ('G', 'H', 'I')
            if (shelltype == 'G') then
                L = 4
            else if (shelltype == 'H') then
                L = 5
            else if (shelltype == 'I') then
                L = 6
            else
                L = -1
            end if
            if (spherical) then
                lang = locase(shelltype(1:1))
                allocate(character(len=3):: comps(2*L+1))
                comps(1) = lang // ' 0'
                fmt = '(a1,a1,i1)'
                do i = 1, L
                    do j = 1, 2
                        write(comps(i*2+j-1), fmt) lang, pm_sign(j:j), i
                    end do
                end do
            else
                allocate(character(len=L):: comps((L+1)*(L+2)/2))
                i = 0
                do ix = 0, L
                    do iy = 0, L-ix
                        iz = L - ix - iy
                        i = i + 1
                        write(comps(i),  '(a,a,a)') repeat('x', ix), &
                            repeat('y', iy), repeat('z', iz)
                    end do
                end do
            end if
        case default
            call RaiseArgError(err, 'shelltype', 'Unsupported shell type')
            return
    end select

end function shell_lmxyz

! ======================================================================

end module gmcd_output
