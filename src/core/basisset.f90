module basisset
    !! Basis sets-related modules
    !!
    !! Procedure related to basis sets and their definition
    use iso_fortran_env, only: real64
    use math, only: build_PascalTriangle, fac => factorial, pi, int_xn_e2ax2, &
        itri_pa, phii_xn_phij
    use exception, only: BaseException, ArgumentError, InitError, RaiseError, &
        RaiseArgError
    use output, only: iu_out, len_int

    implicit none

    integer, parameter :: max_nxyz = 28
    type, public :: PrimitiveFunction
        !! primitive basis function.
        !! alpha and coeffs are coeffs for hybrid functions (ex: SP)
        character(len=2) :: shelltype
        !! Shell type (of the contracted shell)
        integer ::  L
        !! principal angular momentum
        ! character(len=1) :: functype
        ! !! type of the function: S, P...
        integer :: shellid
        !! main shell id
        integer :: ndim
        !! number of components (dimension)
        logical :: pure
        !! true if spherical harmonics, Cartesian otherwise
        logical :: shell_first
        !! true if first primitive of shell
        logical :: shell_last
        !! true if last primitive of shell
        real(real64) :: alpha
        !! alpha: primitive exponent
        real(real64), dimension(:), allocatable :: coeff
        !! Unique normalized contraction coefficients
        !! indexes: 0, -1... : non-normalized original coeffs.
        integer, dimension(:,:), allocatable :: lxyz
        !! Indexes for each unique contraction coefficient
    end type PrimitiveFunction

contains

! ======================================================================



! ======================================================================

subroutine build_bset_DB(n_at, n_shells, pureD, pureF, shell_types, &
                         prim_per_sh, shell_to_at, coef_contr, coef_contrSP, &
                         prim_exp, nprim_per_at, bset, err)
    !! Build basis set database
    !!
    !! Builds a list of the basis functions sorted by atoms and
    !!   primitives.
    integer, intent(in) :: n_at
    !! Number of atoms
    integer, intent(in) :: n_shells
    !! Number of contracted shells
    logical, intent(in) :: pureD
    !! If True, spherical harm. are used for D orbitals, Cart. otherwise
    logical, intent(in) :: pureF
    !! If True, spherical harmonics are used for F and above
    integer, dimension(:), intent(in) :: shell_types
    !! Shell types, as integers
    integer, dimension(:), intent(in) :: prim_per_sh
    !! number of primitives per shell
    integer, dimension(:), intent(in) :: shell_to_at
    !! shell to atom mapping
    real(real64), dimension(:), intent(in) :: coef_contr
    !! contraction coefficients
    real(real64), dimension(:), allocatable, intent(in) :: coef_contrSP
    !! contraction coefficients for shells of type S=P
    real(real64), dimension(:), intent(in) :: prim_exp
    !! primitives exponents
    integer, dimension(:), allocatable, intent(out) :: nprim_per_at
    !! Number of primitive per atom
    type(PrimitiveFunction), dimension(:,:), allocatable, intent(out) :: bset
    !! Basis set DB
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance

    integer :: i, ia, ippa, iprim, iptot, L_ang, ndim, ntot_AOs
    integer :: ish
    integer, dimension(:), allocatable :: shell_per_at
    real(real64) :: c1, c2
    logical :: purefunc
    character(len=2) :: shtype
    class(BaseException), allocatable :: suberr

    err = InitError()

    ! first find the maximum number of primitives / atom
    allocate(nprim_per_at(n_at), shell_per_at(n_at))
    nprim_per_at = 0
    do i = 1, size(shell_to_at)
        ia = shell_to_at(i)
        nprim_per_at(ia) = nprim_per_at(ia) + prim_per_sh(i)
    end do
    allocate(bset(n_at,maxval(nprim_per_at)))
    nprim_per_at = 0
    shell_per_at = 0

    ntot_AOs = 0
    iptot = 0
    do ish = 1, n_shells
        select case(shell_types(ish))
            ! ndim for spherical: 2*L + 1
            ! ndim for cartesian: (L+1)*(L+2)/2
            case (-1)
                shtype = 'SP'
                L_ang = -1
                ndim = 4
                purefunc = .False.
            case (0)
                shtype = 'S'
                L_ang = 0
                ndim = 1
                purefunc = .False.
            case (1)
                shtype = 'P'
                L_ang = 1
                ndim = 3
                purefunc = .False.
            case (-2,2)
                shtype = 'D'
                L_ang = 2
                if (pureD) then
                    ndim = 5
                    purefunc = .True.
                else
                    ndim = 6
                    purefunc = .False.
                end if
            case (-3,3)
                shtype = 'F'
                L_ang = 3
                if (pureF) then
                    ndim = 7
                    purefunc = .True.
                else
                    ndim = 10
                    purefunc = .False.
                end if
            case (-4,4)
                shtype = 'G'
                L_ang = 4
                if (pureF) then
                    ndim = 9
                    purefunc = .True.
                else
                    ndim = 15
                    purefunc = .False.
                end if
            case (-5,5)
                shtype = 'H'
                L_ang = 5
                if (pureF) then
                    ndim = 11
                    purefunc = .True.
                else
                    ndim = 21
                    purefunc = .False.
                end if
            case (-6,6)
                shtype = 'I'
                L_ang = 6
                if (pureF) then
                    ndim = 13
                    purefunc = .True.
                else
                    ndim = 28
                    purefunc = .False.
                end if
            case default
                call RaiseError(err, 'Unrecognized shell type')
                return
        end select
        ntot_AOs = ntot_AOs + ndim
        ia = shell_to_at(ish)
        shell_per_at(ia) = shell_per_at(ia) + 1
        do iprim = 1, prim_per_sh(ish)
            iptot = iptot + 1
            nprim_per_at(ia) = nprim_per_at(ia) + 1
            ippa = nprim_per_at(ia)
            bset(ia,ippa)%shelltype = shtype
            bset(ia,ippa)%L = L_ang
            bset(ia,ippa)%shellid = shell_per_at(ia)
            bset(ia,ippa)%ndim = ndim
            bset(ia,ippa)%pure = purefunc
            bset(ia,ippa)%shell_first = iprim == 1
            bset(ia,ippa)%shell_last = iprim == prim_per_sh(ish)
            bset(ia,ippa)%alpha = prim_exp(iptot)
            c1 = coef_contr(iptot)
            if (allocated(coef_contrSP)) then
                c2 = coef_contrSP(iptot)
            else
                c2 = 0.0_real64
            end if
            call coefs_norm_sh(shtype, c1, bset(ia,ippa)%alpha, c2, &
                               bset(ia,ippa)%coeff, bset(ia,ippa)%lxyz, suberr)
            if (suberr%raised()) then
                select type(suberr)
                    class is (ArgumentError)
                        call RaiseError(err, 'Unrecognized shell type')
                        return
                    class default
                        call RaiseError(err, 'Generic error')
                        return
                end select
            end if
        end do
    end do

end subroutine build_bset_DB

! ======================================================================

subroutine coefs_norm_sh(shtype, coef, alpha, coef2, new_coefs, indexes, err)
    !! Calculate the normalized coefficients for a given primitive shell
    !!
    !! Computes and returned the unique normalized coefficients relevant
    !!   for a given primitive shell
    character(len=*), intent(in) :: shtype
    !! Shell type
    real(real64), intent(in) :: coef
    !! Contracted coefficient
    real(real64), intent(in) :: alpha
    !! Primitive exponent alpha
    real(real64), intent(in) :: coef2
    !! Secondary contracted coefficient, for instance for S=P shell
    real(real64), dimension(:), allocatable, intent(out) :: new_coefs
    !! New, normalized coefficients
    integer, dimension(:,:), allocatable, intent(out) :: indexes
    !! Indexes corresponding to each coefficients
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance

    real(real64), parameter :: f2 = 2.0_real64, sq2pi = sqrt(f2*pi)
    integer :: i, ix, iy, iz, L
    real(real64) :: cnorm, f2sqal

    err = InitError()

    f2sqal = f2*sqrt(alpha)
    select case(shtype(1:1))
        case ('S')
            if (shtype(2:2) == 'P') then
                allocate(new_coefs(-1:2))
                allocate(indexes(3,2))
                ! it is faster to compute the result than call the integral
                ! but this is equivalent to:
                ! sqrt(int_xn_e2ax2(2*1, alpha)*int_xn_e2ax2(0, alpha)**2)
                ! 2*1 since only 1 coordinate (L=1)
                cnorm = sqrt(f2sqal**5/(sq2pi**3))
                new_coefs(2) = coef2*cnorm
                new_coefs(-1) = coef2
                indexes(:,2) = [0, 0, 1]
            else
                allocate(new_coefs(0:1))
                allocate(indexes(3,1))
            end if
            ! Equivalent to:
            ! sqrt(int_xn_e2ax2(0, alpha)**3)  (L=0 for S)
            cnorm = sqrt((f2sqal/sq2pi)**3)
            new_coefs(0) = coef
            new_coefs(1) = coef*cnorm
            indexes(:,1) = [0, 0, 0]
        case ('P')
            allocate(new_coefs(0:1))
            allocate(indexes(3,1))
            ! sqrt(int_xn_e2ax2(2*1, alpha)*int_xn_e2ax2(0, alpha)**2)
            ! 2*1 since only 1 coordinate (L=1)
            cnorm = sqrt(f2sqal**5/(sq2pi**3))
            new_coefs(0) = coef
            new_coefs(1) = coef*cnorm
            indexes(:,1) = [0, 0, 1]
        case ('D')
            allocate(new_coefs(0:2))
            allocate(indexes(3,2))
            ! L = 2, 2 coordinates -> 2 cases: ix=2 / ix=1,iy=1
            ! 2 cases:
            ! - sqrt(int_xn_e2ax2(2*2, alpha)*int_xn_e2ax2(0, alpha)**2)
            ! - sqrt(int_xn_e2ax2(2*1, alpha)**2*int_xn_e2ax2(0, alpha))
            cnorm = sqrt(f2sqal**7/(sq2pi**3))
            new_coefs(0) = coef
            new_coefs(2) = coef*cnorm
            new_coefs(1) = new_coefs(2)/sqrt(3.0_real64)
            indexes(:,2) = [0, 1, 1]
            indexes(:,1) = [0, 0, 2]
        case ('F')
            allocate(new_coefs(0:3))
            allocate(indexes(3,3))
            ! L = 3, 3 coordinates -> 3 cases: ix=3/ix=2,iy=1/ix=1,iy=1,iz=1
            ! 3 cases:
            ! - sqrt(int_xn_e2ax2(2*3, alpha)*int_xn_e2ax2(0, alpha)**2)
            ! - sqrt(int_xn_e2ax2(2*2, alpha)*int_xn_e2ax2(1*2, alpha)
            !     *int_xn_e2ax2(0, alpha))
            ! - sqrt(int_xn_e2ax2(2*1, alpha)**3
            cnorm = sqrt(f2sqal**7/(sq2pi**3))
            new_coefs(0) = coef
            new_coefs(3) = coef*cnorm
            new_coefs(2) = new_coefs(1)/sqrt(3.0_real64)
            new_coefs(1) = new_coefs(2)/sqrt(15.0_real64)
            indexes(:,3) = [1, 1, 1]
            indexes(:,2) = [0, 1, 2]
            indexes(:,1) = [0, 0, 3]
        case ('G')
            allocate(new_coefs(0:15))
            allocate(indexes(3,15))
            new_coefs(0) = coef
            L = 4
            i = 0
            do ix = 0, L
                do iy = 0, L-ix
                    do iz = 0, L-ix-iy
                        if (ix+iy+iz == L) then
                            i = i + 1
                            new_coefs(i) = coef/sqrt(&
                                int_xn_e2ax2(2*ix, alpha) &
                                *int_xn_e2ax2(2*iy, alpha) &
                                *int_xn_e2ax2(2*iz, alpha) &
                                )
                            indexes(:, i) = [ix, iy, iz]
                        end if
                    end do
                end do
            end do
        case ('H')
            allocate(new_coefs(0:21))
            allocate(indexes(3,21))
            new_coefs(0) = coef
            L = 5
            i = 0
            do ix = 0, L
                do iy = 0, L-ix
                    do iz = 0, L-ix-iy
                        if (ix+iy+iz == L) then
                            i = i + 1
                            new_coefs(i) = coef/sqrt(&
                                int_xn_e2ax2(2*ix, alpha) &
                                *int_xn_e2ax2(2*iy, alpha) &
                                *int_xn_e2ax2(2*iz, alpha) &
                                )
                            indexes(:, i) = [ix, iy, iz]
                        end if
                    end do
                end do
            end do
        case ('I')
            allocate(new_coefs(0:28))
            allocate(indexes(3,28))
            new_coefs(0) = coef
            L = 6
            i = 0
            do ix = 0, L
                do iy = 0, L-ix
                    do iz = 0, L-ix-iy
                        if (ix+iy+iz == L) then
                            i = i + 1
                            new_coefs(i) = coef/sqrt(&
                                int_xn_e2ax2(2*ix, alpha) &
                                *int_xn_e2ax2(2*iy, alpha) &
                                *int_xn_e2ax2(2*iz, alpha) &
                                )
                            indexes(:, i) = [ix, iy, iz]
                        end if
                    end do
                end do
            end do
        case default
            call RaiseArgError(err, 'shtype', 'Unsupported shell type')
            return
    end select

end subroutine coefs_norm_sh

! ======================================================================

subroutine coef_transfo_P2C(L_ang, Ncart, Npure, coef2P, coef2C)
    !! Generate coefficients for conversion between pure and Cart. bases
    !!
    !! Generates the coefficients to transform between normalized Cart.
    !!   Gaussians and normalized spherical harmonics (pure).
    integer, intent(in) :: L_ang
    !! Angular momentum to transform
    integer, intent(in) :: Ncart
    !! Number of Cartesian components
    integer, intent(in) :: Npure
    !! Number of pure/spherical components
    real(real64), dimension(Ncart,Npure), intent(out) :: coef2P
    !! Transformation coefficients to pure
    real(real64), dimension(Ncart,Npure), intent(out) :: coef2C
    !! Transformation coefficients to Cartesian

    integer :: i1, i2, k, Lx, Lx1, Lx2, Lxyz, Ly, Ly1, Ly2, Lz, Lz1, &
        Lz2, M
    real(real64) :: a1, a2, s
    ! Form coef2P
    Lxyz = 0
    do Lx = 0, L_ang
        do Ly = 0, (L_ang-Lx)
            Lxyz = Lxyz + 1
            coef2P(Lxyz,1) = coef_C2P(L_ang, 0, Lx, Ly)
            do M = 1, L_ang
                coef2p(Lxyz,2*M) = coef_C2P(L_ang, M, Lx, Ly)
                coef2P(Lxyz,2*M+1) = coef_C2P(L_ang, -M ,Lx, Ly)
            end do
        end do
    end do

    ! Form coef2C
    coef2C = 0.0_real64
    i1 = 0
    do Lx1 = 0, L_ang
        do Ly1 = 0, (L_ang-Lx1)
            i1 = i1 + 1
            Lz1 = L_ang - Lx1 - Ly1
            a1 = sqrt(real(fac(Lx1)*fac(Ly1)*fac(Lz1), kind=real64) &
                      /(fac(2*Lx1)*fac(2*Ly1)*fac(2*Lz1)))
            i2 = 0
            do Lx2 = 0, L_ang
                do Ly2 = 0, (L_ang-Lx2)
                    i2 = i2 + 1
                    Lz2 = L_ang - Lx2 - Ly2
                    a2 = sqrt(real(fac(Lx2)*fac(Ly2)*fac(Lz2), kind=real64) &
                              /(fac(2*Lx2)*fac(2*Ly2)*fac(2*Lz2)))
                    Lx = Lx1 + Lx2
                    Ly = Ly1 + Ly2
                    Lz = Lz1 + Lz2
                    if(mod(Lx,2) == 0 .and. mod(Ly,2) == 0 &
                       .and. mod(Lz,2) == 0) then
                        s = a1*a2*fac(Lx)*fac(Ly)*fac(Lz) &
                            / (fac(Lx/2)*fac(Ly/2)*fac(Lz/2))
                        do k = 1, Npure
                            coef2C(i1,k) = coef2C(i1,k) + s*coef2P(i2,k)
                        end do
                    end if
                end do
            end do
        end do
    end do
end subroutine coef_transfo_P2C

! ======================================================================

function coef_C2P(L, M, Lx, Ly) result(coef)
    !! Calculate the coefficient from Cartesian to pure
    !!
    !! Calculates the coefficient for the conversion from the normalized
    !!   Cartesian component Lx,Ly,Lz=L-Lx-Ly to the spherical harmonics
    !!   component L,M
    real(real64) :: coef
    integer, intent(in) :: L
    !! Angular moment
    integer, intent(in) :: M
    !! M component of the spherical harmonics
    integer, intent(in) :: Lx
    !! Lx component of the Cartesian function
    integer, intent(in) :: Ly
    !! Ly component of the Cartesian function

    integer :: i, j, k, Lx2k, Lz, Ma
    real(real64) :: a, b, c, d, e, g

    Lz = L - Lx - Ly
    Ma = Abs(M)
    if (Ma > L .or. Lz < 0) then
        coef = 0.0_real64
        return
    end if
    a = sqrt(real(fac(2*Lx)*fac(2*Ly)*fac(2*Lz)*fac(L), kind=real64) &
             / (fac(2*L)*fac(Lx)*fac(Ly)*fac(Lz)))
    b = sqrt(real(fac(L-Ma), kind=real64)/fac(L+Ma))/((2**L)*fac(L))
    j = (L-Ma-Lz)/2
    g = 0.0_real64
    if (2*j == L-Ma-Lz) then
        do i = 0, (L-Ma)/2
            c = 0.0_real64
            if (j >= 0 .and. j <= i) then
                c = (real(fac(L), kind=real64)/(fac(i)*fac(L-i))) &
                    *(real(fac(2*L-2*i)*((-1)**i), kind=real64)/fac(L-Ma-2*i)) &
                    *(real(fac(i), kind=real64)/(fac(j)*fac(i-j)))
                do k = 0, j
                    d = 0.0_real64
                    Lx2k = Lx-2*k
                    if (Lx2k >= 0 .and. Lx2k <= Ma) then
                        d = (real(fac(j), kind=real64)/(fac(k)*fac(j-k))) &
                            *(real(fac(Ma), kind=real64) &
                                /(fac(Lx2k)*fac(Ma-Lx2k)))
                        e = 0.0_real64
                        if (M==0 .and. mod(Lx,2) == 0) e = (-1)**(-Lx2k/2)
                        if (M > 0 .and. mod(abs(Ma-Lx),2) == 0) &
                            e = Sqrt(2.0_real64)*(-1)**((Ma-Lx2k)/2)
                        if (M < 0 .and. mod(abs(Ma-Lx),2) == 1) &
                            e = Sqrt(2.0_real64)*(-1)**((Ma-Lx2k)/2)
                        g = g + c*d*e
                    end if
                end do
            end if
        end do
    end if
    coef = a*b*g

    return
end function coef_C2P

! ======================================================================

subroutine fix_norm_AOs(iout, n_at, n_ao, at_crd, nprim_per_at, bsetDB, err, &
                        debug)
    !! Correct basis sets coefficients to normalize the AOs.
    !!
    !! Checks if atomic orbitals are normalized and otherwise correct
    !!   the coefficients to ensure the normalization.
    integer, intent(in) :: iout
    !! Unit for output
    integer, intent(in) :: n_at
    !! Number of atoms
    integer, intent(in) :: n_ao
    !! Number of atomic orbitals
    real(real64), dimension(:,:), intent(in) :: at_crd
    !! Atomic coordinates (in au)
    integer, dimension(:), intent(in) :: nprim_per_at
    !! Number of basis primitives per atom
    type(PrimitiveFunction), dimension(:,:), intent(inout) :: bsetDB
    !! Basis set DB
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance
    logical, intent(in), optional :: debug
    !! Enable debugging printing

    integer, parameter :: &
        dimPa = 24  ! maximum used dimension for Pascal's triangle
    real(real64), parameter :: &
        tol_expm = 100.0_real64, &  ! Max. value of x for e^-x to be relevant
        sqpi3 = sqrt(pi**3)
    integer :: ndi, ndj
    integer :: i, idi, idj, iprim, j, jj, jprim
    integer :: ia, ia0, ja, ja0
    integer, dimension(:,:), allocatable :: ldi, ldj
    real(real64) :: a, ai, aj, cijer, cjer, ebase, r2ij, x
    real(real64), dimension(3) :: ri, rj, ovi
    real(real64), dimension(max_nxyz,max_nxyz) :: ovlp_ij
    real(real64), dimension(max_nxyz,max_nxyz), target :: allones
    real(real64), dimension(:), allocatable :: ao_norms, ci, cj
    real(real64), dimension(:,:), allocatable :: ao_ovlp
    real(real64), dimension(:,:), allocatable, target :: c2p_D, c2p_F, c2p_G, &
        c2p_H, c2p_I
    real(real64), dimension(:,:), pointer :: c2pi, c2pj
    logical :: dbg_print = .false.
    class(BaseException), allocatable :: suberr

    err = InitError()

    if (present(debug)) dbg_print = debug

    if (.not.allocated(itri_pa)) then
        call build_PascalTriangle(dimPa)
    else
        if (size(itri_pa, 1) < 24) then
            deallocate(itri_pa)
            call build_PascalTriangle(dimPa)
        end if
    end if
    ! build an array of ones as reference for the conversion Cart -> Pure if
    !   primitive already in Cartesian.

    if (dbg_print) then
        1000 format(1x,i0,' orbitals read on ',i0,' atoms')
        write(iout, 1000) n_ao, n_at
    end if

    allocate(ao_norms(n_ao), ao_ovlp(max_nxyz,n_ao))
    ao_norms = 0.0_real64
    ao_ovlp = 0.0_real64

    ia0 = 0
    do ia = 1, n_at
        ri = at_crd(:,ia)
        do iprim = 1, nprim_per_at(ia)
            ai = bsetDB(ia,iprim)%alpha
            call set_primC_comp(bsetDB(ia,iprim), ndi, ldi, ci, suberr)
            if (suberr%raised()) then
                select type(suberr)
                    class is (ArgumentError)
                        call RaiseError(err, &
                                        'Unrecognized shell type for iprim')
                        return
                    class default
                        call RaiseError(err, 'Generic error')
                        return
                end select
            end if
            if (bsetDB(ia,iprim)%pure) then
                select case (bsetDB(ia,iprim)%shelltype)
                    case ('D')
                        if (.not.allocated(c2p_D)) &
                            c2p_D = transfo_cart2pure(bsetDB(ia,iprim)%L)
                        c2pi => c2p_D
                    case ('F')
                        if (.not.allocated(c2p_F)) &
                            c2p_F = transfo_cart2pure(bsetDB(ia,iprim)%L)
                        c2pi => c2p_F
                    case ('G')
                        if (.not.allocated(c2p_G)) &
                            c2p_G = transfo_cart2pure(bsetDB(ia,iprim)%L)
                        c2pi => c2p_G
                    case ('H')
                        if (.not.allocated(c2p_H)) &
                            c2p_H = transfo_cart2pure(bsetDB(ia,iprim)%L)
                        c2pi => c2p_H
                    case ('I')
                        if (.not.allocated(c2p_I)) &
                            c2p_I = transfo_cart2pure(bsetDB(ia,iprim)%L)
                        c2pi => c2p_I
                    case default
                        c2pi => allones
                end select
            else
                c2pi => allones
            end if
            ja0 = 0
            do ja = 1, n_at
                rj = at_crd(:,ja)
                r2ij = sum((rj-ri)**2)
                do jprim = 1, nprim_per_at(ja)
                    aj = bsetDB(ja,jprim)%alpha
                    call set_primC_comp(bsetDB(ja,jprim), ndj, ldj, cj, suberr)
                    if (suberr%raised()) then
                        select type(suberr)
                            class is (ArgumentError)
                                call RaiseError(&
                                    err, 'Unrecognized shell type for jprim')
                                return
                            class default
                                call RaiseError(err, 'Generic error')
                                return
                        end select
                    end if
                    if (bsetDB(ja,jprim)%pure) then
                        select case (bsetDB(ja,jprim)%shelltype)
                            case ('D')
                                if (.not.allocated(c2p_D)) &
                                    c2p_D = transfo_cart2pure(bsetDB(ja,jprim)%L)
                                c2pj => c2p_D
                            case ('F')
                                if (.not.allocated(c2p_F)) &
                                    c2p_F = transfo_cart2pure(bsetDB(ja,jprim)%L)
                                c2pj => c2p_F
                            case ('G')
                                if (.not.allocated(c2p_G)) &
                                    c2p_G = transfo_cart2pure(bsetDB(ja,jprim)%L)
                                c2pj => c2p_G
                            case ('H')
                                if (.not.allocated(c2p_H)) &
                                    c2p_H = transfo_cart2pure(bsetDB(ja,jprim)%L)
                                c2pj => c2p_H
                            case ('I')
                                if (.not.allocated(c2p_I)) &
                                    c2p_I = transfo_cart2pure(bsetDB(ja,jprim)%L)
                                c2pj => c2p_I
                            case default
                                c2pj => allones
                        end select
                    else
                        c2pj => allones
                    end if
                    a = 1.0_real64/(ai+aj)
                    x = ai*aj*r2ij*a
                    if (x < tol_expm) then
                        ebase = sqpi3*exp(-x)
                    else
                        if (bsetDB(ja,jprim)%shell_last) &
                            ja0 = ja0 + bsetDB(ja,jprim)%ndim
                        cycle
                    end if

                    ! Now loop over orbitals to which primitives contribute
                    do idj = 1, ndj
                        cjer = ebase*cj(idj)
                        do idi = 1, ndi
                            cijer = cjer*ci(idi)
                            ovi = phii_xn_phij(.true., ldi(:,idi), ri, ai, &
                                               ldj(:,idj), rj, aj, 0)
                            ovlp_ij(idi,idj) = product(ovi)*cijer
                        end do
                    end do

                    ! Check now if we need to convert to pure basis
                    if (bsetDB(ia,iprim)%pure .or. bsetDB(ja,jprim)%pure) then
                        do i = 1, bsetDB(ia,iprim)%ndim
                            do j = 1, bsetDB(ja,jprim)%ndim
                                jj = ja0 + j
                                do idi = 1, ndi
                                    x = c2pi(i,idi)
                                    do idj = 1, ndj
                                        ao_ovlp(i,jj) = ao_ovlp(i,jj) &
                                            + x*ovlp_ij(idi,idj)*c2pj(j,idj)
                                    end do
                                end do
                            end do
                        end do
                    else
                        do idi = 1, ndi
                            do idj = 1, ndj
                                ao_ovlp(idi,ja0+idj) = ao_ovlp(idi,ja0+idj) &
                                    + ovlp_ij(idi,idj)
                            end do
                        end do
                    end if
                    if (bsetDB(ja,jprim)%shell_last) &
                        ja0 = ja0 + bsetDB(ja,jprim)%ndim
                end do
            end do
            if (bsetDB(ia,iprim)%shell_last) then
                do i = 1, bsetDB(ia,iprim)%ndim
                    ao_norms(ia0+i) = ao_ovlp(i,ia0+i)
                end do
                ia0 = ia0 + bsetDB(ia,iprim)%ndim
                ao_ovlp = 0.0_real64
            end if
        end do
    end do

    if (dbg_print) then
        write(iout, '(a)') 'atomic orbital norms before normalization:'
        write(iout, '(6f12.6)') (ao_norms(i),i=1,n_ao)
    end if

    ! Now normalize
    ia0 = 1
    do ia = 1, n_at
        ri = at_crd(:,ia)
        do iprim = 1, nprim_per_at(ia)
            a = 1.0_real64/sqrt(ao_norms(ia0))
            bsetDB(ia,iprim)%coeff = bsetDB(ia,iprim)%coeff*a
            if (bsetDB(ia,iprim)%shell_last) ia0 = ia0 + bsetDB(ia,iprim)%ndim
        end do
    end do

    return
end subroutine fix_norm_AOs


! ======================================================================

subroutine set_primC_comp(bfunc, ndimC, lxyz, coefs, err)
    !! Set component for a Cartesian primitive
    !!
    !! Builds the list of Cartesian components and the associated
    !!   coefficients for a primitive function based on the information
    !!   in `bfunc`.
    type(PrimitiveFunction), intent(in) :: bfunc
    !! Basis set function
    integer, intent(out) :: ndimC
    !! True number of dimension
    integer, dimension(:,:), allocatable, intent(out) :: lxyz
    !! Number of occurrence of each x,y,z coordinate for each dimension
    real(real64), dimension(:), allocatable, intent(out) :: coefs
    !! Contraction coefficient for each dimension
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance

    err = InitError()

    select case (bfunc%shelltype)
        case ('S')
            ndimC = 1
            allocate(lxyz(3,ndimC), coefs(ndimC))
            lxyz = 0
            coefs(1) = bfunc%coeff(1)
        case ('P')
            ndimC = 3
            allocate(lxyz(3,ndimC), coefs(ndimC))
            lxyz(:,1) = [1, 0, 0]
            lxyz(:,2) = [0, 1, 0]
            lxyz(:,3) = [0, 0, 1]
            coefs = bfunc%coeff(1)
        case ('SP')
            ndimC = 4
            allocate(lxyz(3,ndimC), coefs(ndimC))
            lxyz(:,1) = [0, 0, 0]
            lxyz(:,2) = [1, 0, 0]
            lxyz(:,3) = [0, 1, 0]
            lxyz(:,4) = [0, 0, 1]
            coefs(1) = bfunc%coeff(1)
            coefs(2:ndimC) = bfunc%coeff(2)
        case ('D')
            ndimC = 6
            allocate(lxyz(3,ndimC), coefs(ndimC))
            lxyz(:,1) = [2, 0, 0]
            lxyz(:,2) = [0, 2, 0]
            lxyz(:,3) = [0, 0, 2]
            lxyz(:,4) = [1, 1, 0]
            lxyz(:,5) = [1, 0, 1]
            lxyz(:,6) = [0, 1, 1]
            coefs(1:3) = bfunc%coeff(1)
            coefs(4:6) = bfunc%coeff(2)
        case ('F')
            ndimC = 10
            allocate(lxyz(3,ndimC), coefs(ndimC))
            lxyz(:, 1) = [3, 0, 0]
            lxyz(:, 2) = [0, 3, 0]
            lxyz(:, 3) = [0, 0, 3]
            lxyz(:, 4) = [1, 2, 0]
            lxyz(:, 5) = [2, 1, 0]
            lxyz(:, 6) = [2, 0, 1]
            lxyz(:, 7) = [1, 0, 2]
            lxyz(:, 8) = [0, 1, 2]
            lxyz(:, 9) = [0, 2, 1]
            lxyz(:,10) = [1, 1, 1]
            coefs( 1: 3) = bfunc%coeff(1)
            coefs( 4: 9) = bfunc%coeff(2)
            coefs(10:10) = bfunc%coeff(3)
        case ('G')
            ndimC = 15
            allocate(lxyz(3,ndimC), coefs(ndimC))
            lxyz = bfunc%lxyz
            coefs = bfunc%coeff(1:)
        case ('H')
            ndimC = 21
            allocate(lxyz(3,ndimC), coefs(ndimC))
            lxyz = bfunc%lxyz
            coefs = bfunc%coeff(1:)
        case ('I')
            ndimC = 28
            allocate(lxyz(3,ndimC), coefs(ndimC))
            lxyz = bfunc%lxyz
            coefs = bfunc%coeff(1:)
        case default
            call RaiseArgError(err, 'bfunc%shelltype', &
                               'Unsupported shell type')
            return
    end select
end subroutine set_primC_comp

! ======================================================================

function transfo_cart2pure(L_ang) result(convmat)
    !! Provides the transformation matrix from Cartesian to pure
    !!
    !! Provides the conversion matrix to convert from a basis set in
    !!   Cartesian coordinates to a one with spherical harmonics
    !! The conversion table is taken from SOS:
    !!
    !! ## For D orbitals
    !!
    !! Conversion table:
    !!
    !!    |          |    1  |   2  |    3  |   4  |   5  |   6 |
    !!    |          |   xx  |  yy  |   zz  |  xy  |  xz  |  yz |
    !! ---|----------|-------|------|-------|------|------|-----|-----
    !!  1 | 3z^2-R^2 |   -1  |  -1  |    2  |   0  |   0  |   0 | *s1
    !!  2 |    xz    |    0  |   0  |    0  |   0  |   1  |   0 |
    !!  3 |    yz    |    0  |   0  |    0  |   0  |   0  |   1 |
    !!  4 |  x^2-y^2 |    1  |  -1  |    0  |   0  |   0  |   0 | *s5
    !!  5 |    xy    |    0  |   0  |    0  |   1  |   0  |   0 |
    !!
    !! s1 = 1/2 ; s5 = sqrt(3)/2
    !!
    !! ## For F orbitals
    !!
    !!       |     spherical  |  cartesian
    !!  -----|----------------|------------
    !!     1 |  5z^3-3R^2z    |     xxx
    !!     2 |  x(5z^2-R^2)   |     yyy
    !!     3 |  y(5z^2-R^2)   |     zzz
    !!     4 |  xyz           |     xxy
    !!     5 |  z(x^2-y^2)    |     xxz
    !!     6 |  y(y^2-3x^2)   |     yyx
    !!     7 |  x(x^2-3y^2)   |     yyz
    !!     8 |  -             |     zzx
    !!     9 |  -             |     zzy
    !!    10 |  -             |     xyz
    !!
    !! Conversion table
    !!
    !!    |  1  |   2 |    3 |   4 |   5 |   6 |   7 |   8 |   9 |  10 |
    !!    | xxx | yyy |  zzz | xxy | xxz | yyx | yyz | zzx | zzy | xyz |
    !! ---|-----|-----|------|-----|-----|-----|-----|-----|-----|-----|-----
    !!  1 |   0 |   0 | 2*v5 |   0 |  -3 |   0 |  -3 |   0 |   0 |   0 | *s1
    !!  2 | -v5 |   0 |    0 |   0 |   0 |  -1 |   0 |   4 |   0 |   0 | *s2
    !!  3 |   0 | -v5 |    0 |  -1 |   0 |   0 |   0 |   0 |   4 |   0 | *s3
    !!  4 |   0 |   0 |    0 |   0 |   1 |   0 |  -1 |   0 |   0 |   0 | *s4
    !!  5 |   0 |   0 |    0 |   0 |   0 |   0 |   0 |   0 |   0 |   1 | *s5
    !!  6 |  v5 |   0 |    0 |   0 |   0 |  -3 |   0 |   0 |   0 |   0 | *s6
    !!  7 |   0 | -v5 |    0 |   3 |   0 |   0 |   0 |   0 |   0 |   0 | *s7
    !!
    !! @Warning: The function does not check explicitly if the shell is
    !!           supported. Calling procedure should check if the array is
    !!           allocated on exit.
    integer, intent(in) :: L_ang
    !! Angular momentum for the shell of interest
    real(real64), dimension(:,:), allocatable :: convmat
    !! Conversion matrix of dimension (npure,ncart)

    integer :: ncart, npure
    real(real64), parameter :: f0=0.0_real64, f1=1.0_real64, f2=2.0_real64, &
        f3=3.0_real64, f4=4.0_real64, pt5=0.5_real64, sq2=sqrt(2.0_real64), &
        sq3=sqrt(3.0_real64), sq5=sqrt(5.0_real64), sq10=sqrt(10.0_real64)
    real(real64) :: s1, s2, s3, s4, s5, s6, s7
    real(real64), dimension(max_nxyz,max_nxyz) :: tempmat

    select case (L_ang)
        case (0)
            allocate(convmat(1,1))
            convmat = f1
        case (1)
            allocate(convmat(3,3))
            convmat = f0
            convmat(1,1) = f1
            convmat(2,2) = f1
            convmat(3,3) = f1
        case (2)
            allocate(convmat(5,6))
            convmat = f0
            ! Unchanged components
            convmat(2,5) = f1
            convmat(3,6) = f1
            convmat(5,4) = f1
            ! |normalized z^2-r^2> = s1*(two*|norm.zz>-|norm.xx>-|norm.yy>)
            s1 = pt5
            convmat(1,1) = -s1
            convmat(1,2) = -s1
            convmat(1,3) = 2.0_real64*s1
            ! |normalized x^2-y^2> = s5*(|norm.xx>-|norm.yy>)
            s5 = pt5*sq3
            convmat(4,1) = s5
            convmat(4,2) = -s5
        case (3)
            allocate(convmat(7,10))
            convmat = f0
            ! normalization coefficients
            s1 = pt5/sq5
            s2 = pt5*sq3/sq10
            s3 = pt5*sq3/sq10
            s4 = pt5*sq3
            s5 = f1
            s6 = pt5/sq2
            s7 = pt5/sq2
            ! normalized conversion matrix
            convmat(1, 3) =  s1 * f2 * sq5
            convmat(1, 5) = -s1 * f3
            convmat(1, 7) = -s1 * f3
            convmat(2, 1) = -s2 * sq5
            convmat(2, 6) = -s2
            convmat(2, 8) =  s2 * f4
            convmat(3, 2) = -s3 * sq5
            convmat(3, 4) = -s3
            convmat(3, 9) =  s3 * f4
            convmat(4, 5) =  s4
            convmat(4, 7) = -s4
            convmat(5,10) =  s5
            convmat(6, 1) =  s6 * sq5
            convmat(6, 6) = -s6 * f3
            convmat(7, 2) = -s7 * sq5
            convmat(7, 4) =  s7 * f3
        case (4:)
            ncart = (L_ang+1)*(L_ang+2)/2
            npure = 2*L_ang + 1
            allocate(convmat(npure,ncart))
            call coef_transfo_P2C(L_ang, ncart, npure, convmat, tempmat)
    end select

end function transfo_cart2pure

! ======================================================================

subroutine chk_bset_redundancy(n_ao, ovint_i_j, thresh)
    !! Check basis set redundancy
    !!
    !! Checks any redundancy within basis set functions.
    !!
    integer, intent(in) :: n_ao
    !! Number of atomic orbitals
    real(real64), dimension(n_ao,n_ao), intent(in) :: ovint_i_j
    !! Overlap integrals between atomic orbitals, < i | j >
    real(real64), intent(in), optional :: thresh
    !! Threshold to consider redundancy (redundancy if norm below).

    integer :: iao, jao, kao, lao
    real(real64) :: norm, ovlp, thresh0
    real(real64), dimension(n_ao,n_ao) :: trial
    logical, dimension(n_ao) :: is_red
    character(len=60) :: fmt_red

    if (present(thresh)) then
        thresh0 = thresh
        if (thresh0 < epsilon(thresh0)) then
            write(iu_out, '(a)') &
                'WARNING: Threshold for redundancy check too low'
            write(iu_out, '(9x,a,e12.4)') 'Resetting to:', epsilon(thresh0)
            thresh0 = epsilon(thresh0)
        end if
    else
        thresh0 = 1.0e-6_real64
    end if

    write(fmt_red, '("(""Orbital num. "",i",i0,","" - Rest norm:"",f12.8,&
        &:,"" ["",a,""]"")")') len_int(n_ao)

    ! Initialize trial matrix to check overlap
    trial = 0.0_real64
    do iao = 1, n_ao
        trial(iao,iao) = 1.0_real64
    end do
    is_red = .False.

    do iao = 1, n_ao
        do jao = 1, iao-1
            if (.not.is_red(jao)) then
                ovlp = 0.0_real64
                do kao = 1, iao
                    do lao = 1, iao
                        ovlp = ovlp + &
                            trial(iao,kao)*ovint_i_j(kao,lao)*trial(jao,lao)
                    end do
                end do
                trial(iao,:iao) = trial(iao,:iao) - ovlp*trial(jao,:iao)
            end if
        end do
        norm = 0.0
        do kao = 1, iao
            do lao = 1, iao
                norm = norm + trial(iao,kao)*ovint_i_j(kao,lao)*trial(iao,lao)
            end do
        end do
        if (norm > thresh0) then
            write(iu_out, fmt_red) iao, norm
        else
            write(iu_out, fmt_red) iao, norm, 'redundant'
            is_red(iao) = .True.
        end if
        trial(iao,:iao) = trial(iao,:iao)/sqrt(norm)
    end do
            
end subroutine chk_bset_redundancy

! ======================================================================

end module basisset
