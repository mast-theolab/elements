module basisset
    !! Basis sets-related modules
    !!
    !! Procedure related to basis sets and their definition
    use numeric, only: realwp, f0, f1, f2, f3, f4, f5, f10, fhalf, pi
    use math, only: build_PascalTriangle, fac => factorial, int_xn_e2ax2, &
        itri_pa, phii_xn_phij
    use exception, only: BaseException, ArgumentError, InitError, RaiseError, &
        RaiseArgError
    use output, only: iu_out, len_int
    use datatypes, only: BasisSetDB, PrimitiveFunction

    implicit none

    private
    public :: build_bset_DB, chk_bset_redundancy, coef_C2P, coef_transfo_P2C, &
        coefs_norm_sh, set_primC_comp, transfo_cart2pure
    public :: convert_pure2cart, fix_norm_AOs, len_shells_on_atom, &
        num_cart_AOs, num_shells_on_atom
    integer, parameter, public :: max_nxyz = 28

! ----------------------------------------------------------------------

    interface convert_pure2cart
        !! Convert a quantity expressed in pure to Cartesian.
        !!
        !! convert_pure2cart_matrix: convert a matrix
        !! convert_pure2cart_bset: convert basis set information

        module subroutine convert_pure2cart_bsetBF( &
                bsetBF, nprim_per_atom, bsetBF_cart)
            !! Convert basis set information from pure to Cartesian.
            !!
            !! Converts the information for pure basis functions to Cartesian
            !! functions.
            !!
            !! This version expects basis set information as separate arguments.

            type(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
            !! Basis set's basis function information (pure).
            integer, dimension(:), intent(in) :: nprim_per_atom
            !! Number of basis primitives per atom.
            type(PrimitiveFunction), dimension(:,:), intent(out) :: bsetBF_cart
            !! Basis set's basis function information (Cartesian).

        end subroutine convert_pure2cart_bsetBF

        module subroutine convert_pure2cart_bsetDB(bsetDB, bsetBF_cart)
            !! Convert basis set information from pure to Cartesian.
            !!
            !! Converts the information for pure basis functions to Cartesian
            !! functions.
            !!
            !! This version expects a basis set database as argument.

            type(BasisSetDB), intent(in) :: bsetDB
            !! Basis set database (pure).
            type(PrimitiveFunction), dimension(:,:), intent(out) :: bsetBF_cart
            !! Basis set's basis function information (Cartesian).

        end subroutine convert_pure2cart_bsetDB

        module subroutine convert_pure2cart_matrix_bsetBF( &
                bsetBF, nprim_per_atom, matrix_pure, matrix_cart)
            !! Convert matrix from pure to Cartesian basis.
            !!
            !! Converts a matrix from pure to Cartesian basis.
            !!
            !! This version expects basis set information as separate arguments.

            type(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
            !! Basis set's basis function information (pure).
            integer, dimension(:), intent(in) :: nprim_per_atom
            !! Number of basis primitives per atom.
            real(realwp), dimension(:,:), intent(in) :: matrix_pure
            !! Input matrix, in pure basis.
            real(realwp), dimension(:,:), intent(out) :: matrix_cart
            !! Output matrix, in Cartesian basis.

        end subroutine convert_pure2cart_matrix_bsetBF

        module subroutine convert_pure2cart_matrix_bsetDB( &
                bsetDB, matrix_pure, matrix_cart)
            !! Convert matrix from pure to Cartesian basis.
            !!
            !! Converts a matrix from pure to Cartesian basis.
            !!
            !! This version expects basis set information as separate arguments.
            !!
            !! This version expects a basis set database as argument.

            type(BasisSetDB), intent(in) :: bsetDB
            !! Basis set database (pure).
            real(realwp), dimension(:,:), intent(in) :: matrix_pure
            !! Input matrix, in pure basis.
            real(realwp), dimension(:,:), intent(out) :: matrix_cart
            !! Output matrix, in Cartesian basis.

        end subroutine convert_pure2cart_matrix_bsetDB

    end interface convert_pure2cart

! ----------------------------------------------------------------------

    interface fix_norm_AOs
        !! Fix normalization of atomic orbitals
        module subroutine fix_norm_AOs_bsetBF( &
                iout, n_at, n_ao, at_crd, nprim_per_at, bsetBF, err, debug)
            !! Correct basis sets coefficients to normalize the AOs.
            !!
            !! Checks if atomic orbitals are normalized and otherwise correct
            !! the coefficients to ensure the normalization.
            !!
            !! This version expects basis set information as separate arguments.

            integer, intent(in) :: iout
            !! Unit for output.
            integer, intent(in) :: n_at
            !! Number of atoms.
            integer, intent(in) :: n_ao
            !! Number of atomic orbitals.
            real(realwp), dimension(:,:), intent(in) :: at_crd
            !! Atomic coordinates (in au).
            integer, dimension(:), intent(in) :: nprim_per_at
            !! Number of basis primitives per atom.
            type(PrimitiveFunction), dimension(:,:), intent(inout) :: bsetBF
            !! Basis set's basis function information.
            class(BaseException), allocatable, intent(out) :: err
            !! Error instance.
            logical, intent(in), optional :: debug
            !! Enable debugging printing.

        end subroutine fix_norm_AOs_bsetBF

        module subroutine fix_norm_AOs_bsetDB( &
                iout, n_at, n_ao, at_crd, bsetDB, err, debug)
            !! Correct basis sets coefficients to normalize the AOs.
            !!
            !! Checks if atomic orbitals are normalized and otherwise correct
            !! the coefficients to ensure the normalization.
            !!
            !! This version expects a basis set database as argument.

            ! Arguments
            integer, intent(in) :: iout
            !! Unit for output.
            integer, intent(in) :: n_at
            !! Number of atoms.
            integer, intent(in) :: n_ao
            !! Number of atomic orbitals.
            real(realwp), dimension(:,:), intent(in) :: at_crd
            !! Atomic coordinates (in au).
            type(BasisSetDB), intent(inout) :: bsetDB
            !! Basis set's basis function information.
            class(BaseException), allocatable, intent(out) :: err
            !! Error instance.
            logical, intent(in), optional :: debug
            !! Enable debugging printing.

        end subroutine fix_norm_AOs_bsetDB

    end interface fix_norm_AOs

! ----------------------------------------------------------------------

    interface len_shells_on_atom
        !! Get length of shells centered on chosen atom.
        module function len_shells_on_atom_bsetBF( &
                bsetBF, nprim_per_atom, num_shells, ia) &
                result(len_shells_on_atom)
            !! Length of each shell on atom `ia`.
            !!
            !! Calculates the length of each shell centered on atom given in
            !! input.
            !!
            !! This version expects basis set information as separate arguments.

            type(PrimitiveFunction), dimension(:), intent(in) :: bsetBF
            !! Basis set's basis function information.
            integer, dimension(:), intent(in) :: nprim_per_atom
            !! Number of basis primitives per atom.
            integer, intent(in) :: num_shells
            !! Number of shells
            integer, intent(in) :: ia
            !! Atom index.
            integer, dimension(num_shells) :: len_shells_on_atom
            !! Size of each shell on atom `ia`.

        end function len_shells_on_atom_bsetBF

        module function len_shells_on_atom_bsetDB(bsetDB, num_shells, ia) &
                result(len_shells_on_atom)
            !! Length of each shell on atom `ia`.
            !!
            !! Calculates the length of each shell centered on atom given in
            !! input.
            !!
            !! This version expects a basis set database as argument.

            type(BasisSetDB), intent(in) :: bsetDB
            !! Basis set database.
            integer, intent(in) :: num_shells
            !! Number of shells
            integer, intent(in) :: ia
            !! Atom index.
            integer, dimension(num_shells) :: len_shells_on_atom
            !! Size of each shell on atom `ia`.

        end function len_shells_on_atom_bsetDB
    end interface len_shells_on_atom

! ----------------------------------------------------------------------

    interface num_cart_AOs
        !! Number of Cartesian atomic orbitals based on pure.
        module function num_cart_AOs_bsetBF(n_ao, bsetBF, nprim_per_atom) &
                result(num_cart_AOs)
            !! Number of Cartesian atomic orbitals.
            !!
            !! Computes the number of Cartesian-type atomic orbitals.
            !!
            !! This version expects basis set information as separate arguments.

            integer, intent(in) :: n_ao
            !! Number of atomic orbitals.
            type(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
            !! Basis set's basis function information (pure).
            integer, dimension(:), intent(in) :: nprim_per_atom
            !! Number of basis primitives per atom.
            integer :: num_cart_AOs
            !! Number of Cartesian atomic orbitals.

        end function num_cart_AOs_bsetBF

        ! ----------------------------------------------------------------------

        module function num_cart_AOs_bsetDB(n_ao, bsetDB) result(num_cart_AOs)
            !! Number of Cartesian atomic orbitals.
            !!
            !! Computes the number of Cartesian-type atomic orbitals.
            !!
            !! This version expects a basis set database as argument.

            integer, intent(in) :: n_ao
            !! Number of atomic orbitals.
            type(BasisSetDB), intent(in) :: bsetDB
            !! Basis set database (pure).
            integer :: num_cart_AOs
            !! Number of Cartesian atomic orbitals.

        end function num_cart_AOs_bsetDB

    end interface num_cart_AOs

! ----------------------------------------------------------------------

    interface num_shells_on_atom
        !! Get the number of shells centered on chosen atom.
        module function num_shells_on_atom_bsetBF(bsetBF, nprim_per_atom, ia) &
                result(num_shells_on_atom)
            !! Number of shells on atom `ia`.
            !!
            !! Calculates the number of shells centered on atom given in input.
            !!
            !! This version expects basis set information as separate arguments.

            type(PrimitiveFunction), dimension(:), intent(in) :: bsetBF
            !! Basis set's basis function information.
            integer, dimension(:), intent(in) :: nprim_per_atom
            !! Number of basis primitives per atom.
            integer, intent(in) :: ia
            !! Atom index.
            integer :: num_shells_on_atom
            !! Number of shells on atom of interest.

        end function num_shells_on_atom_bsetBF

        module function num_shells_on_atom_bsetDB(bsetDB, ia) &
                result(num_shells_on_atom)
            !! Number of shells on atom `ia`.
            !!
            !! Calculates the number of shells centered on atom given in input.
            !!
            !! This version expects a basis set database as argument.

            type(BasisSetDB), intent(in) :: bsetDB
            !! Basis set database.
            integer, intent(in) :: ia
            !! Atom index.
            integer :: num_shells_on_atom
            !! Number of shells on atom of interest.

        end function num_shells_on_atom_bsetDB
    end interface num_shells_on_atom


contains

! ======================================================================

subroutine build_bset_DB(n_at, n_shells, pureD, pureF, shell_types, &
                         prim_per_sh, shell_to_at, coef_contr, coef_contrSP, &
                         prim_exp, nprim_per_at, bset, err)
    !! Build basis set database.
    !!
    !! Builds a list of the basis functions sorted by atoms and
    !! primitives.

    ! Arguments
    integer, intent(in) :: n_at
    !! Number of atoms.
    integer, intent(in) :: n_shells
    !! Number of contracted shells.
    logical, intent(in) :: pureD
    !! If True, spherical harm. are used for D orbitals, Cart. otherwise.
    logical, intent(in) :: pureF
    !! If True, spherical harmonics are used for F and above.
    integer, dimension(:), intent(in) :: shell_types
    !! Shell types, as integers.
    integer, dimension(:), intent(in) :: prim_per_sh
    !! number of primitives per shell.
    integer, dimension(:), intent(in) :: shell_to_at
    !! shell to atom mapping.
    real(realwp), dimension(:), intent(in) :: coef_contr
    !! contraction coefficients.
    real(realwp), dimension(:), allocatable, intent(in) :: coef_contrSP
    !! contraction coefficients for shells of type S=P.
    real(realwp), dimension(:), intent(in) :: prim_exp
    !! primitives exponents.
    integer, dimension(:), allocatable, intent(out) :: nprim_per_at
    !! Number of primitive per atom.
    type(PrimitiveFunction), dimension(:,:), allocatable, intent(out) :: bset
    !! Basis set database.
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance.

    ! Local
    integer :: i, ia, ippa, iprim, iptot, L_ang, ndim, ntot_AOs
    integer :: ish
    integer, dimension(:), allocatable :: shell_per_at
    real(realwp) :: c1, c2
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
                c2 = f0
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

subroutine chk_bset_redundancy(n_ao, ovint_i_j, thresh)
    !! Check basis set redundancy.
    !!
    !! Checks any redundancy within basis set functions.

    ! Arguments
    integer, intent(in) :: n_ao
    !! Number of atomic orbitals.
    real(realwp), dimension(n_ao,n_ao), intent(in) :: ovint_i_j
    !! Overlap integrals between atomic orbitals, < i | j >.
    real(realwp), intent(in), optional :: thresh
    !! Threshold to consider redundancy (redundancy if norm below).

    ! Local
    integer :: iao, jao, kao, lao
    real(realwp) :: norm, ovlp, thresh0
    real(realwp), dimension(n_ao,n_ao) :: trial
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
        thresh0 = 1.0e-6_realwp
    end if

    write(fmt_red, '("(""Orbital num. "",i",i0,","" - Rest norm:"",f12.8,&
        &:,"" ["",a,""]"")")') len_int(n_ao)

    ! Initialize trial matrix to check overlap
    trial = f0
    do iao = 1, n_ao
        trial(iao,iao) = f1
    end do
    is_red = .false.

    do iao = 1, n_ao
        do jao = 1, iao-1
            if (.not.is_red(jao)) then
                ovlp = f0
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

function coef_C2P(L, M, Lx, Ly) result(coef)
    !! Calculate the coefficient from Cartesian to pure.
    !!
    !! Calculates the coefficient for the conversion from the normalized
    !! Cartesian component Lx,Ly,Lz=L-Lx-Ly to the spherical harmonics
    !! component L,M.

    ! Arguments
    real(realwp) :: coef
    integer, intent(in) :: L
    !! Angular moment
    integer, intent(in) :: M
    !! M component of the spherical harmonics
    integer, intent(in) :: Lx
    !! Lx component of the Cartesian function
    integer, intent(in) :: Ly
    !! Ly component of the Cartesian function

    ! Local
    integer :: i, j, k, Lx2k, Lz, Ma
    real(realwp) :: a, b, c, d, e, g

    Lz = L - Lx - Ly
    Ma = Abs(M)
    if (Ma > L .or. Lz < 0) then
        coef = f0
        return
    end if
    a = sqrt(real(fac(2*Lx)*fac(2*Ly)*fac(2*Lz)*fac(L), kind=realwp) &
             / (fac(2*L)*fac(Lx)*fac(Ly)*fac(Lz)))
    b = sqrt(real(fac(L-Ma), kind=realwp)/fac(L+Ma))/((2**L)*fac(L))
    j = (L-Ma-Lz)/2
    g = f0
    if (2*j == L-Ma-Lz) then
        do i = 0, (L-Ma)/2
            c = f0
            if (j >= 0 .and. j <= i) then
                c = (real(fac(L), kind=realwp)/(fac(i)*fac(L-i))) &
                    *(real(fac(2*L-2*i)*((-1)**i), kind=realwp)/fac(L-Ma-2*i)) &
                    *(real(fac(i), kind=realwp)/(fac(j)*fac(i-j)))
                do k = 0, j
                    d = f0
                    Lx2k = Lx-2*k
                    if (Lx2k >= 0 .and. Lx2k <= Ma) then
                        d = (real(fac(j), kind=realwp)/(fac(k)*fac(j-k))) &
                            *(real(fac(Ma), kind=realwp) &
                                /(fac(Lx2k)*fac(Ma-Lx2k)))
                        e = f0
                        if (M==0 .and. mod(Lx,2) == 0) e = (-1)**(-Lx2k/2)
                        if (M > 0 .and. mod(abs(Ma-Lx),2) == 0) &
                            e = Sqrt(f2)*(-1)**((Ma-Lx2k)/2)
                        if (M < 0 .and. mod(abs(Ma-Lx),2) == 1) &
                            e = Sqrt(f2)*(-1)**((Ma-Lx2k)/2)
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

subroutine coefs_norm_sh(shtype, coef, alpha, coef2, new_coefs, indexes, err)
    !! Calculate the normalized coefficients for a given primitive shell.
    !!
    !! Computes and returned the unique normalized coefficients relevant
    !! for a given primitive shell.

    ! Arguments
    character(len=*), intent(in) :: shtype
    !! Shell type
    real(realwp), intent(in) :: coef
    !! Contracted coefficient
    real(realwp), intent(in) :: alpha
    !! Primitive exponent alpha
    real(realwp), intent(in) :: coef2
    !! Secondary contracted coefficient, for instance for S=P shell
    real(realwp), dimension(:), allocatable, intent(out) :: new_coefs
    !! New, normalized coefficients
    integer, dimension(:,:), allocatable, intent(out) :: indexes
    !! Indexes corresponding to each coefficients
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance

    ! Local
    real(realwp), parameter :: sq2pi = sqrt(f2*pi)
    integer :: i, ix, iy, iz, L
    real(realwp) :: cnorm, f2sqal

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
            new_coefs(1) = new_coefs(2)/sqrt(f3)
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
            cnorm = sqrt(f2sqal**9/(sq2pi**3))
            new_coefs(0) = coef
            new_coefs(3) = coef*cnorm
            new_coefs(2) = new_coefs(3)/sqrt(f3)
            new_coefs(1) = new_coefs(3)/sqrt(15.0_realwp)
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
    !! Generate coefficients for conversion between pure and Cart. bases.
    !!
    !! Generates the coefficients to transform between normalized
    !! Cartesian and normalized spherical harmonics (pure).

    ! Arguments
    integer, intent(in) :: L_ang
    !! Angular momentum to transform
    integer, intent(in) :: Ncart
    !! Number of Cartesian components
    integer, intent(in) :: Npure
    !! Number of pure/spherical components
    real(realwp), dimension(Ncart,Npure), intent(out) :: coef2P
    !! Transformation coefficients to pure
    real(realwp), dimension(Ncart,Npure), intent(out) :: coef2C
    !! Transformation coefficients to Cartesian

    ! Local
    integer :: i1, i2, k, Lx, Lx1, Lx2, Lxyz, Ly, Ly1, Ly2, Lz, Lz1, &
        Lz2, M
    real(realwp) :: a1, a2, s

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
    coef2C = f0
    i1 = 0
    do Lx1 = 0, L_ang
        do Ly1 = 0, (L_ang-Lx1)
            i1 = i1 + 1
            Lz1 = L_ang - Lx1 - Ly1
            a1 = sqrt(real(fac(Lx1)*fac(Ly1)*fac(Lz1), kind=realwp) &
                      /(fac(2*Lx1)*fac(2*Ly1)*fac(2*Lz1)))
            i2 = 0
            do Lx2 = 0, L_ang
                do Ly2 = 0, (L_ang-Lx2)
                    i2 = i2 + 1
                    Lz2 = L_ang - Lx2 - Ly2
                    a2 = sqrt(real(fac(Lx2)*fac(Ly2)*fac(Lz2), kind=realwp) &
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

subroutine set_primC_comp(bfunc, ndimC, lxyz, coefs, err)
    !! Set component for a Cartesian primitive.
    !!
    !! Builds the list of Cartesian components and the associated
    !! coefficients for a primitive function based on the information
    !! in `bfunc`.

    ! Arguments
    type(PrimitiveFunction), intent(in) :: bfunc
    !! Basis set function.
    integer, intent(out) :: ndimC
    !! True number of dimension.
    integer, dimension(:,:), allocatable, intent(out) :: lxyz
    !! Number of occurrence of each x,y,z coordinate for each dimension.
    real(realwp), dimension(:), allocatable, intent(out) :: coefs
    !! Contraction coefficient for each dimension.
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
            ! lxyz(:, 1) = [3, 0, 0]
            ! lxyz(:, 2) = [0, 3, 0]
            ! lxyz(:, 3) = [0, 0, 3]
            ! lxyz(:, 4) = [1, 2, 0]
            ! lxyz(:, 5) = [2, 1, 0]
            ! lxyz(:, 6) = [2, 0, 1]
            ! lxyz(:, 7) = [1, 0, 2]
            ! lxyz(:, 8) = [0, 1, 2]
            ! lxyz(:, 9) = [0, 2, 1]
            ! lxyz(:,10) = [1, 1, 1]
            lxyz(:, 1) = [3, 0, 0]
            lxyz(:, 2) = [0, 3, 0]
            lxyz(:, 3) = [0, 0, 3]
            lxyz(:, 4) = [2, 1, 0]
            lxyz(:, 5) = [2, 0, 1]
            lxyz(:, 6) = [1, 2, 0]
            lxyz(:, 7) = [0, 2, 1]
            lxyz(:, 8) = [1, 0, 2]
            lxyz(:, 9) = [0, 1, 2]
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
    !! Cartesian coordinates to a one with spherical harmonics.
    !!
    !! The conversion table is taken from SOS.
    !!
    !! **For *f* orbitals**
    !!
    !! Conversion table:
    !!
    !! |    |          |    1  |   2  |    3  |   4  |   5  |   6 |      |
    !! | ---|----------|-------|------|-------|------|------|-----|----- |
    !! |    |          |   xx  |  yy  |   zz  |  xy  |  xz  |  yz |      |
    !! | ---|----------|-------|------|-------|------|------|-----|----- |
    !! |  1 | 3z^2-R^2 |   -1  |  -1  |    2  |   0  |   0  |   0 | *s1  |
    !! |  2 |    xz    |    0  |   0  |    0  |   0  |   1  |   0 |      |
    !! |  3 |    yz    |    0  |   0  |    0  |   0  |   0  |   1 |      |
    !! |  4 |  x^2-y^2 |    1  |  -1  |    0  |   0  |   0  |   0 | *s5  |
    !! |  5 |    xy    |    0  |   0  |    0  |   1  |   0  |   0 |      |
    !!
    !! s1 = 1/2 ; s5 = sqrt(3)/2
    !!
    !! **For *f* orbitals**
    !!
    !! |      |     spherical  |  cartesian  |
    !! | -----|----------------|------------ |
    !! |    1 |  5z^3-3R^2z    |     xxx     |
    !! |    2 |  x(5z^2-R^2)   |     yyy     |
    !! |    3 |  y(5z^2-R^2)   |     zzz     |
    !! |    4 |  xyz           |     xxy     |
    !! |    5 |  z(x^2-y^2)    |     xxz     |
    !! |    6 |  y(y^2-3x^2)   |     yyx     |
    !! |    7 |  x(x^2-3y^2)   |     yyz     |
    !! |    8 |  -             |     zzx     |
    !! |    9 |  -             |     zzy     |
    !! |   10 |  -             |     xyz     |
    !!
    !! Conversion table
    !!
    !! |    |  1  |   2 |    3 |   4 |   5 |   6 |   7 |   8 |   9 |  10 |      |
    !! | ---|-----|-----|------|-----|-----|-----|-----|-----|-----|-----|----- |
    !! |    | xxx | yyy |  zzz | xxy | xxz | yyx | yyz | zzx | zzy | xyz |      |
    !! | ---|-----|-----|------|-----|-----|-----|-----|-----|-----|-----|----- |
    !! |  1 |   0 |   0 | 2*v5 |   0 |  -3 |   0 |  -3 |   0 |   0 |   0 | *s1  |
    !! |  2 | -v5 |   0 |    0 |   0 |   0 |  -1 |   0 |   4 |   0 |   0 | *s2  |
    !! |  3 |   0 | -v5 |    0 |  -1 |   0 |   0 |   0 |   0 |   4 |   0 | *s3  |
    !! |  4 |   0 |   0 |    0 |   0 |   1 |   0 |  -1 |   0 |   0 |   0 | *s4  |
    !! |  5 |   0 |   0 |    0 |   0 |   0 |   0 |   0 |   0 |   0 |   1 | *s5  |
    !! |  6 |  v5 |   0 |    0 |   0 |   0 |  -3 |   0 |   0 |   0 |   0 | *s6  |
    !! |  7 |   0 | -v5 |    0 |   3 |   0 |   0 |   0 |   0 |   0 |   0 | *s7  |
    !!
    !! @warning
    !! The function does not check explicitly if the shell is supported.
    !! Calling procedure should check if the array is allocated on exit.
    !! @endwarning

    ! Arguments
    integer, intent(in) :: L_ang
    !! Angular momentum for the shell of interest
    real(realwp), dimension(:,:), allocatable :: convmat
    !! Conversion matrix of dimension (npure,ncart)

    ! Local
    integer :: ncart, npure
    real(realwp), parameter :: sq2=sqrt(f2), sq3=sqrt(f3), sq5=sqrt(f5), &
        sq10=sqrt(f10)
    real(realwp) :: s1, s2, s3, s4, s5, s6, s7
    real(realwp), dimension(max_nxyz,max_nxyz) :: tempmat

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
        case (-1)
            allocate(convmat(4,4))
            convmat = f0
            convmat(1,1) = f1
            convmat(2,2) = f1
            convmat(3,3) = f1
            convmat(4,4) = f1
        case (2)
            allocate(convmat(5,6))
            convmat = f0
            ! Unchanged components
            convmat(2,5) = f1
            convmat(3,6) = f1
            convmat(5,4) = f1
            ! |normalized z^2-r^2> = s1*(two*|norm.zz>-|norm.xx>-|norm.yy>)
            s1 = fhalf
            convmat(1,1) = -s1
            convmat(1,2) = -s1
            convmat(1,3) = f2*s1
            ! |normalized x^2-y^2> = s5*(|norm.xx>-|norm.yy>)
            s5 = fhalf*sq3
            convmat(4,1) = s5
            convmat(4,2) = -s5
        case (3)
            allocate(convmat(7,10))
            convmat = f0
            ! normalization coefficients
            s1 = fhalf/sq5
            s2 = fhalf*sq3/sq10
            s3 = fhalf*sq3/sq10
            s4 = fhalf*sq3
            s5 = f1
            s6 = fhalf/sq2
            s7 = fhalf/sq2
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
        case default
            write(iu_out, '("Unsupported angular momentum: ",i0)') L_ang
            error stop 1
    end select

end function transfo_cart2pure

! ======================================================================

end module basisset
