submodule (basisset) basisset_purecart
    !! Submodule containing the definition of procedures related to the
    !! conversion between Cartesian and pure basis sets.

contains

! ======================================================================

module procedure convert_pure2cart_bsetBF

    integer :: center, iprim

    bsetBF_cart = bsetBF
    do center = 1, size(bsetBF, 1)
        do iprim = 1, nprim_per_atom(center)

            ! bsetBF_cart(center,iprim)%l = bsetBF(center,iprim)%l
            ! bsetBF_cart(center,iprim)%shelltype = bsetBF(center,iprim)%shelltype
            ! bsetBF_cart(center,iprim)%shellid = bsetBF(center,iprim)%shellid
            bsetBF_cart(center,iprim)%pure  = .false.

            ! bsetBF_cart(center,iprim)%shell_first = bsetBF(center,iprim)%shell_first
            ! bsetBF_cart(center,iprim)%shell_last = bsetBF(center,iprim)%shell_last
            ! bsetBF_cart(center,iprim)%alpha = bsetBF(center,iprim)%alpha
            ! bsetBF_cart(center,iprim)%coeff = bsetBF(center,iprim)%coeff

            if (bsetBF(center,iprim)%shelltype == 'D') then
                bsetBF_cart(center,iprim)%ndim = 6
            else if (bsetBF(center,iprim)%shelltype == 'F') then
                bsetBF_cart(center,iprim)%ndim = 10
            else
                bsetBF_cart(center,iprim)%ndim = bsetBF(center,iprim)%ndim
            endif
        enddo
    enddo

end procedure convert_pure2cart_bsetBF

! ======================================================================

module procedure convert_pure2cart_bsetDB

    ! Local
    integer :: center, iprim

    bsetBF_cart = bsetDB%info
    do center = 1, size(bsetDB%info, 1)
        do iprim = 1, bsetDB%nprim_per_at(center)

            ! bsetBF_cart(center,iprim)%l = bsetDB%info(center,iprim)%l
            ! bsetBF_cart(center,iprim)%shelltype = bsetDB%info(center,iprim)%shelltype
            ! bsetBF_cart(center,iprim)%shellid = bsetDB%info(center,iprim)%shellid
            bsetBF_cart(center,iprim)%pure  = .false.

            ! bsetBF_cart(center,iprim)%shell_first = bsetDB%info(center,iprim)%shell_first
            ! bsetBF_cart(center,iprim)%shell_last = bsetDB%info(center,iprim)%shell_last
            ! bsetBF_cart(center,iprim)%alpha = bsetDB%info(center,iprim)%alpha
            ! bsetBF_cart(center,iprim)%coeff = bsetDB%info(center,iprim)%coeff

            if (bsetDB%info(center,iprim)%shelltype == 'D') then
                bsetBF_cart(center,iprim)%ndim = 6
            else if (bsetDB%info(center,iprim)%shelltype == 'F') then
                bsetBF_cart(center,iprim)%ndim = 10
            else
                bsetBF_cart(center,iprim)%ndim = bsetDB%info(center,iprim)%ndim
            endif
        enddo
    enddo

end procedure convert_pure2cart_bsetDB

! ======================================================================

module procedure convert_pure2cart_matrix_bsetBF

    ! Local
    integer :: center_mu, l_mu, len_mu_C, id_mu_C, id_mu_P, &
        len_mu_P, offset_mu, shell_mu, shell_per_atom_mu
    integer :: center_nu, l_nu, len_nu_C, id_nu_C, id_nu_P, &
        len_nu_P, offset_nu, shell_nu, shell_per_atom_nu
    integer :: sp_i, sp_j
    integer, dimension(:), allocatable :: len_shell_per_atom_mu, len_shell_per_atom_nu
    real(realwp), dimension(:,:), allocatable :: mu_t, nu_t
    ! real(realwp), dimension(size(matrix_cart, 1), size(matrix_pure, 1)) :: tmp
    real(realwp), dimension(:,:), allocatable :: tmp

    allocate(tmp(size(matrix_cart, 1), size(matrix_pure, 1)))

    id_mu_P = 1
    id_mu_C = 1

    do center_mu = 1, size(bsetBF, 1)
        offset_mu = 1
        shell_per_atom_mu = num_shells_on_atom(bsetBF(center_mu,:), &
            nprim_per_atom, center_mu)
        allocate(len_shell_per_atom_mu(shell_per_atom_mu))
        len_shell_per_atom_mu = len_shells_on_atom( &
            bsetBF(center_mu, :), nprim_per_atom, shell_per_atom_mu, center_mu)

        do shell_mu = 1, shell_per_atom_mu
            l_mu = bsetBF(center_mu, offset_mu)%l
            if (l_mu == -1) then
                len_mu_C = 4
                len_mu_P = 4
            else
                len_mu_C = (l_mu + 1) * (l_mu + 2) / 2
                len_mu_P = 2 * l_mu + 1
            endif

            mu_t = transfo_cart2pure(l_mu)
            offset_mu = offset_mu + len_shell_per_atom_mu(shell_mu)
            id_nu_P = 1
            id_nu_C = 1

            do center_nu = 1, size(bsetBF, 1)
                offset_nu = 1
                shell_per_atom_nu = num_shells_on_atom( &
                    bsetBF(center_nu, :), nprim_per_atom, center_nu)
                allocate(len_shell_per_atom_nu(shell_per_atom_nu))
                len_shell_per_atom_nu = len_shells_on_atom( &
                    bsetBF(center_nu, :), nprim_per_atom, shell_per_atom_nu, &
                    center_nu)

                do shell_nu = 1, shell_per_atom_nu
                    l_nu = bsetBF(center_nu, offset_nu)%l
                    if (l_nu == -1) then
                        len_nu_C = 4
                        len_nu_P = 4
                    else
                        len_nu_C = (l_nu + 1) * (l_nu + 2) / 2
                        len_nu_P = 2 * l_nu + 1
                    endif
                    nu_t = transfo_cart2pure(l_nu)
                    call xgemm("T", "N", len_nu_C, len_mu_P, len_nu_P, f1, &
                        nu_t, len_nu_P, &
                        matrix_pure(id_nu_P : id_nu_P + len_nu_P - 1, &
                                    id_mu_P : id_mu_P + len_mu_P - 1), &
                        len_nu_P, f0, tmp, len_nu_C)
                    call xgemm("N", "N", len_nu_C, len_mu_C, len_mu_P, f1, &
                        tmp, len_nu_C, mu_t, len_mu_P, f0, &
                        matrix_cart(id_nu_C : id_nu_C + len_nu_C - 1, &
                                    id_mu_C : id_mu_C + len_mu_C - 1), &
                        len_nu_C)
                    offset_nu = offset_nu + len_shell_per_atom_nu(shell_nu)
                    id_nu_P = id_nu_P + len_nu_P
                    id_nu_C = id_nu_C + len_nu_C
                end do
                deallocate(len_shell_per_atom_nu)
            end do
            id_mu_P = id_mu_P + len_mu_P
            id_mu_C = id_mu_C + len_mu_C
        end do
        deallocate(len_shell_per_atom_mu)
    end do

    deallocate(tmp)

end procedure convert_pure2cart_matrix_bsetBF

! ======================================================================

module procedure convert_pure2cart_matrix_bsetDB

    ! Local
    integer :: center_mu, l_mu, len_mu_C, id_mu_C, id_mu_P, &
        len_mu_P, offset_mu, shell_mu, shell_per_atom_mu
    integer :: center_nu, l_nu, len_nu_C, id_nu_C, id_nu_P, &
        len_nu_P, offset_nu, shell_nu, shell_per_atom_nu
    integer :: sp_i, sp_j
    integer, dimension(:), allocatable :: len_shell_per_atom_mu, len_shell_per_atom_nu
    real(realwp), dimension(:,:), allocatable :: mu_t, nu_t
    ! real(realwp), dimension(size(matrix_cart, 1), size(matrix_pure, 1)) :: tmp
    real(realwp), dimension(:,:), allocatable :: tmp

    allocate(tmp(size(matrix_cart, 1), size(matrix_pure, 1)))

    id_mu_P = 1
    id_mu_C = 1

    do center_mu = 1, size(bsetDB%info, 1)
        offset_mu = 1
        shell_per_atom_mu = num_shells_on_atom(bsetDB, center_mu)
        allocate(len_shell_per_atom_mu(shell_per_atom_mu))
        len_shell_per_atom_mu = len_shells_on_atom(bsetDB, shell_per_atom_mu, &
            center_mu)

        do shell_mu = 1, shell_per_atom_mu
            l_mu = bsetDB%info(center_mu, offset_mu)%l
            if (l_mu == -1) then
                len_mu_C = 4
                len_mu_P = 4
            else
                len_mu_C = (l_mu + 1) * (l_mu + 2) / 2
                len_mu_P = 2 * l_mu + 1
            endif

            mu_t = transfo_cart2pure(l_mu)
            offset_mu = offset_mu + len_shell_per_atom_mu(shell_mu)
            id_nu_P = 1
            id_nu_C = 1

            do center_nu = 1, size(bsetDB%info, 1)
                offset_nu = 1
                shell_per_atom_nu = num_shells_on_atom(bsetDB, center_nu)
                allocate(len_shell_per_atom_nu(shell_per_atom_nu))
                len_shell_per_atom_nu = len_shells_on_atom(bsetDB, &
                    shell_per_atom_nu, center_nu)

                do shell_nu = 1, shell_per_atom_nu
                    l_nu = bsetDB%info(center_nu, offset_nu)%l
                    if (l_nu == -1) then
                        len_nu_C = 4
                        len_nu_P = 4
                    else
                        len_nu_C = (l_nu + 1) * (l_nu + 2) / 2
                        len_nu_P = 2 * l_nu + 1
                    endif
                    nu_t = transfo_cart2pure(l_nu)
                    call xgemm("T", "N", len_nu_C, len_mu_P, len_nu_P, f1, &
                        nu_t, len_nu_P, &
                        matrix_pure(id_nu_P : id_nu_P + len_nu_P - 1, &
                                    id_mu_P : id_mu_P + len_mu_P - 1), &
                        len_nu_P, f0, tmp, len_nu_C)
                    call xgemm("N", "N", len_nu_C, len_mu_C, len_mu_P, f1, &
                        tmp, len_nu_C, mu_t, len_mu_P, f0, &
                        matrix_cart(id_nu_C : id_nu_C + len_nu_C - 1, &
                                    id_mu_C : id_mu_C + len_mu_C - 1), &
                        len_nu_C)
                    offset_nu = offset_nu + len_shell_per_atom_nu(shell_nu)
                    id_nu_P = id_nu_P + len_nu_P
                    id_nu_C = id_nu_C + len_nu_C
                end do
                deallocate(len_shell_per_atom_nu)
            end do
            id_mu_P = id_mu_P + len_mu_P
            id_mu_C = id_mu_C + len_mu_C
        end do
        deallocate(len_shell_per_atom_mu)
    end do

    deallocate(tmp)

end procedure convert_pure2cart_matrix_bsetDB

! ======================================================================

module procedure fix_norm_AOs_bsetBF

    integer, parameter :: &
        dimPa = 24  ! maximum used dimension for Pascal's triangle
    real(realwp), parameter :: &
        tol_expm = 100.0_realwp, &  ! Max. value of x for e^-x to be relevant
        sqpi3 = sqrt(pi**3)
    integer :: ndi, ndj
    integer :: i, idi, idj, iprim, j, jj, jprim
    integer :: ia, ia0, ja, ja0
    integer, dimension(:,:), allocatable :: ldi, ldj
    real(realwp) :: a, ai, aj, cijer, cjer, ebase, r2ij, x
    real(realwp), dimension(3) :: ri, rj, ovi
    real(realwp), dimension(max_nxyz,max_nxyz) :: ovlp_ij
    real(realwp), dimension(max_nxyz,max_nxyz), target :: allones
    real(realwp), dimension(:), allocatable :: ao_norms, ci, cj
    real(realwp), dimension(:,:), allocatable :: ao_ovlp
    real(realwp), dimension(:,:), allocatable, target :: c2p_D, c2p_F, c2p_G, &
        c2p_H, c2p_I
    real(realwp), dimension(:,:), pointer :: c2pi, c2pj
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
    ao_norms = f0
    ao_ovlp = f0

    ia0 = 0
    do ia = 1, n_at
        ri = at_crd(:,ia)
        do iprim = 1, nprim_per_at(ia)
            ai = bsetBF(ia,iprim)%alpha
            call set_primC_comp(bsetBF(ia,iprim), ndi, ldi, ci, suberr)
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
            if (bsetBF(ia,iprim)%pure) then
                select case (bsetBF(ia,iprim)%shelltype)
                    case ('D')
                        if (.not.allocated(c2p_D)) &
                            c2p_D = transfo_cart2pure(bsetBF(ia,iprim)%L)
                        c2pi => c2p_D
                    case ('F')
                        if (.not.allocated(c2p_F)) &
                            c2p_F = transfo_cart2pure(bsetBF(ia,iprim)%L)
                        c2pi => c2p_F
                    case ('G')
                        if (.not.allocated(c2p_G)) &
                            c2p_G = transfo_cart2pure(bsetBF(ia,iprim)%L)
                        c2pi => c2p_G
                    case ('H')
                        if (.not.allocated(c2p_H)) &
                            c2p_H = transfo_cart2pure(bsetBF(ia,iprim)%L)
                        c2pi => c2p_H
                    case ('I')
                        if (.not.allocated(c2p_I)) &
                            c2p_I = transfo_cart2pure(bsetBF(ia,iprim)%L)
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
                    aj = bsetBF(ja,jprim)%alpha
                    call set_primC_comp(bsetBF(ja,jprim), ndj, ldj, cj, suberr)
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
                    if (bsetBF(ja,jprim)%pure) then
                        select case (bsetBF(ja,jprim)%shelltype)
                            case ('D')
                                if (.not.allocated(c2p_D)) &
                                    c2p_D = transfo_cart2pure(bsetBF(ja,jprim)%L)
                                c2pj => c2p_D
                            case ('F')
                                if (.not.allocated(c2p_F)) &
                                    c2p_F = transfo_cart2pure(bsetBF(ja,jprim)%L)
                                c2pj => c2p_F
                            case ('G')
                                if (.not.allocated(c2p_G)) &
                                    c2p_G = transfo_cart2pure(bsetBF(ja,jprim)%L)
                                c2pj => c2p_G
                            case ('H')
                                if (.not.allocated(c2p_H)) &
                                    c2p_H = transfo_cart2pure(bsetBF(ja,jprim)%L)
                                c2pj => c2p_H
                            case ('I')
                                if (.not.allocated(c2p_I)) &
                                    c2p_I = transfo_cart2pure(bsetBF(ja,jprim)%L)
                                c2pj => c2p_I
                            case default
                                c2pj => allones
                        end select
                    else
                        c2pj => allones
                    end if
                    a = f1 / (ai+aj)
                    x = ai*aj*r2ij*a
                    if (x < tol_expm) then
                        ebase = sqpi3*exp(-x)
                    else
                        if (bsetBF(ja,jprim)%shell_last) &
                            ja0 = ja0 + bsetBF(ja,jprim)%ndim
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
                    if (bsetBF(ia,iprim)%pure .or. bsetBF(ja,jprim)%pure) then
                        do i = 1, bsetBF(ia,iprim)%ndim
                            do j = 1, bsetBF(ja,jprim)%ndim
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
                    if (bsetBF(ja,jprim)%shell_last) &
                        ja0 = ja0 + bsetBF(ja,jprim)%ndim
                end do
            end do
            if (bsetBF(ia,iprim)%shell_last) then
                do i = 1, bsetBF(ia,iprim)%ndim
                    ao_norms(ia0+i) = ao_ovlp(i,ia0+i)
                end do
                ia0 = ia0 + bsetBF(ia,iprim)%ndim
                ao_ovlp = f0
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
            a = f1 / sqrt(ao_norms(ia0))
            bsetBF(ia,iprim)%coeff = bsetBF(ia,iprim)%coeff*a
            if (bsetBF(ia,iprim)%shell_last) ia0 = ia0 + bsetBF(ia,iprim)%ndim
        end do
    end do

    return
end procedure fix_norm_AOs_bsetBF

! ======================================================================

module procedure fix_norm_AOs_bsetDB

    integer, parameter :: &
        dimPa = 24  ! maximum used dimension for Pascal's triangle
    real(realwp), parameter :: &
        tol_expm = 100.0_realwp, &  ! Max. value of x for e^-x to be relevant
        sqpi3 = sqrt(pi**3)
    integer :: ndi, ndj
    integer :: i, idi, idj, iprim, j, jj, jprim
    integer :: ia, ia0, ja, ja0
    integer, dimension(:,:), allocatable :: ldi, ldj
    real(realwp) :: a, ai, aj, cijer, cjer, ebase, r2ij, x
    real(realwp), dimension(3) :: ri, rj, ovi
    real(realwp), dimension(max_nxyz,max_nxyz) :: ovlp_ij
    real(realwp), dimension(max_nxyz,max_nxyz), target :: allones
    real(realwp), dimension(:), allocatable :: ao_norms, ci, cj
    real(realwp), dimension(:,:), allocatable :: ao_ovlp
    real(realwp), dimension(:,:), allocatable, target :: c2p_D, c2p_F, c2p_G, &
        c2p_H, c2p_I
    real(realwp), dimension(:,:), pointer :: c2pi, c2pj
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
    ao_norms = f0
    ao_ovlp = f0

    ia0 = 0
    do ia = 1, n_at
        ri = at_crd(:,ia)
        do iprim = 1, bsetDB%nprim_per_at(ia)
            ai = bsetDB%info(ia,iprim)%alpha
            call set_primC_comp(bsetDB%info(ia,iprim), ndi, ldi, ci, suberr)
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
            if (bsetDB%info(ia,iprim)%pure) then
                select case (bsetDB%info(ia,iprim)%shelltype)
                    case ('D')
                        if (.not.allocated(c2p_D)) &
                            c2p_D = transfo_cart2pure(bsetDB%info(ia,iprim)%L)
                        c2pi => c2p_D
                    case ('F')
                        if (.not.allocated(c2p_F)) &
                            c2p_F = transfo_cart2pure(bsetDB%info(ia,iprim)%L)
                        c2pi => c2p_F
                    case ('G')
                        if (.not.allocated(c2p_G)) &
                            c2p_G = transfo_cart2pure(bsetDB%info(ia,iprim)%L)
                        c2pi => c2p_G
                    case ('H')
                        if (.not.allocated(c2p_H)) &
                            c2p_H = transfo_cart2pure(bsetDB%info(ia,iprim)%L)
                        c2pi => c2p_H
                    case ('I')
                        if (.not.allocated(c2p_I)) &
                            c2p_I = transfo_cart2pure(bsetDB%info(ia,iprim)%L)
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
                do jprim = 1, bsetDB%nprim_per_at(ja)
                    aj = bsetDB%info(ja,jprim)%alpha
                    call set_primC_comp(bsetDB%info(ja,jprim), ndj, ldj, cj, &
                        suberr)
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
                    if (bsetDB%info(ja,jprim)%pure) then
                        select case (bsetDB%info(ja,jprim)%shelltype)
                            case ('D')
                                if (.not.allocated(c2p_D)) &
                                    c2p_D = transfo_cart2pure( &
                                        bsetDB%info(ja,jprim)%L)
                                c2pj => c2p_D
                            case ('F')
                                if (.not.allocated(c2p_F)) &
                                    c2p_F = transfo_cart2pure( &
                                        bsetDB%info(ja,jprim)%L)
                                c2pj => c2p_F
                            case ('G')
                                if (.not.allocated(c2p_G)) &
                                    c2p_G = transfo_cart2pure( &
                                        bsetDB%info(ja,jprim)%L)
                                c2pj => c2p_G
                            case ('H')
                                if (.not.allocated(c2p_H)) &
                                    c2p_H = transfo_cart2pure( &
                                        bsetDB%info(ja,jprim)%L)
                                c2pj => c2p_H
                            case ('I')
                                if (.not.allocated(c2p_I)) &
                                    c2p_I = transfo_cart2pure( &
                                        bsetDB%info(ja,jprim)%L)
                                c2pj => c2p_I
                            case default
                                c2pj => allones
                        end select
                    else
                        c2pj => allones
                    end if
                    a = f1 / (ai+aj)
                    x = ai*aj*r2ij*a
                    if (x < tol_expm) then
                        ebase = sqpi3*exp(-x)
                    else
                        if (bsetDB%info(ja,jprim)%shell_last) &
                            ja0 = ja0 + bsetDB%info(ja,jprim)%ndim
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
                    if (bsetDB%info(ia,iprim)%pure .or. &
                        bsetDB%info(ja,jprim)%pure) then
                        do i = 1, bsetDB%info(ia,iprim)%ndim
                            do j = 1, bsetDB%info(ja,jprim)%ndim
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
                    if (bsetDB%info(ja,jprim)%shell_last) &
                        ja0 = ja0 + bsetDB%info(ja,jprim)%ndim
                end do
            end do
            if (bsetDB%info(ia,iprim)%shell_last) then
                do i = 1, bsetDB%info(ia,iprim)%ndim
                    ao_norms(ia0+i) = ao_ovlp(i,ia0+i)
                end do
                ia0 = ia0 + bsetDB%info(ia,iprim)%ndim
                ao_ovlp = f0
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
        do iprim = 1, bsetDB%nprim_per_at(ia)
            a = f1 / sqrt(ao_norms(ia0))
            bsetDB%info(ia,iprim)%coeff = bsetDB%info(ia,iprim)%coeff*a
            if (bsetDB%info(ia,iprim)%shell_last) &
                ia0 = ia0 + bsetDB%info(ia,iprim)%ndim
        end do
    end do

    return
end procedure fix_norm_AOs_bsetDB

! ======================================================================

module procedure len_shells_on_atom_bsetBF

    integer :: i, counter, offset

    counter = 1
    offset = 1
    do i = 1, nprim_per_atom(ia)
        if (.not.bsetBF(i)%shell_last) then
           counter = counter + 1
        else
           len_shells_on_atom(offset) = counter
           counter = 1
           offset = offset + 1
        endif
    enddo

end procedure len_shells_on_atom_bsetBF

! ======================================================================

module procedure len_shells_on_atom_bsetDB

    integer :: i, counter, offset

    counter = 1
    offset = 1
    do i = 1, bsetDB%nprim_per_at(ia)
        if (.not.bsetDB%info(ia,i)%shell_last) then
           counter = counter + 1
        else
           len_shells_on_atom(offset) = counter
           counter = 1
           offset = offset + 1
        endif
    enddo

end procedure len_shells_on_atom_bsetDB

! ======================================================================

module procedure num_cart_AOs_bsetBF

    integer :: n_ao_cartesian, n_d, n_f

    ! Local
    integer :: i, center

    num_cart_AOs = 0
    n_d = 0
    n_f = 0

    do center = 1, size(bsetBF, 1)
        do i = 1, nprim_per_atom(center)
            if (bsetBF(center, i)%shell_first) then
                if (bsetBF(center, i)%shelltype == 'D' ) then
                    n_d = n_d + 1
                else if (bsetBF(center, i)%shelltype == 'F' ) then
                    n_f = n_f + 1
                endif
            endif
        enddo
    enddo
    num_cart_AOs = n_ao + n_d + n_f * 3

end procedure num_cart_AOs_bsetBF

! ======================================================================

module procedure num_cart_AOs_bsetDB
    
    integer :: i, center, n_d, n_f

    num_cart_AOs = 0
    n_d = 0
    n_f = 0

    do center = 1, size(bsetDB%info, 1)
        do i = 1, bsetDB%nprim_per_at(center)
            if (bsetDB%info(center, i)%shell_first) then
                if (bsetDB%info(center, i)%shelltype == 'D' ) then
                    n_d = n_d + 1
                else if (bsetDB%info(center, i)%shelltype == 'F' ) then
                    n_f = n_f + 1
                endif
            endif
        enddo
    enddo
    num_cart_AOs = n_ao + n_d + n_f * 3

end procedure num_cart_AOs_bsetDB

! ======================================================================

module procedure num_shells_on_atom_bsetBF
    
    integer :: i

    num_shells_on_atom = 0

    do i = 1, nprim_per_atom(ia)
        if (bsetBF(i)%shell_first) &
            num_shells_on_atom = num_shells_on_atom + 1
    enddo

end procedure num_shells_on_atom_bsetBF

! ======================================================================

module procedure num_shells_on_atom_bsetDB

    integer :: i

    num_shells_on_atom = 0

    do i = 1, bsetDB%nprim_per_at(ia)
        if (bsetDB%info(ia,i)%shell_first) &
            num_shells_on_atom = num_shells_on_atom + 1
    enddo

end procedure num_shells_on_atom_bsetDB

! ======================================================================

end submodule basisset_purecart