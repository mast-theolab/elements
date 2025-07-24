module orbital
    !! Orbitals (atomic/molecular)-related module.
    !!
    !! The module defines procedures related to atomic and molecular
    !! orbitals.

    use basisset, only: get_cart_L_der_sh_at, get_cart_L_norms_sh, &
        get_cart_L_sh_at, get_cart_r_der_sh_at, get_cart_r_sh_at
    use datatypes, only: BasisSetDB, MoleculeDB, PrimitiveFunction
    use exception, only: runstat
    use numeric, only: realwp, f0, f10m10

    implicit none

    real(realwp), parameter, private :: thresh_at_center = f10m10

    interface eval_AOs_nabla_chi_at
        !! Evaluate AOs and first derivatives at position.
        !!
        !! Evaluate atomic orbitals (chi) and their first derivative (nabla
        !! chi) at a chosen position.
        module procedure eval_AOs_nabla_chi_at_arr, eval_AOs_nabla_chi_at_db, &
            eval_AOs_nabla_chi_at_dim
    end interface eval_AOs_nabla_chi_at

    interface eval_AOs_chi_at
        !! Evaluate AOs and first derivatives at position.
        !!
        !! Evaluate atomic orbitals (chi) and their first derivative (nabla
        !! chi) at a chosen position.
        module procedure eval_AOs_chi_at_arr, eval_AOs_chi_at_db, &
            eval_AOs_chi_at_dim
    end interface eval_AOs_chi_at

    interface convert_AO2MO
        module procedure convert_AO2MO_1, convert_AO2MO_N
    end interface convert_AO2MO

contains

! ======================================================================

subroutine convert_AO2MO_1(n_ao, n_mo, c_ia, q_ao, q_mo, tmp_arr)
    !! Convert a scalar quantity from atomic to molecular orbitals
    !!
    !! Takes a quantity in atomic orbitals (`q_ao`) to molecular orbitals
    !!   (`q_mo`).
    !! Note: The quantity must be scalar (see convert_AO2MO_N otherwise)
    integer, intent(in) :: n_ao
    !! Number of atomic orbitals
    integer, intent(in) :: n_mo
    !! Number of molecular orbitals
    real(realwp), dimension(:,:), intent(in) :: c_ia
    !! Coefficients of MOs in AOs basis (geometry: n_mo, n_ao)
    real(realwp), dimension(:,:), intent(in) :: q_ao
    !! Quantity in AO basis
    real(realwp), dimension(:,:), intent(out) :: q_mo
    !! Quantity in MO basis
    real(realwp), dimension(:,:) :: tmp_arr
    !! Temporary array

    integer :: a, b, i, j

    !$omp parallel do collapse(2)
    do j = 1, n_mo
        do a = 1, n_ao
            tmp_arr(a,j) = f0
            do b = 1, n_ao
                tmp_arr(a,j) = tmp_arr(a,j) + c_ia(j,b)*q_ao(a,b)
            end do
        end do
    end do
    !$omp end parallel do

    !$omp parallel do collapse(2)
    do i = 1, n_mo
        do j = 1, n_mo
            q_mo(i,j) = f0
            do a = 1, n_ao
                q_mo(i,j) = q_mo(i,j) + c_ia(i,a)*tmp_arr(a,j)
            end do
        end do
    end do
    !$omp end parallel do

end subroutine convert_AO2MO_1

! ======================================================================

subroutine convert_AO2MO_N(n_ao, n_mo, c_ia, q_ao, q_mo, tmp_arr)
    !! Convert a vector from atomic to molecular orbitals
    !!
    !! Takes a quantity in atomic orbitals (`q_ao`) to molecular orbitals
    !!   (`q_mo`).
    !! Note: The quantity is expected to be a vector.
    integer, intent(in) :: n_ao
    !! Number of atomic orbitals
    integer, intent(in) :: n_mo
    !! Number of molecular orbitals
    real(realwp), dimension(:,:), intent(in) :: c_ia
    !! Coefficients of MOs in AOs basis (geometry: n_mo, n_ao)
    real(realwp), dimension(:,:,:), intent(in) :: q_ao
    !! Quantity in AO basis
    real(realwp), dimension(:,:,:), intent(out) :: q_mo
    !! Quantity in MO basis
    real(realwp), dimension(:,:,:) :: tmp_arr
    !! Temporary array

    integer :: a, b, i, j

    !$omp parallel do collapse(2)
    do j = 1, n_mo
        do a = 1, n_ao
            tmp_arr(:,a,j) = f0
            do b = 1, n_ao
                tmp_arr(:,a,j) = tmp_arr(:,a,j) + c_ia(j,b)*q_ao(:,a,b)
            end do
        end do
    end do
    !$omp end parallel do

    !$omp parallel do collapse(2)
    do i = 1, n_mo
        do j = 1, n_mo
            q_mo(:,i,j) = f0
            do a = 1, n_ao
                q_mo(:,i,j) = q_mo(:,i,j) + c_ia(i,a)*tmp_arr(:,a,j)
            end do
        end do
    end do
    !$omp end parallel do

end subroutine convert_AO2MO_N

! ======================================================================

subroutine eval_AOs_chi_at_arr(at_crd, nprim_per_at, bsetBF, x, y, z, chi_at, &
                               prt_warn)
    !! Evaluate atomic orbitals at position.
    !!
    !! Evaluate atomic orbitals (chi) at a chosen position.
    !!
    !! @warning
    !! Only Cartesian basis sets are supported for now.
    !! @endwarning
    !!
    !! @note "version"
    !! This version uses arrays as arguments, except for bsetBF, which
    !! is a database of basis set functions.
    !! @endnote
    real(realwp), dimension(:,:), intent(in) :: at_crd
        !! Atomic coordinates.
    integer, dimension(:), intent(in) :: nprim_per_at
        !! Number of primitives on each atomic center.
    class(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
        !! Basis set's basis function information.
    real(realwp), intent(in) :: x, y, z
        !! Cartesian components of the point of interest.
    real(realwp), dimension(:), intent(inout) :: chi_at
        !! Atomic orbitals (chi) evaluated at chosen point.
    logical, intent(in), optional :: prt_warn
        !! Print warning messages if risk of divergences.

    integer :: beg_bset, beg_sh, end_bset, end_sh, i_sh, ia, n_prim, &
        n_prim_at, n_prim_sh
    real(realwp) :: x_rel, y_rel, z_rel
    logical :: warn
    character(len=256) :: msg

    if (present(prt_warn)) then
        warn = prt_warn
    else
        warn = .false.
    end if

    chi_at = f0
    beg_sh = 1
    do ia = 1, size(at_crd, 2)
        x_rel = (x - at_crd(1,ia))
        y_rel = (y - at_crd(2,ia))
        z_rel = (z - at_crd(3,ia))

        if (sum(abs([x_rel, y_rel, z_rel])) < thresh_at_center .and. warn) then
            write(msg, '("Chosen point close to atom center ",i0,&
                &". Divergence may occur!")') ia
            call runstat%raise_warning(trim(msg))
        end if

        n_prim_at = nprim_per_at(ia)

        i_sh = 1
        n_prim_sh = 1
        n_prim = 0
        beg_bset = 1
        do while (n_prim < n_prim_at)
            do while (.not.bsetBF(ia,i_sh)%shell_last)
                i_sh = i_sh + 1
                n_prim_sh = n_prim_sh + 1
            end do
            n_prim = n_prim + n_prim_sh
            end_bset = beg_bset + n_prim_sh - 1
            end_sh = beg_sh + bsetBF(ia,i_sh)%ndim - 1

            chi_at(beg_sh:end_sh) = get_AOs_sh_at( &
                bsetBF(ia,i_sh)%ndim, bsetBF(ia,beg_bset:end_bset), &
                x_rel, y_rel, z_rel)

            beg_sh = beg_sh + bsetBF(ia,i_sh)%ndim

            i_sh = i_sh + 1
            beg_bset = i_sh
            n_prim_sh = 1
        end do
        end_sh = end_sh + 1

    end do

end subroutine eval_AOs_chi_at_arr

! ======================================================================

subroutine eval_AOs_chi_at_db(molDB, bsetDB, x, y, z, chi_at, prt_warn)
    !! Evaluate atomic orbitals at position.
    !!
    !! Evaluate atomic orbitals (chi) at a chosen position.
    !!
    !! @warning
    !! Only Cartesian basis sets are supported for now.
    !! @endwarning
    !!
    !! @note "version"
    !! This version takes general databases as arguments.
    !! @endnote
    class(MoleculeDB), intent(in) :: molDB
        !! Molecule database.
    class(BasisSetDB), intent(in) :: bsetDB
        !! Basis set database
    real(realwp), intent(in) :: x, y, z
        !! Cartesian components of the point of interest.
    real(realwp), dimension(:), intent(inout) :: chi_at
        !! Atomic orbitals (chi) evaluated at chosen point.
    logical, intent(in), optional :: prt_warn
        !! Print warning messages if risk of divergences.

    integer :: beg_bset, beg_sh, end_bset, end_sh, i_sh, ia, n_prim, &
        n_prim_at, n_prim_sh
    real(realwp) :: x_rel, y_rel, z_rel
    logical :: warn
    character(len=256) :: msg

    if (present(prt_warn)) then
        warn = prt_warn
    else
        warn = .false.
    end if

    chi_at = f0
    beg_sh = 1
    do ia = 1, molDB%n_at
        x_rel = (x - molDB%at_crd(1,ia))
        y_rel = (y - molDB%at_crd(2,ia))
        z_rel = (z - molDB%at_crd(3,ia))

        if (sum(abs([x_rel, y_rel, z_rel])) < thresh_at_center .and. warn) then
            write(msg, '("Chosen point close to atom center ",i0,&
                &". Divergence may occur!")') ia
            call runstat%raise_warning(trim(msg))
        end if

        n_prim_at = bsetDB%nprim_per_at(ia)

        i_sh = 1
        n_prim_sh = 1
        n_prim = 0
        beg_bset = 1
        do while (n_prim < n_prim_at)
            do while (.not.bsetDB%info(ia,i_sh)%shell_last)
                i_sh = i_sh + 1
                n_prim_sh = n_prim_sh + 1
            end do
            n_prim = n_prim + n_prim_sh
            end_bset = beg_bset + n_prim_sh - 1
            end_sh = beg_sh + bsetDB%info(ia, i_sh)%ndim - 1

            chi_at(beg_sh:end_sh) = get_AOs_sh_at( &
                bsetDB%info(ia,i_sh)%ndim, bsetDB%info(ia,beg_bset:end_bset), &
                x_rel, y_rel, z_rel)

            beg_sh = beg_sh + bsetDB%info(ia, i_sh)%ndim

            i_sh = i_sh + 1
            beg_bset = i_sh
            n_prim_sh = 1
        end do
        end_sh = end_sh + 1

    end do

end subroutine eval_AOs_chi_at_db

! ======================================================================

subroutine eval_AOs_chi_at_dim(n_at, at_crd, nprim_per_at, bsetBF, x, y, z, &
                               chi_at, prt_warn)
    !! Evaluate atomic orbitals at position.
    !!
    !! Evaluate atomic orbitals (chi) at a chosen position.
    !!
    !! @warning
    !! Only Cartesian basis sets are supported for now.
    !! @endwarning
    !!
    !! @note "version"
    !! This version uses arrays as arguments, except for bsetBF, which
    !! is a database of basis set functions.
    !! Relevant dimensions are also provided
    !! @endnote
    integer, intent(in) :: n_at
        !! Number of atoms.
    real(realwp), dimension(:,:), intent(in) :: at_crd
        !! Atomic coordinates.
    integer, dimension(:), intent(in) :: nprim_per_at
        !! Number of primitives on each atomic center.
    class(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
        !! Basis set's basis function information.
    real(realwp), intent(in) :: x, y, z
        !! Cartesian components of the point of interest.
    real(realwp), dimension(:), intent(inout) :: chi_at
        !! Atomic orbitals (chi) evaluated at chosen point.
    logical, intent(in), optional :: prt_warn
        !! Print warning messages if risk of divergences.

    integer :: beg_bset, beg_sh, end_bset, end_sh, i_sh, ia, n_prim, &
        n_prim_at, n_prim_sh
    real(realwp) :: x_rel, y_rel, z_rel
    logical :: warn
    character(len=256) :: msg

    if (present(prt_warn)) then
        warn = prt_warn
    else
        warn = .false.
    end if

    chi_at = f0
    beg_sh = 1
    do ia = 1, n_at
        x_rel = (x - at_crd(1,ia))
        y_rel = (y - at_crd(2,ia))
        z_rel = (z - at_crd(3,ia))

        if (sum(abs([x_rel, y_rel, z_rel])) < thresh_at_center .and. warn) then
            write(msg, '("Chosen point close to atom center ",i0,&
                &". Divergence may occur!")') ia
            call runstat%raise_warning(trim(msg))
        end if

        n_prim_at = nprim_per_at(ia)

        i_sh = 1
        n_prim_sh = 1
        n_prim = 0
        beg_bset = 1
        do while (n_prim < n_prim_at)
            do while (.not.bsetBF(ia,i_sh)%shell_last)
                i_sh = i_sh + 1
                n_prim_sh = n_prim_sh + 1
            end do
            n_prim = n_prim + n_prim_sh
            end_bset = beg_bset + n_prim_sh - 1
            end_sh = beg_sh + bsetBF(ia,i_sh)%ndim - 1

            chi_at(beg_sh:end_sh) = get_AOs_sh_at( &
                bsetBF(ia,i_sh)%ndim, bsetBF(ia,beg_bset:end_bset), &
                x_rel, y_rel, z_rel)

            beg_sh = beg_sh + bsetBF(ia,i_sh)%ndim

            i_sh = i_sh + 1
            beg_bset = i_sh
            n_prim_sh = 1
        end do
        end_sh = end_sh + 1

    end do

end subroutine eval_AOs_chi_at_dim

! ======================================================================

subroutine eval_AOs_nabla_chi_at_arr(at_crd, nprim_per_at, bsetBF, x, y, z, &
                                     chi_at, d1_chi_at, prt_warn)
    !! Evaluate AOs and first derivatives at position.
    !!
    !! Evaluate atomic orbitals (chi) and their first derivative (nabla
    !! chi) at a chosen position.
    !!
    !! @warning
    !! Only Cartesian basis sets are supported for now.
    !! @endwarning
    !!
    !! @note "version"
    !! This version uses arrays as arguments, except for bsetBF, which
    !! is a database of basis set functions.
    !! @endnote
    real(realwp), dimension(:,:), intent(in) :: at_crd
        !! Atomic coordinates.
    integer, dimension(:), intent(in) :: nprim_per_at
        !! Number of primitives on each atomic center.
    class(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
        !! Basis set's basis function information.
    real(realwp), intent(in) :: x, y, z
        !! Cartesian components of the point of interest.
    real(realwp), dimension(:), intent(inout) :: chi_at
        !! Atomic orbitals (chi) evaluated at chosen point.
    real(realwp), dimension(:,:), intent(inout) :: d1_chi_at
        !! First derivatives of AOs (nabla chi) evaluated at chosen point.
    logical, intent(in), optional :: prt_warn
        !! Print warning messages if risk of divergences.

    integer :: beg_bset, beg_sh, end_bset, end_sh, i_sh, ia, n_prim, &
        n_prim_at, n_prim_sh
    real(realwp) :: x_rel, y_rel, z_rel
    logical :: warn
    character(len=256) :: msg

    if (present(prt_warn)) then
        warn = prt_warn
    else
        warn = .false.
    end if

    chi_at = f0
    d1_chi_at = f0
    beg_sh = 1
    do ia = 1, size(at_crd, 2)
        x_rel = (x - at_crd(1,ia))
        y_rel = (y - at_crd(2,ia))
        z_rel = (z - at_crd(3,ia))

        if (sum(abs([x_rel, y_rel, z_rel])) < thresh_at_center .and. warn) then
            write(msg, '("Chosen point close to atom center ",i0,&
                &". Divergence may occur!")') ia
            call runstat%raise_warning(trim(msg))
        end if

        n_prim_at = nprim_per_at(ia)

        i_sh = 1
        n_prim_sh = 1
        n_prim = 0
        beg_bset = 1
        do while (n_prim < n_prim_at)
            do while (.not.bsetBF(ia,i_sh)%shell_last)
                i_sh = i_sh + 1
                n_prim_sh = n_prim_sh + 1
            end do
            n_prim = n_prim + n_prim_sh
            end_bset = beg_bset + n_prim_sh - 1
            end_sh = beg_sh + bsetBF(ia,i_sh)%ndim - 1

            call get_AOs_d1_sh_at(bsetBF(ia,i_sh)%ndim, &
                                  bsetBF(ia,beg_bset:end_bset), &
                                  x_rel, y_rel, z_rel, &
                                  chi_at(beg_sh:end_sh), &
                                  d1_chi_at(:,beg_sh:end_sh))

            beg_sh = beg_sh + bsetBF(ia,i_sh)%ndim

            i_sh = i_sh + 1
            beg_bset = i_sh
            n_prim_sh = 1
        end do
        end_sh = end_sh + 1

    end do

end subroutine eval_AOs_nabla_chi_at_arr

! ======================================================================

subroutine eval_AOs_nabla_chi_at_db(molDB, bsetDB, x, y, z, chi_at, &
                                    d1_chi_at, prt_warn)
    !! Evaluate AOs and first derivatives at position.
    !!
    !! Evaluate atomic orbitals (chi) and their first derivative (nabla
    !! chi) at a chosen position.
    !!
    !! @warning
    !! Only Cartesian basis sets are supported for now.
    !! @endwarning
    !!
    !! @note "version"
    !! This version takes general databases as arguments.
    !! @endnote
    class(MoleculeDB), intent(in) :: molDB
        !! Molecule database.
    class(BasisSetDB), intent(in) :: bsetDB
        !! Basis set database
    real(realwp), intent(in) :: x, y, z
        !! Cartesian components of the point of interest.
    real(realwp), dimension(:), intent(inout) :: chi_at
        !! Atomic orbitals (chi) evaluated at chosen point.
    real(realwp), dimension(:,:), intent(inout) :: d1_chi_at
        !! First derivatives of AOs (nabla chi) evaluated at chosen point.
    logical, intent(in), optional :: prt_warn
        !! Print warning messages if risk of divergences.

    integer :: beg_bset, beg_sh, end_bset, end_sh, i_sh, ia, n_prim, &
        n_prim_at, n_prim_sh
    real(realwp) :: x_rel, y_rel, z_rel
    logical :: warn
    character(len=256) :: msg

    if (present(prt_warn)) then
        warn = prt_warn
    else
        warn = .false.
    end if

    chi_at = f0
    d1_chi_at = f0
    beg_sh = 1
    do ia = 1, molDB%n_at
        x_rel = (x - molDB%at_crd(1,ia))
        y_rel = (y - molDB%at_crd(2,ia))
        z_rel = (z - molDB%at_crd(3,ia))

        if (sum(abs([x_rel, y_rel, z_rel])) < thresh_at_center .and. warn) then
            write(msg, '("Chosen point close to atom center ",i0,&
                &". Divergence may occur!")') ia
            call runstat%raise_warning(trim(msg))
        end if

        n_prim_at = bsetDB%nprim_per_at(ia)

        i_sh = 1
        n_prim_sh = 1
        n_prim = 0
        beg_bset = 1
        do while (n_prim < n_prim_at)
            do while (.not.bsetDB%info(ia,i_sh)%shell_last)
                i_sh = i_sh + 1
                n_prim_sh = n_prim_sh + 1
            end do
            n_prim = n_prim + n_prim_sh
            end_bset = beg_bset + n_prim_sh - 1
            end_sh = beg_sh + bsetDB%info(ia, i_sh)%ndim - 1

            call get_AOs_d1_sh_at(bsetDB%info(ia,i_sh)%ndim, &
                                  bsetDB%info(ia,beg_bset:end_bset), &
                                  x_rel, y_rel, z_rel, &
                                  chi_at(beg_sh:end_sh), &
                                  d1_chi_at(:,beg_sh:end_sh))

            beg_sh = beg_sh + bsetDB%info(ia,i_sh)%ndim

            i_sh = i_sh + 1
            beg_bset = i_sh
            n_prim_sh = 1
        end do
        end_sh = end_sh + 1

    end do

end subroutine eval_AOs_nabla_chi_at_db

! ======================================================================

subroutine eval_AOs_nabla_chi_at_dim(n_at, at_crd, nprim_per_at, bsetBF, x, &
                                     y, z, chi_at, d1_chi_at, prt_warn)
    !! Evaluate AOs and first derivatives at position.
    !!
    !! Evaluate atomic orbitals (chi) and their first derivative (nabla
    !! chi) at a chosen position.
    !!
    !! @warning
    !! Only Cartesian basis sets are supported for now.
    !! @endwarning
    !!
    !! @note "version"
    !! This version uses arrays as arguments, except for bsetBF, which
    !! is a database of basis set functions.
    !! Relevant dimensions are also provided
    !! @endnote
    integer, intent(in) :: n_at
        !! Number of atoms.
    real(realwp), dimension(3,n_at), intent(in) :: at_crd
        !! Atomic coordinates.
    integer, dimension(:), intent(in) :: nprim_per_at
        !! Number of primitives on each atomic center.
    class(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
        !! Basis set's basis function information.
    real(realwp), intent(in) :: x, y, z
        !! Cartesian components of the point of interest.
    real(realwp), dimension(:), intent(inout) :: chi_at
        !! Atomic orbitals (chi) evaluated at chosen point.
    real(realwp), dimension(:,:), intent(inout) :: d1_chi_at
        !! First derivatives of AOs (nabla chi) evaluated at chosen point.
    logical, intent(in), optional :: prt_warn
        !! Print warning messages if risk of divergences.

    integer :: beg_bset, beg_sh, end_bset, end_sh, i_sh, ia, n_prim, &
        n_prim_at, n_prim_sh
    real(realwp) :: x_rel, y_rel, z_rel
    logical :: warn
    character(len=256) :: msg

    if (present(prt_warn)) then
        warn = prt_warn
    else
        warn = .false.
    end if

    chi_at = f0
    d1_chi_at = f0
    beg_sh = 1
    do ia = 1, n_at
        x_rel = (x - at_crd(1,ia))
        y_rel = (y - at_crd(2,ia))
        z_rel = (z - at_crd(3,ia))

        if (sum(abs([x_rel, y_rel, z_rel])) < thresh_at_center .and. warn) then
            write(msg, '("Chosen point close to atom center ",i0,&
                &". Divergence may occur!")') ia
            call runstat%raise_warning(trim(msg))
        end if

        n_prim_at = nprim_per_at(ia)

        i_sh = 1
        n_prim_sh = 1
        n_prim = 0
        beg_bset = 1
        do while (n_prim < n_prim_at)
            do while (.not.bsetBF(ia,i_sh)%shell_last)
                i_sh = i_sh + 1
                n_prim_sh = n_prim_sh + 1
            end do
            n_prim = n_prim + n_prim_sh
            end_bset = beg_bset + n_prim_sh - 1
            end_sh = beg_sh + bsetBF(ia,i_sh)%ndim - 1

            call get_AOs_d1_sh_at(bsetBF(ia,i_sh)%ndim, &
                                  bsetBF(ia,beg_bset:end_bset), &
                                  x_rel, y_rel, z_rel, &
                                  chi_at(beg_sh:end_sh), &
                                  d1_chi_at(:,beg_sh:end_sh))

            beg_sh = beg_sh + bsetBF(ia,i_sh)%ndim

            i_sh = i_sh + 1
            beg_bset = i_sh
            n_prim_sh = 1
        end do
        end_sh = end_sh + 1

    end do

end subroutine eval_AOs_nabla_chi_at_dim

! ======================================================================

subroutine get_AOs_d1_sh_at(len_sh, bset_sh, x, y, z, chi_sh_at, &
                            d1_chi_sh_at)
    !! Get AOs and first derivatives for a shell at position.
    !!
    !! Computes and returns the atomic orbitals (chi) and their
    !! first derivative (nabla chi) for a given shell at a chosen
    !! position.
    integer, intent(in) :: len_sh
        !! Length of the shell.
    type(PrimitiveFunction), dimension(:), intent(in) :: bset_sh
        !! Basis set functions database for a given shell.
    real(realwp), intent(in) :: x, y, z
        !! Cartesian components of the point of interest.
    real(realwp), dimension(len_sh), intent(out) :: chi_sh_at
        !! Atomic orbitals for the shell at chosen point.
    real(realwp), dimension(3,len_sh), intent(out) :: d1_chi_sh_at
        !! First derivatives of AOs for the shell at chosen point.

    integer :: idx_coeff, L_ang, len_sub_sh, n_prim_sh
    real(realwp) :: rad
    real(realwp), dimension(3) :: d1_rad
    real(realwp), dimension(len_sh) :: ang
    real(realwp), dimension(len_sh, 3) :: d1_ang
    real(realwp), dimension(len_sh) :: anorm
    real(realwp), dimension(:), allocatable :: anorm_SP

    idx_coeff = 0
    n_prim_sh = size(bset_sh)

    if (bset_sh(1)%shelltype == 'SP') then
        ! Treat S part
        L_ang = 0
        len_sub_sh = 1
        rad = get_cart_r_sh_at(L_ang, idx_coeff, n_prim_sh, bset_sh, x, y, z)
        d1_rad = get_cart_r_der_sh_at(L_ang, idx_coeff, n_prim_sh, bset_sh, &
                                      x, y, z)

        allocate(anorm_SP(len_sub_sh))
        anorm_SP = get_cart_L_norms_sh(L_ang, len_sub_sh)

        ang(1:1) = get_cart_L_sh_at(L_ang, len_sub_sh, x, y, z)
        d1_ang(1:1,:) = get_cart_L_der_sh_at(L_ang, len_sub_sh, x, y, z)
        chi_sh_at(1:1) = rad * ang(1:1) / anorm_SP
        d1_chi_sh_at(1,1:1) = &
            (rad * d1_ang(1:1,1) + d1_rad(1) * ang(1:1)) / anorm_SP
        d1_chi_sh_at(2,1:1) = &
            (rad * d1_ang(1:1,2) + d1_rad(2) * ang(1:1)) / anorm_SP
        d1_chi_sh_at(3,1:1) = &
            (rad * d1_ang(1:1,3) + d1_rad(3) * ang(1:1)) / anorm_SP
        deallocate(anorm_SP)

        ! P part
        L_ang = 1
        idx_coeff = -1
        len_sub_sh = 3

        rad = get_cart_r_sh_at(L_ang, idx_coeff, n_prim_sh, bset_sh, x, y, z)
        d1_rad = get_cart_r_der_sh_at(L_ang, idx_coeff, n_prim_sh, bset_sh, &
                                      x, y, z)

        allocate(anorm_SP(len_sub_sh))
        anorm_SP = get_cart_L_norms_sh(L_ang, len_sub_sh)

        ang(2:4) = get_cart_L_sh_at(L_ang, len_sub_sh, x, y, z)
        d1_ang(2:4,:) = get_cart_L_der_sh_at(L_ang, len_sub_sh, x, y, z)
        chi_sh_at(2:4) = rad * ang(2:4) / anorm_SP
        d1_chi_sh_at(1,2:4) = &
            (rad * d1_ang(2:4,1) + d1_rad(1) * ang(2:4)) / anorm_SP
        d1_chi_sh_at(2,2:4) = &
            (rad * d1_ang(2:4,2) + d1_rad(2) * ang(2:4)) / anorm_SP
        d1_chi_sh_at(3,2:4) = &
            (rad * d1_ang(2:4,3) + d1_rad(3) * ang(2:4)) / anorm_SP
        deallocate(anorm_SP)

    else
        rad = get_cart_r_sh_at(bset_sh(1)%L, idx_coeff, n_prim_sh, bset_sh, &
                               x, y, z)
        d1_rad = get_cart_r_der_sh_at(bset_sh(1)%L, idx_coeff, n_prim_sh, &
                                      bset_sh, x, y, z)
        anorm = get_cart_L_norms_sh(bset_sh(1)%L, len_sh)
        ang = get_cart_L_sh_at(bset_sh(1)%L, len_sh, x, y, z)
        d1_ang = get_cart_L_der_sh_at(bset_sh(1)%L, len_sh, x, y, z)
        chi_sh_at = rad * ang / anorm
        d1_chi_sh_at(1,:) = (rad * d1_ang(:,1) + d1_rad(1) * ang) / anorm
        d1_chi_sh_at(2,:) = (rad * d1_ang(:,2) + d1_rad(2) * ang) / anorm
        d1_chi_sh_at(3,:) = (rad * d1_ang(:,3) + d1_rad(3) * ang ) / anorm

    end if

end subroutine get_AOs_d1_sh_at

! ======================================================================

function get_AOs_sh_at(len_sh, bset_sh, x, y, z) result(chi_sh_at)
    !! Get atomic orbitals for a shell at position.
    !!
    !! Computes and returns the atomic orbitals (chi) for a given &
    !! shell at a chosen position.
    integer, intent(in) :: len_sh
        !! Length of the shell.
    type(PrimitiveFunction), dimension(:), intent(in) :: bset_sh
        !! Basis set functions database for a given shell.
    real(realwp), intent(in) :: x, y, z
        !! Cartesian components of the point of interest.
    real(realwp), dimension(len_sh) :: chi_sh_at
        !! Atomic orbitals for the shell at chosen point.

    integer :: idx_coeff, L_ang, len_sub_sh, nprim_sh
    real(realwp) :: rad
    real(realwp), dimension(len_sh) :: ang
    real(realwp), dimension(len_sh) :: anorm
    real(realwp), dimension(:), allocatable :: anorm_SP

    chi_sh_at = f0
    nprim_sh = size(bset_sh)
    idx_coeff = 0

    if (bset_sh(1)%shelltype == 'SP') then
        ! Treat S part
        L_ang = 0
        len_sub_sh = 1
        rad = get_cart_r_sh_at(L_ang, idx_coeff, nprim_sh, bset_sh, x, y, z)
        ang(1:1) = get_cart_L_sh_at(L_ang, len_sub_sh, x, y, z)

        allocate(anorm_SP(len_sub_sh))
        anorm_SP = get_cart_L_norms_sh(L_ang, len_sub_sh)
        chi_sh_at(1:1) =  rad * ang(1:1) / anorm_SP
        deallocate(anorm_SP)

        ! P part
        L_ang = 1
        idx_coeff = -1
        len_sub_sh = 3
        rad = get_cart_r_sh_at(L_ang, idx_coeff, nprim_sh, bset_sh, x, y, z)
        ang(2:4) = get_cart_L_sh_at(L_ang, len_sub_sh, x, y, z)

        allocate(anorm_SP(len_sub_sh))
        anorm_SP = get_cart_L_norms_sh(L_ang, len_sub_sh)
        chi_sh_at(2:4) =  rad * ang(2:4) / anorm_SP
        deallocate(anorm_SP)

    else
        rad = get_cart_r_sh_at(bset_sh(1)%L, idx_coeff, nprim_sh, bset_sh, &
                               x, y, z)
        ang = get_cart_L_sh_at(bset_sh(1)%L, len_sh, x, y, z)
        anorm = get_cart_L_norms_sh(bset_sh(1)%L, len_sh)
        chi_sh_at = rad * ang / anorm
    end if

end function get_AOs_sh_at

! ======================================================================

end module orbital