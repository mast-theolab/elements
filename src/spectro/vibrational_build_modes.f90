submodule (vibrational) vib_build_modes

    implicit none

contains

! ======================================================================

module procedure build_modes_arr_lt

    integer :: ia, ioff, ja, n_at, n_at3, n_at3tt
    real(realwp) :: sqmas_i, sqmas_j
    real(realwp), dimension(:,:), allocatable :: F_mweigh
    logical :: do_weigh

    n_at = size(at_crd, 2)
    n_at3 = 3*n_at
    n_at3tt = n_at3*(n_at3+1)/2
    if (size(F_cart) /= n_at3tt .or. size(at_mass) /= n_at) then
        call run%error%raise_argerror('size', &
            'inconsistency in input arrays', source='build_modes')
        return
    end if

    if (present(is_weighted)) then
        do_weigh = .not.is_weighted
    else
        do_weigh = .true.
    end if
    allocate(F_mweigh(n_at3,n_at3))
    if (do_weigh) then
        do ia = 1, n_at3
            sqmas_i = f1/sqrt(at_mass((ia+2)/3))
            ioff = ia*(ia-1)/2
            do ja = 1, ia-1
                sqmas_j = f1/sqrt(at_mass((ja+2)/3))
                F_mweigh(ja,ia) = F_cart(ioff+ja)*sqmas_i*sqmas_j
                F_mweigh(ia,ja) = F_mweigh(ja,ia)
            end do
            F_mweigh(ia,ia) = F_cart(ioff+ia)*sqmas_i**2
        end do
    else
        do ia = 1, n_at3
            ioff = ia*(ia-1)/2
            do ja = 1, ia-1
                F_mweigh(ja,ia) = F_cart(ioff+ja)
                F_mweigh(ia,ja) = F_mweigh(ja,ia)
            end do
            F_mweigh(ia,ia) = F_cart(ioff+ia)
        end do
    end if

    call build_modes_algo(n_at, F_mweigh, at_crd, at_mass, remove_rottrans, &
                          n_vib, L_mat, freq, L_mweigh, red_mass)
    deallocate(F_mweigh)

end procedure build_modes_arr_lt

! ======================================================================

module procedure build_modes_arr_sq

    integer :: ia, ja, n_at, n_at3
    real(realwp) :: sqmas_i, sqmas_j
    real(realwp), dimension(:,:), allocatable :: F_mweigh
    logical :: do_weigh

    n_at = size(at_crd, 2)
    n_at3 = 3*n_at
    if (size(F_mweigh, 1) /= n_at3 .or. size(at_mass, 1) /= n_at) then
        call run%error%raise_argerror('size', &
            'inconsistency in input arrays', source='build_modes')
        return
    end if

    if (present(is_weighted)) then
        do_weigh = .not.is_weighted
    else
        do_weigh = .true.
    end if
    allocate(F_mweigh(n_at3,n_at3))
    if (do_weigh) then
        do ia = 1, n_at3
            sqmas_i = f1/sqrt(at_mass((ia+2)/3))
            do ja = 1, ia-1
                sqmas_j = f1/sqrt(at_mass((ja+2)/3))
                F_mweigh(ja,ia) = F_cart(ja,ia)*sqmas_i*sqmas_j
                F_mweigh(ia,ja) = F_mweigh(ja,ia)
            end do
            F_mweigh(ia,ia) = F_cart(ia,ia)*sqmas_i**2
        end do
    else
        F_mweigh = F_cart
    end if

    call build_modes_algo(n_at, F_mweigh, at_crd, at_mass, remove_rottrans, &
                          n_vib, L_mat, freq, L_mweigh, red_mass)
    deallocate(F_mweigh)

end procedure build_modes_arr_sq

! ======================================================================

module procedure build_modes_db_lt

    integer :: ia, ioff, ja, n_at3, n_at3tt, nvib
    real(realwp) :: sqmas_i, sqmas_j
    real(realwp), dimension(:,:), allocatable :: F_mweigh
    logical :: do_freq, do_Lmat, do_Lwgt, do_nvib, do_rmas, do_weigh

    n_at3 = 3*molDB%n_at
    n_at3tt = n_at3*(n_at3+1)/2
    if (size(F_cart) /= n_at3tt) then
        call run%error%raise_argerror('size', &
            'inconsistency in input arrays', source='build_modes')
        return
    end if

    if (present(is_weighted)) then
        do_weigh = .not.is_weighted
    else
        do_weigh = .true.
    end if
    allocate(F_mweigh(n_at3,n_at3))
    if (do_weigh) then
        do ia = 1, n_at3
            sqmas_i = f1/sqrt(molDB%at_mas((ia+2)/3))
            ioff = ia*(ia-1)/2
            do ja = 1, ia-1
                sqmas_j = f1/sqrt(molDB%at_mas((ja+2)/3))
                F_mweigh(ja,ia) = F_cart(ioff+ja)*sqmas_i*sqmas_j
                F_mweigh(ia,ja) = F_mweigh(ja,ia)
            end do
            F_mweigh(ia,ia) = F_cart(ioff+ia)*sqmas_i**2
        end do
    else
        do ia = 1, n_at3
            ioff = ia*(ia-1)/2
            do ja = 1, ia-1
                F_mweigh(ja,ia) = F_cart(ioff+ja)
                F_mweigh(ia,ja) = F_mweigh(ja,ia)
            end do
            F_mweigh(ia,ia) = F_cart(ioff+ia)
        end do
    end if

    ! Deallocate the content of `vib` if already loaded
    vibDB%n_vib = 0
    if (vibDB%loaded) then
        if (allocated(vibDB%freq)) deallocate(vibDB%freq)
        if (allocated(vibDB%L_mwg)) deallocate(vibDB%L_mwg)
        if (allocated(vibDB%L_mat)) deallocate(vibDB%L_mat)
        if (allocated(vibDB%red_mass)) deallocate(vibDB%red_mass)
    end if

    ! Check which information to provide
    if (present(set_nvib)) then
        do_nvib = set_nvib
    else
        do_nvib = .true.
    end if
    if (present(set_Lmat)) then
        do_Lmat = set_Lmat
    else
        do_Lmat = .true.
    end if
    if (present(set_freq)) then
        do_freq = set_freq
    else
        do_freq = .true.
    end if
    if (present(set_Lmweigh)) then
        do_Lwgt = set_Lmweigh
    else
        do_Lwgt = .true.
    end if
    if (present(set_redmas)) then
        do_rmas = set_redmas
    else
        do_rmas = .true.
    end if

    ! now do allocation
    if (do_Lmat) allocate(vibDB%L_mat(n_at3,n_at3))
    if (do_freq) then
        allocate(vibDB%freq(n_at3))
        allocate(vibDB%red_freq(n_at3))
    end if
    if (do_Lwgt) allocate(vibDB%L_mwg(n_at3,n_at3))
    if (do_rmas) allocate(vibDB%red_mass(n_at3))
    nvib = 0

    if (do_Lmat.and.do_freq.and.do_Lwgt.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, freq=vibDB%freq, &
            L_mweigh=vibDB%L_mwg, red_mass=vibDB%red_mass)
    else if (do_Lmat.and.do_freq.and.do_Lwgt) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, freq=vibDB%freq, &
            L_mweigh=vibDB%L_mwg)
    else if (do_Lmat.and.do_freq.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, freq=vibDB%freq, &
            red_mass=vibDB%red_mass)
    else if (do_Lmat.and.do_Lwgt.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, L_mweigh=vibDB%L_mwg, &
            red_mass=vibDB%red_mass)
    else if (do_freq.and.do_Lwgt.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, freq=vibDB%freq, L_mweigh=vibDB%L_mwg, &
            red_mass=vibDB%red_mass)
    else if (do_Lmat.and.do_freq) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, freq=vibDB%freq)
    else if (do_Lmat.and.do_Lwgt) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, L_mweigh=vibDB%L_mwg)
    else if (do_Lmat.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, &
            red_mass=vibDB%red_mass)
    else if (do_freq.and.do_Lwgt) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, freq=vibDB%freq, L_mweigh=vibDB%L_mwg)
    else if (do_freq.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, freq=vibDB%freq, red_mass=vibDB%red_mass)
    else if (do_Lwgt.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mweigh=vibDB%L_mwg, &
            red_mass=vibDB%red_mass)
    else if (do_Lmat) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat)
    else if (do_freq) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, freq=vibDB%freq)
    else if (do_Lwgt) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mweigh=vibDB%L_mwg)
    else if (do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, red_mass=vibDB%red_mass)
    end if
    deallocate(F_mweigh)
    if (run%error%is_ok()) then
        if (do_Lmat) vibDB%L_mat = vibDB%L_mat(:n_at3,:nvib)
        if (do_freq) then
            vibDB%freq = vibDB%freq(:nvib)
            vibDB%red_freq = phys_conv%au2cm1(vibDB%freq(:nvib), .true.)
        end if
        if (do_Lwgt) vibDB%L_mwg = vibDB%L_mwg(:n_at3,:nvib)
        if (do_rmas) vibDB%red_mass = vibDB%red_mass(:nvib)
        if (do_nvib) vibDB%n_vib = nvib
        vibDB%loaded = .true.
    end if

end procedure build_modes_db_lt

! ======================================================================

module procedure build_modes_db_sq

    integer :: ia, ja, n_at3, nvib
    real(realwp) :: sqmas_i, sqmas_j
    real(realwp), dimension(:,:), allocatable :: F_mweigh
    logical :: do_freq, do_Lmat, do_Lwgt, do_nvib, do_rmas, do_weigh

    n_at3 = 3*molDB%n_at
    if (size(F_cart, 1) /= n_at3) then
        call run%error%raise_argerror('size', &
            'inconsistency in input arrays', source='build_modes')
        return
    end if

    if (present(is_weighted)) then
        do_weigh = .not.is_weighted
    else
        do_weigh = .true.
    end if
    allocate(F_mweigh(n_at3,n_at3))
    if (do_weigh) then
        do ia = 1, n_at3
            sqmas_i = f1/sqrt(molDB%at_mas((ia+2)/3))
            do ja = 1, ia-1
                sqmas_j = f1/sqrt(molDB%at_mas((ja+2)/3))
                F_mweigh(ja,ia) = F_cart(ja,ia)*sqmas_i*sqmas_j
                F_mweigh(ia,ja) = F_mweigh(ja,ia)
            end do
            F_mweigh(ia,ia) = F_cart(ia,ia)*sqmas_i**2
        end do
    else
        F_mweigh = F_cart
    end if

    ! Deallocate the content of `vib` if already loaded
    vibDB%n_vib = 0
    if (vibDB%loaded) then
        if (allocated(vibDB%freq)) deallocate(vibDB%freq)
        if (allocated(vibDB%L_mwg)) deallocate(vibDB%L_mwg)
        if (allocated(vibDB%L_mat)) deallocate(vibDB%L_mat)
        if (allocated(vibDB%red_mass)) deallocate(vibDB%red_mass)
    end if

    ! Check which information to provide
    if (present(set_nvib)) then
        do_nvib = set_nvib
    else
        do_nvib = .true.
    end if
    if (present(set_Lmat)) then
        do_Lmat = set_Lmat
    else
        do_Lmat = .true.
    end if
    if (present(set_freq)) then
        do_freq = set_freq
    else
        do_freq = .true.
    end if
    if (present(set_Lmweigh)) then
        do_Lwgt = set_Lmweigh
    else
        do_Lwgt = .true.
    end if
    if (present(set_redmas)) then
        do_rmas = set_redmas
    else
        do_rmas = .true.
    end if

    ! now do allocation
    if (do_Lmat) allocate(vibDB%L_mat(n_at3,n_at3))
    if (do_freq) then
        allocate(vibDB%freq(n_at3))
        allocate(vibDB%red_freq(n_at3))
    end if
    if (do_Lwgt) allocate(vibDB%L_mwg(n_at3,n_at3))
    if (do_rmas) allocate(vibDB%red_mass(n_at3))
    nvib = 0

    if (do_Lmat.and.do_freq.and.do_Lwgt.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, freq=vibDB%freq, &
            L_mweigh=vibDB%L_mwg, red_mass=vibDB%red_mass)
    else if (do_Lmat.and.do_freq.and.do_Lwgt) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, freq=vibDB%freq, &
            L_mweigh=vibDB%L_mwg)
    else if (do_Lmat.and.do_freq.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, freq=vibDB%freq, &
            red_mass=vibDB%red_mass)
    else if (do_Lmat.and.do_Lwgt.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, L_mweigh=vibDB%L_mwg, &
            red_mass=vibDB%red_mass)
    else if (do_freq.and.do_Lwgt.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, freq=vibDB%freq, L_mweigh=vibDB%L_mwg, &
            red_mass=vibDB%red_mass)
    else if (do_Lmat.and.do_freq) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, freq=vibDB%freq)
    else if (do_Lmat.and.do_Lwgt) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, L_mweigh=vibDB%L_mwg)
    else if (do_Lmat.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat, &
            red_mass=vibDB%red_mass)
    else if (do_freq.and.do_Lwgt) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, freq=vibDB%freq, L_mweigh=vibDB%L_mwg)
    else if (do_freq.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, freq=vibDB%freq, red_mass=vibDB%red_mass)
    else if (do_Lwgt.and.do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mweigh=vibDB%L_mwg, &
            red_mass=vibDB%red_mass)
    else if (do_Lmat) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mat=vibDB%L_mat)
    else if (do_freq) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, freq=vibDB%freq)
    else if (do_Lwgt) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, L_mweigh=vibDB%L_mwg)
    else if (do_rmas) then
        call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
            remove_rottrans, n_vib=nvib, red_mass=vibDB%red_mass)
    end if
    deallocate(F_mweigh)
    if (run%error%is_ok()) then
        if (do_Lmat) vibDB%L_mat = vibDB%L_mat(:n_at3,:nvib)
        if (do_freq) then
            vibDB%freq = vibDB%freq(:nvib)
            vibDB%red_freq = phys_conv%au2cm1(vibDB%freq(:nvib), .true.)
        end if
        if (do_Lwgt) vibDB%L_mwg = vibDB%L_mwg(:n_at3,:nvib)
        if (do_rmas) vibDB%red_mass = vibDB%red_mass(:nvib)
        if (do_nvib) vibDB%n_vib = nvib
        vibDB%loaded = .true.
    end if

end procedure build_modes_db_sq

! ======================================================================

module procedure build_modes_dim_lt

    integer :: ia, ioff, ja, n_at3, n_at3tt
    real(realwp) :: sqmas_i, sqmas_j
    real(realwp), dimension(:,:), allocatable :: F_mweigh
    logical :: do_weigh

    n_at3 = 3*n_at
    n_at3tt = n_at3*(n_at3+1)/2
    if (size(F_cart) /= n_at3tt .or. size(at_mass, 1) /= n_at .or. &
            size(at_crd, 2) /= n_at) then
        call run%error%raise_argerror('size', &
            'inconsistency in input arrays', source='build_modes')
        return
    end if

    if (present(is_weighted)) then
        do_weigh = .not.is_weighted
    else
        do_weigh = .true.
    end if
    allocate(F_mweigh(n_at3,n_at3))
    if (do_weigh) then
        do ia = 1, n_at3
            sqmas_i = f1/sqrt(at_mass((ia+2)/3))
            ioff = ia*(ia-1)/2
            do ja = 1, ia-1
                sqmas_j = f1/sqrt(at_mass((ja+2)/3))
                F_mweigh(ja,ia) = F_cart(ioff+ja)*sqmas_i*sqmas_j
                F_mweigh(ia,ja) = F_mweigh(ja,ia)
            end do
            F_mweigh(ia,ia) = F_cart(ioff+ia)*sqmas_i**2
        end do
    else
        do ia = 1, n_at3
            ioff = ia*(ia-1)/2
            do ja = 1, ia-1
                F_mweigh(ja,ia) = F_cart(ioff+ja)
                F_mweigh(ia,ja) = F_mweigh(ja,ia)
            end do
            F_mweigh(ia,ia) = F_cart(ioff+ia)
        end do
    end if

    call build_modes_algo(n_at, F_mweigh, at_crd, at_mass, remove_rottrans, &
                          n_vib, L_mat, freq, L_mweigh, red_mass)
    deallocate(F_mweigh)

end procedure build_modes_dim_lt

! ======================================================================

module procedure build_modes_dim_sq

    integer :: ia, ja, n_at3
    real(realwp) :: sqmas_i, sqmas_j
    real(realwp), dimension(:,:), allocatable :: F_mweigh
    logical :: do_weigh

    n_at3 = 3*n_at
    if (size(F_cart, 1) /= n_at3 .or. size(at_mass, 1) /= n_at .or. &
            size(at_crd, 2) /= n_at) then
        call run%error%raise_argerror('size', &
            'inconsistency in input arrays', source='build_modes')
        return
    end if

    if (present(is_weighted)) then
        do_weigh = .not.is_weighted
    else
        do_weigh = .true.
    end if
    allocate(F_mweigh(n_at3,n_at3))
    if (do_weigh) then
        do ia = 1, n_at3
            sqmas_i = f1/sqrt(at_mass((ia+2)/3))
            do ja = 1, ia-1
                sqmas_j = f1/sqrt(at_mass((ja+2)/3))
                F_mweigh(ja,ia) = F_cart(ja,ia)*sqmas_i*sqmas_j
                F_mweigh(ia,ja) = F_mweigh(ja,ia)
            end do
            F_mweigh(ia,ia) = F_cart(ia,ia)*sqmas_i**2
        end do
    else
        F_mweigh = F_cart
    end if

    call build_modes_algo(n_at, F_mweigh, at_crd, at_mass, remove_rottrans, &
                          n_vib, L_mat, freq, L_mweigh, red_mass)
    deallocate(F_mweigh)

end procedure build_modes_dim_sq

! ======================================================================

module procedure build_modes_moldb_lt

    integer :: ia, ioff, ja, n_at3, n_at3tt
    real(realwp) :: sqmas_i, sqmas_j
    real(realwp), dimension(:,:), allocatable :: F_mweigh
    logical :: do_weigh

    n_at3 = 3*molDB%n_at
    n_at3tt = n_at3*(n_at3+1)/2
    if (size(F_cart) /= n_at3tt) then
        call run%error%raise_argerror('size', &
            'inconsistency in input arrays', source='build_modes')
        return
    end if

    if (present(is_weighted)) then
        do_weigh = .not.is_weighted
    else
        do_weigh = .true.
    end if
    allocate(F_mweigh(n_at3,n_at3))
    if (do_weigh) then
        do ia = 1, n_at3
            sqmas_i = f1/sqrt(molDB%at_mas((ia+2)/3))
            ioff = ia*(ia-1)/2
            do ja = 1, ia-1
                sqmas_j = f1/sqrt(molDB%at_mas((ja+2)/3))
                F_mweigh(ja,ia) = F_cart(ioff+ja)*sqmas_i*sqmas_j
                F_mweigh(ia,ja) = F_mweigh(ja,ia)
            end do
            F_mweigh(ia,ia) = F_cart(ioff+ia)*sqmas_i**2
        end do
    else
        do ia = 1, n_at3
            ioff = ia*(ia-1)/2
            do ja = 1, ia-1
                F_mweigh(ja,ia) = F_cart(ioff+ja)
                F_mweigh(ia,ja) = F_mweigh(ja,ia)
            end do
            F_mweigh(ia,ia) = F_cart(ioff+ia)
        end do
    end if

    call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
                          remove_rottrans, n_vib, L_mat, freq, L_mweigh, &
                          red_mass)
    deallocate(F_mweigh)

end procedure build_modes_moldb_lt

! ======================================================================

module procedure build_modes_moldb_sq

    integer :: ia, ja, n_at3
    real(realwp) :: sqmas_i, sqmas_j
    real(realwp), dimension(:,:), allocatable :: F_mweigh
    logical :: do_weigh

    n_at3 = 3*molDB%n_at
    if (size(F_cart, 1) /= n_at3) then
        call run%error%raise_argerror('size', &
            'inconsistency in input arrays', source='build_modes')
        return
    end if

    if (present(is_weighted)) then
        do_weigh = .not.is_weighted
    else
        do_weigh = .true.
    end if
    allocate(F_mweigh(n_at3,n_at3))
    if (do_weigh) then
        do ia = 1, n_at3
            sqmas_i = f1/sqrt(molDB%at_mas((ia+2)/3))
            do ja = 1, ia-1
                sqmas_j = f1/sqrt(molDB%at_mas((ja+2)/3))
                F_mweigh(ja,ia) = F_cart(ja,ia)*sqmas_i*sqmas_j
                F_mweigh(ia,ja) = F_mweigh(ja,ia)
            end do
            F_mweigh(ia,ia) = F_cart(ia,ia)*sqmas_i**2
        end do
    else
        F_mweigh = F_cart
    end if

    call build_modes_algo(molDB%n_at, F_mweigh, molDB%at_crd, molDB%at_mas, &
                          remove_rottrans, n_vib, L_mat, freq, L_mweigh, &
                          red_mass)
    deallocate(F_mweigh)

end procedure build_modes_moldb_sq

! ======================================================================

subroutine build_modes_algo(n_at, F_mweigh, at_crd, at_mass, remove_rottrans, &
                            n_vib, L_mat, freq, L_mweigh, red_mass)
    !! Build the normal modes from the Cartesian force constant matrix.
    !!
    !! This is the core algorithm to build the normal modes
    !! (eigenvectors and eigenvalues) by diagonalizing the mass-weighted
    !! force constants matrix in `F_mweigh`.
    !! Quantities to return are simply provided as arguments, and must
    !! be allocated before.
    !! The removal of residual rotational and translational components
    !! is done as described in the white paper of Gaussian.
    !!
    !! @note
    !! The algorithm assumes that basic tests on data consistency have
    !! been done before.
    !! @endnote
    integer, intent(in) :: n_at
    !! Number of atoms.
    real(realwp), dimension(3*n_at,3*n_at), intent(in) :: F_mweigh
    !! Square force constant matrix, mass_weighted
    real(realwp), dimension(3,n_at), intent(in) :: at_crd
    !! Atomic coordinates.
    real(realwp), dimension(n_at), intent(in) :: at_mass
    !! Atomic masses.
    logical, intent(in), optional :: remove_rottrans
    !! Remove residual rotations and translations.
    integer, intent(inout), optional :: n_vib
    !! Number of normal modes.
    !! If set in input and positive, the routine checks that its results are
    !!   consistent.
    real(realwp), dimension(:,:), intent(out), optional :: L_mat
    !! Dimensionless eigenvectors matrix.
    real(realwp), dimension(:), intent(out), optional :: freq
    !! Wavenumbers (in cm-1).
    !! Note that imaginary frequencies are reported as negative.
    real(realwp), dimension(:,:), intent(out), optional :: L_mweigh
    !! Mass-weighted eigenvectors matrix.
    !! Typically used for conversion from Cart to mass-weighted Q.
    real(realwp), dimension(:), intent(out), optional :: red_mass
    !! Reduced mass of each vibration.

    integer :: i, ia, info, ioff, ja, lwork, n, n_at3, n_trro, n_vib1
    real(realwp) :: sqmas_i, thresh, x
    real(realwp), dimension(6) :: trro_norm
    real(realwp), dimension(3,3) :: rot_mat, tensor
    real(realwp), dimension(:), allocatable :: eval, work, work2
    real(realwp), dimension(:,:), allocatable :: crd_new, evec, evec_new, F_int
    logical :: proj_rt
    logical, dimension(:), allocatable :: mask
    character(len=80) :: msg

    if (.not.present(L_mat) .and. .not.present(L_mweigh) .and. &
        .not.present(freq)) then
        ! No quantity asked, nothing to do.
        return
    end if

    n_at3 = 3*n_at

    ! Diagonalize the matrix to get first estimate of eval/evec
    allocate(evec(n_at3,n_at3), source=F_mweigh)
    allocate(eval(n_at3))
    call xsyev('V', 'L', n_at3, evec, n_at3, eval, eval, -1, info)
    if (info /= 0) then
        write(msg, '("error code INFO=",i0)') info
        call run%error%raise_error('calc', 'gen', &
            'failed to build normal modes', &
            details='failed to get optimal work size from xsyev', &
            extra=trim(msg))
        return
    end if
    ! we want work to have at least n_at3 size for later operations.
    lwork = max(int(eval(1)), n_at3)
    allocate(work(lwork))
    call xsyev('V', 'L', n_at3, evec, n_at3, eval, work, lwork, info)
    if (info /= 0) then
        write(msg, '("error code INFO=",i0)') info
        call run%error%raise_error('calc', 'singularity', &
            'failed to build normal modes', &
            details='failed to diagonalize mass-weighted force constants &
                &matrix', &
            extra=trim(msg))
        return
    end if

    ! Check if we will project out the rotations/translations
    if (present(remove_rottrans)) then
        proj_rt = remove_rottrans
    else
        proj_rt = .true.
    end if
    if (proj_rt) then
        n = n_at3
    else
        n = 6
    end if

    ! Construct the pure rotations and translations
    ! ---------------------------------------------
    allocate(evec_new(n,n_at3), crd_new(3,n_at))
    call Eckart_orient(n_at, at_crd, at_mass, rot_mat=rot_mat, new_crd=crd_new)
    rot_mat = transpose(rot_mat)
    evec_new = f0
    do ia = 1, n_at
        sqmas_i = sqrt(at_mass(ia))
        ioff = 3*(ia-1)
        evec_new(ioff+1:ioff+3,1) = matmul([sqmas_i, f0, f0], rot_mat)
        evec_new(ioff+1:ioff+3,2) = matmul([f0, sqmas_i, f0], rot_mat)
        evec_new(ioff+1:ioff+3,3) = matmul([f0, f0, sqmas_i], rot_mat)
        tensor = sqmas_i*(crd_new(:,ia).x.rot_mat)
        evec_new(ioff+1:ioff+3,4) = tensor(1,:)
        evec_new(ioff+1:ioff+3,5) = tensor(2,:)
        evec_new(ioff+1:ioff+3,6) = tensor(3,:)
    end do
    deallocate(crd_new)

    trro_norm = sum(evec_new(:,:6)**2, dim=1)
    n_trro = count(trro_norm > 1.0e-9_realwp)
    if (n_trro < 5) then
        call run%error%raise_error('calc', 'gen', &
            'failed to build normal modes', &
            details='unable to identify rotations/translations')
        return
    end if
    n_vib1 = n_at3 - n_trro
    if (present(n_vib)) then
        if (n_vib > 0 .and. n_vib1 /= n_vib) then
            call run%error%raise_error('calc', 'gen', &
                'failed to build normal modes', &
                details='mismatch between computed and input number of modes')
            return
        end if
        if (n_vib <= 0) n_vib = n_vib1
    end if
    ! Normalize the translation/rotation vectors
    do i = 1, n_trro
        evec_new(:,i) = evec_new(:,i)/trro_norm(i)
    end do

    ! Identification of rotation/translations in original eigenvectors
    ! ----------------------------------------------------------------
    do ia = 1, n_at3
        work(ia) = f0
        do i = 1, n_trro
            x = f0
            do ja = 1, n_at3
                x = x + evec(ja,ia)*evec_new(ja,i)
            end do
            work(ia) = work(ia) + x**2
        end do
    end do
    allocate(work2(n_at3))
    work2 = work(:n_at3)
    call xlasrt('D', n_at3, work2, info)
    if (info /= 0) then
        write(msg, '("error code from sorting INFO=",i0)') info
        call run%error%raise_error('calc', 'gen', &
            'failed to build normal modes', &
            details='unable to identify the rotations/translations', &
            extra=trim(msg))
        return
    end if
    ! thresh is set slightly higher to ensure mask is properly set
    thresh = work2(n_trro) - epsilon(thresh)
    allocate(mask(n_at3))
    mask = work(:n_at3) < thresh

    ! Build corrected normal modes
    ! ----------------------------
    if (proj_rt) then
        ! We now build the full matrix, copying the missing vectors in evec_new
        i = n_trro
        do ia = 1, n_at3
            if (mask(ia)) then
                i = i + 1
                evec_new(:,i) = evec(:,ia)
            end if
        end do

        ! Gram-Schmidt orthogonalization through QR factorization
        call xgeqrf(n_at3, n_at3, evec_new, n_at3, work2, work, lwork, info)
        if (info /= 0) then
            write(msg, '("error code from mode orthogonalization INFO=",i0)') &
                info
            call run%error%raise_error('calc', 'gen', &
                'failed to build normal modes', &
                details='failed to orthogonalize the normal modes', &
                extra=trim(msg))
            return
        end if
        call xorgqr(n_at3, n_at3, n_at3, evec_new, n_at3, work2, work, lwork, &
                    info)
        if (info /= 0) then
            write(msg, &
                 '("error code when building orthgonal modes INFO=",i0)') info
            call run%error%raise_error('calc', 'gen', &
                'failed to build normal modes', &
                details='failed to orthogonalize the normal modes', &
                extra=trim(msg))
            return
        end if

        ! Project mass-weighted force const. matrix onto new basis set
        ! and recomput the new eigenvectors
        allocate(F_int(n_at3,n_at3))
        call xgemm('T', 'N', n_at3, n_at3, n_at3, f1, evec_new, n_at3, &
                   F_mweigh, n_at3, f0, evec, n_at3)
        call xgemm('N', 'N', n_at3, n_at3, n_at3, f1, evec, n_at3, &
                   evec_new, n_at3, f0, F_int, n_at3)
        call xsyev('V', 'L', n_at3, F_int, n_at3, eval, work, lwork, info)
        call xgemm('N', 'N', n_at3, n_at3, n_at3, f1, evec_new, n_at3, &
                   F_int, n_at3, f0, evec, n_at3)
        deallocate(F_int)
    end if
    deallocate(work2)

    ! Finalization: Copy necessary arrays
    ! -----------------------------------
    if (present(freq)) then
        i = 0
        do ia = 1, n_at3
            if (mask(ia)) then
                i = i + 1
                freq(i) = sign(sqrt(phys_conv%dE_au2cm(abs(eval(ia)), 2)), &
                               eval(ia))
            end if
        end do
    end if
    if (present(L_mat)) then
        i = 0
        do ia = 1, n_at3
            if (mask(ia)) then
                i = i + 1
                L_mat(:,i) = evec(:,ia)
            end if
        end do
    end if
    if (present(red_mass) .or. present(L_mweigh)) then
        do ia = 1, n_at3
            x = f0
            do ja = 1, n_at3
                evec(ja,ia) = evec(ja,ia) / sqrt(at_mass((ja+2)/3))
                x = x + evec(ja,ia)**2
            end do
            work(ia) = f1/x
        end do
        if (present(red_mass)) then
            i = 0
            do ia = 1, n_at3
                if (mask(ia)) then
                    i = i + 1
                    red_mass(i) = work(ia)
                end if
            end do
        end if
        if (present(L_mweigh)) then
            i = 0
            do ia = 1, n_at3
                if (mask(ia)) then
                    i = i + 1
                    L_mweigh(:,i) = evec(:,ia)*sqrt(work(ia))
                end if
            end do
        end if
    end if

    deallocate(eval, evec, mask, work)

end subroutine build_modes_algo

! ======================================================================

end submodule vib_build_modes