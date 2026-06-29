submodule (cubegen) cubegen_grid

    use atominfo, only: atdata
    use datatypes, only: MoleculeDB, PrimitiveFunction
    use numeric, only: realwp, f0, f2, f10, f10m4, f10m10, is_close_to
    use run_env, only: run

    integer, dimension(3), parameter :: &
        npts_cubic_grid_vmin = [21, 21, 21], &
        npts_cubic_grid_vlo = [41, 41, 41], &
        npts_cubic_grid_low = [61, 61, 61], &
        npts_cubic_grid_med = [81, 81, 81], &
        npts_cubic_grid_hig = [101, 101, 101], &
        npts_cubic_grid_vhi = [121, 121, 121], &
        npts_cubic_grid_def = npts_cubic_grid_med
    real(realwp), parameter :: grid_nucleus_tol = 1.0e-6_realwp

    interface calc_grid_bounds_dens
        !! Calculate the cube grid bounds based on electronic density.
        !!
        !! Compute the grid bounds (min and max) to enclose the electronic
        !! density.
        !!
        !! @warning
        !! The basis set is assumed to be Cartesian, no check is done here.
        !! @endwarning
        module procedure calc_grid_bounds_dens_arr, calc_grid_bounds_dens_db, &
            calc_grid_bounds_dens_ful
    end interface calc_grid_bounds_dens

    interface calc_grid_bounds_tcd
        !! Calculate the cube grid bounds compatible with TCD calculations.
        !!
        !! Compute the grid bounds (min and max) to enclose the transition
        !! current density.
        !!
        !! @warning
        !! The basis set is assumed to be Cartesian, no check is done here.
        !! @endwarning
        module procedure calc_grid_bounds_tcd_arr, calc_grid_bounds_tcd_db, &
            calc_grid_bounds_tcd_ful
    end interface calc_grid_bounds_tcd

    interface calc_grid_bounds_tdd
        !! Calculate the cube grid bounds compatible with TDD calculations.
        !!
        !! Compute the grid bounds (min and max) to enclose the transition
        !! dipole density.
        !!
        !! @warning
        !! The basis set is assumed to be Cartesian, no check is done here.
        !! @endwarning
        module procedure calc_grid_bounds_tdd_arr, calc_grid_bounds_tdd_db, &
            calc_grid_bounds_tdd_ful
    end interface calc_grid_bounds_tdd

contains

! ======================================================================

subroutine calc_grid_bounds_dens_arr(n_aos, at_num, at_crd, nprim_per_at, &
                                     bsetBF, e_dens, min_grid, max_grid)
    !! @note "version"
    !! This version expects arrays in arguments, except for bsetBF,
    !! which is a database of basis set functions.
    !! Dimensions easily recovered from arrays are recomputed.
    !! @endnote
    integer, intent(in) :: n_aos
        !! Number of atomic orbitals
    integer, dimension(:), intent(in) :: at_num
        !! Atomic numbers.
    real(realwp), dimension(:,:), intent(in) :: at_crd
        !! Atomic coordinates.
    integer, dimension(:), intent(in) :: nprim_per_at
        !! Number of primitives on each atomic center.
    class(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
        !! Basis set's basis function information.
    real(realwp), dimension(:,:), intent(in) :: e_dens
        !! Electronic density.
    real(realwp), dimension(3), intent(out) :: min_grid
        !! Grid lower bound (minimum).
    real(realwp), dimension(3), intent(out) :: max_grid
        !! Grid upper bound (maximum).

    if (size(at_crd, 2) /= size(at_num)) then
        call run%error%raise_argerror('size', &
            'size inconsistency between at_crd and at_num', &
            source='calc_grid_bounds_dens_arr')
        return
    end if

    call calc_grid_bounds_dens_ful( &
        size(at_num), n_aos, at_num, at_crd, nprim_per_at, bsetBF, e_dens, &
        min_grid, max_grid)

end subroutine calc_grid_bounds_dens_arr

! ======================================================================

subroutine calc_grid_bounds_dens_db(molDB, bsetDB, e_dens, min_grid, max_grid)
    !! @note "version"
    !! This version expects databases as arguments for core quantities.
    !! @endnote
    class(MoleculeDB), intent(in) :: molDB
        !! Molecule database.
    class(BasisSetDB), intent(in) :: bsetDB
        !! Basis set database.
    real(realwp), dimension(:,:), intent(in) :: e_dens
        !! Electronic density.
    real(realwp), dimension(3), intent(out) :: min_grid
        !! Grid lower bound (minimum).
    real(realwp), dimension(3), intent(out) :: max_grid
        !! Grid upper bound (maximum).

    call calc_grid_bounds_dens_ful( &
        molDB%n_at, bsetDB%n_basis, molDB%at_num, molDB%at_crd, &
        bsetDB%nprim_per_at, bsetDB%info, e_dens, min_grid, max_grid)
    ! real(realwp), parameter :: density_threshold = f10m10
    ! real(realwp), parameter :: buffer_factor = (f2 + f10) / f10
    ! real(realwp) :: tmp_dens
    ! real(realwp), dimension(6, 3) :: vertices
    ! real(realwp), dimension(:), allocatable :: aos_at_point
    ! real(realwp), dimension(:), allocatable :: tmp
    ! integer :: i, ia, ix
    ! logical :: density_ok

    ! min_grid = molDB%at_crd(:,1) - atdata(molDB%at_num(1))%rcov(1)
    ! max_grid = molDB%at_crd(:,1) + atdata(molDB%at_num(1))%rcov(1)

    ! do ia = 2, molDB%n_at
    !     do ix = 1, 3
    !         min_grid(ix) = min( &
    !             min_grid(ix), &
    !             molDB%at_crd(ix,ia) - atdata(molDB%at_num(ia))%rcov(1))
    !         max_grid(ix) = max( &
    !             max_grid(ix), &
    !             molDB%at_crd(ix,ia) + atdata(molDB%at_num(ia))%rcov(1))
    !     end do
    ! end do

    ! tmp_dens = f0

    ! ! Define the vertices of the box
    ! vertices(1,:) = [min_grid(1), min_grid(2), min_grid(3)]
    ! vertices(2,:) = [max_grid(1), min_grid(2), min_grid(3)]
    ! vertices(3,:) = [max_grid(1), max_grid(2), min_grid(3)]
    ! vertices(4,:) = [max_grid(1), min_grid(2), max_grid(3)]
    ! vertices(5,:) = [min_grid(1), max_grid(2), min_grid(3)]
    ! vertices(6,:) = [max_grid(1), max_grid(2), max_grid(3)]

    ! allocate(aos_at_point(bsetDB%n_basis))
    ! allocate(tmp(bsetDB%n_basis))

    ! do while (.true.)

    !     density_ok = .true.

    !     aos_at_point = f0
    !     tmp = f0
    !     tmp_dens = f0

    !     ! Evaluate density for each vertex of the box

    !     !$omp parallel do private(aos_at_point, tmp, tmp_dens) shared(density_ok)
    !     do i = 1, 6

    !         call eval_AOs_chi_at(molDB, bsetDB, vertices(i,1), vertices(i,2), &
    !                              vertices(i,3), aos_at_point)

    !         call xgemv('N', bsetDB%n_basis, bsetDB%n_basis, f1, e_dens, &
    !                    bsetDB%n_basis, aos_at_point, 1, f0, tmp, 1)

    !         ! Calculate density for this vertex
    !         tmp_dens = xdot(bsetDB%n_basis, aos_at_point, 1, tmp, 1)

    !         ! Check if the density at this vertex is below the threshold
    !         if (abs(tmp_dens) > density_threshold) then
    !             !$omp atomic write
    !             density_ok = .false.
    !             vertices(i,:) = vertices(i,:) * buffer_factor
    !         end if

    !     end do
    !     !$omp end parallel do

    !     ! If density is okay at all vertices, stop the loop
    !     if (density_ok) exit

    ! end do

    ! deallocate(tmp)
    ! deallocate(aos_at_point)

    ! min_grid(1) = minval(vertices(:, 1))
    ! min_grid(2) = minval(vertices(:, 2))
    ! min_grid(3) = minval(vertices(:, 3))

    ! max_grid(1) = maxval(vertices(:, 1))
    ! max_grid(2) = maxval(vertices(:, 2))
    ! max_grid(3) = maxval(vertices(:, 3))

end subroutine calc_grid_bounds_dens_db

! ======================================================================

subroutine calc_grid_bounds_dens_ful(n_at, n_aos, at_num, at_crd, &
                                     nprim_per_at, bsetBF, e_dens, &
                                     min_grid, max_grid)
    !! @note "version"
    !! Full version, expecting dimensions and arrays specification,
    !! except for bsetBF, which is a database of basis set functions.
    !! @endnote
    integer, intent(in) :: n_at
        !! Number of atoms
    integer, intent(in) :: n_aos
        !! Number of atomic orbitals
    integer, dimension(n_at), intent(in) :: at_num
        !! Atomic numbers.
    real(realwp), dimension(3,n_at), intent(in) :: at_crd
        !! Atomic coordinates.
    integer, dimension(:), intent(in) :: nprim_per_at
        !! Number of primitives on each atomic center.
    class(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
        !! Basis set's basis function information.
    real(realwp), dimension(:,:), intent(in) :: e_dens
        !! Electronic density.

    real(realwp), dimension(3), intent(out) :: min_grid
        !! Grid lower bound (minimum).
    real(realwp), dimension(3), intent(out) :: max_grid
        !! Grid upper bound (maximum).

    integer :: i, ia, iter, ix
    real(realwp), parameter :: density_threshold = f10m10
    real(realwp), parameter :: buffer_factor = (f2 + f10) / f10
    real(realwp) :: tmp_dens
    real(realwp), dimension(6, 3) :: vertices
    real(realwp), dimension(:), allocatable :: aos_at_point
    real(realwp), dimension(:), allocatable :: tmp
    logical, dimension(6) :: density_ok

    min_grid = at_crd(:,1) - atdata(at_num(1))%rcov(1)
    max_grid = at_crd(:,1) + atdata(at_num(1))%rcov(1)

    do ia = 2, n_at
        do ix = 1, 3
            min_grid(ix) = min( &
                min_grid(ix), &
                at_crd(ix,ia) - atdata(at_num(ia))%rcov(1))
            max_grid(ix) = max( &
                max_grid(ix), &
                at_crd(ix,ia) + atdata(at_num(ia))%rcov(1))
        end do
    end do

    tmp_dens = f0

    ! Define the vertices of the box
    vertices(1,:) = [min_grid(1), min_grid(2), min_grid(3)]
    vertices(2,:) = [max_grid(1), min_grid(2), min_grid(3)]
    vertices(3,:) = [max_grid(1), max_grid(2), min_grid(3)]
    vertices(4,:) = [max_grid(1), min_grid(2), max_grid(3)]
    vertices(5,:) = [min_grid(1), max_grid(2), min_grid(3)]
    vertices(6,:) = [max_grid(1), max_grid(2), max_grid(3)]

    allocate(aos_at_point(n_aos), tmp(n_aos))

    iter = 0
    do while (.true.)

        density_ok = .true.
        iter = iter + 1
        aos_at_point = f0
        tmp = f0

        ! Evaluate density for each vertex of the box
        !$omp parallel do private(i, aos_at_point, tmp, tmp_dens)
        do i = 1, 6
            call eval_AOs_chi_at(n_at, at_crd, nprim_per_at, bsetBF, &
                                 vertices(i,1), vertices(i,2), vertices(i,3), &
                                 aos_at_point)

            call xgemv('N', n_aos, n_aos, f1, e_dens, n_aos, aos_at_point, 1, &
                       f0, tmp, 1)

            ! Calculate density for this vertex
            tmp_dens = xdot(n_aos, aos_at_point, 1, tmp, 1)

            ! Check if the density at this vertex is below the threshold
            density_ok(i) = abs(tmp_dens) <= density_threshold
            if (.not.density_ok(i)) &
                vertices(i,:) = vertices(i,:) * buffer_factor

        end do
        !$omp end parallel do

        ! If density is okay at all vertices, stop the loop
        if (all(density_ok)) exit

    end do

    deallocate(aos_at_point, tmp)

    min_grid(1) = minval(vertices(:, 1))
    min_grid(2) = minval(vertices(:, 2))
    min_grid(3) = minval(vertices(:, 3))

    max_grid(1) = maxval(vertices(:, 1))
    max_grid(2) = maxval(vertices(:, 2))
    max_grid(3) = maxval(vertices(:, 3))

end subroutine calc_grid_bounds_dens_ful

! ======================================================================

subroutine calc_grid_bounds_tcd_arr(n_aos, at_num, at_crd, nprim_per_at, &
                                    bsetBF, e_trans_dens, min_grid, max_grid)
    !! @note "version"
    !! This version expects arrays in arguments, except for bsetBF,
    !! which is a database of basis set functions.
    !! Dimensions easily recovered from arrays are recomputed.
    !! @endnote
    integer, intent(in) :: n_aos
        !! Number of atomic orbitals
    integer, dimension(:), intent(in) :: at_num
        !! Atomic numbers.
    real(realwp), dimension(:,:), intent(in) :: at_crd
        !! Atomic coordinates.
    integer, dimension(:), intent(in) :: nprim_per_at
        !! Number of primitives on each atomic center.
    class(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
        !! Basis set's basis function information.
    real(realwp), dimension(:,:), intent(in) :: e_trans_dens
        !! Electronic transition density.
    real(realwp), dimension(3), intent(out) :: min_grid
        !! Grid lower bound (minimum).
    real(realwp), dimension(3), intent(out) :: max_grid
        !! Grid upper bound (maximum).

    if (size(at_crd, 2) /= size(at_num)) then
        call run%error%raise_argerror('size', &
            'size inconsistency between at_crd and at_num', &
            source='calc_grid_bounds_tcd_arr')
        return
    end if

    call calc_grid_bounds_tcd_ful( &
        size(at_num), n_aos, at_num, at_crd, nprim_per_at, bsetBF, &
        e_trans_dens, min_grid, max_grid)

end subroutine calc_grid_bounds_tcd_arr

! ======================================================================

subroutine calc_grid_bounds_tcd_db(molDB, bsetDB, e_trans_dens, min_grid, &
                                   max_grid)
    !! @note "version"
    !! This version expects databases as arguments for core quantities.
    !! @endnote
    class(MoleculeDB), intent(in) :: molDB
        !! Molecule Database.
    class(BasisSetDB), intent(in) :: bsetDB
        !! Basis set database.
    real(realwp), dimension(:,:), intent(in) :: e_trans_dens
        !! Electronic transition density.
    real(realwp), dimension(3), intent(out) :: min_grid
        !! Grid lower bound (minimum).
    real(realwp), dimension(3), intent(out) :: max_grid
        !! Grid upper bound (maximum).

    call calc_grid_bounds_tcd_ful( &
        molDB%n_at, bsetDB%n_basis, molDB%at_num, molDB%at_crd, &
        bsetDB%nprim_per_at, bsetDB%info, e_trans_dens, min_grid, max_grid)
    ! real(realwp), parameter :: density_threshold = f10m10
    ! real(realwp), dimension(:), allocatable :: aos_at_point
    ! real(realwp), dimension(:,:), allocatable :: aos_at_point_1der
    ! real(realwp), dimension(:), allocatable :: tmp
    ! real(realwp), dimension(3) :: tmp_dens
    ! real(realwp), parameter :: buffer_factor = (f2 + f10) / f10
    ! integer :: i, ia, ix
    ! real(realwp), dimension(6, 3) :: vertices
    ! logical :: density_ok

    ! min_grid = molDB%at_crd(:,1) - atdata(molDB%at_num(1))%rcov(1)
    ! max_grid = molDB%at_crd(:,1) + atdata(molDB%at_num(1))%rcov(1)

    ! do ia = 2, molDB%n_at
    !     do ix = 1, 3
    !         min_grid(ix) = min( &
    !             min_grid(ix), &
    !             molDB%at_crd(ix,ia) - atdata(molDB%at_num(ia))%rcov(1))
    !         max_grid(ix) = max( &
    !             max_grid(ix), &
    !             molDB%at_crd(ix,ia) + atdata(molDB%at_num(ia))%rcov(1))
    !     end do
    ! end do

    ! tmp_dens = f0

    ! ! Define the vertices of the box
    ! vertices(1, :) = [min_grid(1), min_grid(2), min_grid(3)]
    ! vertices(2, :) = [max_grid(1), min_grid(2), min_grid(3)]
    ! vertices(3, :) = [max_grid(1), max_grid(2), min_grid(3)]
    ! vertices(4, :) = [max_grid(1), min_grid(2), max_grid(3)]
    ! vertices(5, :) = [min_grid(1), max_grid(2), min_grid(3)]
    ! vertices(6, :) = [max_grid(1), max_grid(2), max_grid(3)]

    ! allocate(aos_at_point(bsetDB%n_basis))
    ! allocate(aos_at_point_1der(3, bsetDB%n_basis))
    ! allocate(tmp(bsetDB%n_basis))

    ! do while (.true.)

    !     density_ok = .true.

    !     aos_at_point = f0
    !     aos_at_point_1der = f0
    !     tmp = f0
    !     tmp_dens = f0

    !     ! Evaluate density for each vertex of the box

    !     !!$omp parallel do private(aos_at_point, aos_at_point_1der, tmp, tmp_dens) shared(density_ok)
    !     do i = 1, 6
    !         call eval_AOs_nabla_chi_at(molDB, bsetDB, vertices(i,1), &
    !                                    vertices(i,2), vertices(i,3), &
    !                                    aos_at_point, aos_at_point_1der)

    !         call xgemv('N', bsetDB%n_basis, bsetDB%n_basis, f1, e_trans_dens, &
    !                    bsetDB%n_basis, aos_at_point, 1, f0, tmp, 1)

    !         ! Calculate density for this vertex
    !         do ix = 1, 3
    !             tmp_dens(ix) = - xdot( &
    !                 bsetDB%n_basis, aos_at_point_1der(ix,:), 1, tmp, 1)
    !         enddo

    !         ! Check if the density at this vertex is below the threshold
    !         if (sum(abs(tmp_dens)) > density_threshold) then
    !             !$omp atomic write
    !             density_ok = .false.
    !             vertices(i,:) = vertices(i,:) * buffer_factor
    !         end if

    !     end do
    !     !!$omp end parallel do

    !     ! If density is okay at all vertices, stop the loop
    !     if (density_ok) exit

    ! end do

    ! deallocate(tmp)
    ! deallocate(aos_at_point_1der)
    ! deallocate(aos_at_point)

    ! min_grid(1) = minval(vertices(:, 1))
    ! min_grid(2) = minval(vertices(:, 2))
    ! min_grid(3) = minval(vertices(:, 3))

    ! max_grid(1) = maxval(vertices(:, 1))
    ! max_grid(2) = maxval(vertices(:, 2))
    ! max_grid(3) = maxval(vertices(:, 3))

end subroutine calc_grid_bounds_tcd_db

! ======================================================================

subroutine calc_grid_bounds_tcd_ful(n_at, n_aos, at_num, at_crd, &
                                    nprim_per_at, bsetBF, e_trans_dens, &
                                    min_grid, max_grid)
    !! @note "version"
    !! Full version, expecting dimensions and arrays specification,
    !! except for bsetBF, which is a database of basis set functions.
    !! @endnote
    integer, intent(in) :: n_at
        !! Number of atoms
    integer, intent(in) :: n_aos
        !! Number of atomic orbitals
    integer, dimension(n_at), intent(in) :: at_num
        !! Atomic numbers.
    real(realwp), dimension(3,n_at), intent(in) :: at_crd
        !! Atomic coordinates.
    integer, dimension(:), intent(in) :: nprim_per_at
        !! Number of primitives on each atomic center.
    class(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
        !! Basis set's basis function information.
    real(realwp), dimension(:,:), intent(in) :: e_trans_dens
        !! Electronic transition density.
    real(realwp), dimension(3), intent(out) :: min_grid
        !! Grid lower bound (minimum).
    real(realwp), dimension(3), intent(out) :: max_grid
        !! Grid upper bound (maximum).

    integer :: i, ia, ix
    real(realwp), parameter :: density_threshold = f10m10
    real(realwp), parameter :: buffer_factor = (f2 + f10) / f10
    real(realwp), dimension(3) :: tmp_dens
    real(realwp), dimension(6, 3) :: vertices
    real(realwp), dimension(:), allocatable :: tmp
    real(realwp), dimension(:), allocatable :: aos_at_point
    real(realwp), dimension(:,:), allocatable :: aos_at_point_1der
    logical, dimension(6) :: density_ok

    min_grid = at_crd(:,1) - atdata(at_num(1))%rcov(1)
    max_grid = at_crd(:,1) + atdata(at_num(1))%rcov(1)

    do ia = 2, n_at
        do ix = 1, 3
            min_grid(ix) = min( &
                min_grid(ix), &
                at_crd(ix,ia) - atdata(at_num(ia))%rcov(1))
            max_grid(ix) = max( &
                max_grid(ix), &
                at_crd(ix,ia) + atdata(at_num(ia))%rcov(1))
        end do
    end do

    ! Define the vertices of the box
    vertices(1,:) = [min_grid(1), min_grid(2), min_grid(3)]
    vertices(2,:) = [max_grid(1), min_grid(2), min_grid(3)]
    vertices(3,:) = [max_grid(1), max_grid(2), min_grid(3)]
    vertices(4,:) = [max_grid(1), min_grid(2), max_grid(3)]
    vertices(5,:) = [min_grid(1), max_grid(2), min_grid(3)]
    vertices(6,:) = [max_grid(1), max_grid(2), max_grid(3)]

    allocate(aos_at_point(n_aos), aos_at_point_1der(3,n_aos), tmp(n_aos))

    do while (.true.)

        density_ok = .true.

        aos_at_point = f0
        aos_at_point_1der = f0
        tmp = f0
        tmp_dens = f0

        ! Evaluate density for each vertex of the box

        !!$omp parallel do private(aos_at_point, aos_at_point_1der, tmp, tmp_dens) shared(density_ok)
        do i = 1, 6
            call eval_AOs_nabla_chi_at(n_at, at_crd, nprim_per_at, bsetBF, &
                                       vertices(i,1), vertices(i,2), &
                                       vertices(i,3), aos_at_point, &
                                       aos_at_point_1der)

            call xgemv('N', n_aos, n_aos, f1, e_trans_dens, n_aos, &
                       aos_at_point, 1, f0, tmp, 1)

            ! Calculate density for this vertex
            do ix = 1, 3
                tmp_dens(ix) = - xdot(n_aos, aos_at_point_1der(ix,:), 1, &
                                      tmp, 1)
            enddo

            ! Check if the density at this vertex is below the threshold
            density_ok(i) = sum(abs(tmp_dens)) <= density_threshold
            if (.not.density_ok(i)) &
                vertices(i,:) = vertices(i,:) * buffer_factor

        end do
        !!$omp end parallel do

        ! If density is okay at all vertices, stop the loop
        if (all(density_ok)) exit

    end do

    deallocate(tmp, aos_at_point_1der,aos_at_point)

    min_grid(1) = minval(vertices(:,1))
    min_grid(2) = minval(vertices(:,2))
    min_grid(3) = minval(vertices(:,3))

    max_grid(1) = maxval(vertices(:,1))
    max_grid(2) = maxval(vertices(:,2))
    max_grid(3) = maxval(vertices(:,3))

end subroutine calc_grid_bounds_tcd_ful

! ======================================================================

subroutine calc_grid_bounds_tdd_arr(n_aos, at_num, at_crd, nprim_per_at, &
                                    bsetBF, e_trans_dens, min_grid, max_grid)
    !! @note "version"
    !! This version expects arrays in arguments, except for bsetBF,
    !! which is a database of basis set functions.
    !! Dimensions easily recovered from arrays are recomputed.
    !! @endnote
    integer, intent(in) :: n_aos
        !! Number of atomic orbitals
    integer, dimension(:), intent(in) :: at_num
        !! Atomic numbers.
    real(realwp), dimension(:,:), intent(in) :: at_crd
        !! Atomic coordinates.
    integer, dimension(:), intent(in) :: nprim_per_at
        !! Number of primitives on each atomic center.
    class(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
        !! Basis set's basis function information.
    real(realwp), dimension(:,:), intent(in) :: e_trans_dens
        !! Electronic transition density.
    real(realwp), dimension(3), intent(out) :: min_grid
        !! Grid lower bound (minimum).
    real(realwp), dimension(3), intent(out) :: max_grid
        !! Grid upper bound (maximum).

    if (size(at_crd, 2) /= size(at_num)) then
        call run%error%raise_argerror('size', &
            'size inconsistency between at_crd and at_num', &
            source='calc_grid_bounds_tdd_arr')
        return
    end if

    call calc_grid_bounds_tdd_ful( &
        size(at_num), n_aos, at_num, at_crd, nprim_per_at, bsetBF, &
        e_trans_dens, min_grid, max_grid)

end subroutine calc_grid_bounds_tdd_arr

! ======================================================================

subroutine calc_grid_bounds_tdd_db(molDB, bsetDB, e_trans_dens, min_grid, &
                                   max_grid)
    !! @note "version"
    !! This version expects databases as arguments for core quantities.
    !! @endnote
    class(MoleculeDB), intent(in) :: molDB
        !! Molecule Database.
    class(BasisSetDB), intent(in) :: bsetDB
        !! Basis set database.
    real(realwp), dimension(:,:), intent(in) :: e_trans_dens
        !! Electronic transition density.
    real(realwp), dimension(3), intent(out) :: min_grid
        !! Grid lower bound (minimum).
    real(realwp), dimension(3), intent(out) :: max_grid
        !! Grid upper bound (maximum).

    call calc_grid_bounds_tdd_ful( &
        molDB%n_at, bsetDB%n_basis, molDB%at_num, molDB%at_crd, &
        bsetDB%nprim_per_at, bsetDB%info, e_trans_dens, min_grid, max_grid)

end subroutine calc_grid_bounds_tdd_db

! ======================================================================

subroutine calc_grid_bounds_tdd_ful(n_at, n_aos, at_num, at_crd, &
                                    nprim_per_at, bsetBF, e_trans_dens, &
                                    min_grid, max_grid)
    !! @note "version"
    !! Full version, expecting dimensions and arrays specification,
    !! except for bsetBF, which is a database of basis set functions.
    !! @endnote
    integer, intent(in) :: n_at
        !! Number of atoms
    integer, intent(in) :: n_aos
        !! Number of atomic orbitals
    integer, dimension(n_at), intent(in) :: at_num
        !! Atomic numbers.
    real(realwp), dimension(3,n_at), intent(in) :: at_crd
        !! Atomic coordinates.
    integer, dimension(:), intent(in) :: nprim_per_at
        !! Number of primitives on each atomic center.
    class(PrimitiveFunction), dimension(:,:), intent(in) :: bsetBF
        !! Basis set's basis function information.
    real(realwp), dimension(:,:), intent(in) :: e_trans_dens
        !! Electronic transition density.
    real(realwp), dimension(3), intent(out) :: min_grid
        !! Grid lower bound (minimum).
    real(realwp), dimension(3), intent(out) :: max_grid
        !! Grid upper bound (maximum).

    integer :: i, ia, ix
    real(realwp), parameter :: density_threshold = f10m10
    real(realwp), parameter :: buffer_factor = (f2 + f10) / f10
    real(realwp) :: rho_at
    real(realwp), dimension(3) :: tdd_at
    real(realwp), dimension(6, 3) :: vertices
    real(realwp), dimension(:), allocatable :: tmp
    real(realwp), dimension(:), allocatable :: aos_at_point
    logical, dimension(6) :: density_ok

    min_grid = at_crd(:,1) - atdata(at_num(1))%rcov(1)
    max_grid = at_crd(:,1) + atdata(at_num(1))%rcov(1)

    do ia = 2, n_at
        do ix = 1, 3
            min_grid(ix) = min( &
                min_grid(ix), &
                at_crd(ix,ia) - atdata(at_num(ia))%rcov(1))
            max_grid(ix) = max( &
                max_grid(ix), &
                at_crd(ix,ia) + atdata(at_num(ia))%rcov(1))
        end do
    end do

    ! Define the vertices of the box
    vertices(1,:) = [min_grid(1), min_grid(2), min_grid(3)]
    vertices(2,:) = [max_grid(1), min_grid(2), min_grid(3)]
    vertices(3,:) = [max_grid(1), max_grid(2), min_grid(3)]
    vertices(4,:) = [max_grid(1), min_grid(2), max_grid(3)]
    vertices(5,:) = [min_grid(1), max_grid(2), min_grid(3)]
    vertices(6,:) = [max_grid(1), max_grid(2), max_grid(3)]

    allocate(aos_at_point(n_aos), tmp(n_aos))

    do while (.true.)

        density_ok = .true.
        aos_at_point = f0
        tmp = f0
        rho_at = f0
        tdd_at = f0

        ! Evaluate transition dipole density for each vertex of the box
        do i = 1, 6
            call eval_AOs_chi_at(n_at, at_crd, nprim_per_at, bsetBF, &
                                 vertices(i,1), vertices(i,2), vertices(i,3), &
                                 aos_at_point)

            call xgemv('N', n_aos, n_aos, f1, e_trans_dens, n_aos, &
                       aos_at_point, 1, f0, tmp, 1)

            rho_at = xdot(n_aos, aos_at_point, 1, tmp, 1)
            tdd_at = -vertices(i,:)*rho_at

            density_ok(i) = sum(abs(tdd_at)) <= density_threshold
            if (.not.density_ok(i)) &
                vertices(i,:) = vertices(i,:) * buffer_factor

        end do

        ! If density is okay at all vertices, stop the loop
        if (all(density_ok)) exit

    end do

    deallocate(tmp, aos_at_point)

    min_grid(1) = minval(vertices(:,1))
    min_grid(2) = minval(vertices(:,2))
    min_grid(3) = minval(vertices(:,3))

    max_grid(1) = maxval(vertices(:,1))
    max_grid(2) = maxval(vertices(:,2))
    max_grid(3) = maxval(vertices(:,3))

end subroutine calc_grid_bounds_tdd_ful

! ======================================================================

module procedure init_gridDB_params_int

    integer :: key

    ! Set shape of grid, currently only cubic supported.
    gridDB%shape = 'cubic'

    ! Set number of points based on sparsity definition
    select case (sparsity)
    case (0)
        gridDB%n_points = npts_cubic_grid_def
    case (1)
        gridDB%n_points = npts_cubic_grid_low
    case (2)
        gridDB%n_points = npts_cubic_grid_med
    case (3)
        gridDB%n_points = npts_cubic_grid_hig
    case (4)
        gridDB%n_points = npts_cubic_grid_vhi
    case (5)
        gridDB%n_points = npts_cubic_grid_vlo
    case (10)
        gridDB%n_points = npts_cubic_grid_vmin
    case (-1)
        if (.not.present(n_points)) then
            call gridDB%error%raise_argerror('missing', &
                'Missing number of points to set grid')
            return
        end if
        gridDB%n_points = n_points
    case default
        call gridDB%error%raise_argerror('val', &
            'Unrecognized grid sparsity level')
        return
    end select

    ! Compute the grid
    key = test_cubetype_key(type)
    select case (key)
    case (1)
        if (.not.present(bsetDB)) then
            call gridDB%error%raise_argerror('missing', &
                'Missing basis set to compute grid bounds')
            return
        else if (.not.bsetDB%is_cart()) then
            call gridDB%error%raise_argerror('missing', &
                'Cartesian basis set needed to set up grid')
            return
        end if
        if (.not.present(e_dens)) then
            call gridDB%error%raise_argerror('missing', &
                'Missing electronic density to compute grid bounds')
            return
        end if
        call calc_grid_bounds_dens(molDB, bsetDB, e_dens, gridDB%min, &
                                   gridDB%max)
    case (2:4)
        if (.not.present(bsetDB)) then
            call gridDB%error%raise_argerror('missing', &
                'Missing basis set to compute grid bounds')
            return
        else if (.not.bsetDB%is_cart()) then
            call gridDB%error%raise_argerror('missing', &
                'Cartesian basis set needed to set up grid')
            return
        end if
        if (.not.present(e_trans_dens)) then
            call gridDB%error%raise_argerror('missing', &
                'Missing electronic transition density to compute grid bounds')
            return
        end if
        call calc_grid_bounds_tcd(molDB, bsetDB, e_trans_dens, gridDB%min, &
                                  gridDB%max)
    case (5:6)
        if (.not.present(bsetDB)) then
            call gridDB%error%raise_argerror('missing', &
                'Missing basis set to compute grid bounds')
            return
        else if (.not.bsetDB%is_cart()) then
            call gridDB%error%raise_argerror('missing', &
                'Cartesian basis set needed to set up grid')
            return
        end if
        if (.not.present(e_trans_dens)) then
            call gridDB%error%raise_argerror('missing', &
                'Missing electronic transition density to compute grid bounds')
            return
        end if
        call calc_grid_bounds_tdd(molDB, bsetDB, e_trans_dens, gridDB%min, &
                                  gridDB%max)
    case default
        call gridDB%error%raise_argerror('missing', &
            'Unsupported type of cube to compute grid boundaries.')
        return
    end select
    gridDB%step_size = (gridDB%max - gridDB%min) &
        / real(gridDB%n_points-1, realwp)
    if (present(molDB)) call warning_grid_points_near_nuclei(gridDB, molDB)

end procedure init_gridDB_params_int

! ======================================================================

subroutine warning_grid_points_near_nuclei(gridDB, molDB)
    class(CubegridDB), intent(in) :: gridDB
    class(MoleculeDB), intent(in) :: molDB

    integer :: iat
    integer, dimension(3) :: ipoint
    real(realwp) :: dist
    real(realwp), dimension(3) :: grid_point
    character(len=256) :: msg

    do iat = 1, molDB%n_at
        ipoint = nint((molDB%at_crd(:,iat) - gridDB%min) / gridDB%step_size) + 1
        if (any(ipoint < 1) .or. any(ipoint > gridDB%n_points)) cycle

        grid_point = gridDB%min + real(ipoint - 1, realwp) * gridDB%step_size
        dist = norm2(grid_point - molDB%at_crd(:,iat))
        if (dist <= grid_nucleus_tol) then
            write(msg, '("Grid point very close to nuclear position. Atom: ",i0,&
                &" Grid index: ",3(i0,1x),"Distance: ",es13.6)') &
                iat, ipoint, dist
            call run%error%raise_warning('calc', 'singularity', trim(msg))
        end if
    end do
end subroutine warning_grid_points_near_nuclei

! ======================================================================

module procedure init_gridDB_params_lab

    integer :: lvl_sparse

    lvl_sparse = sparsity_lab2int(sparsity)
    if (lvl_sparse == -2) then
        call gridDB%error%raise_argerror('value', &
            'Unrecognized grid sparsity level')
        return
    end if

    call init_gridDB_params_int(gridDB, type, lvl_sparse, n_points, &
                                molDB, bsetDB, e_dens, e_trans_dens)

end procedure init_gridDB_params_lab

! ======================================================================

module procedure read_gridDB_params

    integer, parameter :: MAX_KEYS = 5, MAX_ALIASES = 2
    integer :: i, ialias, iend, ikey, ios, ios_read, iu
    logical :: exists, found, unknown_ok
    character(len=:), allocatable :: trueline, vals
    character(len=256) :: line, msg
    ! Keys parameters and associated flags
    character(len=*), dimension(MAX_ALIASES, MAX_KEYS), parameter :: &
        keys = reshape([ &
            'grid type   ', '            ', &
            'grid origin ', 'grid minimum', &
            'n points    ', '            ', &
            'steps       ', 'step        ', &
            'grid maximum', '            ' ], [MAX_ALIASES, MAX_KEYS])
        !! Keys in the file to parse.
        !! Note: grid min, max, npoints and steps must be together to facilitate
        !! test on all data present.
    logical, dimension(MAX_KEYS) :: keys_found = [(.false., i = 1, MAX_KEYS)]
    integer, dimension(MAX_KEYS), parameter :: vals_num = [1, 3, -3, -3, 3]
        !! number of expected values. <0 if 1 also accepted (then duplicated)
    integer, dimension(MAX_KEYS), parameter :: vals_type = [0, 2, 1, 2, 2]
        !! expected type of the values: 0: string, 1: integer, 2: real
    integer :: nvals
    character(len=256), dimension(1) :: vals_c
    integer, dimension(3) :: vals_i
    real(realwp), dimension(3) :: vals_r

    inquire(file=file_grid, exist=exists)
    if (.not.exists) then
        write(msg, '("File ",a," does not exist")') trim(file_grid)
        call gridDB%error%raise_error('file', 'missing', trim(msg))
        return
    end if

    open(newunit=iu, file=file_grid, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        write(msg, '("Error while opening file ",a,".")') trim(file_grid)
        call gridDB%error%raise_error('file', 'open', trim(msg))
        return
    end if

    if (present(ignore_unknown)) then
        unknown_ok = ignore_unknown
    else
        unknown_ok = .false.
    end if

    read(iu, '(a)', iostat=ios) line
    do while (ios == 0)
        ! Check if comment character present
        iend = index(line, '#') - 1
        if (iend < 0) iend = len_trim(line)
        trueline = locase(trim(adjustl(line(:iend))))
        if (trueline /= ' ') then
            found = .false.
            keyloop: do ikey = 1, MAX_KEYS
                do ialias = 1, MAX_ALIASES
                    if (keys(ialias,ikey) /= ' ' &
                        .and. index(trueline, trim(keys(ialias,ikey))) == 1) then
                        if (keys_found(ikey)) then
                            write(msg, &
                                '("Duplicate/equivalent key found: ",a,".")') &
                                trim(keys(ialias,ikey))
                            call gridDB%error%raise_error('key', 'duplicate', &
                                trim(msg))
                            return
                        end if
                        keys_found(ikey) = .true.
                        found = .true.
                        ! Build string containing the values
                        vals = trim(adjustl( &
                            line(len_trim(keys(ialias,ikey))+1:)))
                        ! We accept ":", "=" as separators between keys, values
                        ! Only 1 separator accepted, we simply "cancel" the
                        ! character.
                        if (vals(1:1) == ':' .or. vals(1:1) == '=') &
                            vals(1:1) = ' '
                        exit keyloop
                    end if
                end do
            end do keyloop
            if (.not.found) then
                if (unknown_ok) then
                    read(iu, '(a)', iostat=ios) line
                    cycle
                end if
                write(msg, '("Unable to parse line """,a,""".")') trim(line)
                call gridDB%error%raise_error('data', 'struct', trim(msg))
                return
            end if
            ! key found, now parsing
            select case (vals_type(ikey))
            case (0)
                nvals = abs(vals_num(ikey))
                read(vals, *, iostat=ios_read) (vals_c(i), i=1, nvals)
                if (ios_read < 0) then
                    ! ios_read: end-of... reached, so too many args asked
                    ! check if alternative (1) number requested.
                    if (vals_num(ikey) < 0) then
                        nvals = 1
                        read(vals, *, iostat=ios_read) vals_c(1)
                        if (ios_read /= 0) then
                            write(msg, '("Failed to parse value for ",a,".")')&
                                trim(keys(ialias,ikey))
                            call gridDB%error%raise_error('data', 'struct', &
                                                          trim(msg))
                            return
                        end if
                    else
                        write(msg, &
                              '("Incorrect number of values provided for "&
                              &,a,". ",i0," expected.")') &
                                trim(keys(ialias,ikey)), nvals
                        call gridDB%error%raise_error('data', 'num', &
                                                      trim(msg))
                        return
                    end if
                end if
            case (1)
                nvals = abs(vals_num(ikey))
                read(vals, *, iostat=ios_read) (vals_i(i), i=1, nvals)
                if (ios_read < 0) then
                    ! ios_read: end-of... reached, so too many args asked
                    ! check if alternative (1) number requested.
                    if (vals_num(ikey) < 0) then
                        nvals = 1
                        read(vals, *, iostat=ios_read) vals_i(1)
                        if (ios_read /= 0) then
                            write(msg, '("Failed to parse value for ",a,".")')&
                                trim(keys(ialias,ikey))
                            call gridDB%error%raise_error('data', 'struct', &
                                                          trim(msg))
                            return
                        end if
                    else
                        write(msg, &
                                '("Incorrect number of values provided for "&
                                &,a,". ",i0," expected.")') &
                            trim(keys(ialias,ikey)), nvals
                        call gridDB%error%raise_error('data', 'num', trim(msg))
                        return
                    end if
                else if (ios_read > 0) then
                    write(msg, '("Integer values expected for ",a,".")') &
                        trim(keys(ialias,ikey))
                    call gridDB%error%raise_error('data', 'type', trim(msg))
                    return
                end if
            case (2)
                nvals = abs(vals_num(ikey))
                read(vals, *, iostat=ios_read) (vals_r(i), i=1, nvals)
                if (ios_read < 0) then
                    ! ios_read: end-of... reached, so too many args asked
                    ! check if alternative (1) number requested.
                    if (vals_num(ikey) < 0) then
                        nvals = 1
                        read(vals, *, iostat=ios_read) vals_r(1)
                        if (ios_read /= 0) then
                            write(msg, '("Failed to parse value for ",a,".")')&
                                trim(keys(ialias,ikey))
                            call gridDB%error%raise_error('data', 'struct', &
                                                          trim(msg))
                            return
                        end if
                    else
                        write(msg, &
                                '("Incorrect number of values provided for "&
                                &,a,". ",i0," expected.")') &
                            trim(keys(ialias,ikey)), nvals
                        call gridDB%error%raise_error('data', 'num', trim(msg))
                        return
                    end if
                else if (ios_read > 0) then
                    write(msg, '("Real values expected for ",a,".")') &
                        trim(keys(ialias,ikey))
                    call gridDB%error%raise_error('data', 'type', trim(msg))
                    return
                end if
            end select
            ! Parsing finished, now collect data.
            ! This part depends on the quantity to assign, so
            ! case-by-case assignment.
            select case (ikey)
            case(1)  ! Grid type
                gridDB%shape = trim(vals_c(1))
            case(2)  ! Grid minimum coordinates
                if (nvals == 3) then
                    gridDB%min = vals_r(:3)
                else
                    gridDB%min = [vals_r(1), vals_r(1), vals_r(1)]
                end if
            case(3)  ! Grid number of points
                if (nvals == 3) then
                    gridDB%n_points = vals_i(:3)
                else
                    gridDB%n_points = [vals_i(1), vals_i(1), vals_i(1)]
                end if
                if (any(gridDB%n_points <= 0)) then
                    call gridDB%error%raise_error('value', 'wrong', &
                        'Number of grid points cannot be negative.')
                    return
                end if
            case(4)  ! Grid step
                if (nvals == 3) then
                    gridDB%step_size = vals_r(:3)
                else
                    gridDB%step_size = [vals_r(1), vals_r(1), vals_r(1)]
                end if
                if (any(gridDB%step_size <= f0)) then
                    call gridDB%error%raise_error('value', 'wrong',  &
                        'Only positive steps supported.')
                    return
                end if
            case(5)  ! Grid maximum coordinates
                if (nvals == 3) then
                    gridDB%min = vals_r(:3)
                else
                    gridDB%min = [vals_r(1), vals_r(1), vals_r(1)]
                end if
            end select
        end if
        read(iu, '(a)', iostat=ios) line
    end do

    ! Now let us check if some information was missing or over parametrization
    if (all(keys_found(2:5))) then
        ! Over-parametrization, check if reasonable
        vals_r = (gridDB%max-gridDB%min)  / real(gridDB%n_points-1, realwp)
        if (.not.all(is_close_to(vals_r, gridDB%step_size, f10m4))) then
            call gridDB%error%raise_error('data', 'inconsitent', &
                "Inconsistency in grid parameters given in input.")
            return
        end if
    else if (count(.not.(keys_found(2:5))) > 1) then
        call gridDB%error%raise_error('data', 'missing', &
                "Missing parameters for the grid specifications.")
        return
    else if (.not.keys_found(2)) then
        gridDB%min = gridDB%max - gridDB%n_points * gridDB%step_size
    else if (.not.keys_found(5)) then
        gridDB%max = gridDB%min + gridDB%n_points * gridDB%step_size
    else if (.not.keys_found(3)) then
        gridDB%n_points = nint((gridDB%max-gridDB%min)/gridDB%step_size) + 1
        ! Check if step is consistent
        vals_r = (gridDB%max-gridDB%min)  / real(gridDB%n_points-1, realwp)
        if (.not.all(is_close_to(vals_r, gridDB%step_size, f10m4))) then
            call gridDB%error%raise_error('data', 'inconsistent', &
                "Step size incompatible with grid size.")
            return
        end if
    else if (.not.keys_found(4)) then
        gridDB%step_size = (gridDB%max - gridDB%min) &
            / real(gridDB%n_points-1, realwp)
    end if
    if (.not. keys_found(1)) gridDB%shape = 'cubic'

end procedure read_gridDB_params

! ======================================================================

integer function sparsity_lab2int(label) result(level)
    !! Convert the a sparsity level given as keyword to integer.
    !! The function returns -1 if the label is not supported.
    character(len=*), intent(in) :: label
        !! Sparsity label.

    select case (locase(trim(label)))
    case ('default')
        level = 0
    case ('low')
        level = 1
    case ('medium')
        level = 2
    case ('high')
        level = 3
    case ('veryhigh')
        level = 4
    case ('verylow')
        level = 5
    case ('scarce')
        level = 10
    case ('input', 'read')
        level = -1
    case default
        level = -2
    end select

end function sparsity_lab2int

! ======================================================================

end submodule cubegen_grid
