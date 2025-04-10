module geometry
    !! This module defines procedures related to structures.
    !! This includes definition of geometrical centers, rotation,
    !! superposition...
    use numeric, only: f0, f1, f2, realwp, small
    use lapack_drv, only: xsyev
    use output, only: iu_out, write_err
    use datatypes, only: MoleculeDB
    use exception, only: runstat

    implicit none

    private
    public :: center_of_mass, Eckart_orient, inertia_moments, superpose

    interface center_of_mass
        module procedure center_of_mass_arr, &
            center_of_mass_db, &
            center_of_mass_dim
    end interface center_of_mass

    interface Eckart_orient
        module procedure Eckart_orient_arr, &
            Eckart_orient_arr_upd, &
            Eckart_orient_db, &
            Eckart_orient_db_upd, &
            Eckart_orient_dim, &
            Eckart_orient_dim_upd
    end interface Eckart_orient

    interface inertia_moments
        module procedure inertia_moments_arr, &
            inertia_moments_db, &
            inertia_moments_dim
    end interface inertia_moments

    interface superpose
        module procedure superpose_arr, &
            superpose_db, &
            superpose_dim
    end interface superpose

contains

! ======================================================================

function center_of_mass_arr(at_crd, at_mass) result(com)
    !! Compute the center of mass.
    !!
    !! Computes the center of mass of a given structure.
    !!
    !! @note "Version"
    !! This version takes data arrays as elements, recovering dimensions
    !! from their shape.
    !! Information that can be easily recovered, like the number of
    !! atoms from the array sizes, is recomputed.
    !! @endnote
    real(realwp), dimension(:,:), intent(in) :: at_crd
    !! Atomic coordinates.
    real(realwp), dimension(:), intent(in) :: at_mass
    !! Atomic masses.
    real(realwp), dimension(3) :: com
    !! Center of mass.

    integer :: ia
    real(realwp) :: tot_mass
    character(len=256) :: msg

    if (size(at_crd, 2) /= size(at_mass)) then
        write(msg, '(i0," atoms from masses vs ",i0," from coordinates")') &
            size(at_mass), size(at_crd, 2)
        call runstat%raise_error('Unable to find the center of mass', &
            details='Inconsistency between atomic masses and coordinates', &
            extra=trim(msg), source='center_of_mass', cat='dev')
        return
    end if

    com = f0
    tot_mass = f0
    !$omp parallel do reduction(+:com,tot_mass)
    do ia = 1, size(at_mass)
        com = com + at_mass(ia)*at_crd(:,ia)
        tot_mass = tot_mass + at_mass(ia)
    end do
    !$omp end parallel do
    com = com / tot_mass

end function center_of_mass_arr


! ======================================================================

function center_of_mass_db(molDB) result(com)
    !! Compute the center of mass.
    !!
    !! Computes the center of mass of a given structure.
    !!
    !! @note "Version"
    !! This version takes a moleculeDB object as argument.
    !! @endnote
    class(MoleculeDB), intent(in) :: molDB
    !! Molecule database.
    real(realwp), dimension(3) :: com
    !! Center of mass.

    integer :: ia
    real(realwp) :: tot_mass

    com = f0
    tot_mass = f0
    !$omp parallel do reduction(+:com,tot_mass)
    do ia = 1, molDB%n_at
        com = com + molDB%at_mas(ia)*molDB%at_crd(:,ia)
        tot_mass = tot_mass + molDB%at_mas(ia)
    end do
    !$omp end parallel do
    com = com / tot_mass

end function center_of_mass_db

! ======================================================================

function center_of_mass_dim(n_at, at_crd, at_mass) result(com)
    !! Compute the center of mass.
    !!
    !! Computes the center of mass of a given structure.
    !!
    !! @note "Version"
    !! This version takes data arrays and explicit size specifications.
    !! No check on arrays size is done besides the explicit
    !! specification.
    !! @endnote
    integer, intent(in) :: n_at
    !! Number of atoms.
    real(realwp), dimension(3,n_at), intent(in) :: at_crd
    !! Atomic coordinates.
    real(realwp), dimension(n_at), intent(in) :: at_mass
    !! Atomic masses.
    real(realwp), dimension(3) :: com
    !! Center of mass.

    integer :: ia
    real(realwp) :: tot_mass

    com = f0
    tot_mass = f0
    !$omp parallel do reduction(+:com,tot_mass)
    do ia = 1, n_at
        com = com + at_mass(ia)*at_crd(:,ia)
        tot_mass = tot_mass + at_mass(ia)
    end do
    !$omp end parallel do
    com = com / tot_mass

end function center_of_mass_dim

! ======================================================================

subroutine Eckart_orient_arr(at_crd, at_mass, rot_mat, p_mom, new_crd, &
                             trans_vec)
    !! Get Eckart orientation parameters.
    !!
    !! Builds transformation information to Eckart orientation
    !!
    !! * rot_mat : rotation matrix from input to Eckart orientation.
    !! * p_mom : principal moments of inertia (in Eckart orientation).
    !! * new_crd : new geometry, in Eckart orientation
    !! * trans_vec : translation vector to center of mass.
    !!
    !! The output quantities are all optional, it is possible to choose
    !! which quantity to get in return.
    !!
    !! `at_crd` and `new_crd` are expected different because of the
    !! strict intent, but the code itself is safe for overlapping data.
    !!
    !! @note "Version"
    !! This version takes data arrays as elements, recovering dimensions
    !! from their shape.
    !! Information that can be easily recovered, like the number of
    !! atoms from the array sizes, is recomputed.
    !! @endnote
    real(realwp), dimension(:,:), intent(in) :: at_crd
    !! Atomic coordinates.
    real(realwp), dimension(:), intent(in) :: at_mass
    !! Atomic masses.
    real(realwp), dimension(3,3), intent(out), optional :: rot_mat
    !! Rotation matrix to Eckart orientation.
    real(realwp), dimension(3), intent(out), optional :: p_mom
    !! Principal moments of inertia.
    real(realwp), dimension(:,:), intent(inout), optional :: new_crd
    !! Coordinates after orientation.
    real(realwp), dimension(3), intent(out), optional :: trans_vec
    !! Translation vector to Eckart orientation (center of mass).

    integer :: ia, info, n_at
    real(realwp) :: x
    real(realwp), dimension(3) :: com, eval
    real(realwp), dimension(9) :: scratch
    real(realwp), dimension(3,3) :: evec, tensor
    real(realwp), dimension(3,size(at_crd, 2)) :: crd
    character(len=256) :: msg


    if (size(at_crd, 2) /= size(at_mass)) then
        write(msg, '(i0," atoms from masses vs ",i0," from coordinates")') &
            size(at_mass), size(at_crd, 2)
        call runstat%raise_error('Unable to build Eckart orientation', &
            details='Inconsistency between atomic masses and coordinates', &
            extra=trim(msg), source='Eckart_orient', cat='dev')
        return
    end if

    n_at = size(at_mass)
    com = center_of_mass_dim(n_at, at_crd, at_mass)
    do ia = 1, n_at
        crd(:,ia) = at_crd(:,ia) - com
    end do
    tensor = inertia_moments_dim(n_at, crd, at_mass)

    if (abs(tensor(1,2)) + abs(tensor(1,3)) + abs(tensor(2,3)) < small) then
        eval(1) = tensor(1,1)
        eval(2) = tensor(2,2)
        eval(3) = tensor(3,3)
        evec = f0
        evec(1,1) = f1
        evec(2,2) = f1
        evec(3,3) = f1
    else
        call xsyev('V', 'L', 3, tensor, 3, eval, scratch, 9, info)
        if (info /= 0) then
            write(msg, '("xSyEV failed in row: ",i0)') info
            call runstat%raise_error('Unable to build Eckart orientation', &
                details='Failure to diagonalize the inertia tensor matrix', &
                extra=trim(msg))
            return
        end if
        evec = tensor
    end if
    ! Check that the chirality is preserved
    x = evec(1,1)*(evec(2,2)*evec(3,3) - evec(3,2)*evec(2,3)) &
        + evec(1,2)*(evec(2,3)*evec(3,1) - evec(2,1)*evec(3,3)) &
        + evec(1,3)*(evec(2,1)*evec(3,2) - evec(2,2)*evec(3,1))
    if (x < f0) evec(:,1) = -evec(:,1)

    if (present(rot_mat)) rot_mat = evec
    if (present(p_mom)) p_mom = eval
    if (present(new_crd)) new_crd = matmul(transpose(evec), crd)
    if (present(trans_vec)) trans_vec = -com

end subroutine Eckart_orient_arr

! ======================================================================

subroutine Eckart_orient_arr_upd(at_crd, at_mass, update, rot_mat, p_mom, &
                                 trans_vec)
    !! Get Eckart orientation parameters.
    !!
    !! Builds transformation information to Eckart orientation
    !!
    !! * rot_mat : rotation matrix from input to Eckart orientation.
    !! * p_mom : principal moments of inertia (in Eckart orientation).
    !! * new_crd : new geometry, in Eckart orientation
    !! * trans_vec : translation vector to center of mass.
    !!
    !! The output quantities are all optional, it is possible to choose
    !! which quantity to get in return.
    !!
    !! `at_crd` is rewritten if `update` is true.
    !!
    !! @note "Version"
    !! This version takes data arrays as elements, recovering dimensions
    !! from their shape.
    !! Information that can be easily recovered, like the number of
    !! atoms from the array sizes, is recomputed.
    !! @endnote
    real(realwp), dimension(:,:), intent(inout) :: at_crd
    !! Atomic coordinates.
    real(realwp), dimension(:), intent(in) :: at_mass
    !! Atomic masses.
    logical, intent(in) :: update
    !! Update `at_crd` with the new coordinates.
    real(realwp), dimension(3,3), intent(out), optional :: rot_mat
    !! Rotation matrix to Eckart orientation.
    real(realwp), dimension(3), intent(out), optional :: p_mom
    !! Principal moments of inertia.
    real(realwp), dimension(3), intent(out), optional :: trans_vec
    !! Translation vector to Eckart orientation (center of mass).

    integer :: ia, info, n_at
    real(realwp) :: x
    real(realwp), dimension(3) :: com, eval
    real(realwp), dimension(9) :: scratch
    real(realwp), dimension(3,3) :: evec, tensor
    real(realwp), dimension(3,size(at_crd, 2)) :: crd
    character(len=256) :: msg


    if (size(at_crd, 2) /= size(at_mass)) then
        write(msg, '(i0," atoms from masses vs ",i0," from coordinates")') &
            size(at_mass), size(at_crd, 2)
        call runstat%raise_error('Unable to build Eckart orientation', &
            details='Inconsistency between atomic masses and coordinates', &
            extra=trim(msg), source='Eckart_orient', cat='dev')
        return
    end if

    n_at = size(at_mass)
    com = center_of_mass_dim(n_at, at_crd, at_mass)
    do ia = 1, n_at
        crd(:,ia) = at_crd(:,ia) - com
    end do
    tensor = inertia_moments_dim(n_at, crd, at_mass)

    if (abs(tensor(1,2)) + abs(tensor(1,3)) + abs(tensor(2,3)) < small) then
        eval(1) = tensor(1,1)
        eval(2) = tensor(2,2)
        eval(3) = tensor(3,3)
        evec = f0
        evec(1,1) = f1
        evec(2,2) = f1
        evec(3,3) = f1
    else
        call xsyev('V', 'L', 3, tensor, 3, eval, scratch, 9, info)
        if (info /= 0) then
            write(msg, '("xSyEV failed in row: ",i0)') info
            call runstat%raise_error('Unable to build Eckart orientation', &
                details='Failure to diagonalize the inertia tensor matrix', &
                extra=trim(msg))
            return
        end if
        evec = tensor
    end if
    ! Check that the chirality is preserved
    x = evec(1,1)*(evec(2,2)*evec(3,3) - evec(3,2)*evec(2,3)) &
        + evec(1,2)*(evec(2,3)*evec(3,1) - evec(2,1)*evec(3,3)) &
        + evec(1,3)*(evec(2,1)*evec(3,2) - evec(2,2)*evec(3,1))
    if (x < f0) evec(:,1) = -evec(:,1)

    if (present(rot_mat)) rot_mat = evec
    if (present(p_mom)) p_mom = eval
    if (update) at_crd = matmul(transpose(evec), crd)
    if (present(trans_vec)) trans_vec = -com

end subroutine Eckart_orient_arr_upd

! ======================================================================

subroutine Eckart_orient_db(molDB, rot_mat, p_mom, new_crd, new_molDB, &
                            trans_vec)
    !! Get Eckart orientation parameters.
    !!
    !! Builds transformation information to Eckart orientation
    !!
    !! * rot_mat : rotation matrix from input to Eckart orientation.
    !! * p_mom : principal moments of inertia (in Eckart orientation).
    !! * new_crd : new geometry, in Eckart orientation, as an array
    !! * new_molDB : new molecule database, with geom. in Eckart orient.
    !! * trans_vec : translation vector to center of mass.
    !!
    !! The output quantities are all optional, it is possible to choose
    !! which quantity to get in return.
    !!
    !! `molDB` and `new_molDB` are expected different because of the
    !! strict intent, but the code itself is safe for overlapping data.
    !!
    !! @note "Version"
    !! This version takes a moleculeDB object as argument.
    !! @endnote
    class(MoleculeDB), intent(in) :: molDB
    !! Molecule database.
    real(realwp), dimension(3,3), intent(out), optional :: rot_mat
    !! Rotation matrix to Eckart orientation.
    real(realwp), dimension(3), intent(out), optional :: p_mom
    !! Principal moments of inertia.
    real(realwp), dimension(:,:), intent(out), optional :: new_crd
    !! Coordinates after orientation.
    class(MoleculeDB), intent(out), optional :: new_molDB
    !! Molecule database.
    real(realwp), dimension(3), intent(out), optional :: trans_vec
    !! Translation vector to Eckart orientation (center of mass).

    integer :: ia, info
    real(realwp) :: x
    real(realwp), dimension(3) :: com, eval
    real(realwp), dimension(9) :: scratch
    real(realwp), dimension(3,3) :: evec, tensor
    real(realwp), dimension(3,molDB%n_at) :: crd
    character(len=256) :: msg

    com = center_of_mass_db(molDB)
    do ia = 1, molDB%n_at
        crd(:,ia) = molDB%at_crd(:,ia) - com
    end do
    tensor = inertia_moments_dim(molDB%n_at, crd, molDB%at_mas)

    if (abs(tensor(1,2)) + abs(tensor(1,3)) + abs(tensor(2,3)) < small) then
        eval(1) = tensor(1,1)
        eval(2) = tensor(2,2)
        eval(3) = tensor(3,3)
        evec = f0
        evec(1,1) = f1
        evec(2,2) = f1
        evec(3,3) = f1
    else
        call xsyev('V', 'L', 3, tensor, 3, eval, scratch, 9, info)
        if (info /= 0) then
            write(msg, '("xSyEV failed in row: ",i0)') info
            call runstat%raise_error('Unable to build Eckart orientation', &
                details='Failure to diagonalize the inertia tensor matrix', &
                extra=trim(msg))
            return
        end if
        evec = tensor
    end if

    ! Check that the chirality is preserved
    x = evec(1,1)*(evec(2,2)*evec(3,3) - evec(3,2)*evec(2,3)) &
        + evec(1,2)*(evec(2,3)*evec(3,1) - evec(2,1)*evec(3,3)) &
        + evec(1,3)*(evec(2,1)*evec(3,2) - evec(2,2)*evec(3,1))
    if (x < f0) evec(:,1) = -evec(:,1)

    if (present(rot_mat)) rot_mat = evec
    if (present(p_mom)) p_mom = eval
    if (present(new_crd)) new_crd = matmul(transpose(evec), crd)
    if (present(new_molDB)) then
        call molDB%copy_to(new_molDB)
        new_molDB%at_crd = matmul(transpose(evec), crd)
    end if
    if (present(trans_vec)) trans_vec = -com

end subroutine Eckart_orient_db

! ======================================================================

subroutine Eckart_orient_db_upd(molDB, update, rot_mat, p_mom, new_crd, &
                                trans_vec)
    !! Get Eckart orientation parameters.
    !!
    !! Builds transformation information to Eckart orientation.
    !!
    !! * rot_mat : rotation matrix from input to Eckart orientation.
    !! * p_mom : principal moments of inertia (in Eckart orientation).
    !! * trans_vec : translation vector to center of mass.
    !!
    !! The output quantities are all optional, it is possible to choose
    !! which quantity to get in return.
    !!
    !! `molDB` is rewritten if `update` is true.
    !!
    !! @note "Version"
    !! This version takes a moleculeDB object as argument.
    !! @endnote
    class(MoleculeDB), intent(inout) :: molDB
    !! Molecule database.
    logical, intent(in) :: update
    !! Update `at_crd` with the new coordinates.
    real(realwp), dimension(3,3), intent(out), optional :: rot_mat
    !! Rotation matrix to Eckart orientation.
    real(realwp), dimension(3), intent(out), optional :: p_mom
    !! Principal moments of inertia.
    real(realwp), dimension(:,:), intent(out), optional :: new_crd
    !! Coordinates after orientation.
    real(realwp), dimension(3), intent(out), optional :: trans_vec
    !! Translation vector to Eckart orientation (center of mass).

    integer :: ia, info
    real(realwp) :: x
    real(realwp), dimension(3) :: com, eval
    real(realwp), dimension(9) :: scratch
    real(realwp), dimension(3,3) :: evec, tensor
    real(realwp), dimension(3,molDB%n_at) :: crd
    character(len=256) :: msg

    com = center_of_mass_db(molDB)
    do ia = 1, molDB%n_at
        crd(:,ia) = molDB%at_crd(:,ia) - com
    end do
    tensor = inertia_moments_dim(molDB%n_at, crd, molDB%at_mas)

    if (abs(tensor(1,2)) + abs(tensor(1,3)) + abs(tensor(2,3)) < small) then
        eval(1) = tensor(1,1)
        eval(2) = tensor(2,2)
        eval(3) = tensor(3,3)
        evec = f0
        evec(1,1) = f1
        evec(2,2) = f1
        evec(3,3) = f1
    else
        call xsyev('V', 'L', 3, tensor, 3, eval, scratch, 9, info)
        if (info /= 0) then
            write(msg, '("xSyEV failed in row: ",i0)') info
            call runstat%raise_error('Unable to build Eckart orientation', &
                details='Failure to diagonalize the inertia tensor matrix', &
                extra=trim(msg))
            return
        end if
        evec = tensor
    end if

    ! Check that the chirality is preserved
    x = evec(1,1)*(evec(2,2)*evec(3,3) - evec(3,2)*evec(2,3)) &
        + evec(1,2)*(evec(2,3)*evec(3,1) - evec(2,1)*evec(3,3)) &
        + evec(1,3)*(evec(2,1)*evec(3,2) - evec(2,2)*evec(3,1))
    if (x < f0) evec(:,1) = -evec(:,1)

    if (present(rot_mat)) rot_mat = evec
    if (present(p_mom)) p_mom = eval
    if (present(new_crd)) new_crd = matmul(transpose(evec), crd)
    if (update) molDB%at_crd = matmul(transpose(evec), crd)
    if (present(trans_vec)) trans_vec = -com

end subroutine Eckart_orient_db_upd

! ======================================================================

subroutine Eckart_orient_dim(n_at, at_crd, at_mass, rot_mat, p_mom, new_crd, &
                             trans_vec)
    !! Get Eckart orientation parameters.
    !!
    !! Builds transformation information to Eckart orientation
    !!
    !! * rot_mat : rotation matrix from input to Eckart orientation.
    !! * p_mom : principal moments of inertia (in Eckart orientation).
    !! * new_crd : new geometry, in Eckart orientation
    !! * trans_vec : translation vector to center of mass.
    !!
    !! The output quantities are all optional, it is possible to choose
    !! which quantity to get in return.
    !!
    !!
    !! `at_crd` and `new_crd` are expected different because of the
    !! strict intent, but the code itself is safe for overlapping data.
    !!
    !! @note "Version"
    !! This version takes data arrays and explicit size specifications.
    !! No check on arrays size is done besides the explicit
    !! specification.
    !! @endnote
    integer, intent(in) :: n_at
    !! Number of atoms.
    real(realwp), dimension(3,n_at), intent(in) :: at_crd
    !! Atomic coordinates.
    real(realwp), dimension(n_at), intent(in) :: at_mass
    !! Atomic masses.
    real(realwp), dimension(3,3), intent(out), optional :: rot_mat
    !! Rotation matrix to Eckart orientation.
    real(realwp), dimension(3), intent(out), optional :: p_mom
    !! Principal moments of inertia.
    real(realwp), dimension(:,:), intent(out), optional :: new_crd
    !! Coordinates after orientation.
    real(realwp), dimension(3), intent(out), optional :: trans_vec
    !! Translation vector to Eckart orientation (center of mass).

    integer :: ia, info
    real(realwp) :: x
    real(realwp), dimension(3) :: com, eval
    real(realwp), dimension(9) :: scratch
    real(realwp), dimension(3,3) :: evec, tensor
    real(realwp), dimension(3,n_at) :: crd
    character(len=256) :: msg

    com = center_of_mass_dim(n_at, at_crd, at_mass)
    do ia = 1, n_at
        crd(:,ia) = at_crd(:,ia) - com
    end do
    tensor = inertia_moments_dim(n_at, crd, at_mass)

    if (abs(tensor(1,2)) + abs(tensor(1,3)) + abs(tensor(2,3)) < small) then
        eval(1) = tensor(1,1)
        eval(2) = tensor(2,2)
        eval(3) = tensor(3,3)
        evec = f0
        evec(1,1) = f1
        evec(2,2) = f1
        evec(3,3) = f1
    else
        call xsyev('V', 'L', 3, tensor, 3, eval, scratch, 9, info)
        if (info /= 0) then
            write(msg, '("xSyEV failed in row: ",i0)') info
            call runstat%raise_error('Unable to build Eckart orientation', &
                details='Failure to diagonalize the inertia tensor matrix', &
                extra=trim(msg))
            return
        end if
        evec = tensor
    end if

    ! Check that the chirality is preserved
    x = evec(1,1)*(evec(2,2)*evec(3,3) - evec(3,2)*evec(2,3)) &
        + evec(1,2)*(evec(2,3)*evec(3,1) - evec(2,1)*evec(3,3)) &
        + evec(1,3)*(evec(2,1)*evec(3,2) - evec(2,2)*evec(3,1))
    if (x < f0) evec(:,1) = -evec(:,1)

    if (present(rot_mat)) rot_mat = evec
    if (present(p_mom)) p_mom = eval
    if (present(new_crd)) new_crd = matmul(transpose(evec), crd)
    if (present(trans_vec)) trans_vec = -com

end subroutine Eckart_orient_dim

! ======================================================================

subroutine Eckart_orient_dim_upd(n_at, at_crd, at_mass, update, rot_mat, &
                                 p_mom, trans_vec)
    !! Get Eckart orientation parameters.
    !!
    !! Builds transformation information to Eckart orientation
    !!
    !! * rot_mat : rotation matrix from input to Eckart orientation.
    !! * p_mom : principal moments of inertia (in Eckart orientation).
    !! * new_crd : new geometry, in Eckart orientation
    !! * trans_vec : translation vector to center of mass.
    !!
    !! The output quantities are all optional, it is possible to choose
    !! which quantity to get in return.
    !!
    !! `at_crd` is rewritten if `update` is true.
    !!
    !! @note "Version"
    !! This version takes data arrays and explicit size specifications.
    !! No check on arrays size is done besides the explicit
    !! specification.
    !! @endnote
    integer, intent(in) :: n_at
    !! Number of atoms.
    real(realwp), dimension(3,n_at), intent(inout) :: at_crd
    !! Atomic coordinates.
    real(realwp), dimension(n_at), intent(in) :: at_mass
    !! Atomic masses.
    logical, intent(in) :: update
    !! Update `at_crd` with the new coordinates.
    real(realwp), dimension(3,3), intent(out), optional :: rot_mat
    !! Rotation matrix to Eckart orientation.
    real(realwp), dimension(3), intent(out), optional :: p_mom
    !! Principal moments of inertia.
    real(realwp), dimension(3), intent(out), optional :: trans_vec
    !! Translation vector to Eckart orientation (center of mass).

    integer :: ia, info
    real(realwp) :: x
    real(realwp), dimension(3) :: com, eval
    real(realwp), dimension(9) :: scratch
    real(realwp), dimension(3,3) :: evec, tensor
    real(realwp), dimension(3,n_at) :: crd
    character(len=256) :: msg

    com = center_of_mass_dim(n_at, at_crd, at_mass)
    do ia = 1, n_at
        crd(:,ia) = at_crd(:,ia) - com
    end do
    tensor = inertia_moments_dim(n_at, crd, at_mass)

    if (abs(tensor(1,2)) + abs(tensor(1,3)) + abs(tensor(2,3)) < small) then
        eval(1) = tensor(1,1)
        eval(2) = tensor(2,2)
        eval(3) = tensor(3,3)
        evec = f0
        evec(1,1) = f1
        evec(2,2) = f1
        evec(3,3) = f1
    else
        call xsyev('V', 'L', 3, tensor, 3, eval, scratch, 9, info)
        if (info /= 0) then
            write(msg, '("xSyEV failed in row: ",i0)') info
            call runstat%raise_error('Unable to build Eckart orientation', &
                details='Failure to diagonalize the inertia tensor matrix', &
                extra=trim(msg))
            return
        end if
        evec = tensor
    end if

    ! Check that the chirality is preserved
    x = evec(1,1)*(evec(2,2)*evec(3,3) - evec(3,2)*evec(2,3)) &
        + evec(1,2)*(evec(2,3)*evec(3,1) - evec(2,1)*evec(3,3)) &
        + evec(1,3)*(evec(2,1)*evec(3,2) - evec(2,2)*evec(3,1))
    if (x < f0) evec(:,1) = -evec(:,1)

    if (present(rot_mat)) rot_mat = evec
    if (present(p_mom)) p_mom = eval
    if (update) at_crd = matmul(transpose(evec), crd)
    if (present(trans_vec)) trans_vec = -com

end subroutine Eckart_orient_dim_upd

! ======================================================================

function inertia_moments_arr(at_crd, at_mass) result(tensor)
    !! Compute the inertia moments tensor.
    !!
    !! Computes the inertia moments tensor.
    !!
    !! @note "Version"
    !! This version takes data arrays as elements, recovering dimensions
    !! from their shape.
    !! Information that can be easily recovered, like the number of
    !! atoms from the array sizes, is recomputed.
    !! @endnote
    real(realwp), dimension(:,:), intent(in) :: at_crd
    !! Atomic coordinates.
    real(realwp), dimension(:), intent(in) :: at_mass
    !! Atomic masses.
    real(realwp), dimension(3,3) :: tensor
    !! Center of mass.

    integer :: ia
    real(realwp) :: atm
    character(len=256) :: msg

    if (size(at_crd, 2) /= size(at_mass)) then
        write(msg, '(i0," atoms from masses vs ",i0," from coordinates")') &
            size(at_mass), size(at_crd, 2)
        call runstat%raise_error('Unable to build inertia tensor', &
            details='Inconsistency between atomic masses and coordinates', &
            extra=trim(msg), source='inertia_moments', cat='dev')
        return
    end if

    ! Definition of the tensor of moments of inertia
    tensor = f0
    !$omp parallel do default(shared) private(atm) reduction(+:tensor)
    do ia = 1, size(at_mass)
        atm = at_mass(ia)
        tensor(1,1) = tensor(1,1) + atm*(at_crd(2,ia)**2 + at_crd(3,ia)**2)
        tensor(2,2) = tensor(2,2) + atm*(at_crd(1,ia)**2 + at_crd(3,ia)**2)
        tensor(3,3) = tensor(3,3) + atm*(at_crd(1,ia)**2 + at_crd(2,ia)**2)
        tensor(1,2) = tensor(1,2) - atm*(at_crd(1,ia) * at_crd(2,ia))
        tensor(1,3) = tensor(1,3) - atm*(at_crd(1,ia) * at_crd(3,ia))
        tensor(2,3) = tensor(2,3) - atm*(at_crd(2,ia) * at_crd(3,ia))
    end do
    !$omp end parallel do
    tensor(2,1) = tensor(1,2)
    tensor(3,1) = tensor(1,3)
    tensor(3,2) = tensor(2,3)

end function inertia_moments_arr

! ======================================================================

function inertia_moments_db(molDB) result(tensor)
    !! Compute the inertia moments tensor.
    !!
    !! Computes the inertia moments tensor.
    !!
    !! @note "Version"
    !! This version takes a moleculeDB object as argument.
    !! @endnote
    class(MoleculeDB), intent(in) :: molDB
    !! Molecule database.
    real(realwp), dimension(3,3) :: tensor
    !! Center of mass.

    integer :: ia
    real(realwp) :: atm

    ! Definition of the tensor of moments of inertia
    tensor = f0
    !$omp parallel do default(shared) private(atm) reduction(+:tensor)
    do ia = 1, molDB%n_at
        associate (xyz => molDB%at_crd)
        atm = molDB%at_mas(ia)
        tensor(1,1) = tensor(1,1) + atm*(xyz(2,ia)**2 + xyz(3,ia)**2)
        tensor(2,2) = tensor(2,2) + atm*(xyz(1,ia)**2 + xyz(3,ia)**2)
        tensor(3,3) = tensor(3,3) + atm*(xyz(1,ia)**2 + xyz(2,ia)**2)
        tensor(1,2) = tensor(1,2) - atm*(xyz(1,ia) * xyz(2,ia))
        tensor(1,3) = tensor(1,3) - atm*(xyz(1,ia) * xyz(3,ia))
        tensor(2,3) = tensor(2,3) - atm*(xyz(2,ia) * xyz(3,ia))
        end associate
    end do
    !$omp end parallel do
    tensor(2,1) = tensor(1,2)
    tensor(3,1) = tensor(1,3)
    tensor(3,2) = tensor(2,3)

end function inertia_moments_db

! ======================================================================

function inertia_moments_dim(n_at, at_crd, at_mass) result(tensor)
    !! Compute the inertia moments tensor.
    !!
    !! Computes the inertia moments tensor.
    !!
    !! @note "Version"
    !! This version takes data arrays and explicit size specifications.
    !! No check on arrays size is done besides the explicit
    !! specification.
    !! @endnote
    integer, intent(in) :: n_at
    !! Number of atoms.
    real(realwp), dimension(3,n_at), intent(in) :: at_crd
    !! Atomic coordinates.
    real(realwp), dimension(n_at), intent(in) :: at_mass
    !! Atomic masses.
    real(realwp), dimension(3,3) :: tensor
    !! Center of mass.

    integer :: ia
    real(realwp) :: atm

    ! Definition of the tensor of moments of inertia
    tensor = f0
    !$omp parallel do default(shared) private(atm) reduction(+:tensor)
    do ia = 1, n_at
        atm = at_mass(ia)
        tensor(1,1) = tensor(1,1) + atm*(at_crd(2,ia)**2 + at_crd(3,ia)**2)
        tensor(2,2) = tensor(2,2) + atm*(at_crd(1,ia)**2 + at_crd(3,ia)**2)
        tensor(3,3) = tensor(3,3) + atm*(at_crd(1,ia)**2 + at_crd(2,ia)**2)
        tensor(1,2) = tensor(1,2) - atm*(at_crd(1,ia) * at_crd(2,ia))
        tensor(1,3) = tensor(1,3) - atm*(at_crd(1,ia) * at_crd(3,ia))
        tensor(2,3) = tensor(2,3) - atm*(at_crd(2,ia) * at_crd(3,ia))
    end do
    !$omp end parallel do
    tensor(2,1) = tensor(1,2)
    tensor(3,1) = tensor(1,3)
    tensor(3,2) = tensor(2,3)

end function inertia_moments_dim

! ======================================================================

subroutine superpose_arr(at_crd_new, at_crd_ref, weights, mask, rot_mat, &
                         trans_vec)
    !! Superpose a geometrical structure onto another.
    !!
    !! The routine tries to find the best overlap between a reference
    !! structure and a structure to superpose.  The coordinates can be
    !! weighted, each atom having its own weight specification.
    !! A mask can be provided so the overlap is only done on the
    !! selected atoms.
    !!
    !! The routine use the quaternion method, as documented in:
    !!
    !! 1. G.R. Kneller, Mol. Sim. 7, 113-119 (1991)
    !! 2. G.R. Kneller, J. Chim. Phys. 88, 2709-2715 (1991)
    !!
    !! @warning "Limitations"
    !! `mask`, `weights` and the atomic coordinates must have the same
    !! "leading" dimensions.
    !! @endwarning
    !!
    !! @note "Version"
    !! This version takes data arrays as elements, recovering dimensions
    !! from their shape.
    !! Information that can be easily recovered, like the number of
    !! atoms from the array sizes, is recomputed.
    !! @endnote
    real(realwp), dimension(:,:), intent(inout) :: at_crd_new
    !! Atomic coordinates to transform, rotated on output.
    real(realwp), dimension(:,:), intent(in) :: at_crd_ref
    !! Atomic coordinates (reference).
    real(realwp), dimension(:), intent(in), optional :: weights
    !! Atomic weights (e.g., masses).
    logical, dimension(:), intent(in), optional :: mask
    !! Mask, specifying which atom to include in superposition
    real(realwp), dimension(3,3), intent(out), optional :: rot_mat
    !! Rotation matrix from original to final orientation.
    real(realwp), dimension(3), intent(out), optional :: trans_vec
    !! Translation vector from original to final position.

    integer :: ia, info, n_at
    real(realwp) :: fac, q0, q0q0, q0q1, q0q2, q0q3, q1, q1q1, q1q2, q1q3, &
        q2, q2q2, q2q3, q3, q3q3, xnew, xref, xy, xz, ynew, yref, yx, yz, &
        znew, zref, zx, zy
    real(realwp), dimension(3) :: com_new, com_ref
    real(realwp), dimension(4) :: qeval
    real(realwp), dimension(16) :: scratch
    real(realwp), dimension(3,3) :: rmat
    real(realwp), dimension(4,4) :: qmat, qevec
    character(len=256) :: msg

    real(realwp), dimension(:), allocatable :: at_mass
    logical, dimension(:), allocatable :: at_mask

    n_at = size(at_crd_ref, 2)
    if (n_at /= size(at_crd_new, 2)) then
        write(msg, '(i0," atoms from reference vs ",i0,&
                   &" from new coordinates")') n_at, size(at_crd_new, 2)
        call runstat%raise_error('Unable to superpose structures', &
            details='Inconsistency between atomic coordinates', &
                   extra=trim(msg), source='superpose', cat='dev')
        return
    end if

    ! Check optional input arrays and set alternatives if needed
    if (present(weights)) then
        if (size(weights) /= n_at) then
            write(msg, '("weights provided for ",i0," atoms while ",i0,&
                       &" expected")')  n_at, size(weights)
            call runstat%raise_error('Unable to superpose structures', &
                details='Inconsistency between weights and coordinates', &
                extra=trim(msg), source='superpose', cat='dev')
            return
        end if
        at_mass = weights
    else
        allocate(at_mass(n_at))
        at_mass = f1
    end if

    if (present(mask)) then
        if (size(mask) /= n_at) then
            write(msg, '("Mask provided for ",i0," atoms while ",i0,&
                       &" expected")')  n_at, size(mask)
            call runstat%raise_error('Unable to superpose structures', &
                details='Inconsistency between mask and coordinates', &
                extra=trim(msg), source='superpose', cat='dev')
            return
        end if
        at_mask = mask
    else
        allocate(at_mask(n_at))
        at_mask = .true.
    end if

    ! Main algorithm
    ! -- 1. get center of mass
    com_ref = center_of_mass_dim(n_at, at_crd_ref, at_mass)
    com_new = center_of_mass_dim(n_at, at_crd_new, at_mass)
    ! -- 2. build the quaternion matrix
    qmat = f0
    do ia = 1, n_at
        if (.not.at_mask(ia)) cycle
        fac = at_mass(ia)

        xref = at_crd_ref(1,ia) - com_ref(1)
        yref = at_crd_ref(2,ia) - com_ref(2)
        zref = at_crd_ref(3,ia) - com_ref(3)
        xnew = at_crd_new(1,ia) - com_new(1)
        ynew = at_crd_new(2,ia) - com_new(2)
        znew = at_crd_new(3,ia) - com_new(3)

        xy = xnew * yref
        xz = xnew * zref
        yx = ynew * xref
        yz = ynew * zref
        zx = znew * xref
        zy = znew * yref
        q0 = xnew**2 + ynew**2 + znew**2 + xref**2 + yref**2 + zref**2
        q1 = f2 * xnew * xref
        q2 = f2 * ynew * yref
        q3 = f2 * znew * zref
        qmat(1,1) = qmat(1,1) + fac * (q0 - q1 - q2 - q3)
        qmat(2,2) = qmat(2,2) + fac * (q0 - q1 + q2 + q3)
        qmat(3,3) = qmat(3,3) + fac * (q0 + q1 - q2 + q3)
        qmat(4,4) = qmat(4,4) + fac * (q0 + q1 + q2 - q3)
        qmat(2,1) = qmat(2,1) + f2 * fac * (yz - zy)
        qmat(3,1) = qmat(3,1) + f2 * fac * (zx - xz)
        qmat(4,1) = qmat(4,1) + f2 * fac * (xy - yx)
        qmat(3,2) = qmat(3,2) - f2 * fac * (xy + yx)
        qmat(4,2) = qmat(4,2) - f2 * fac * (xz + zx)
        qmat(4,3) = qmat(4,3) - f2 * fac * (yz + zy)
    end do

    ! -- 3. Diagonalize
    call xsyev('V', 'L', 4, qmat, 4, qeval, scratch, 16, info)
    if (info /= 0) then
        write(msg, '("xSyEV failed in row: ",i0)') info
        call runstat%raise_error('Unable to superpose structures', &
            details='Failure to diagonalize the quaternion matrix', &
            extra=trim(msg))
        return
    end if
    qevec = qmat
    ! We take the eigenvector corresponding to smalles eigenvalue
    ! By construction, it is the first one
    q0 = qevec(1,1)
    q1 = qevec(2,1)
    q2 = qevec(3,1)
    q3 = qevec(4,1)
    ! -- 4. Build rotation matrix
    q0q0 = q0**2
    q1q1 = q1**2
    q2q2 = q2**2
    q3q3 = q3**2
    rmat(1,1) = q0q0 + q1q1 - q2q2 - q3q3
    rmat(2,2) = q0q0 - q1q1 + q2q2 - q3q3
    rmat(3,3) = q0q0 - q1q1 - q2q2 + q3q3
    q0q1 = q0 * q1
    q0q2 = q0 * q2
    q0q3 = q0 * q3
    q1q2 = q1 * q2
    q1q3 = q1 * q3
    q2q3 = q2 * q3
    rmat(1,2) = f2 * (q1q2 - q0q3)
    rmat(2,1) = f2 * (q1q2 + q0q3)
    rmat(1,3) = f2 * (q1q3 + q0q2)
    rmat(3,1) = f2 * (q1q3 - q0q2)
    rmat(2,3) = f2 * (q2q3 - q0q1)
    rmat(3,2) = f2 * (q2q3 + q0q1)

    do ia = 1, n_at
        xnew = at_crd_new(1,ia) - com_new(1)
        ynew = at_crd_new(2,ia) - com_new(2)
        znew = at_crd_new(3,ia) - com_new(3)
        at_crd_new(1,ia) = rmat(1,1)*xnew + rmat(2,1)*ynew + rmat(3,1)*znew
        at_crd_new(2,ia) = rmat(1,2)*xnew + rmat(2,2)*ynew + rmat(3,2)*znew
        at_crd_new(3,ia) = rmat(1,3)*xnew + rmat(2,3)*ynew + rmat(3,3)*znew
        at_crd_new(:,ia) = at_crd_new(:,ia) + com_ref
    end do

    if (present(rot_mat)) rot_mat = rmat
    if (present(trans_vec)) trans_vec = com_ref - com_new

end subroutine superpose_arr

! ======================================================================

subroutine superpose_db(molDB, at_crd_ref, use_mass, weights, mask, &
                        rot_mat, trans_vec)
    !! Superpose a geometrical structure onto another.
    !!
    !! The routine tries to find the best overlap between a reference
    !! structure and a structure to superpose.  The coordinates can be
    !! weighted, each atom having its own weight specification.
    !! A mask can be provided so the overlap is only done on the
    !! selected atoms.
    !!
    !! The routine use the quaternion method, as documented in:
    !!
    !! 1. G.R. Kneller, Mol. Sim. 7, 113-119 (1991)
    !! 2. G.R. Kneller, J. Chim. Phys. 88, 2709-2715 (1991)
    !!
    !! @warning "Limitations"
    !! `mask`, `weights` and the atomic coordinates must have the same
    !! "leading" dimensions.
    !! @endwarning
    !!
    !! @note "Version"
    !! This version takes a moleculeDB object as argument.
    !! @endnote
    class(MoleculeDB), intent(inout) :: molDB
    !! Molecule database.
    real(realwp), dimension(3,molDB%n_at), intent(in) :: at_crd_ref
    !! Atomic coordinates (reference).
    logical, intent(in), optional :: use_mass
    !! Use atomic masses stored in `molDB` as weights.
    real(realwp), dimension(molDB%n_at), intent(in), optional :: weights
    !! Atomic weights (e.g., masses).
    logical, dimension(molDB%n_at), intent(in), optional :: mask
    !! Mask, specifying which atom to include in superposition
    real(realwp), dimension(3,3), intent(out), optional :: rot_mat
    !! Rotation matrix from original to final orientation.
    real(realwp), dimension(3), intent(out), optional :: trans_vec
    !! Translation vector from original to final position.

    integer :: ia, info
    real(realwp) :: fac, q0, q0q0, q0q1, q0q2, q0q3, q1, q1q1, q1q2, q1q3, &
        q2, q2q2, q2q3, q3, q3q3, xnew, xref, xy, xz, ynew, yref, yx, yz, &
        znew, zref, zx, zy
    real(realwp), dimension(3) :: com_new, com_ref
    real(realwp), dimension(4) :: qeval
    real(realwp), dimension(16) :: scratch
    real(realwp), dimension(3,3) :: rmat
    real(realwp), dimension(4,4) :: qmat, qevec
    character(len=256) :: msg

    real(realwp), dimension(molDB%n_at) :: at_mass
    logical, dimension(molDB%n_at) :: at_mask

    ! Check optional input arrays and set alternatives if needed
    if (present(weights)) then
        if (present(use_mass)) then
            write(msg, '("Weight specifications given while requesting use &
                       &of internal atomic masses")')
            call runstat%raise_error('Unable to superpose structures', &
                details='Confusing request on weight specifications', &
                extra=trim(msg), cat='dev', source='superpose')
            return
        end if
        at_mass = weights
    else if (present(use_mass)) then
        if (use_mass) then
            at_mass = molDB%at_mas
        else
            at_mass = f1
        end if
    else
        at_mass = molDB%at_mas
    end if

    if (present(mask)) then
        at_mask = mask
    else
        at_mask = .true.
    end if

    ! Main algorithm
    ! -- 1. get center of mass
    com_ref = center_of_mass_dim(molDB%n_at, at_crd_ref, at_mass)
    com_new = center_of_mass_dim(molDB%n_at, molDB%at_crd, at_mass)
    ! -- 2. build the quaternion matrix
    qmat = f0
    do ia = 1, molDB%n_at
        if (.not.at_mask(ia)) cycle
        fac = at_mass(ia)

        xref = at_crd_ref(1,ia) - com_ref(1)
        yref = at_crd_ref(2,ia) - com_ref(2)
        zref = at_crd_ref(3,ia) - com_ref(3)
        xnew = molDB%at_crd(1,ia) - com_new(1)
        ynew = molDB%at_crd(2,ia) - com_new(2)
        znew = molDB%at_crd(3,ia) - com_new(3)

        xy = xnew * yref
        xz = xnew * zref
        yx = ynew * xref
        yz = ynew * zref
        zx = znew * xref
        zy = znew * yref
        q0 = xnew**2 + ynew**2 + znew**2 + xref**2 + yref**2 + zref**2
        q1 = f2 * xnew * xref
        q2 = f2 * ynew * yref
        q3 = f2 * znew * zref
        qmat(1,1) = qmat(1,1) + fac * (q0 - q1 - q2 - q3)
        qmat(2,2) = qmat(2,2) + fac * (q0 - q1 + q2 + q3)
        qmat(3,3) = qmat(3,3) + fac * (q0 + q1 - q2 + q3)
        qmat(4,4) = qmat(4,4) + fac * (q0 + q1 + q2 - q3)
        qmat(2,1) = qmat(2,1) + f2 * fac * (yz - zy)
        qmat(3,1) = qmat(3,1) + f2 * fac * (zx - xz)
        qmat(4,1) = qmat(4,1) + f2 * fac * (xy - yx)
        qmat(3,2) = qmat(3,2) - f2 * fac * (xy + yx)
        qmat(4,2) = qmat(4,2) - f2 * fac * (xz + zx)
        qmat(4,3) = qmat(4,3) - f2 * fac * (yz + zy)
    end do

    ! -- 3. Diagonalize
    call xsyev('V', 'L', 4, qmat, 4, qeval, scratch, 16, info)
    if (info /= 0) then
        write(msg, '("xSyEV failed in row: ",i0)') info
        call runstat%raise_error('Unable to superpose structures', &
            details='Failure to diagonalize the quaternion matrix', &
            extra=trim(msg))
        return
    end if
    qevec = qmat
    ! We take the eigenvector corresponding to smalles eigenvalue
    ! By construction, it is the first one
    q0 = qevec(1,1)
    q1 = qevec(2,1)
    q2 = qevec(3,1)
    q3 = qevec(4,1)
    ! -- 4. Build rotation matrix
    q0q0 = q0**2
    q1q1 = q1**2
    q2q2 = q2**2
    q3q3 = q3**2
    rmat(1,1) = q0q0 + q1q1 - q2q2 - q3q3
    rmat(2,2) = q0q0 - q1q1 + q2q2 - q3q3
    rmat(3,3) = q0q0 - q1q1 - q2q2 + q3q3
    q0q1 = q0 * q1
    q0q2 = q0 * q2
    q0q3 = q0 * q3
    q1q2 = q1 * q2
    q1q3 = q1 * q3
    q2q3 = q2 * q3
    rmat(1,2) = f2 * (q1q2 - q0q3)
    rmat(2,1) = f2 * (q1q2 + q0q3)
    rmat(1,3) = f2 * (q1q3 + q0q2)
    rmat(3,1) = f2 * (q1q3 - q0q2)
    rmat(2,3) = f2 * (q2q3 - q0q1)
    rmat(3,2) = f2 * (q2q3 + q0q1)

    do ia = 1, molDB%n_at
        xnew = molDB%at_crd(1,ia) - com_new(1)
        ynew = molDB%at_crd(2,ia) - com_new(2)
        znew = molDB%at_crd(3,ia) - com_new(3)
        molDB%at_crd(1,ia) = rmat(1,1)*xnew + rmat(2,1)*ynew + rmat(3,1)*znew
        molDB%at_crd(2,ia) = rmat(1,2)*xnew + rmat(2,2)*ynew + rmat(3,2)*znew
        molDB%at_crd(3,ia) = rmat(1,3)*xnew + rmat(2,3)*ynew + rmat(3,3)*znew
        molDB%at_crd(:,ia) = molDB%at_crd(:,ia) + com_ref
    end do

    if (present(rot_mat)) rot_mat = rmat
    if (present(trans_vec)) trans_vec = com_ref - com_new

end subroutine superpose_db

! ======================================================================

subroutine superpose_dim(n_at, at_crd_new, at_crd_ref, weights, mask, &
                         rot_mat, trans_vec)
    !! Superpose a geometrical structure onto another.
    !!
    !! The routine tries to find the best overlap between a reference
    !! structure and a structure to superpose.  The coordinates can be
    !! weighted, each atom having its own weight specification.
    !! A mask can be provided so the overlap is only done on the
    !! selected atoms.
    !!
    !! The routine use the quaternion method, as documented in:
    !!
    !! 1. G.R. Kneller, Mol. Sim. 7, 113-119 (1991)
    !! 2. G.R. Kneller, J. Chim. Phys. 88, 2709-2715 (1991)
    !!
    !! @warning "Limitations"
    !! `mask`, `weights` and the atomic coordinates must have the same
    !! "leading" dimensions.
    !! @endwarning
    !!
    !! @note "Version"
    !! This version takes data arrays and explicit size specifications.
    !! No check on arrays size is done besides the explicit
    !! specification.
    !! @endnote
    integer, intent(in) :: n_at
    !! Number of atoms.
    real(realwp), dimension(3,n_at), intent(inout) :: at_crd_new
    !! Atomic coordinates to transform, rotated on output.
    real(realwp), dimension(3,n_at), intent(in) :: at_crd_ref
    !! Atomic coordinates (reference).
    real(realwp), dimension(n_at), intent(in), optional :: weights
    !! Atomic weights (e.g., masses).
    logical, dimension(n_at), intent(in), optional :: mask
    !! Mask, specifying which atom to include in superposition
    real(realwp), dimension(3,3), intent(out), optional :: rot_mat
    !! Rotation matrix from original to final orientation.
    real(realwp), dimension(3), intent(out), optional :: trans_vec
    !! Translation vector from original to final position.

    integer :: ia, info
    real(realwp) :: fac, q0, q0q0, q0q1, q0q2, q0q3, q1, q1q1, q1q2, q1q3, &
        q2, q2q2, q2q3, q3, q3q3, xnew, xref, xy, xz, ynew, yref, yx, yz, &
        znew, zref, zx, zy
    real(realwp), dimension(3) :: com_new, com_ref
    real(realwp), dimension(4) :: qeval
    real(realwp), dimension(16) :: scratch
    real(realwp), dimension(3,3) :: rmat
    real(realwp), dimension(4,4) :: qmat, qevec
    character(len=256) :: msg

    real(realwp), dimension(n_at) :: at_mass
    logical, dimension(n_at) :: at_mask

    ! Check optional input arrays and set alternatives if needed
    if (present(weights)) then
        at_mass = weights
    else
        at_mass = f1
    end if

    if (present(mask)) then
        at_mask = mask
    else
        at_mask = .true.
    end if

    ! Main algorithm
    ! -- 1. get center of mass
    com_ref = center_of_mass_dim(n_at, at_crd_ref, at_mass)
    com_new = center_of_mass_dim(n_at, at_crd_new, at_mass)
    ! -- 2. build the quaternion matrix
    qmat = f0
    do ia = 1, n_at
        if (.not.at_mask(ia)) cycle
        fac = at_mass(ia)

        xref = at_crd_ref(1,ia) - com_ref(1)
        yref = at_crd_ref(2,ia) - com_ref(2)
        zref = at_crd_ref(3,ia) - com_ref(3)
        xnew = at_crd_new(1,ia) - com_new(1)
        ynew = at_crd_new(2,ia) - com_new(2)
        znew = at_crd_new(3,ia) - com_new(3)

        xy = xnew * yref
        xz = xnew * zref
        yx = ynew * xref
        yz = ynew * zref
        zx = znew * xref
        zy = znew * yref
        q0 = xnew**2 + ynew**2 + znew**2 + xref**2 + yref**2 + zref**2
        q1 = f2 * xnew * xref
        q2 = f2 * ynew * yref
        q3 = f2 * znew * zref
        qmat(1,1) = qmat(1,1) + fac * (q0 - q1 - q2 - q3)
        qmat(2,2) = qmat(2,2) + fac * (q0 - q1 + q2 + q3)
        qmat(3,3) = qmat(3,3) + fac * (q0 + q1 - q2 + q3)
        qmat(4,4) = qmat(4,4) + fac * (q0 + q1 + q2 - q3)
        qmat(2,1) = qmat(2,1) + f2 * fac * (yz - zy)
        qmat(3,1) = qmat(3,1) + f2 * fac * (zx - xz)
        qmat(4,1) = qmat(4,1) + f2 * fac * (xy - yx)
        qmat(3,2) = qmat(3,2) - f2 * fac * (xy + yx)
        qmat(4,2) = qmat(4,2) - f2 * fac * (xz + zx)
        qmat(4,3) = qmat(4,3) - f2 * fac * (yz + zy)
    end do

    ! -- 3. Diagonalize
    call xsyev('V', 'L', 4, qmat, 4, qeval, scratch, 16, info)
    if (info /= 0) then
        write(msg, '("xSyEV failed in row: ",i0)') info
        call runstat%raise_error('Unable to superpose structures', &
            details='Failure to diagonalize the quaternion matrix', &
            extra=trim(msg))
        return
    end if
    qevec = qmat
    ! We take the eigenvector corresponding to smalles eigenvalue
    ! By construction, it is the first one
    q0 = qevec(1,1)
    q1 = qevec(2,1)
    q2 = qevec(3,1)
    q3 = qevec(4,1)
    ! -- 4. Build rotation matrix
    q0q0 = q0**2
    q1q1 = q1**2
    q2q2 = q2**2
    q3q3 = q3**2
    rmat(1,1) = q0q0 + q1q1 - q2q2 - q3q3
    rmat(2,2) = q0q0 - q1q1 + q2q2 - q3q3
    rmat(3,3) = q0q0 - q1q1 - q2q2 + q3q3
    q0q1 = q0 * q1
    q0q2 = q0 * q2
    q0q3 = q0 * q3
    q1q2 = q1 * q2
    q1q3 = q1 * q3
    q2q3 = q2 * q3
    rmat(1,2) = f2 * (q1q2 - q0q3)
    rmat(2,1) = f2 * (q1q2 + q0q3)
    rmat(1,3) = f2 * (q1q3 + q0q2)
    rmat(3,1) = f2 * (q1q3 - q0q2)
    rmat(2,3) = f2 * (q2q3 - q0q1)
    rmat(3,2) = f2 * (q2q3 + q0q1)

    do ia = 1, n_at
        xnew = at_crd_new(1,ia) - com_new(1)
        ynew = at_crd_new(2,ia) - com_new(2)
        znew = at_crd_new(3,ia) - com_new(3)
        at_crd_new(1,ia) = rmat(1,1)*xnew + rmat(2,1)*ynew + rmat(3,1)*znew
        at_crd_new(2,ia) = rmat(1,2)*xnew + rmat(2,2)*ynew + rmat(3,2)*znew
        at_crd_new(3,ia) = rmat(1,3)*xnew + rmat(2,3)*ynew + rmat(3,3)*znew
        at_crd_new(:,ia) = at_crd_new(:,ia) + com_ref
    end do

    if (present(rot_mat)) rot_mat = rmat
    if (present(trans_vec)) trans_vec = com_ref - com_new

end subroutine superpose_dim

! ======================================================================

end module geometry