module vibrational
    !! A module storing procedures related to vibrations.

    use numeric, only: f0, f1, realwp, small
    use lapack_drv, only: xsyev
    use physics, only: phys_conv
    use geometry, only: Eckart_orient
    use math, only: operator(.x.)
    use datatypes, only: MoleculeDB, VibrationsDB
    use blas_drv, only: xgemm
    use lapack_drv, only: xgeqrf, xlasrt, xorgqr
    use exception, only: runstat

    implicit none

    private
    public :: build_modes, set_orientation

! ----------------------------------------------------------------------

    interface build_modes
        !! Build the normal modes from the Cartesian force constant matrix.
        module subroutine build_modes_arr_lt( &
                F_cart, at_crd, at_mass, is_weighted, remove_rottrans, n_vib, &
                L_mat, freq, L_mweigh, red_mass)
            !! Interface routine to call main algorithm with Cartesian force
            !! constants provided as a rank-1 array containing the
            !! lower-triangular part of the matrix and molecular data as arrays.
            !! The force constants can be mass-weighted or not.
            !!
            !! In this version, the number of atoms is read from the input
            !! quantities.
            real(realwp), dimension(:), intent(in) :: F_cart
            !! Square force constant matrix, mass_weighted or not
            real(realwp), dimension(:,:), intent(in) :: at_crd
            !! Atomic coordinates.
            real(realwp), dimension(:), intent(in) :: at_mass
            !! Atomic masses.
            logical, intent(in), optional :: is_weighted
            !! Cartesian force constants contain mass-weighted quantities.
            logical, intent(in), optional :: remove_rottrans
            !! Remove residual rotations and translations.
            integer, intent(inout), optional :: n_vib
            !! Number of normal modes.
            !! If set in input, the routine checks that its results are consistent.
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

        end subroutine build_modes_arr_lt

        module subroutine build_modes_arr_sq( &
                F_cart, at_crd, at_mass, is_weighted, remove_rottrans, n_vib, &
                L_mat, freq, L_mweigh, red_mass)
            !! Interface routine to call main algorithm with Cartesian force
            !! constants provided as a square matrix and molecular data as arrays.
            !! The force constants can be mass-weighted or not.
            !!
            !! In this version, the number of atoms is read from the input
            !! quantities.
            real(realwp), dimension(:,:), intent(in) :: F_cart
            !! Square force constant matrix, mass_weighted or not
            real(realwp), dimension(:,:), intent(in) :: at_crd
            !! Atomic coordinates.
            real(realwp), dimension(:), intent(in) :: at_mass
            !! Atomic masses.
            logical, intent(in), optional :: is_weighted
            !! Cartesian force constants contain mass-weighted quantities.
            logical, intent(in), optional :: remove_rottrans
            !! Remove residual rotations and translations.
            integer, intent(inout), optional :: n_vib
            !! Number of normal modes.
            !! If set in input, the routine checks that its results are consistent.
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

        end subroutine build_modes_arr_sq

        module subroutine build_modes_dim_lt( &
                n_at, F_cart, at_crd, at_mass, is_weighted, remove_rottrans, &
                n_vib, L_mat, freq, L_mweigh, red_mass)
            !! Interface routine to call main algorithm with Cartesian force
            !! constants provided as a rank-1 array containing the
            !! lower-triangular part of the matrix and molecular data as arrays.
            !! The force constants can be mass-weighted or not.
            integer, intent(in) :: n_at
            !! Number of atoms.
            real(realwp), dimension(:), intent(in) :: F_cart
            !! Square force constant matrix, mass_weighted or not
            real(realwp), dimension(3,n_at), intent(in) :: at_crd
            !! Atomic coordinates.
            real(realwp), dimension(n_at), intent(in) :: at_mass
            !! Atomic masses.
            logical, intent(in), optional :: is_weighted
            !! Cartesian force constants contain mass-weighted quantities.
            logical, intent(in), optional :: remove_rottrans
            !! Remove residual rotations and translations.
            integer, intent(inout), optional :: n_vib
            !! Number of normal modes.
            !! If set in input, the routine checks that its results are consistent.
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

        end subroutine build_modes_dim_lt

        module subroutine build_modes_dim_sq( &
                n_at, F_cart, at_crd, at_mass, is_weighted, remove_rottrans, &
                n_vib, L_mat, freq, L_mweigh, red_mass)
            !! Interface routine to call main algorithm with Cartesian force
            !! constants provided as a square matrix and molecular data as arrays.
            !! The force constants can be mass-weighted or not.
            integer, intent(in) :: n_at
            !! Number of atoms.
            real(realwp), dimension(3*n_at,3*n_at), intent(in) :: F_cart
            !! Square force constant matrix, mass_weighted or not
            real(realwp), dimension(3,n_at), intent(in) :: at_crd
            !! Atomic coordinates.
            real(realwp), dimension(n_at), intent(in) :: at_mass
            !! Atomic masses.
            logical, intent(in), optional :: is_weighted
            !! Cartesian force constants contain mass-weighted quantities.
            logical, intent(in), optional :: remove_rottrans
            !! Remove residual rotations and translations.
            integer, intent(inout), optional :: n_vib
            !! Number of normal modes.
            !! If set in input, the routine checks that its results are consistent.
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

        end subroutine build_modes_dim_sq

        module subroutine build_modes_mol_lt( &
                F_cart, mol, is_weighted, remove_rottrans, n_vib, L_mat, &
                freq, L_mweigh, red_mass)
            !! Interface routine to call main algorithm with Cartesian force
            !! constants provided as a rank-1 array containing the
            !! lower-triangular part of the matrix and molecular data as an
            !! object.
            !! The force constants can be mass-weighted or not.
            real(realwp), dimension(:), intent(in) :: F_cart
            !! Square force constant matrix, mass_weighted or not
            class(MoleculeDB), intent(in) :: mol
            !! Molecule database.
            logical, intent(in), optional :: is_weighted
            !! Cartesian force constants contain mass-weighted quantities.
            logical, intent(in), optional :: remove_rottrans
            !! Remove residual rotations and translations.
            integer, intent(inout), optional :: n_vib
            !! Number of normal modes.
            !! If set in input, the routine checks that its results are consistent.
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

        end subroutine build_modes_mol_lt

        module subroutine build_modes_mol_sq( &
            F_cart, mol, is_weighted, remove_rottrans, n_vib, L_mat, freq, &
            L_mweigh, red_mass)
            !! Interface routine to call main algorithm with Cartesian force
            !! constants provided as a square matrix and molecular data as an
            !! object.
            !! The force constants can be mass-weighted or not.
            real(realwp), dimension(:,:), intent(in) :: F_cart
            !! Square force constant matrix, mass_weighted or not
            class(MoleculeDB), intent(in) :: mol
            !! Molecule database.
            logical, intent(in), optional :: is_weighted
            !! Cartesian force constants contain mass-weighted quantities.
            logical, intent(in), optional :: remove_rottrans
            !! Remove residual rotations and translations.
            integer, intent(inout), optional :: n_vib
            !! Number of normal modes.
            !! If set in input, the routine checks that its results are consistent.
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

        end subroutine build_modes_mol_sq

        module subroutine build_modes_vib_lt( &
                F_cart, mol, vib, is_weighted, remove_rottrans, set_nvib, &
                set_Lmat, set_freq, set_Lmweigh, set_redmas)
            !! Interface routine to call main algorithm with Cartesian force
            !! constants provided as a rank-1 array containing the
            !! lower-triangular part of the matrix and molecular data as an
            !! object.  The data are stored in a vibrational data object.
            !! The force constants can be mass-weighted or not.
            real(realwp), dimension(:), intent(in) :: F_cart
            !! Square force constant matrix, mass_weighted or not
            class(MoleculeDB), intent(in) :: mol
            !! Molecule database.
            type(VibrationsDB), intent(out) :: vib
            !! Vibration database.
            logical, intent(in), optional :: is_weighted
            !! Cartesian force constants contain mass-weighted quantities.
            logical, intent(in), optional :: remove_rottrans
            !! Remove residual rotations and translations.
            logical, intent(in), optional :: set_nvib
            !! Set the number of normal modes.
            logical, intent(in), optional :: set_Lmat
            !! Set dimensionless eigenvectors matrix.
            logical, intent(in), optional :: set_freq
            !! Set the wavenumbers vector (in cm-1).
            !! Imaginary frequencies are reported as negative.
            logical, intent(in), optional :: set_Lmweigh
            !! Set dimensionless eigenvectors matrix.
            logical, intent(in), optional :: set_redmas
            !! Reduced mass of each vibration.

        end subroutine build_modes_vib_lt

        module subroutine build_modes_vib_sq( &
            F_cart, mol, vib, is_weighted, remove_rottrans, set_nvib, &
            set_Lmat, set_freq, set_Lmweigh, set_redmas)
            !! Interface routine to call main algorithm with Cartesian force
            !! constants provided as a square matrix and molecular data as an
            !! object.  The data are stored in a vibrational data object.
            !! The force constants can be mass-weighted or not.
            real(realwp), dimension(:,:), intent(in) :: F_cart
            !! Square force constant matrix, mass_weighted or not
            class(MoleculeDB), intent(in) :: mol
            !! Molecule database.
            type(VibrationsDB), intent(out) :: vib
            !! Vibration database.
            logical, intent(in), optional :: is_weighted
            !! Cartesian force constants contain mass-weighted quantities.
            logical, intent(in), optional :: remove_rottrans
            !! Remove residual rotations and translations.
            logical, intent(in), optional :: set_nvib
            !! Set the number of normal modes.
            logical, intent(in), optional :: set_Lmat
            !! Set dimensionless eigenvectors matrix.
            logical, intent(in), optional :: set_freq
            !! Set the wavenumbers vector (in cm-1).
            !! Imaginary frequencies are reported as negative.
            logical, intent(in), optional :: set_Lmweigh
            !! Set dimensionless eigenvectors matrix.
            logical, intent(in), optional :: set_redmas
            !! Reduced mass of each vibration.

        end subroutine build_modes_vib_sq

    end interface build_modes

! ----------------------------------------------------------------------
    interface set_orientation
        !! Set orientation of normal modes so the largest component is positive
        module procedure &
            set_orientation_arr, set_orientation_dim, set_orientation_vib
    end interface set_orientation

contains

! ======================================================================

subroutine set_orientation_arr(L_mat)
    !! Set the orientation of an array of modes.
    !!
    !! Set the orientation of an array of modes so the largest component
    !! is positive.
    !! If there is more than one component with the largest value, the
    !! first component is used to set the orientation.
    !!
    !! @note "version"
    !! This version takes a single array, finding the dimensions from
    !! its shape.
    !! @endnote
    real(realwp), dimension(:,:), intent(inout) :: L_mat
    !! List of normal modes.

    integer :: i, ia
    real(realwp) :: x
    logical :: is_neg

    do i = 1, size(L_mat, 2)
        x = f0
        is_neg = .false.
        do ia = 1, size(L_mat, 1)
            if (abs(L_mat(ia,i)) > x) then
                is_neg = L_mat(ia,i) < f0
                x = abs(L_mat(ia,i)) + small
            end if
        end do
        if (is_neg) L_mat(:,i) = -L_mat(:,i)
    end do

end subroutine set_orientation_arr

! ======================================================================

subroutine set_orientation_dim(n_at3, n_vib, L_mat)
    !! Set the orientation of an array of modes.
    !!
    !! Set the orientation of an array of modes so the largest component
    !! is positive.
    !! If there is more than one component with the largest value, the
    !! first component is used to set the orientation.
    !!
    !! @note "version"
    !! This version takes the explicit dimension of the array.
    !! @endnote
    integer, intent(in) :: n_at3
    !! Number of atomic coordinates.
    integer, intent(in) :: n_vib
    !! Number of normal modes.
    real(realwp), dimension(n_at3,n_vib), intent(inout) :: L_mat
    !! List of normal modes.

    integer :: i, ia
    real(realwp) :: x
    logical :: is_neg

    do i = 1, n_vib
        x = f0
        is_neg = .false.
        do ia = 1, n_at3
            if (abs(L_mat(ia,i)) > x) then
                is_neg = L_mat(ia,i) < f0
                x = abs(L_mat(ia,i)) + small
            end if
        end do
        if (is_neg) L_mat(:,i) = -L_mat(:,i)
    end do

end subroutine set_orientation_dim

! ======================================================================

subroutine set_orientation_vib(vib)
    !! Set the orientation of an array of modes.
    !!
    !! Set the orientation of an array of modes so the largest component
    !! is positive.
    !! If there is more than one component with the largest value, the
    !! first component is used to set the orientation.
    !!
    !! @note "version"
    !! This version sets the normal modes arrays populated in a
    !! vibrationalDB object.
    !! @endnote
    class(VibrationsDB), intent(inout) :: vib
    !! VibrationsDB instance.

    integer :: i, ia
    real(realwp) :: x
    logical :: is_neg

    if (allocated(vib%L_mat)) then
        do i = 1, vib%n_vib
            x = f0
            is_neg = .false.
            do ia = 1, size(vib%L_mat, 1)
                if (abs(vib%L_mat(ia,i)) > x) then
                    is_neg = vib%L_mat(ia,i) < f0
                    x = abs(vib%L_mat(ia,i)) + small
                end if
            end do
            if (is_neg) vib%L_mat(:,i) = -vib%L_mat(:,i)
        end do
    end if

    if (allocated(vib%L_mwg)) then
        do i = 1, vib%n_vib
            x = f0
            is_neg = .false.
            do ia = 1, size(vib%L_mwg, 1)
                if (abs(vib%L_mwg(ia,i)) > x) then
                    is_neg = vib%L_mwg(ia,i) < f0
                    x = abs(vib%L_mwg(ia,i)) + small
                end if
            end do
            if (is_neg) vib%L_mwg(:,i) = -vib%L_mwg(:,i)
        end do
    end if

end subroutine set_orientation_vib

! ======================================================================


end module vibrational