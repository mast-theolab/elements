module vibrational
    !! A module storing procedures related to vibrations.

    use numeric, only: f0, f1, f10m1, f2, realwp, small
    use lapack_drv, only: xsyev
    use physics, only: phys_conv, boltzmann, slight, planck
    use geometry, only: Eckart_orient
    use math, only: operator(.x.)
    use datatypes, only: MoleculeDB, VibrationsDB
    use blas_drv, only: xgemm
    use lapack_drv, only: xgeqrf, xlasrt, xorgqr
    use exception, only: runstat

    implicit none

    private
    public :: boltz_pop_max_quanta, build_modes, full_boltz_pop, set_orientation

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

        module subroutine build_modes_db_lt( &
                F_cart, molDB, vibDB, is_weighted, remove_rottrans, set_nvib, &
                set_Lmat, set_freq, set_Lmweigh, set_redmas)
            !! Interface routine to call main algorithm with Cartesian force
            !! constants provided as a rank-1 array containing the
            !! lower-triangular part of the matrix and molecular data as an
            !! object.  The data are stored in a vibrational data object.
            !! The force constants can be mass-weighted or not.
            real(realwp), dimension(:), intent(in) :: F_cart
            !! Square force constant matrix, mass_weighted or not
            class(MoleculeDB), intent(in) :: molDB
            !! Molecule database.
            type(VibrationsDB), intent(out) :: vibDB
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

        end subroutine build_modes_db_lt

        module subroutine build_modes_db_sq( &
            F_cart, molDB, vibDB, is_weighted, remove_rottrans, set_nvib, &
            set_Lmat, set_freq, set_Lmweigh, set_redmas)
            !! Interface routine to call main algorithm with Cartesian force
            !! constants provided as a square matrix and molecular data as an
            !! object.  The data are stored in a vibrational data object.
            !! The force constants can be mass-weighted or not.
            real(realwp), dimension(:,:), intent(in) :: F_cart
            !! Square force constant matrix, mass_weighted or not
            class(MoleculeDB), intent(in) :: molDB
            !! Molecule database.
            type(VibrationsDB), intent(out) :: vibDB
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

        end subroutine build_modes_db_sq

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

        module subroutine build_modes_moldb_lt( &
                F_cart, molDB, is_weighted, remove_rottrans, n_vib, L_mat, &
                freq, L_mweigh, red_mass)
            !! Interface routine to call main algorithm with Cartesian force
            !! constants provided as a rank-1 array containing the
            !! lower-triangular part of the matrix and molecular data as an
            !! object.
            !! The force constants can be mass-weighted or not.
            real(realwp), dimension(:), intent(in) :: F_cart
            !! Square force constant matrix, mass_weighted or not
            class(MoleculeDB), intent(in) :: molDB
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

        end subroutine build_modes_moldb_lt

        module subroutine build_modes_moldb_sq( &
            F_cart, molDB, is_weighted, remove_rottrans, n_vib, L_mat, freq, &
            L_mweigh, red_mass)
            !! Interface routine to call main algorithm with Cartesian force
            !! constants provided as a square matrix and molecular data as an
            !! object.
            !! The force constants can be mass-weighted or not.
            real(realwp), dimension(:,:), intent(in) :: F_cart
            !! Square force constant matrix, mass_weighted or not
            class(MoleculeDB), intent(in) :: molDB
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

        end subroutine build_modes_moldb_sq

    end interface build_modes

! ----------------------------------------------------------------------

    interface boltz_pop_max_quanta
        !! Compute the maximum number of quanta for each mode set in `nq_index`.
        module procedure boltz_pop_max_quanta_arr, boltz_pop_max_quanta_dim, &
            boltz_pop_max_quanta_db, boltz_pop_max_quanta_db_dim
    end interface boltz_pop_max_quanta

! ----------------------------------------------------------------------

    interface full_boltz_pop
        !! Compute the total Boltzmann population using analytical sums.
        module procedure full_boltz_pop_arr, full_boltz_pop_db
    end interface full_boltz_pop

! ----------------------------------------------------------------------

    interface set_orientation
        !! Set orientation of normal modes so the largest component is positive
        module procedure &
            set_orientation_arr, set_orientation_dim, set_orientation_db
    end interface set_orientation

contains

! ======================================================================

subroutine boltz_pop_max_quanta_arr(freq, nq_index, nq_max, temperature, &
                                    pop_min, ignore_0, list_states, list_pops)
    !! Compute the maximum number of quanta for each mode set in
    !! `nq_index`.
    !! Optionally, the list of states and their populations can be
    !! returned if `list_states` and `list_pops` are provided.
    !! If `ignore_0` is true, only states with all modes excited with
    !! at least 1 quantum are considered, states with less excited modes
    !! are ignored.
    !!
    !! @warning
    !! The maximum number of quanta reported is for the case where all
    !! chosen modes are excited by at least one quantum.
    !! The list of states does not include higher-energy states whose
    !! population is lower than `pop_min`.
    !! @endwarning
    !!
    !! @note "version"
    !! This version takes data arrays, recovering dimensions from their
    !! shape.
    !! Dimensions are recovered from the arrays.
    !! @endnote
    real(realwp), dimension(:), intent(in) :: freq
    !! Harmonic frequencies (in cm^-1).
    integer, dimension(:), intent(in) :: nq_index
    !! Index of the modes to be excited.  It should not contain null values.
    integer, dimension(:), intent(out) :: nq_max
    !! Maximum number of quanta achievable for each mode.
    real(realwp), intent(in), optional :: temperature
    !! Temperature, in K.
    real(realwp), intent(in), optional :: pop_min
    !! Minimum population of a state compared to ground state for inclusion.
    logical, intent(in), optional :: ignore_0
    !! If true, ignore sub-states with less excited modes in `list_states`.
    integer, dimension(:,:), allocatable, intent(out), optional :: list_states
    !! List of states with population > `pop_min`.
    !! If `ignore_0=.true.`, all sub-states are included, otherwise only states
    !! with same number of quanta.
    real(realwp), dimension(:), allocatable, intent(out), optional :: list_pops
    !! List of populations with respect to vibrational ground state.

    integer :: n_vib, n_modes

    n_vib = size(freq)
    n_modes = size(nq_index)
    if (size(nq_max) < n_modes) then
        call runstat%raise_error('inconsistency in input arrays', cat='dev', &
                                 source='boltz_pop_max_quanta_arr')
        return
    end if

    call boltz_pop_max_quanta_dim(n_vib, n_modes, freq, nq_index, nq_max, &
                                  temperature, pop_min, ignore_0, &
                                  list_states, list_pops)

end subroutine boltz_pop_max_quanta_arr

! ======================================================================

subroutine boltz_pop_max_quanta_db(vibDB, nq_index, nq_max, temperature, &
                                   pop_min, ignore_0, list_states, list_pops)
    !! Compute the maximum number of quanta for each mode set in
    !! `nq_index`.
    !! Optionally, the list of states and their populations can be
    !! returned if `list_states` and `list_pops` are provided.
    !! If `ignore_0` is true, only states with all modes excited with
    !! at least 1 quantum are considered, states with less excited modes
    !! are ignored.
    !!
    !! @warning
    !! The maximum number of quanta reported is for the case where all
    !! chosen modes are excited by at least one quantum.
    !! The list of states does not include higher-energy states whose
    !! population is lower than `pop_min`.
    !! @endwarning
    !!
    !! @note "version
    !! This version takes a VibrationsDB object.
    !! @endnote
    class(VibrationsDB), intent(in) :: vibDB
    !! VibrationsDB instance.
    integer, dimension(:), intent(in) :: nq_index
    !! Index of the modes to be excited.  It should not contain null values.
    integer, dimension(:), intent(out) :: nq_max
    !! Maximum number of quanta achievable for each mode.
    real(realwp), intent(in), optional :: temperature
    !! Temperature, in K.
    real(realwp), intent(in), optional :: pop_min
    !! Minimum population of a state compared to ground state for inclusion.
    logical, intent(in), optional :: ignore_0
    !! If true, ignore sub-states with less excited modes in `list_states`.
    integer, dimension(:,:), allocatable, intent(out), optional :: list_states
    !! List of states with population > `pop_min`.
    !! If `ignore_0=.true.`, all sub-states are included, otherwise only states
    !! with same number of quanta.
    real(realwp), dimension(:), allocatable, intent(out), optional :: list_pops
    !! List of populations with respect to vibrational ground state.

    integer :: n_modes

    n_modes = size(nq_index)
    if (size(nq_max) < n_modes) then
        call runstat%raise_error('inconsistency in input arrays', cat='dev', &
                                 source='boltz_pop_max_quanta_arr')
        return
    end if

    call boltz_pop_max_quanta_dim(vibDB%n_vib, n_modes, vibDB%freq, nq_index, &
                                  nq_max, temperature, pop_min, ignore_0, &
                                  list_states, list_pops)

end subroutine boltz_pop_max_quanta_db

! ======================================================================

subroutine boltz_pop_max_quanta_db_dim(vibDB, n_modes, nq_index, nq_max, &
                                       temperature, pop_min, ignore_0, &
                                       list_states, list_pops)
    !! Compute the maximum number of quanta for each mode set in
    !! `nq_index`.
    !! Optionally, the list of states and their populations can be
    !! returned if `list_states` and `list_pops` are provided.
    !! If `ignore_0` is true, only states with all modes excited with
    !! at least 1 quantum are considered, states with less excited modes
    !! are ignored.
    !!
    !! @warning
    !! The maximum number of quanta reported is for the case where all
    !! chosen modes are excited by at least one quantum.
    !! The list of states does not include higher-energy states whose
    !! population is lower than `pop_min`.
    !! @endwarning
    !!
    !! @note "version
    !! This version takes a VibrationsDB object and expects the number
    !! of excited modes to be provided.
    !! @endnote
    class(VibrationsDB), intent(in) :: vibDB
    !! VibrationsDB instance.
    integer, intent(in) :: n_modes
    !! Number of excited modes
    integer, dimension(n_modes), intent(in) :: nq_index
    !! Index of the modes to be excited.  It should not contain null values.
    integer, dimension(n_modes), intent(out) :: nq_max
    !! Maximum number of quanta achievable for each mode.
    real(realwp), intent(in), optional :: temperature
    !! Temperature, in K.
    real(realwp), intent(in), optional :: pop_min
    !! Minimum population of a state compared to ground state for inclusion.
    logical, intent(in), optional :: ignore_0
    !! If true, ignore sub-states with less excited modes in `list_states`.
    integer, dimension(:,:), allocatable, intent(out), optional :: list_states
    !! List of states with population > `pop_min`.
    !! If `ignore_0=.true.`, all sub-states are included, otherwise only states
    !! with same number of quanta.
    real(realwp), dimension(:), allocatable, intent(out), optional :: list_pops
    !! List of populations with respect to vibrational ground state.

    call boltz_pop_max_quanta_dim(vibDB%n_vib, n_modes, vibDB%freq, nq_index, &
                                  nq_max, temperature, pop_min, ignore_0, &
                                  list_states, list_pops)

end subroutine boltz_pop_max_quanta_db_dim

! ======================================================================

subroutine boltz_pop_max_quanta_dim(n_vib, n_modes, freq, nq_index, nq_max, &
                                    temperature, pop_min, ignore_0, &
                                    list_states, list_pops)
    !! Compute the maximum number of quanta for each mode set in
    !! `nq_index`.
    !! Optionally, the list of states and their populations can be
    !! returned if `list_states` and `list_pops` are provided.
    !! If `ignore_0` is true, only states with all modes excited with
    !! at least 1 quantum are considered, states with less excited modes
    !! are ignored.
    !!
    !! @warning
    !! The maximum number of quanta reported is for the case where all
    !! chosen modes are excited by at least one quantum.
    !! The list of states does not include higher-energy states whose
    !! population is lower than `pop_min`.
    !! @endwarning
    !!
    !! @note "version
    !! This version takes data arrays and dimensions.
    !! @endnote
    integer, intent(in) :: n_vib
    !! Total number of vibrational modes.
    integer, intent(in) :: n_modes
    !! Number of excited modes.
    real(realwp), dimension(n_vib), intent(in) :: freq
    !! Harmonic frequencies (in cm^-1).
    integer, dimension(n_modes), intent(in), target :: nq_index
    !! Index of the modes to be excited.  It should not contain null values.
    integer, dimension(n_modes), intent(out) :: nq_max
    !! Maximum number of quanta achievable for each mode.
    real(realwp), intent(in), optional :: temperature
    !! Temperature, in K.
    real(realwp), intent(in), optional :: pop_min
    !! Minimum population of a state compared to ground state for inclusion.
    logical, intent(in), optional :: ignore_0
    !! If true, ignore sub-states with less excited modes in `list_states`.
    integer, dimension(:,:), allocatable, intent(out), optional :: list_states
    !! List of states with population > `pop_min`.
    !! If `ignore_0=.true.`, all sub-states are included, otherwise only states
    !! with same number of quanta.
    real(realwp), dimension(:), allocatable, intent(out), optional :: list_pops
    !! List of populations with respect to vibrational ground state.

    integer :: i, ibase, i_state, n_states
    integer, dimension(:), allocatable :: nqi
    real(realwp) :: bz_kT, E_base, E_max, E_val, p_min, T_val
    logical :: do_pops, do_states, include_0
    character(len=128) :: msg

    ! Set the work temperature
    if (present(temperature)) then
        if (temperature <= f0) then
            call runstat%raise_error( &
                'Invalid temperature', &
                details='Failed to set populated vibrational states')
            return
        end if
        T_val = temperature
    else
        T_val = real(300, kind=realwp)
    end if

    ! Set the minimum population
    if (present(pop_min)) then
        if (temperature <= f0) then
            call runstat%raise_error( &
                'Invalid value for the minimum population', &
                details='Failed to set populated vibrational states')
            return
        end if
        p_min = pop_min
    else
        p_min = f10m1
    end if

    ! Check if all modes are set in `nq_index`.
    if (any(nq_index <= 0)) then
        call runstat%raise_error( &
            'Invalid list of indexes', &
            'Mode indexes cannot be null or negative')
        return
    else
        nq_max = 0
        ! We should have taken care of null or negative indexes
        ! Now check duplicate.
        allocate(nqi(n_vib))
        nqi = 0
        do i = 1, n_modes
            nqi(nq_index(i)) = nqi(nq_index(i)) + 1
        end do
        if (any(nqi > 1)) then
            call runstat%raise_error( &
                'Invalid list of indexes', &
                'Duplicate indexes in list of indexes')
            return
        end if
        deallocate(nqi)
    end if

    ! Compute maximum energy
    bz_kT = T_val * boltzmann/(planck*slight)
    E_max = -bz_kT*log(p_min)

    ! Compute the basic energy
    E_base = sum(freq(nq_index))

    ! Find the maximum number of states
    ! If E_base > E_max, then no states possible, stopping.
    if (E_base <= E_max) then
        nq_max = int((E_max - E_base)/freq(nq_index)) + 1
    else
        nq_max = 0
        return
    end if

    do_states = present(list_states)
    do_pops = present(list_pops)
    if (do_states .or. do_pops) then
        if (present(ignore_0)) then
            include_0 = .not.ignore_0
        else
            include_0 = .true.
        end if
        if (include_0) then
            n_states = product(nq_max+1)
        else
            n_states = product(nq_max)
        end if
        if (do_states) allocate(list_states(n_modes,n_states))
        if (do_pops) allocate(list_pops(n_states))
        allocate(nqi(n_modes))
        if (include_0) then
            ibase = 0
        else
            ibase = 1
        end if
        nqi = ibase
        i_state = 1
        if (do_states) list_states(:,i_state) = nqi
        if (do_pops) list_pops(i_state) = exp(-sum(nqi*freq(nq_index))/bz_kT)
        i = 1
        state_counter: do
            if (nqi(i) == nq_max(i)) then
                do while (nqi(i) == nq_max(i))
                    i = i + 1
                    if (i > n_modes) exit state_counter
                end do
                nqi(:i-1) = ibase
                nqi(i) = nqi(i) + 1
                i = 1
            else
                nqi(i) = nqi(i) + 1
            end if
            E_val = sum(nqi*freq(nq_index))
            if (E_val <= E_max) then
                i_state = i_state + 1
                if (do_states) list_states(:,i_state) = nqi
                if (do_pops) list_pops(i_state) = exp(-E_val/bz_kT)
            end if
        end do state_counter
        if (i_state > n_states) then
            write(msg, '(i0," states expected but ",i0," built.")') &
                n_states, i_state
            call runstat%raise_error( &
                'Inconsistency in listing of populated states', &
                msg)
            return
        else if(i_state < n_states) then
            list_states = list_states(:,:i_state)
            list_pops = list_pops(:i_state)
        end if
    end if

end subroutine boltz_pop_max_quanta_dim

! ======================================================================

function full_boltz_pop_arr(freq, temperature, incl_ZPVE) result(pop)
    !! Compute the total Boltzmann population at a given temperature.
    !! The function uses analytic sums true within the harmonic
    !! approximation.
    !! The total sum can be computed considering the (electronic)
    !! minimum of the potential energy surface or using the ground
    !! vibrational state as reference.
    !!
    !! Within the harmonic approximation, the total population is given
    !! by,
    !! \[
    !!     Z = \sum_s
    !!         e^{-\sum_{i=1}^N \frac{\hbar \omega_i}{k_B T}
    !!         (n^s_i + \frac{1}{2})}
    !!     = \sum_s \prod_{i=1}^N e^{-\frac{\hbar \omega_i}{k_B T}
    !!                            (n^s_i + \frac{1}{2})}
    !!     = \sum_{n_1=0}^\infty \sum_{n_2=0}^\infty \ldots
    !!     \sum_{n_N=0}^\infty
    !!     \prod_{i=1}^N e^{-\frac{\hbar \omega_i}{k_B T}
    !!                     (n_i + \frac{1}{2})}
    !! \]
    !! which can be simplified to
    !! \[
    !!     Z = \prod_{i=1}^N \sum_{n_i=0}^\infty
    !!      e^{-\frac{\hbar \omega_i}{k_B T}(n_i + \frac{1}{2})}
    !!     = \prod_{i=1}^N Z_i
    !! \]
    !! with,
    !! \[  
    !!     Z_i = \sum_{n_i=0}^\infty
    !!      e^{-\frac{\hbar \omega_i}{k_B T}(n_i + \frac{1}{2})}
    !!     = \frac{e^{-\frac{\hbar \omega_i}{2 k_B T}}}
    !!         {1 - e^{-\frac{\hbar \omega_i}{k_B T}}}
    !! \]
    !!
    !! @note "version"
    !! Array version, dimensions are not necessary since elemental
    !! operations used.
    !! @endnote
    real(realwp), dimension(:), intent(in) :: freq
        !! Harmonic frequencies (in cm^-1).
    real(realwp), intent(in), optional :: temperature
        !! Temperature, in K.
    logical, intent(in), optional :: incl_ZPVE
        !! Include ZPVE in calculation of the total population (default: false)
    real(realwp) :: pop
        !! Total Boltzmann population

    real(realwp) :: bz_kT, fnum
    logical :: add_zpve

    if (temperature < f0) then
        call runstat%raise_error( &
            'Invalid temperature', &
            details='Failed to compute total Boltzmann population')
        return
    end if

    bz_kT = temperature * boltzmann/(planck*slight)

    if (present(incl_ZPVE)) then
        add_zpve = incl_ZPVE
    else
        add_zpve = .false.
    end if

    if (add_zpve) then
        fnum = product(exp(-freq/(f2*bz_kT)))
    else
        fnum = f1
    end if

    pop = product(f1/(f1 - exp(-freq/(bz_kT))))

end function full_boltz_pop_arr

! ======================================================================

function full_boltz_pop_db(vibDB, temperature, incl_ZPVE) result(pop)
    !! Compute the total Boltzmann population at a given temperature.
    !! The function uses analytic sums true within the harmonic
    !! approximation.
    !! The total sum can be computed considering the (electronic)
    !! minimum of the potential energy surface or using the ground
    !! vibrational state as reference.
    !!
    !! See [[full_boltz_pop_arr]] for implementation details.
    !!
    !! @note "version
    !! This version takes a VibrationsDB object.
    !! @endnote
    class(VibrationsDB), intent(in) :: vibDB
        !! VibrationsDB instance.
    real(realwp), intent(in), optional :: temperature
        !! Temperature, in K.
    logical, intent(in), optional :: incl_ZPVE
        !! Include ZPVE in calculation of the total population (default: false)
    real(realwp) :: pop
        !! Total Boltzmann population
    
    pop = full_boltz_pop_arr(vibDB%freq, temperature, incl_ZPVE)

end function full_boltz_pop_db

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

subroutine set_orientation_db(vibDB)
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
    class(VibrationsDB), intent(inout) :: vibDB
    !! VibrationsDB instance.

    integer :: i, ia
    real(realwp) :: x
    logical :: is_neg

    if (allocated(vibDB%L_mat)) then
        do i = 1, vibDB%n_vib
            x = f0
            is_neg = .false.
            do ia = 1, size(vibDB%L_mat, 1)
                if (abs(vibDB%L_mat(ia,i)) > x) then
                    is_neg = vibDB%L_mat(ia,i) < f0
                    x = abs(vibDB%L_mat(ia,i)) + small
                end if
            end do
            if (is_neg) vibDB%L_mat(:,i) = -vibDB%L_mat(:,i)
        end do
    end if

    if (allocated(vibDB%L_mwg)) then
        do i = 1, vibDB%n_vib
            x = f0
            is_neg = .false.
            do ia = 1, size(vibDB%L_mwg, 1)
                if (abs(vibDB%L_mwg(ia,i)) > x) then
                    is_neg = vibDB%L_mwg(ia,i) < f0
                    x = abs(vibDB%L_mwg(ia,i)) + small
                end if
            end do
            if (is_neg) vibDB%L_mwg(:,i) = -vibDB%L_mwg(:,i)
        end do
    end if

end subroutine set_orientation_db

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

end module vibrational
