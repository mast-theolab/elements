submodule (vibronic) vibronic_Duschinsky 
    !! Submodule containing the Duschinsky matrix and shift calculation

    use blas_drv, only: xgemm
    use exception, only: runstat
    use numeric, only: f0, f1
    use physics, only: phys_conv
    use string, only: upcase

contains

! ======================================================================

module procedure Duschinsky_matrix_identity

    ! Local variables
    logical :: is_identity_local
    integer :: i

    is_identity_local = .true.
    if (present(is_identity)) is_identity_local = is_identity

    if (is_identity_local) then
        allocate(Jmat(n_vib, n_vib))

        Jmat = f0
        !$omp parallel do private(i)
        do i = 1, n_vib
            Jmat(i, i) = f1
        end do
        !$omp end parallel do

    else
        call runstat%raise_error('Duschinsky matrix cannot be constructed', &
            details='Only n_vib as parameter but is_identity is false', &
            source='Duschinsky_matrix_identity', cat='dev')
    end if

end procedure Duschinsky_matrix_identity

! ======================================================================

module procedure Duschinsky_matrix_arr

    ! Local variables
    integer :: i
    integer :: local_n_vib, local_n_at
    logical :: local_is_identity

    ! Sanity check
    if (size(L_mat1, 1) /= size(L_mat2, 1) &
        .or. size(L_mat1, 2) /= size(L_mat2, 2)) then
        call runstat%raise_error('L_mat1 and L_mat2 must have the same dimensions', &
            source='Duschinsky_matrix_arr', cat='dev')
    end if
    if (present(at_mass1)) then
        if (size(at_mass1) * 3 /= size(L_mat1, 1) ) then
            call runstat%raise_error('at_mass1 and L_mat1 must have 3*n_at elements', &
                source='Duschinsky_matrix_arr', cat='dev')
        end if
    end if
    if (present(at_mass2)) then
        if (size(at_mass2) * 3 /= size(L_mat2, 1) ) then
            call runstat%raise_error('at_mass2 and L_mat2 must have 3*n_at elements', &
                source='Duschinsky_matrix_arr', cat='dev')
        end if
    end if

    local_is_identity = .false.
    if (present(is_identity)) local_is_identity = is_identity

    local_n_vib = size(L_mat1, 2) 
    local_n_at = size(L_mat1, 1) / 3

    if (local_is_identity) then
        Jmat = Duschinsky_matrix_identity(local_n_vib, .true.)
    elseif (present(Ginv_mat)) then
        Jmat = Duschinsky_matrix_internal(local_n_vib, local_n_at, L_mat1, L_mat2, Ginv_mat)
    elseif (present(at_mass1) .and. present(at_mass2)) then
        Jmat = Duschinsky_matrix_orthogonal_with_mass(local_n_vib, local_n_at, L_mat1, L_mat2, at_mass1, at_mass2)
    else
        Jmat = Duschinsky_matrix_orthogonal(local_n_vib, local_n_at, L_mat1, L_mat2)
    end if

end procedure Duschinsky_matrix_arr

! ======================================================================

module procedure Duschinsky_matrix_db

    ! Local variables
    integer :: local_n_vib, local_n_at
    logical :: local_is_identity

    ! Sanity check
    if (.not.(vib1%loaded)) then
        call runstat%raise_error('vib1 VibrationsDB is not loaded', &
            source='Duschinsky_matrix_db', cat='dev')
    else if (.not.allocated(vib1%L_mat)) then
        call runstat%raise_error('vib1 L_mat is not loaded', &
            source='Duschinsky_matrix_db', cat='dev')
    end if

    if (.not.(vib2%loaded)) then
        call runstat%raise_error('vib2 VibrationsDB is not loaded', &
            source='Duschinsky_matrix_db', cat='dev')
    else if (.not.allocated(vib2%L_mat)) then
        call runstat%raise_error('vib2 L_mat is not loaded', &
            source='Duschinsky_matrix_db', cat='dev')
    end if

    if (present(mol1)) then
        if (.not.(mol1%loaded)) then
            call runstat%raise_error('mol1 MoleculeDB is not loaded', &
                source='Duschinsky_matrix_db', cat='dev')
        else if (.not.allocated(mol1%at_mas)) then
            call runstat%raise_error('mol1 mass is not loaded', &
                source='Duschinsky_matrix_db', cat='dev')
        end if
    end if

    if (present(mol2)) then
        if (.not.(mol2%loaded)) then
            call runstat%raise_error('mol2 MoleculeDB is not loaded', &
                source='Duschinsky_matrix_db', cat='dev')
        else if (.not.allocated(mol2%at_mas)) then
            call runstat%raise_error('mol2 mass is not loaded', &
                source='Duschinsky_matrix_db', cat='dev')
        end if
    end if

    local_is_identity = .false.
    if (present(is_identity)) local_is_identity = is_identity

    local_n_vib = vib1%n_vib
    if (local_n_vib /= vib2%n_vib) then
        call runstat%raise_error('vib1 and vib2 must have the same number of vibrational modes', &
            source='Duschinsky_matrix_db', cat='dev')
    end if

    local_n_at = size(vib1%L_mat, 1) / 3
    if (local_n_at /= size(vib2%L_mat, 1) / 3) then
        call runstat%raise_error('vib1 and vib2 must have the same number of atoms', &
            source='Duschinsky_matrix_db', cat='dev')
    end if

    if (local_is_identity) then
        Jmat = Duschinsky_matrix_identity(local_n_vib, .true.)
    elseif (present(Ginv_mat)) then
        Jmat = Duschinsky_matrix_internal(local_n_vib, local_n_at, vib1%L_mat, vib2%L_mat, Ginv_mat)
    elseif (present(mol1) .and. present(mol2)) then
        Jmat = Duschinsky_matrix_orthogonal_with_mass(local_n_vib, local_n_at, vib1%L_mat, vib2%L_mat, &
            mol1%at_mas, mol2%at_mas)
    else
        Jmat = Duschinsky_matrix_orthogonal(local_n_vib, local_n_at, vib1%L_mat, vib2%L_mat)
    end if

end procedure Duschinsky_matrix_db

! ======================================================================

module function Duschinsky_matrix_orthogonal(n_vib, n_at, L_mat1, L_mat2) result(Jmat) 
    !! Computes the Duschinsky matrix J for orthogonal coordinates
    !!
    !! The Duschinsky matrix can be calculated as:
    !!   J = L_mat1.T * L_mat2
    !!
    integer, intent(in) :: n_vib
    !! Number of vibrational modes
    integer, intent(in) :: n_at
    !! Number of atoms
    real(realwp), dimension(n_at * 3, n_vib), intent(in) :: L_mat1, L_mat2
    !! L matrix for the 1st and 2nd state
    real(realwp), dimension(:, :), allocatable :: Jmat

    allocate(Jmat(n_vib, n_vib))

    call xgemm('T', 'N', n_vib, n_vib, n_at * 3, f1, L_mat1, n_at * 3, &
        L_mat2, n_at * 3, f0, Jmat, n_vib)

end function Duschinsky_matrix_orthogonal

! ======================================================================

module function Duschinsky_matrix_orthogonal_with_mass(n_vib, n_at, L_mat1, L_mat2, at_mass1, at_mass2) result(Jmat) 
    !! Computes the Duschinsky matrix J for orthogonal coordinates with masses
    !!
    !! The Duschinsky matrix can be calculated as:
    !!   J = L_mat1.T * at_mass1^-0.5 * at_mass2^0.5 * L_mat2
    !!
    integer, intent(in) :: n_vib
    !! Number of vibrational modes
    integer, intent(in) :: n_at
    !! Number of atoms
    real(realwp), dimension(3 * n_at, n_vib), intent(in) :: L_mat1, L_mat2
    !! L matrix for the 1st and 2nd state
    real(realwp), dimension(n_at), intent(in) :: at_mass1, at_mass2
    !! Masses of the atoms
    real(realwp), dimension(:, :), allocatable :: Jmat

    ! Local
    real(realwp), dimension(:, :), allocatable :: tmp_mat
    integer :: i

    allocate(tmp_mat(3 * n_at, n_vib))

    !$omp parallel do private(i)
    do i = 1, n_at
        tmp_mat(3*(i-1):3*i,:) = L_mat1(3*(i-1):3*i,:) &
            / sqrt(phys_conv%au2amu(at_mass1(i), .true.)) &
            * sqrt(phys_conv%au2amu(at_mass2(i), .true.))
    end do
    !$omp end parallel do

    allocate(Jmat(n_vib, n_vib))

    call xgemm('T', 'N', n_vib, n_vib, n_at * 3, f1, tmp_mat, n_at * 3, &
        L_mat2, n_at * 3, f0, Jmat, n_vib)

    deallocate(tmp_mat)

end function Duschinsky_matrix_orthogonal_with_mass

! ======================================================================

module function Duschinsky_matrix_internal(n_vib, n_at, L_mat1, L_mat2, Ginv_mat) result(Jmat) 
    !! Computes the Duschinsky matrix J for internal coordinates
    !!
    !! @note "Not implemented"
    integer, intent(in) :: n_vib
    !! Number of vibrational modes
    integer, intent(in) :: n_at
    !! Number of atoms
    real(realwp), dimension(n_at * 3, n_vib), intent(in) :: L_mat1, L_mat2
    !! L matrix for the 1st and 2nd state
    real(realwp), dimension(n_at * 3, n_at * 3), intent(in) :: Ginv_mat
    !! Masses of the atoms
    real(realwp), dimension(:, :), allocatable :: Jmat

    call runstat%raise_error('Not implemented', &
        details='Construction of Duschinsky matrix in internal coordinates is not implemented', &
        source='Duschinsky_matrix_internal', cat='dev')

end function Duschinsky_matrix_internal

! ======================================================================

module procedure Duschinsky_shift_arr

    ! Local variables
    real(realwp), dimension(:, :), allocatable :: tmp_r2
    real(realwp) :: tmp
    integer :: i, j, k
    character(len=2) :: local_mode
    integer :: local_n_vib, local_n_at

    ! Sanity checks for alway present arguments
    if (present(n_vib)) then
        if (size(L_mat1, 2) /= n_vib) then
            call runstat%raise_error('L_mat1 must have n_vib columns', &
                source='Duschinsky_shift_arr', cat='dev')
        end if
    end if
    if (present(n_at)) then
        if (size(L_mat1, 1) /= 3 * n_at) then
            call runstat%raise_error('L_mat1 must have 3*n_at rows', &
                source='Duschinsky_shift_arr', cat='dev')
        end if
        if (size(at_mass1) /= n_at) then
            call runstat%raise_error('at_mass1 must have n_at elements', &
                source='Duschinsky_shift_arr', cat='dev')
        end if
    end if

    local_mode = 'AH'
    if (present(mode)) local_mode = upcase(trim(mode))

    local_n_vib = size(L_mat1, 2)
    local_n_at = size(L_mat1, 1) / 3

    select case (local_mode)
        case ('AS', 'AH')
            ! Sanity checks for optional arguments
            if (.not. present(coord1) .or. .not. present(coord2)) then
                call runstat%raise_error('Provide coord1 and coord2 for AS and AH', &
                    source='Duschinsky_shift_arr', cat='dev')
            end if
            if (local_n_at /= size(coord1, 2) .or. local_n_at /= size(coord2, 2)) then
                call runstat%raise_error('coord1 and coord2 must have n_at columns', &
                    source='Duschinsky_shift_arr', cat='dev')
            end if

            Kvec = Duschinsky_shift_adiabatic(local_n_vib, local_n_at, L_mat1, at_mass1, coord1, coord2)
        case ('VG')
            ! Sanity checks for optional arguments
            if (.not. present(red_freq2) .or. .not. present(grad2)) then
                call runstat%raise_error('Provide red_freq2 and grad2 for VG', &
                    source='Duschinsky_shift_arr', cat='dev')
            end if
            if (local_n_vib /= size(red_freq2)) then
                call runstat%raise_error('red_freq2 must have n_vib elements', &
                    source='Duschinsky_shift_arr', cat='dev')
            end if
            Kvec = Duschinsky_shift_vertical_gradient(local_n_vib, local_n_at, L_mat1, at_mass1, red_freq2, grad2)
        case ('VH')
            ! Sanity checks for optional arguments
            if (.not. present(red_freq2) .or. .not. present(grad2) .or. .not. present(Jmat)) then
                call runstat%raise_error('Provide red_freq2, grad2 and Jmat for VH', &
                    source='Duschinsky_shift_arr', cat='dev')
            end if
            if (local_n_vib /= size(red_freq2)) then
                call runstat%raise_error('red_freq2 must have n_vib elements', &
                    source='Duschinsky_shift_arr', cat='dev')
            end if
            if (size(Jmat, 1) /= local_n_vib .or. size(Jmat, 2) /= local_n_vib) then
                call runstat%raise_error('Jmat must be a square matrix of size n_vib', &
                    source='Duschinsky_shift_arr', cat='dev')
            end if
            Kvec = Duschinsky_shift_vertical_hessian(local_n_vib, local_n_at, L_mat1, at_mass1, red_freq2, grad2, Jmat)
        case default
            call runstat%raise_error('Provide appropriate mode for K-vector calculation', &
                source='Duschinsky_shift_arr', cat='dev')
    end select

end procedure Duschinsky_shift_arr

! ======================================================================

module procedure Duschinsky_shift_db

    ! Local variables
    character(len=2) :: local_mode

    ! Sanity check for always present arguments
    if (.not.(vib1%loaded)) then
        call runstat%raise_error('vib1 VibrationsDB is not loaded', &
            source='Duschinsky_shift_db', cat='dev')
    else if (.not.allocated(vib1%L_mat)) then
        call runstat%raise_error('vib1 L_mat is not loaded', &
            source='Duschinsky_shift_db', cat='dev')
    end if
    if (.not.(mol1%loaded)) then
        call runstat%raise_error('mol1 MoleculeDB is not loaded', &
            source='Duschinsky_shift_db', cat='dev')
    else if (.not.allocated(mol1%at_mas)) then
        call runstat%raise_error('mol1 mass is not loaded', &
            source='Duschinsky_shift_db', cat='dev')
    end if

    local_mode = 'AH'
    if (present(mode)) local_mode = upcase(trim(mode))

    select case (local_mode)
        case ('AS', 'AH')
            ! Sanity checks for optional arguments
            if (.not. present(mol2)) then
                call runstat%raise_error('Provide mol2 for AS and AH', &
                    source='Duschinsky_shift_db', cat='dev')
            end if
            if (.not.(mol2%loaded)) then
                call runstat%raise_error('mol2 MoleculeDB is not loaded', &
                    source='Duschinsky_shift_db', cat='dev')
            else if (.not.allocated(mol2%at_crd)) then
                call runstat%raise_error('mol2 coords is not loaded', &
                    source='Duschinsky_shift_db', cat='dev')
            end if
            Kvec = Duschinsky_shift_adiabatic(vib1%n_vib, mol1%n_at, vib1%L_mat, &
                mol1%at_mas, mol1%at_crd, mol2%at_crd)
        case ('VG')
            ! Sanity checks for optional arguments
            if (.not. present(grad2)) then
                call runstat%raise_error('Provide grad2 for VG', &
                    source='Duschinsky_shift_db', cat='dev')
            end if
            Kvec = Duschinsky_shift_vertical_gradient(vib1%n_vib, mol1%n_at, vib1%L_mat, &
                mol1%at_mas, vib1%red_freq, grad2)
        case ('VH')
            ! Sanity checks for optional arguments
            if (.not. present(vib2) .or. .not. present(grad2) .or. .not. present(Jmat)) then
                call runstat%raise_error('Provide vib2, grad2 and Jmat for VH', &
                    source='Duschinsky_shift_db', cat='dev')
            end if
            if (.not.(vib2%loaded)) then
                call runstat%raise_error('vib2 VibrationsDB is not loaded', &
                    source='Duschinsky_shift_db', cat='dev')
            else if (.not.allocated(vib2%red_freq)) then
                call runstat%raise_error('vib2 red_freq is not loaded', &
                    source='Duschinsky_shift_db', cat='dev')
            end if
            Kvec = Duschinsky_shift_vertical_hessian(vib1%n_vib, mol1%n_at, vib1%L_mat, &
                mol1%at_mas, vib2%red_freq, grad2, Jmat)
        case default
            call runstat%raise_error('Provide appropriate mode for K-vector calculation', &
                source='Duschinsky_shift_db', cat='dev')
    end select

end procedure Duschinsky_shift_db

! ======================================================================

module function Duschinsky_shift_adiabatic(n_vib, n_at, L_mat1, at_mass1, coord1, coord2) result(Kvec) 
    !! Computes the Duschinsky vector K for adiabatic states (AS and AH)
    !!
    !! The Duschinsky vector can be calculated as:
    !!  K = L_mat1 * sqrt(at_mass1) * (coord2 - coord1)
    !! 
    integer, intent(in) :: n_vib
    !! Number of vibrational modes
    integer, intent(in) :: n_at
    !! Number of atoms
    real(realwp), dimension(:, :), intent(in) :: L_mat1
    !! L matrix for the 1st state
    real(realwp), dimension(:), intent(in) :: at_mass1
    !! Masses of the atoms
    real(realwp), dimension(:, :), intent(in) :: coord1, coord2
    !! Coordinates of the 1st and 2nd state
    real(realwp), dimension(:), allocatable :: Kvec
    !! Duschinsky vector

    ! Local variables
    real(realwp) :: tmp
    integer :: i, j, k

    allocate(Kvec(n_vib))
    kvec = f0

    ! Constructing for AS and AH
    do i = 1, n_at
        do k = 1, 3
            tmp = sqrt(phys_conv%au2amu(at_mass1(i), .true.)) * (coord2(k, i) - coord1(k, i))
            !$omp parallel do
            do j = 1, n_vib
                Kvec(j) = Kvec(j) + L_mat1(k + (i - 1) * 3, j) * tmp
            end do
            !$omp end parallel do
        end do
    end do

end function Duschinsky_shift_adiabatic

! ======================================================================

module function Duschinsky_shift_vertical_gradient(n_vib, n_at, L_mat1, at_mass1, red_freq2, grad2) result(Kvec)
    !! Computes the Duschinsky vector K for vertical gradient states (VG)
    !!
    !! The Duschinsky vector can be calculated as:
    !!  K = -L_mat1.T * omega2^(-2) * sqrt(at_mass1) * grad2
    !! 
    integer, intent(in) :: n_vib
    !! Number of vibrational modes
    integer, intent(in) :: n_at
    !! Number of atoms
    real(realwp), dimension(:), intent(in) :: at_mass1
    !! Masses of the atoms
    real(realwp), dimension(:, :), intent(in) :: L_mat1
    !! L matrix for the 1st state
    real(realwp), dimension(:), intent(in) :: red_freq2
    !! Frequencies of the 2nd state
    real(realwp), dimension(:, :), intent(in) :: grad2
    !! Gradient of the 2nd state
    real(realwp), dimension(:), allocatable :: Kvec
    !! Duschinsky vector

    ! Local variables
    real(realwp) :: tmp
    integer :: i, j, k

    allocate(Kvec(n_vib))
    Kvec = f0

    ! Constructing for VG
    do i = 1, n_at
        do k = 1, 3
            tmp = grad2(k, i) / sqrt(phys_conv%au2amu(at_mass1(i), .true.))
            !$omp parallel do private(j)
            do j = 1, n_vib
                Kvec(j) = Kvec(j) &
                        - (f1 / (abs(red_freq2(j))*red_freq2(j))) &
                            * L_mat1(k + (i - 1) * 3, j) * tmp
            end do
            !$omp end parallel do
        end do
    end do

end function Duschinsky_shift_vertical_gradient

! ======================================================================

module function Duschinsky_shift_vertical_hessian(n_vib, n_at, L_mat1, at_mass1, &
        red_freq2, grad2, Jmat) result(Kvec)
    !! Computes the Duschinsky vector K for vertical hessian states (VH)
    !!
    !! The Duschinsky vector can be calculated as:
    !!  K = -L_mat1.T * J * red_freq2^(-2) * J.T * sqrt(at_mass1) * grad2
    !!
    integer, intent(in) :: n_vib
    !! Number of vibrational modes
    integer, intent(in) :: n_at
    !! Number of atoms
    real(realwp), dimension(n_at * 3, n_vib), intent(in) :: L_mat1
    !! L matrix for the 1st state
    real(realwp), dimension(n_at), intent(in) :: at_mass1
    !! Masses of the atoms
    real(realwp), dimension(n_vib), intent(in) :: red_freq2
    !! Reduced red_frequencies of the 2nd state
    real(realwp), dimension(:, :), intent(in) :: grad2, Jmat
    !! Gradient of the 2nd state and Duschinsky matrix

    real(realwp), dimension(:), allocatable :: Kvec
    !! Duschinsky vector

    ! Local variables
    real(realwp), dimension(:, :), allocatable :: tmp1_r2, tmp2_r2
    real(realwp) :: tmp
    integer :: i, j, k

    ! Constructing for VH
    allocate(tmp1_r2(n_vib, n_vib))

    !$omp parallel do private(i)
    do i = 1, n_vib
        tmp1_r2(:, i) = Jmat(:, i) &
            / (abs(red_freq2(i))*red_freq2(i))
    end do
    !$omp end parallel do

    allocate(tmp2_r2(n_vib, n_vib))

    call xgemm('N', 'T', n_vib, n_vib, n_vib, f1, tmp1_r2, n_vib, &
        Jmat, n_vib, f0, tmp2_r2, n_vib)

    deallocate(tmp1_r2)
    allocate(tmp1_r2(3 * n_at, n_vib))

    call xgemm('N', 'N', 3 * n_at, n_vib, n_vib, f1, L_mat1, 3 * n_at, &
        tmp2_r2, n_vib, f0, tmp1_r2, 3 * n_at)

    allocate(Kvec(n_vib))
    Kvec = f0

    do i = 1, n_at
        do k = 1, 3
            tmp = grad2(k, i) / sqrt(phys_conv%au2amu(at_mass1(i), .true.))
            !$omp parallel do private(j)
            do j = 1, n_vib
                Kvec(j) = Kvec(j) - tmp1_r2(k + (i - 1) * 3, j) * tmp
            end do
            !$omp end parallel do
        end do
    end do

end function Duschinsky_shift_vertical_hessian

! ======================================================================

end submodule vibronic_Duschinsky
