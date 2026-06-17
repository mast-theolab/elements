module vibrational_PT2
    !! Module dedicated to vibrational perturbation theory at 2nd order.

    use arrays, only: ij2lin => ij2lin_lt
    use numeric, only: f0, f1, f2, realwp
    use run_env, only: run

    interface calc_en_vib
        !! Calculate the VPT2 vibrational energy
        !!
        !! Considering a state specification as:
        !! non-null modes, number_of_quanta
        !! the routine computes and returns its vibrational energy.
        !!
        !! `e_zpve` corresponds to the zero-point vib. energy.
        !! It must be calculated outside and provided in input,
        !! otherwise the function computes the energy from the
        !! vibrational ground state.
        module procedure calc_en_vib_arr_lt, calc_en_vib_dim_lt, &
            calc_en_vib_arr_sq, calc_en_vib_dim_sq

    end interface calc_en_vib

contains

! ======================================================================

function calc_en_vib_arr_lt(nq_index, nq_value, freq, anh_X_mat, e_zpve) &
        result(energy)
    !! This version expects anh_X_mat to be in linear lower-triangular
    !! form, the dimensions are recovered from the array sizes.
    integer, dimension(:), intent(in) :: nq_index
        !! Index of the excited modes.  It should not contain null values.
    integer, dimension(:), intent(in) :: nq_value
        !! Number of quanta associated to each excited mode.
    real(realwp), dimension(:), intent(in) :: freq
        !! Wavenumbers, in cm^-1^.
    real(realwp), dimension(:), intent(in) :: anh_X_mat
        !! Anharmonic X matrix stored in linear lower-tri. form, in cm^-1^.
    real(realwp), intent(in), optional :: e_zpve
        !! Energy of the zero-point vibrational energy, in cm^-1^.
    real(realwp) :: energy
        !! Resulting vibrational energy.

    integer :: lnq, nvib

    if (size(nq_index) /= size(nq_value)) then
        call run%error%raise_argerror('size', &
            'inconsistency in definition of the state', &
            source='calc_en_vib_arr_lt')
        return
    end if

    lnq = size(nq_index)

    nvib = size(freq)
    if (size(anh_X_mat) /= nvib*(nvib+1)/2) then
        call run%error%raise_argerror('size', &
            'inconsistency between sizes of freq and anh. X matrix', &
            source='calc_en_vib_arr_lt')
        return
    end if

    energy = calc_en_vib_dim_lt(lnq, nvib, nq_index, nq_value, freq, &
                                anh_X_mat, e_zpve)

end function calc_en_vib_arr_lt

! ======================================================================

function calc_en_vib_dim_lt(n_vib, n_modes, nq_index, nq_value, freq, &
                            anh_X_mat, e_zpve) result(energy)
    !! This version expects anh_X_mat to be in linear lower-triangular
    !! form, correct dimensions should be provided.
    integer, intent(in) :: n_vib
        !! Total number of vibrational modes.
    integer, intent(in) :: n_modes
        !! Number of excited modes.
    integer, dimension(n_modes), intent(in) :: nq_index
        !! Index of the excited modes.  It should not contain null values.
    integer, dimension(n_modes), intent(in) :: nq_value
        !! Number of quanta associated to each excited mode.
    real(realwp), dimension(n_vib), intent(in) :: freq
        !! Wavenumbers, in cm^-1^.
    real(realwp), dimension(:), intent(in) :: anh_X_mat
        !! Anharmonic X matrix stored in linear lower-tri. form, in cm^-1^.
    real(realwp), intent(in), optional :: e_zpve
        !! Energy of the zero-point vibrational energy, in cm^-1^.
    real(realwp) :: energy
        !! Resulting vibrational energy.

    integer :: i, i0, j, k, l, ni, nj, nk
    real(realwp) :: x

    if (any(nq_index == 0) .or. any(nq_value == 0)) then
        call run%error%raise_argerror('value', &
           'state specification contains null data', &
           source='calc_en_vib_dim_lt')
        return
    end if

    if (present(e_zpve)) then
        energy = e_zpve
    else
        energy = f0
    end if

    select case(n_modes)
    case(1)  ! 1 excited mode
        ni = nq_value(1)
        i = nq_index(1)
        energy = energy + ni*freq(i) + (ni**2+ni) * anh_X_mat(ij2lin(i,i))
        do j = 1, i-1
            energy = energy + ni*anh_X_mat(ij2lin(i,j))/f2
        end do
        do j = i+1, n_vib
            energy = energy + ni*anh_X_mat(ij2lin(i,j))/f2
        end do

    case(2)  ! 2 excited modes
        ni = nq_value(1)
        nj = nq_value(2)
        i = nq_index(1)
        j = nq_index(2)
        energy = energy + ni*freq(i) + nj*freq(j) &
            + (ni**2 + ni) * anh_X_mat(ij2lin(i,i)) &
            + (nj**2 + nj) * anh_X_mat(ij2lin(j,j)) &
            + (ni*nj + (ni+nj)/f2) * anh_X_mat(ij2lin(i,j))
        do k = 1, n_vib
            if (k /= i .and. k /= j) &
                energy = energy + ni*anh_X_mat(ij2lin(i,k))/f2 &
                    + nj*anh_X_mat(ij2lin(j,k))/f2
        end do

    case(3)  ! 2 excited modes
        ni = nq_value(1)
        nj = nq_value(2)
        nk = nq_value(3)
        i = nq_index(1)
        j = nq_index(2)
        k = nq_index(3)
        energy = energy + ni*freq(i) + nj*freq(j) + nk*freq(k) &
            + (ni**2 + ni) * anh_X_mat(ij2lin(i,i)) &
            + (nj**2 + nj) * anh_X_mat(ij2lin(j,j)) &
            + (nk**2 + nk) * anh_X_mat(ij2lin(k,k)) &
            + (ni*nj + (ni+nj)/f2) * anh_X_mat(ij2lin(i,j)) &
            + (ni*nk + (ni+nk)/f2) * anh_X_mat(ij2lin(i,k)) &
            + (nj*nk + (nj+nk)/f2) * anh_X_mat(ij2lin(j,k))
        do l = 1, n_vib
            if (l /= i .and. l /= j .and. l /= k) &
                energy = energy + ni*anh_X_mat(ij2lin(i,l))/f2 &
                    + nj*anh_X_mat(ij2lin(j,l))/f2 &
                    + nk*anh_X_mat(ij2lin(k,l))/f2
        end do

    case default
        do i0 = 1, n_modes
            ni = nq_value(i0)
            i = nq_index(i0)
            energy = energy + ni*freq(i)
            do j = 1, i0
                energy = energy &
                    + ni*nq_value(j)*anh_X_mat(ij2lin(i,nq_index(j)))
            end do
            do j = 1, n_vib
                if (j == i) then
                    x = f2
                else
                    x = f1
                end if
                energy = energy &
                    + x*ni*anh_X_mat(ij2lin(i,nq_index(j)))/f2
            end do
        end do

    end select

end function calc_en_vib_dim_lt

! ======================================================================

function calc_en_vib_arr_sq(nq_index, nq_value, freq, anh_X_mat, e_zpve) &
        result(energy)
    !! This version expects anh_X_mat to be in square form, fully filled,
    !! the dimensions are recovered from the array sizes.
    integer, dimension(:), intent(in) :: nq_index
        !! Index of the excited modes.  It should not contain null values.
    integer, dimension(:), intent(in) :: nq_value
        !! Number of quanta associated to each excited mode.
    real(realwp), dimension(:), intent(in) :: freq
        !! Wavenumbers, in cm^-1^.
    real(realwp), dimension(:,:), intent(in) :: anh_X_mat
        !! Anharmonic X matrix, in cm^-1^.
    real(realwp), intent(in), optional :: e_zpve
        !! Energy of the zero-point vibrational energy, in cm^-1^.
    real(realwp) :: energy
        !! Resulting vibrational energy.

    integer :: lnq, nvib

    if (size(nq_index) /= size(nq_value)) then
        call run%error%raise_argerror('size', &
            'inconsistency in definition of the state', &
             source='calc_en_vib_arr_sq')
        return
    end if

    lnq = size(nq_index)

    nvib = size(freq)
    if (size(anh_X_mat) /= nvib**2) then
        call run%error%raise_argerror('size', &
            'inconsistency between sizes of freq and anh. X matrix', &
            source='calc_en_vib_arr_sq')
        return
    end if

    energy = calc_en_vib_dim_sq(lnq, nvib, nq_index, nq_value, freq, &
                                anh_X_mat, e_zpve)

end function calc_en_vib_arr_sq

! ======================================================================

function calc_en_vib_dim_sq(n_vib, n_modes, nq_index, nq_value, freq, &
                            anh_X_mat, e_zpve) result(energy)
    !! This version expects anh_X_mat to be in square form, fully filled,
    !! correct dimensions should be provided.
    integer, intent(in) :: n_vib
        !! Total number of vibrational modes.
    integer, intent(in) :: n_modes
        !! Number of excited modes.
    integer, dimension(n_modes), intent(in) :: nq_index
        !! Index of the excited modes.  It should not contain null values.
    integer, dimension(n_modes), intent(in) :: nq_value
        !! Number of quanta associated to each excited mode.
    real(realwp), dimension(n_vib), intent(in) :: freq
        !! Wavenumbers, in cm^-1^.
    real(realwp), dimension(n_vib,n_vib), intent(in) :: anh_X_mat
        !! Anharmonic X matrix, in cm^-1^.
    real(realwp), intent(in), optional :: e_zpve
        !! Energy of the zero-point vibrational energy, in cm^-1^.
    real(realwp) :: energy
        !! Resulting vibrational energy.

    integer :: i, i0, j, k, l, ni, nj, nk
    real(realwp) :: x

    if (any(nq_index == 0) .or. any(nq_value == 0)) then
        call run%error%raise_argerror('value', &
            'state specification contains null data', &
            source='calc_en_vib_dim_lt')
        return
    end if

    if (present(e_zpve)) then
        energy = e_zpve
    else
        energy = f0
    end if

    select case(n_modes)
    case(1)  ! 1 excited mode
        ni = nq_value(1)
        i = nq_index(1)
        energy = energy + ni*freq(i) + (ni**2+ni) * anh_X_mat(i,i)
        do j = 1, i-1
            energy = energy + ni*anh_X_mat(i,j)/f2
        end do
        do j = i+1, n_vib
            energy = energy + ni*anh_X_mat(i,j)/f2
        end do

    case(2)  ! 2 excited modes
        ni = nq_value(1)
        nj = nq_value(2)
        i = nq_index(1)
        j = nq_index(2)
        energy = energy + ni*freq(i) + nj*freq(j) &
            + (ni**2 + ni) * anh_X_mat(i,i) &
            + (nj**2 + nj) * anh_X_mat(j,j) &
            + (ni*nj + (ni+nj)/f2) * anh_X_mat(i,j)
        do k = 1, n_vib
            if (k /= i .and. k /= j) &
                energy = energy + ni*anh_X_mat(i,k)/f2 &
                    + nj*anh_X_mat(j,k)/f2
        end do

    case(3)  ! 2 excited modes
        ni = nq_value(1)
        nj = nq_value(2)
        nk = nq_value(3)
        i = nq_index(1)
        j = nq_index(2)
        k = nq_index(3)
        energy = energy + ni*freq(i) + nj*freq(j) + nk*freq(k) &
            + (ni**2 + ni) * anh_X_mat(i,i) &
            + (nj**2 + nj) * anh_X_mat(j,j) &
            + (nk**2 + nk) * anh_X_mat(k,k) &
            + (ni*nj + (ni+nj)/f2) * anh_X_mat(i,j) &
            + (ni*nk + (ni+nk)/f2) * anh_X_mat(i,k) &
            + (nj*nk + (nj+nk)/f2) * anh_X_mat(j,k)
        do l = 1, n_vib
            if (l /= i .and. l /= j .and. l /= k) &
                energy = energy + ni*anh_X_mat(i,l)/f2 &
                    + nj*anh_X_mat(j,l)/f2 + nk*anh_X_mat(k,l)/f2
        end do

    case default
        do i0 = 1, n_modes
            ni = nq_value(i0)
            i = nq_index(i0)
            energy = energy + ni*freq(i)
            do j = 1, i0
                energy = energy + ni*nq_value(j)*anh_X_mat(i,nq_index(j))
            end do
            do j = 1, n_vib
                if (j == i) then
                    x = f2
                else
                    x = f1
                end if
                energy = energy + x*ni*anh_X_mat(i,nq_index(j))/f2
            end do
        end do

    end select

end function calc_en_vib_dim_sq

! ======================================================================

end module vibrational_PT2