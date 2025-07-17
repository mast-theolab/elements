module orbital
    !! Orbitals (atomic/molecular)-related module.
    !!
    !! The module defines procedures related to atomic and molecular
    !! orbitals.
    
    use numeric, only: realwp, f0

    implicit none

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

end module orbital