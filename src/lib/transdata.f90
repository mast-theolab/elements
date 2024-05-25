module transdata
    use iso_fortran_env, only: int32, int64, real64

    integer :: &
        n_states = 0, &   ! molecular charge
        id_state = -1     ! index of reference state (0: ground)
    integer, dimension(:), allocatable :: &
        ispin_exc         ! spin of electronic excited states
    real(real64), dimension(:,:,:,:), allocatable :: &
        exc_dens, &       ! excited-state densities
        g2e_dens          ! ground-to-excited transition densities
    real(real64), dimension(:,:), allocatable :: &
        g2e_eldip, &      ! ground-to-excited electric dipoles
        g2e_magdip        ! ground-to-excited magnetic dipoles
    real(real64), dimension(:), allocatable :: &
        g2e_energy        ! excitation energies

end module transdata
