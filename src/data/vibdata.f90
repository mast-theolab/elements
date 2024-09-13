module vibdata
    use iso_fortran_env, only: int32, int64, real64

    integer :: &
        n_vib = 0  ! number of normal modes
    real(real64), dimension(:), allocatable :: &
        freq, &   ! Harmonic wavenumbers (cm-1)
        red_mass  ! Reduced masses of the vibrations
    real(real64), dimension(:,:), allocatable :: &
        L_mwg, &   ! Eigenvectors of the mass-weighted force constants
        L_mat      ! Eigenvectors of the Hessian matrix, dimensionless

end module vibdata
