module moldata
    use iso_fortran_env, only: int32, int64, real64
    use basisset, only: PrimitiveFunction

    integer :: &
        charge = 0, &    ! molecular charge
        multip = 0, &    ! molecular multiplicity
        n_ab = 1, &      ! number of unique alpha-beta orbitals (1=close shell)
        n_ao = 0, &      ! number of atomic orbitals
        n_at = 0, &      ! number of atoms
        n_basis = 0, &   ! number of basis functions
        n_el = 0, &      ! number of electrons
        n_mo = 0, &      ! Max. number of molecular orbitals, mainly for storage
        n_shells = 0     ! number of primitive shells
    integer, dimension(2) :: &  ! n_ab elements are expected to be set
        n_mos = [0, 0], &   ! number of alpha/beta molecular orbitals
        n_els = [0, 0]      ! number of alpha/beta electrons
    integer(int32), dimension(:), allocatable :: &
        nprim_per_at ! number of primitive basis funcs / atom
    real(real64), dimension(:,:), allocatable :: &  ! last dim depend on n_ab
        en_mos           ! energies of alpha/beta orbitals
    real(real64), dimension(:), allocatable :: &
        at_chg, &        ! nuclear charges
        at_mas           ! atomic masses
    real(real64), dimension(:,:,:), allocatable :: &  ! last dim depend on n_ab
        coef_mos         ! coefficients of alpha/beta molecular orbitals
    real(real64), dimension(:,:), allocatable :: &
        at_crd           ! atomic coordinates
    logical :: &
        openshell = .True., &  ! if molecule is open shell
        pureD = .True., &      ! pure D basis functions
        pureF = .True.         ! pure F... basis functions
    character(len=2), dimension(:), allocatable :: &
        at_lab           ! atomic labels
    type(PrimitiveFunction), dimension(:,:), allocatable :: &
        bset_info        ! basis set information
    
end module moldata
