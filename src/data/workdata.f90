module workdata
    use iso_fortran_env, only: int32, int64, real64
    use basisset, only: PrimitiveFunction

    type, public :: MoleculeDB
        integer :: &
            charge = 0, &  ! molecular charge
            multip = 0, &  ! molecular multiplicity
            n_at = 0, &    ! number of atoms
            n_el = 0       ! number of electrons
        integer, dimension(:), allocatable :: &
            at_num         ! atomic numbers
        real(real64) :: &
            energy = 0.0_real64  ! total energy (electronic state-dependent)
        real(real64), dimension(:), allocatable :: &
            at_chg, &      ! nuclear charges
            at_mas         ! atomic masses
        real(real64), dimension(:,:), allocatable :: &
            at_crd         ! atomic coordinates
        character(len=2), dimension(:), allocatable :: &
            at_lab         ! atomic labels
        logical :: loaded = .false.
    end type MoleculeDB

    type, public :: BasisSetDB
        integer :: &
            n_basis = 0, &   ! number of basis functions
            n_basok = 0, &   ! number of basis functions actually used
            n_shells = 0     ! number of primitive shells
        integer(int32), dimension(:), allocatable :: &
            nprim_per_at     ! number of primitive basis funcs / atom
        logical :: &
            pureD = .true., &  ! pure D basis functions
            pureF = .true.     ! pure F... basis functions
        type(PrimitiveFunction), dimension(:,:), allocatable :: &
            info ! basis set information
        logical :: loaded = .false.
    end type BasisSetDB

    type, public :: OrbitalsDB
        integer :: &
            n_ab = 1, &  ! num. of unique alpha-beta orbitals (1=closed shell)
            n_ao = 0, &  ! number of atomic orbitals
            n_mo = 0     ! Max. number of molecular orbitals, mainly for storage
        integer, dimension(2) :: &  ! n_ab elements are expected to be set
            n_mos = [0, 0], &   ! number of alpha/beta molecular orbitals
            n_els = [0, 0]      ! number of alpha/beta electrons
        real(real64), dimension(:,:), allocatable :: &  ! last dim depend on n_ab
            en_mos   ! energies of alpha/beta orbitals
        real(real64), dimension(:,:,:), allocatable :: &  ! last dim depend on n_ab
            coef_mos  ! coefficients of alpha/beta molecular orbitals
        logical :: &
            openshell = .false.  ! if molecule is open shell
        logical :: loaded = .false.
    end type OrbitalsDB

    type, public :: VibrationsDB
        integer :: &
            n_vib = 0  ! number of normal modes
        real(real64), dimension(:), allocatable :: &
            freq, &   ! Harmonic wavenumbers (cm-1)
            red_mass  ! Reduced masses of the vibrations
        real(real64), dimension(:,:), allocatable :: &
            L_mwg, &   ! Eigenvectors of the mass-weighted force constants
            L_mat      ! Eigenvectors of the Hessian matrix, dimensionless
        logical :: loaded = .false.
    end type VibrationsDB

    type, public :: ExcitationDB
        integer :: &
            n_states = 0, &   ! molecular charge
            id_state = -1     ! index of reference state (0: ground)
        integer, dimension(:), allocatable :: &
            ispin_exc         ! spin of electronic excited states
        real(real64) :: &
            gs_energy = 0.0_real64  ! ground-state energy
        real(real64), dimension(:,:,:,:), allocatable :: &
            exc_dens, &       ! excited-state densities
            g2e_dens          ! ground-to-excited transition densities
        real(real64), dimension(:,:), allocatable :: &
            g2e_eldip, &      ! ground-to-excited electric dipoles
            g2e_magdip        ! ground-to-excited magnetic dipoles
        real(real64), dimension(:), allocatable :: &
            exc_energy        ! excited-state energies
            g2e_energy        ! excitation energies
        logical :: &
            dens_loaded = .false., &  ! Density information has been loaded
            prop_loaded = .false.     ! Properties have been loaded
    end type ExcitationDB

end module workdata
