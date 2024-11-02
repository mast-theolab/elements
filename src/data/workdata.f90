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
            exc_energy, &     ! excited-state energies
            g2e_energy        ! excitation energies
        logical :: &
            dens_loaded = .false., &  ! Density information has been loaded
            prop_loaded = .false.     ! Properties have been loaded
    end type ExcitationDB

    type PropertyDB
        integer :: id
        !! Internal identification of property.
        integer :: order
        !! Derivative order (`0` if reference value of property).
        character(len=:), allocatable :: label
        !! Short label for the property (for internal use/tests).
        character(len=:), allocatable :: name
        !! Full name of the quantity, for printing.
        character(len=:), allocatable :: unit
        !! Unit of the quantity.
        integer, dimension(:), allocatable :: pdim
        !! Dimension/shape of the property itself.
        integer, dimension(:), allocatable :: shape
        !! Actual shape of the extracted quantity (e.g., for reshaping).
        integer, dimension(:), allocatable :: dim_shape
        !! Internal shape of each dimension:
        !! 0. Unknown/not relevant
        !! 1. Linear storage
        !! 2: lower-triangular part of symmetric 2D matrix.
        !! 3: lower-triangular part of symmetric 3D matrix.
        !! -2: lower-triangular part of antisymm. 2D matrix.
        !! -3: lower-triangular part of antisymm. 3D matrix.
        real(real64), dimension(:), allocatable :: data
        !! Extracted data.
        integer, dimension(2) :: states = [0, 0]
        !! Electronic states involved (same for state-specific)
        !! Last value can be -1 to include all final states.
        logical :: loaded = .false.
        !! Data have been loaded
        logical :: known = .false.
        !! Property is known/supported
        integer :: istat = 0
        !! Store a simple status code
        !! -2. problem to read file
        !! -1. corrupted/unrecognized data file.
        !!  0. no error, all operations proceeded properly.
        !!  1. property not supported by program.
        !!  2. property not found.
        !! 10. unspecified error.
        !! 99. inconsistency in query, e.g., end state < start state
    end type PropertyDB

end module workdata
