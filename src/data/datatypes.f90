module datatypes
    use numeric, only: realwp, f0
    use string, only: locase
    use exception, only: BaseException, InitError, RaiseArgError

    implicit none
    private

    type, public :: AtomDB
        character(len=2) :: symbol           ! atomic symbol
        character(len=40) :: name            ! atom name
        integer :: number                    ! atomic number
        real(realwp) :: mass                 ! atomic mass (in amu)
        real(realwp), dimension(3) :: rcov   ! covalent radii
        ! Covalent radii are provided for single, double and triple bonds
        real(realwp), dimension(4) :: rvdw   ! Van der Waals radii
        real(realwp) :: rvis                 ! Visible radius, for display
        real(realwp) :: rcovG                ! Covalent rad. used by Gaussian
        integer, dimension(3) :: color       ! Color as RGB vector.
    contains
        procedure :: get_rcov => get_atom_rcov
        procedure :: get_rvdw => get_atom_rvdw
    end type AtomDB

    type, public :: BasisSetDB
        integer :: &
            n_basis = 0, &   ! number of basis functions
            n_basok = 0, &   ! number of basis functions actually used
            n_shells = 0     ! number of primitive shells
        integer, dimension(:), allocatable :: &
            nprim_per_at     ! number of primitive basis funcs / atom
        logical :: &
            pureD = .true., &  ! pure D basis functions
            pureF = .true.     ! pure F... basis functions
        type(PrimitiveFunction), dimension(:,:), allocatable :: &
            info ! basis set information
        logical :: loaded = .false.
    end type BasisSetDB

    type, public :: ExcitationDB
        integer :: &
            n_states = 0, &   ! molecular charge
            id_state = -1     ! index of reference state (0: ground)
        integer, dimension(:), allocatable :: &
            ispin_exc         ! spin of electronic excited states
        real(realwp) :: &
            gs_energy = f0  ! ground-state energy
        real(realwp), dimension(:,:,:,:), allocatable :: &
            exc_dens, &       ! excited-state densities
            g2e_dens          ! ground-to-excited transition densities
        real(realwp), dimension(:,:), allocatable :: &
            g2e_eldip, &      ! ground-to-excited electric dipoles
            g2e_magdip        ! ground-to-excited magnetic dipoles
        real(realwp), dimension(:), allocatable :: &
            exc_energy, &     ! excited-state energies
            g2e_energy        ! excitation energies
        logical :: &
            dens_loaded = .false., &  ! Density information has been loaded
            prop_loaded = .false.     ! Properties have been loaded
    end type ExcitationDB

    type, public :: MoleculeDB
        integer :: &
            charge = 0, &  ! molecular charge
            multip = 0, &  ! molecular multiplicity
            n_at = 0, &    ! number of atoms
            n_el = 0       ! number of electrons
        integer, dimension(:), allocatable :: &
            at_num         ! atomic numbers
        real(realwp) :: &
            energy = f0  ! total energy (electronic state-dependent)
        real(realwp), dimension(:), allocatable :: &
            at_chg, &      ! nuclear charges
            at_mas         ! atomic masses
        real(realwp), dimension(:,:), allocatable :: &
            at_crd         ! atomic coordinates
        real(realwp), dimension(:,:,:), allocatable :: &
            el_dens        ! electronic density
        character(len=2), dimension(:), allocatable :: &
            at_lab         ! atomic labels
        logical :: &
            loaded = .false., &
            dens_loaded = .false.   ! electronic density has been loaded.
    end type MoleculeDB

    type, public :: OrbitalsDB
        integer :: &
            n_ab = 1, &  ! num. of unique alpha-beta orbitals (1=closed shell)
            n_ao = 0, &  ! number of atomic orbitals
            n_mo = 0     ! Max. number of molecular orbitals, mainly for storage
        integer, dimension(2) :: &  ! n_ab elements are expected to be set
            n_mos = [0, 0], &   ! number of alpha/beta molecular orbitals
            n_els = [0, 0]      ! number of alpha/beta electrons
        real(realwp), dimension(:,:), allocatable :: &  ! last dim depend on n_ab
            en_mos   ! energies of alpha/beta orbitals
        real(realwp), dimension(:,:,:), allocatable :: &  ! last dim depend on n_ab
            coef_mos  ! coefficients of alpha/beta molecular orbitals
        logical :: &
            openshell = .false.  ! if molecule is open shell
        logical :: loaded = .false.
    end type OrbitalsDB

    type, public :: PrimitiveFunction
        !! primitive basis function.
        !! alpha and coeffs are coeffs for hybrid functions (ex: SP)
        character(len=2) :: shelltype
        !! Shell type (of the contracted shell)
        integer ::  L
        !! principal angular momentum
        ! character(len=1) :: functype
        ! !! type of the function: S, P...
        integer :: shellid
        !! main shell id
        integer :: ndim
        !! number of components (dimension)
        logical :: pure
        !! true if spherical harmonics, Cartesian otherwise
        logical :: shell_first
        !! true if first primitive of shell
        logical :: shell_last
        !! true if last primitive of shell
        real(realwp) :: alpha
        !! alpha: primitive exponent
        real(realwp), dimension(:), allocatable :: coeff
        !! Unique normalized contraction coefficients
        !! indexes: 0, -1... : non-normalized original coeffs.
        integer, dimension(:,:), allocatable :: lxyz
        !! Indexes for each unique contraction coefficient
    end type PrimitiveFunction

    type, public :: PropertyDB
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
        real(realwp), dimension(:), allocatable :: data
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

    type, public :: VibrationsDB
        integer :: &
            n_vib = 0  ! number of normal modes
        real(realwp), dimension(:), allocatable :: &
            freq, &   ! Harmonic wavenumbers (cm-1)
            red_mass  ! Reduced masses of the vibrations
        real(realwp), dimension(:,:), allocatable :: &
            L_mwg, &   ! Eigenvectors of the mass-weighted force constants
            L_mat      ! Eigenvectors of the Hessian matrix, dimensionless
        logical :: loaded = .false.
    end type VibrationsDB

    contains

    ! ======================================================================
    
    function get_atom_rcov(this, what, err) result(rcov)
        !! Return the value of the covalent radius
        !!
        !! Returns a single value for the covalent radius based on the
        !! choice by the user.
        !!
        !! Available values are:
        !!
        !! - single: covalent radius for single bond (default)
        !! - double: covalent radius for double bond
        !! - triple: covalent radius for triple bond
        !! - gaussian: covalent radius based on Gaussian values
        !!
        !! Note: if err is not provided and the keyword is incorrect,
        !!       the function simply returns 0.0 to not block operations
        !!       during execution.
        class(AtomDB), intent(in) :: this
        !! Instance of the AtomDB type.
        character(len=*), intent(in), optional :: what
        !! Which type of covalent radius to return
        class(BaseException), intent(out), allocatable, optional :: err
        !! Error instance
        real(realwp) :: rcov
        !! Resulting covalent radius.
    
        character(len=:), allocatable :: key
    
        if (present(err)) err = InitError()
    
        if (present(what)) then
            key = locase(trim(what))
        else
            key = 'single'
        end if
    
        select case(key)
        case('single')
            rcov = this%rcov(1)
        case('double')
            rcov = this%rcov(2)
        case('triple')
            rcov = this%rcov(3)
        case('gaussian', 'gxx')
            rcov = this%rcovG
        case default
            if (present(err)) then
                call RaiseArgError(err, 'Unrecognized category of covalent bond')
            else
                rcov = 0.0_realwp
            end if
        end select
    
    end function get_atom_rcov
    
    ! ======================================================================
    
    function get_atom_rvdw(this, db, err) result(rvdw)
        !! Return a value for the van der Waals radius
        !!
        !! Returns a value for the van der Waals radius among available
        !! data.
        !!
        !! Available values are:
        !!
        !! - bondi64:
        !!   A. Bondi, J. Phys. Chem. A 1964 (68) 441
        !!   https://doi.org/10.1021/j100785a001
        !! - truhlar09:
        !!   M. Mantina, A.C. Chamberlin, R. Valero, C.J. Cramer,
        !!   D.G. Truhlar, J. Phys. Chem. A 2009 (113) 5809.
        !!   https://doi.org/10.1021/jp8111556
        !!   Extension of Bondi's set with some additional atoms from
        !!   main group, but some values from Bondi were ignored.
        !! - alvarez13
        !!   S. Alvarez, Dalt. Trans. 2013 (42) 8617
        !!   https://dx.doi.org/10.1039/c3dt50599e
        !!   Statistical analysis from Cambridge Structure Database
        !! - rahm16
        !!   M. Rahm, R. Hoffmann, N.W. Ashcroft,
        !!   Chem. Eur. J. 2016 (22) 14625.
        !!   https://doi.org/10.1002/chem.201602949
        !!   Radii built by considering a threshold in density of
        !!   0.001 e.bohr^-3, computed at the PBE0/ANO-RCC level
        !! - truhlar_ext
        !!   Extended basis set considering values of bondi64 if not
        !!   provided in original work of Trular and coworkers.
        !!
        !! Note: if err is not provided and the keyword is incorrect,
        !!       the function simply returns 0.0 to not block operations
        !!       during execution.
        class(AtomDB), intent(in) :: this
        !! Instance of the AtomDB type.
        character(len=*), intent(in), optional :: db
        !! Which type of covalent radius to return
        class(BaseException), intent(out), allocatable, optional :: err
        !! Error instance
        real(realwp) :: rvdw
        !! Resulting van der Waals radius.
    
        character(len=:), allocatable :: key
    
        if (present(err)) err = InitError()
    
        if (present(db)) then
            key = locase(trim(db))
        else
            key = 'truhlar_ext'
        end if
    
        select case(key)
        case('bondi', 'bondi64')
            rvdw = this%rvdw(1)
        case('truhlar', 'truhlar09')
            rvdw = this%rvdw(2)
        case('truhlar_ext', 'hybrid')
            if (this%rvdw(2) < epsilon(rvdw)) then
                rvdw = this%rvdw(1)
            else
                rvdw = this%rvdw(2)
            end if
        case('alvarez', 'alvarez13')
            rvdw = this%rvdw(3)
        case('rahm', 'rahm16')
            rvdw = this%rvdw(4)
        case default
            if (present(err)) then
                call RaiseArgError(err, 'Unrecognized source for vdW radii')
            else
                rvdw = 0.0_realwp
            end if
        end select
    
    end function get_atom_rvdw
    
    ! ======================================================================
    
    end module datatypes
