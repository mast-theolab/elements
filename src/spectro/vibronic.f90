module vibronic

    use numeric, only: f0, f1, realwp
    use datatypes, only: MoleculeDB, VibrationsDB
    use physics, only: phys_conv
    use run_env, only: run

    implicit none

    private

    public :: Duschinsky_matrix, Duschinsky_shift, extrapolate_geom

    interface Duschinsky_matrix
        !! Computes the Duschinsky matrix J
        !!
        !! The Duschinsky matrix can be calculated as:
        !!
        !!  * VH/AH: J = L_mat1.T * L_mat2
        !!  * VG/AS: J = Identity
        !!
        !!  The type of calculation is determined by the presence of the arguments

        module function Duschinsky_matrix_identity(n_vib, is_identity) result(Jmat)
            !! The Duschinsky matrix is set to identity
            integer, intent(in) :: n_vib
                !! Number of vibrational modes
            logical, optional, intent(in) :: is_identity
                !! If true, the Duschinsky matrix is set to identity
            real(realwp), dimension(:, :), allocatable :: Jmat
                !! Duschinsky matrix
        end function Duschinsky_matrix_identity

        module function Duschinsky_matrix_arr(L_mat1, L_mat2, at_mass1, at_mass2, &
                                       Ginv_mat, is_identity) result(Jmat)
            !! This version takes arrays as arguments.
            real(realwp), dimension(:, :), intent(in) :: L_mat1, L_mat2
                !! L matrix for the 1st and 2nd states
            real(realwp), dimension(:), optional, intent(in) :: at_mass1, at_mass2
                !! Masses of the atoms
            logical, optional, intent(in) :: is_identity
                !! If true, the Duschinsky matrix is set to identity
            real(realwp), dimension(:, :), optional, intent(in) :: Ginv_mat
                !! Inverse of the Wilson B matrix
            real(realwp), dimension(:, :), allocatable :: Jmat
                !! Duschinsky matrix
        end function Duschinsky_matrix_arr

        module function Duschinsky_matrix_db(vib1, vib2, mol1, mol2, Ginv_mat, is_identity) result(Jmat)
            !! This version takes VibrationsDB and MoleculeDB objects as arguments.
            use datatypes, only: MoleculeDB, VibrationsDB

            class(VibrationsDB), intent(in) :: vib1, vib2
            !! Number of vibrational modes
            class(MoleculeDB), optional, intent(in) :: mol1, mol2
            !! Molecule database.
            real(realwp), optional, dimension(:, :), intent(in) :: Ginv_mat
            !! Inverse of the Wilson B matrix
            logical, optional, intent(in) :: is_identity
            !! If true, the Duschinsky matrix is set to identity
            real(realwp), dimension(:, :), allocatable :: Jmat
            !! Duschinsky matrix
        end function Duschinsky_matrix_db

    end interface Duschinsky_matrix

! ----------------------------------------------------------------------

    interface Duschinsky_shift
        !! Computes the Duschinsky vector K
        !!
        !! The Duschinsky vector can be calculated as:
        !!
        !! * AS/AH: K = L_mat1 * sqrt(at_mass1) * (R2 - R1)
        !! * VG: K = -L_mat1.T * omega2^(-2) * sqrt(at_mass1) * gradient
        !! * VH: K = -L_mat1.T * J * omega2^(-2) * J.T * sqrt(at_mass1) * grad
        !!
        !!  The type of calculation is determined by the mode argument

        module function Duschinsky_shift_arr(L_mat1, at_mass1, mode, n_vib, n_at, coord1, coord2, &
                              red_freq2, grad2, Jmat) result(Kvec)
            !! This version takes arrays as arguments.
            real(realwp), dimension(:, :), intent(in) :: L_mat1
                !! L matrix for the 1st state
            real(realwp), dimension(:), intent(in) :: at_mass1
                !! Masses of the atoms
            character(len=*), optional, intent(in) :: mode
                !! Type of calculation (AS, AH, VG, VH)
            integer, optional, intent(in) :: n_vib
                !! Number of vibrational modes
            integer, optional, intent(in) :: n_at
                !! Number of atoms
            real(realwp), dimension(:, :), optional, intent(in) :: coord1
                !! Coordinates of the 1st state
            real(realwp), dimension(:, :), optional, intent(in) :: coord2
                !! Coordinates of the 2nd state
            real(realwp), dimension(:), optional, intent(in) :: red_freq2
                !! Frequencies of the 2nd state
            real(realwp), dimension(:, :), optional, intent(in) :: grad2, Jmat
                !! Gradient of the 2nd state and Duschinsky matrix
            real(realwp), dimension(:), allocatable :: Kvec
                !! Duschinsky vector
        end function Duschinsky_shift_arr

        module function Duschinsky_shift_db(vib1, mol1, mol2, mode, vib2, grad2, Jmat) result(Kvec)
            !! This version takes a DB objects as argument.
            class(VibrationsDB), intent(in) :: vib1
                !! Vibrations database for the 1st state
            class(MoleculeDB), intent(in) :: mol1
                !! Molecule database for the 1st state
            character(len=*), optional, intent(in) :: mode
                !! Type of calculation (AS, AH, VG, VH)
            class(MoleculeDB), optional, intent(in) :: mol2
                !! Molecule database for the 2nd state
            class(VibrationsDB), optional, intent(in) :: vib2
                !! Vibrations database for the 2nd state
            real(realwp), dimension(:, :), optional, intent(in) :: grad2
                !! Gradient of the 2nd state
            real(realwp), dimension(:, :), optional, intent(in) :: Jmat
                !! Duschinsky matrix
            real(realwp), dimension(:), allocatable :: Kvec
                !! Duschinsky vector
        end function Duschinsky_shift_db

    end interface Duschinsky_shift

! ----------------------------------------------------------------------

    interface extrapolate_geom
        !! Extrapolate geometry from shift vector.
        !!
        !! Considering a reference, vertical geometry, the routine
        !! extrapolates a new equilibrium geometry based on the
        !! Duschinsky vector.
        module procedure extrapolate_geom_arr, extrapolate_geom_dim, &
            extrapolate_geom_db
    end interface extrapolate_geom

contains

! ======================================================================

function extrapolate_geom_arr(coord_ref, at_mass, L_mat, dusch_vec) &
        result(coord)
    !! Extrapolate the geometry using data arrays as arguments,
    !! recovering dimensions from their shape.
    !! Information that can be easily recovered, like the number of
    !! atoms from the array sizes, is recomputed.
    real(realwp), dimension(:,:), intent(in) :: coord_ref
        !! Coordinates of the reference point, in au.
    real(realwp), dimension(:), intent(in) :: at_mass
        !! Atomic masses, in u.
    real(realwp), dimension(:,:), intent(in) :: L_mat
        !! Dimensionless Hessian eigenvectors matrix.
    real(realwp), dimension(:), intent(in) :: dusch_vec
        !! Duschinsky shift vector, in au.
    real(realwp), dimension(:,:), allocatable :: coord
        !! Extrapolated coordinates, in au.

    integer :: n_atoms, n_vib

    n_atoms = size(at_mass)
    n_vib = size(dusch_vec)
    if (size(coord_ref, 2) /= n_atoms) then
        call run%error%raise_argerror('size', &
            'unable to build extrapolated geometry', &
            details='inconsistent size between coordinates and masses', &
            source='extrapolate_geom_arr')
        return
    end if

    if (size(L_mat, 2) /= n_vib .and. size(L_mat, 1) /= 3*n_atoms) then
        call run%error%raise_argerror('size', &
            'unable to build extrapolated geometry', &
            details='inconsistent size between Lmat and coords/shift vec', &
            source='extrapolate_geom_arr')
        return
    end if

    coord = extrapolate_geom_dim(n_atoms, n_vib, coord_ref, at_mass, L_mat, &
                                 dusch_vec)

end function extrapolate_geom_arr

! ======================================================================

function extrapolate_geom_db(molDB, vibDB, dusch_vec) result(coord)
    !! Extrapolate the geometry using data data objects as arguments.
    class(MoleculeDB), intent(in) :: molDB
        !! Molecule database.
    type(VibrationsDB), intent(in) :: vibDB
        !! Vibration database.
    real(realwp), dimension(:), intent(in) :: dusch_vec
        !! Duschinsky shift vector, in au.
    real(realwp), dimension(:,:), allocatable :: coord
        !! Extrapolated coordinates, in au.

    if (.not.molDB%loaded) then
        call run%error%raise_argerror('missing', &
            'unable to build extrapolated geometry', &
            details='molecular database is not loaded', &
            source='extrapolate_geom_db')
        return
    end if

    if (.not.vibDB%loaded .or. .not.allocated(vibDB%L_mat)) then
        call run%error%raise_argerror('missing', &
            'unable to build extrapolated geometry', &
            details='missing Hessian eigenvectors matrix', &
            source='extrapolate_geom_db')
        return
    end if

    if (size(vibDB%L_mat, 2) /= size(dusch_vec)) then
        call run%error%raise_argerror('size', &
            'unable to build extrapolated geometry', &
            details='inconsistent size between Lmat and shift vec', &
            source='extrapolate_geom_db')
        return
    end if

    coord = extrapolate_geom_dim(molDB%n_at, vibDB%n_vib, molDB%at_crd, &
                                 molDB%at_mas, vibDB%L_mat, dusch_vec)

end function extrapolate_geom_db

! ======================================================================

function extrapolate_geom_dim(n_atoms, n_vib, coord_ref, at_mass, L_mat, &
                              dusch_vec) result(coord)
    !! Extrapolate the geometry using data arrays as arguments.
    !! No check on arrays size is done besides the explicit
    !! specification.
    integer, intent(in) :: n_atoms
        !! Number of atoms.
    integer, intent(in) :: n_vib
        !! Number of vibrations.
    real(realwp), dimension(3,n_atoms), intent(in) :: coord_ref
        !! Coordinates of the reference point, in au.
    real(realwp), dimension(n_atoms), intent(in) :: at_mass
        !! Atomic masses, in u.
    real(realwp), dimension(3*n_atoms,n_vib), intent(in) :: L_mat
        !! Dimensionless Hessian eigenvectors matrix.
    real(realwp), dimension(n_vib), intent(in) :: dusch_vec
        !! Duschinsky shift vector, in au.
    real(realwp), dimension(:,:), allocatable :: coord
        !! Extrapolated coordinates, in au.

    integer :: ia, ix
    real(realwp) :: weight

    allocate(coord(3,n_atoms))

    do ia = 1, n_atoms
        weight = f1/sqrt(phys_conv%au2amu(at_mass(ia), reverse=.true.))
        do ix = 1, 3
            coord(ix,ia) = coord_ref(ix,ia) &
                + sum(L_mat((ia-1)+ix,:)*dusch_vec(:))*weight
        end do
    end do

end function extrapolate_geom_dim

! ======================================================================


end module vibronic
