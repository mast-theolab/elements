module vibronic

    use numeric, only: realwp
    use datatypes, only: MoleculeDB, VibrationsDB

    implicit none

    private

    public :: Duschinsky_matrix, Duschinsky_shift

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

contains

end module vibronic
