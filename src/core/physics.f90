module physics
    !! Physics-related constants and functions
    !!
    !! Provides common phyical constants and subroutines.
    !! Content:
    !! - PhysFact : derived-type with physical conversion factors
    !! - PhysConst : derived-typw with physical constants 
    use iso_fortran_env, only: real32, real64

    implicit none
    
    private

    real(real64), parameter, public :: &
        planck = 6.62606896e-34_real64, &       ! Planck (J.s)
        avogadro = 6.02214179e23_real64, &      ! Avogadro number (mol^-1)
        slight = 2.99792458e10_real64, &        ! Speed of light (cm/s)
        boltzmann = 1.3806504e-23_real64, &     ! Boltzmann (J/K)
        fine_struct = 1._real64/137.035999679_real64, &  ! Fine structure (no unit)
        mol_vol = 22.413996e-3_real64, &        ! Molar volume of ideal gas (m^3@273.15K)
        e_mag_mom = -928.476377e-26_real64, &   ! Electron Magnetic Moment (J/Tesla)
        p_rest_mass = 1.672621637e-27_real64, & ! Proton rest mass (kg)
        g_factor = 2.0023193043622_real64, &    ! Free electron g-factor (no unit)
        bohr_radius = 0.52917720859_real64, &   ! bohr radius in Ang
        u_at_mass = 1.660538782e-27_real64, &   ! u in kg
        e_charge = 1.602176487e-19_real64, &    ! electron charge in coulomb
        calorie = 4.184_real64, &               ! 1 calorie in joule
        E_hartree = 4.35974394e-18_real64       ! 1 hartree in joule

    type, public :: PhysFact
        contains
            procedure, nopass, private :: &
                s_conv_bohr_to_Ang, d_conv_bohr_to_Ang, &
                s_conv_Ang_to_bohr, d_conv_Ang_to_bohr, &
                s_conv_amu_to_kg, d_conv_amu_to_kg, &
                s_conv_e_to_C, d_conv_e_to_C, &
                s_conv_cal_to_J, d_conv_cal_to_J, &
                s_conv_hartree_to_J, d_conv_hartree_to_J
            generic, public :: bohr2Ang => s_conv_bohr_to_Ang, d_conv_bohr_to_Ang
            generic, public :: Ang2bohr => s_conv_Ang_to_bohr, d_conv_Ang_to_bohr
            generic, public :: amu2kg => s_conv_amu_to_kg, d_conv_amu_to_kg
            generic, public :: e2C => s_conv_e_to_C, d_conv_e_to_C
            generic, public :: cal2J => s_conv_cal_to_J, d_conv_cal_to_J
            generic, public :: Eh2J => s_conv_hartree_to_J, d_conv_hartree_to_J
    end type PhysFact

contains

! ======================================================================

elemental real(real32) function s_conv_bohr_to_Ang(x)
    real(real32), intent(in) :: x
    s_conv_bohr_to_Ang = x * real(bohr_radius, kind=real32)
end function s_conv_bohr_to_Ang

! ======================================================================

elemental real(real64) function d_conv_bohr_to_Ang(x)
    real(real64), intent(in) :: x
    d_conv_bohr_to_Ang = x * bohr_radius
end function d_conv_bohr_to_Ang

! ======================================================================

elemental real(real32) function s_conv_Ang_to_bohr(x)
    real(real32), intent(in) :: x
    s_conv_Ang_to_bohr = x / real(bohr_radius, kind=real32)
end function s_conv_Ang_to_bohr

! ======================================================================

elemental real(real64) function d_conv_Ang_to_bohr(x)
    real(real64), intent(in) :: x
    d_conv_Ang_to_bohr = x / bohr_radius
end function d_conv_Ang_to_bohr

! ======================================================================

elemental real(real32) function s_conv_amu_to_kg(x)
    real(real32), intent(in) :: x
    s_conv_amu_to_kg = x * real(u_at_mass, kind=real32)
end function s_conv_amu_to_kg

! ======================================================================

elemental real(real64) function d_conv_amu_to_kg(x)
    real(real64), intent(in) :: x
    d_conv_amu_to_kg = x * u_at_mass
end function d_conv_amu_to_kg

! ======================================================================

elemental real(real32) function s_conv_e_to_C(x)
    real(real32), intent(in) :: x
    s_conv_e_to_C = x * real(e_charge, kind=real32)
end function s_conv_e_to_C

! ======================================================================

elemental real(real64) function d_conv_e_to_C(x)
    real(real64), intent(in) :: x
    d_conv_e_to_C = x * e_charge
end function d_conv_e_to_C

! ======================================================================

elemental real(real32) function s_conv_cal_to_J(x)
    real(real32), intent(in) :: x
    s_conv_cal_to_J = x * real(calorie, kind=real32)
end function s_conv_cal_to_J

! ======================================================================

elemental real(real64) function d_conv_cal_to_J(x)
    real(real64), intent(in) :: x
    d_conv_cal_to_J = x * calorie
end function d_conv_cal_to_J

! ======================================================================

elemental real(real32) function s_conv_hartree_to_J(x)
    real(real32), intent(in) :: x
    s_conv_hartree_to_J = x * real(E_hartree, kind=real32)
end function s_conv_hartree_to_J

! ======================================================================

elemental real(real64) function d_conv_hartree_to_J(x)
    real(real64), intent(in) :: x
    d_conv_hartree_to_J = x * E_hartree
end function d_conv_hartree_to_J

! ======================================================================

end module physics
