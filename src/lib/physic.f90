module physic
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
        finestruct = 1._real64/137.035999679_real64, &  ! Fine structure (no unit)
        molvol = 22.413996e-3_real64, &         ! Molar volume of ideal gas (m^3@273.15K)
        emagmom = -928.476377e-26_real64, &     ! Electron Magnetic Moment (J/Tesla)
        prestmass = 1.672621637e-27_real64, &   ! Proton rest mass (kg)
        gfactor = 2.0023193043622_real64        ! Free electron g-factor (no unit)

    type, public :: PhysFact
        private
        real(real64) :: &
            bohr = 0.52917720859_real64, &       ! bohr radius in Ang
            amu = 1.660538782e-27_real64, &      ! u in kg
            echarge = 1.602176487e-19_real64, &  ! electron charge in coulomb
            cal = 4.184_real64, &                ! 1 calorie in joule
            hartree = 4.35974394e-18_real64      ! 1 hartree in joule
        contains
            private
            procedure :: &
                s_conv_bohr_to_Ang, d_conv_bohr_to_Ang, &
                s_conv_amu_to_kg, d_conv_amu_to_kg, &
                s_conv_e_to_C, d_conv_e_to_C, &
                s_conv_cal_to_J, d_conv_cal_to_J, &
                s_conv_hartree_to_J, d_conv_hartree_to_J
            generic, public :: bohr2Ang => s_conv_bohr_to_Ang, d_conv_bohr_to_Ang
            generic, public :: amu2kg => s_conv_amu_to_kg, d_conv_amu_to_kg
            generic, public :: e2C => s_conv_e_to_C, d_conv_e_to_C
            generic, public :: cal2J => s_conv_cal_to_J, d_conv_cal_to_J
            generic, public :: Eh2J => s_conv_hartree_to_J, d_conv_hartree_to_J
    end type PhysFact

contains

! ======================================================================

elemental real(real32) function s_conv_bohr_to_Ang(this, x)
    class(PhysFact), intent(in) :: this
    real(real32), intent(in) :: x
    s_conv_bohr_to_Ang = x * real(this%bohr, kind=real32)
end function s_conv_bohr_to_Ang

! ======================================================================

elemental real(real64) function d_conv_bohr_to_Ang(this, x)
    class(PhysFact), intent(in) :: this
    real(real64), intent(in) :: x
    d_conv_bohr_to_Ang = x * this%bohr
end function d_conv_bohr_to_Ang

! ======================================================================

elemental real(real32) function s_conv_amu_to_kg(this, x)
    class(PhysFact), intent(in) :: this
    real(real32), intent(in) :: x
    s_conv_amu_to_kg = x * real(this%amu, kind=real32)
end function s_conv_amu_to_kg

! ======================================================================

elemental real(real64) function d_conv_amu_to_kg(this, x)
    class(PhysFact), intent(in) :: this
    real(real64), intent(in) :: x
    d_conv_amu_to_kg = x * this%amu
end function d_conv_amu_to_kg

! ======================================================================

elemental real(real32) function s_conv_e_to_C(this, x)
    class(PhysFact), intent(in) :: this
    real(real32), intent(in) :: x
    s_conv_e_to_C = x * real(this%echarge, kind=real32)
end function s_conv_e_to_C

! ======================================================================

elemental real(real64) function d_conv_e_to_C(this, x)
    class(PhysFact), intent(in) :: this
    real(real64), intent(in) :: x
    d_conv_e_to_C = x * this%echarge
end function d_conv_e_to_C

! ======================================================================

elemental real(real32) function s_conv_cal_to_J(this, x)
    class(PhysFact), intent(in) :: this
    real(real32), intent(in) :: x
    s_conv_cal_to_J = x * real(this%cal, kind=real32)
end function s_conv_cal_to_J

! ======================================================================

elemental real(real64) function d_conv_cal_to_J(this, x)
    class(PhysFact), intent(in) :: this
    real(real64), intent(in) :: x
    d_conv_cal_to_J = x * this%cal
end function d_conv_cal_to_J

! ======================================================================

elemental real(real32) function s_conv_hartree_to_J(this, x)
    class(PhysFact), intent(in) :: this
    real(real32), intent(in) :: x
    s_conv_hartree_to_J = x * real(this%hartree, kind=real32)
end function s_conv_hartree_to_J

! ======================================================================

elemental real(real64) function d_conv_hartree_to_J(this, x)
    class(PhysFact), intent(in) :: this
    real(real64), intent(in) :: x
    d_conv_hartree_to_J = x * this%hartree
end function d_conv_hartree_to_J

! ======================================================================

end module physic
