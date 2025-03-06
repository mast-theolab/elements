module physics
    !! Physics-related constants and functions
    !!
    !! Provides common phyical constants and subroutines.
    !! Content:
    !!
    !! - PhysFact : derived-type with physical conversion factors
    !! - phys_conv: Instance of PhysFact
    !!
    !! Most functions in phys_conv can operate in reverse to convert
    !! backward.
    use iso_fortran_env, only: real32, real64

    implicit none

    private

    real(real64), parameter, public :: &
        planck = 6.62606896e-34_real64, &
        !! Planck constant (J.s)
        avogadro = 6.02214179e23_real64, &
        !! Avogadro number (mol^-1)
        slight = 2.99792458e10_real64, &
        !! Speed of light (cm/s)
        boltzmann = 1.3806504e-23_real64, &
        !! Boltzmann (J/K)
        fine_struct = 1._real64/137.035999679_real64, &
        !! Fine structure (no unit)
        mol_vol = 22.413996e-3_real64, &
        !! Molar volume of ideal gas (m^3@273.15K)
        e_mag_mom = -928.476377e-26_real64, &
        !! Electron Magnetic Moment (J/Tesla)
        p_rest_mass = 1.672621637e-27_real64, &
        !! Proton rest mass (kg)
        g_factor = 2.0023193043622_real64, &
        !! Free electron g-factor (no unit)
        bohr_radius = 0.52917720859_real64, &
        ! bohr radius in Ang
        u_at_mass = 1.660538782e-27_real64, &
        !! u in kg
        e_charge = 1.602176487e-19_real64, &
        !! electron charge in coulomb
        calorie = 4.184_real64, &
        !! 1 calorie in joule
        E_hartree = 4.35974394e-18_real64, &
        !! 1 hartree in joule
        e_mass = E_hartree*1.0e4_real64/(slight*fine_struct)**2
        !! electron mass (atomic units of mass)

    real(real64), parameter, private :: &
        pi = 4.0_real64*atan(1.0_real64), &
        hbar_amu_Ang_s1 = planck*1.0e20_real64/(2.0_real64*pi*u_at_mass), &
        !! Reduced Planck constant in amu.Ang^2.s^-1
        hc = planck*slight, &
        !! h*c, in J.cm
        factG = 2.0_real64*pi*slight/hbar_amu_Ang_s1
        !! 2 pi c / hbar in cm.amu^-1.Ang^-2

    type, public :: PhysFact
        contains
            procedure, nopass, private :: &
                s_conv_bohr_to_Ang, d_conv_bohr_to_Ang, &
                s_conv_amu_to_kg, d_conv_amu_to_kg, &
                s_conv_e_to_C, d_conv_e_to_C, &
                s_conv_cal_to_J, d_conv_cal_to_J, &
                s_conv_hartree_to_J, d_conv_hartree_to_J, &
                s_conv_au_freq_to_cm1, d_conv_au_freq_to_cm1, &
                s_conv_dE_au_to_cgs, d_conv_dE_au_to_cgs
            generic, public :: bohr2Ang => s_conv_bohr_to_Ang, d_conv_bohr_to_Ang
            generic, public :: amu2kg => s_conv_amu_to_kg, d_conv_amu_to_kg
            generic, public :: e2C => s_conv_e_to_C, d_conv_e_to_C
            generic, public :: cal2J => s_conv_cal_to_J, d_conv_cal_to_J
            generic, public :: Eh2J => s_conv_hartree_to_J, d_conv_hartree_to_J
            generic, public :: au2cm1 => s_conv_au_freq_to_cm1, &
                d_conv_au_freq_to_cm1
            generic, public :: dE_au2cm => s_conv_dE_au_to_cgs, &
            d_conv_dE_au_to_cgs
    end type PhysFact

    type(PhysFact), public :: phys_conv

contains

! ======================================================================

elemental real(real32) function s_conv_bohr_to_Ang(x, reverse)
    real(real32), intent(in) :: x
    logical, intent(in), optional :: reverse
    logical :: do_reverse
    if (present(reverse)) then
        do_reverse = reverse
    else
        do_reverse = .false.
    end if
    if (do_reverse) then
        s_conv_bohr_to_Ang = x / real(bohr_radius, kind=real32)
    else
        s_conv_bohr_to_Ang = x * real(bohr_radius, kind=real32)
    end if
end function s_conv_bohr_to_Ang

! ======================================================================

elemental real(real64) function d_conv_bohr_to_Ang(x, reverse)
    real(real64), intent(in) :: x
    logical, intent(in), optional :: reverse
    logical :: do_reverse
    if (present(reverse)) then
        do_reverse = reverse
    else
        do_reverse = .false.
    end if
    if (do_reverse) then
        d_conv_bohr_to_Ang = x / bohr_radius
    else
        d_conv_bohr_to_Ang = x * bohr_radius
    end if
end function d_conv_bohr_to_Ang

! ======================================================================

elemental real(real32) function s_conv_amu_to_kg(x, reverse)
    real(real32), intent(in) :: x
    logical, intent(in), optional :: reverse
    logical :: do_reverse
    if (present(reverse)) then
        do_reverse = reverse
    else
        do_reverse = .false.
    end if
    if (do_reverse) then
        s_conv_amu_to_kg = x / real(u_at_mass, kind=real32)
    else
        s_conv_amu_to_kg = x * real(u_at_mass, kind=real32)
    end if
end function s_conv_amu_to_kg

! ======================================================================

elemental real(real64) function d_conv_amu_to_kg(x, reverse)
    real(real64), intent(in) :: x
    logical, intent(in), optional :: reverse
    logical :: do_reverse
    if (present(reverse)) then
        do_reverse = reverse
    else
        do_reverse = .false.
    end if
    if (do_reverse) then
        d_conv_amu_to_kg = x / u_at_mass
    else
        d_conv_amu_to_kg = x * u_at_mass
    end if
end function d_conv_amu_to_kg

! ======================================================================

elemental real(real32) function s_conv_e_to_C(x, reverse)
    real(real32), intent(in) :: x
    logical, intent(in), optional :: reverse
    logical :: do_reverse
    if (present(reverse)) then
        do_reverse = reverse
    else
        do_reverse = .false.
    end if
    if (do_reverse) then
        s_conv_e_to_C = x / real(e_charge, kind=real32)
    else
        s_conv_e_to_C = x * real(e_charge, kind=real32)
    end if
end function s_conv_e_to_C

! ======================================================================

elemental real(real64) function d_conv_e_to_C(x, reverse)
    real(real64), intent(in) :: x
    logical, intent(in), optional :: reverse
    logical :: do_reverse
    if (present(reverse)) then
        do_reverse = reverse
    else
        do_reverse = .false.
    end if
    if (do_reverse) then
        d_conv_e_to_C = x / e_charge
    else
        d_conv_e_to_C = x * e_charge
    end if
end function d_conv_e_to_C

! ======================================================================

elemental real(real32) function s_conv_cal_to_J(x, reverse)
    real(real32), intent(in) :: x
    logical, intent(in), optional :: reverse
    logical :: do_reverse
    if (present(reverse)) then
        do_reverse = reverse
    else
        do_reverse = .false.
    end if
    if (do_reverse) then
        s_conv_cal_to_J = x / real(calorie, kind=real32)
    else
        s_conv_cal_to_J = x * real(calorie, kind=real32)
    end if
end function s_conv_cal_to_J

! ======================================================================

elemental real(real64) function d_conv_cal_to_J(x, reverse)
    real(real64), intent(in) :: x
    logical, intent(in), optional :: reverse
    logical :: do_reverse
    if (present(reverse)) then
        do_reverse = reverse
    else
        do_reverse = .false.
    end if
    if (do_reverse) then
        d_conv_cal_to_J = x / calorie
    else
        d_conv_cal_to_J = x * calorie
    end if
end function d_conv_cal_to_J

! ======================================================================

elemental real(real32) function s_conv_au_freq_to_cm1(x, reverse)
    !! Convert frequency in atomic unit to energy in cm-1.
    real(real32), intent(in) :: x
    logical, intent(in), optional :: reverse
    logical :: do_reverse
    if (present(reverse)) then
        do_reverse = reverse
    else
        do_reverse = .false.
    end if
    if (do_reverse) then
        s_conv_au_freq_to_cm1 = x/real(E_hartree/(planck*slight), kind=real32)
    else
        s_conv_au_freq_to_cm1 = x*real(E_hartree/(planck*slight), kind=real32)
    end if
end function s_conv_au_freq_to_cm1

! ======================================================================

elemental real(real64) function d_conv_au_freq_to_cm1(x, reverse)
    !! Convert frequency in atomic unit to energy in cm-1.
    real(real64), intent(in) :: x
    logical, intent(in), optional :: reverse
    logical :: do_reverse
    if (present(reverse)) then
        do_reverse = reverse
    else
        do_reverse = .false.
    end if
    if (do_reverse) then
        d_conv_au_freq_to_cm1 = x / (E_hartree / (planck * slight))
    else
        d_conv_au_freq_to_cm1 = x * E_hartree / (planck * slight)
    end if
end function d_conv_au_freq_to_cm1

! ======================================================================

elemental real(real32) function s_conv_hartree_to_J(x, reverse)
    real(real32), intent(in) :: x
    logical, intent(in), optional :: reverse
    logical :: do_reverse
    if (present(reverse)) then
        do_reverse = reverse
    else
        do_reverse = .false.
    end if
    if (do_reverse) then
        s_conv_hartree_to_J = x / real(E_hartree, kind=real32)
    else
        s_conv_hartree_to_J = x * real(E_hartree, kind=real32)
    end if
end function s_conv_hartree_to_J

! ======================================================================

elemental real(real64) function d_conv_hartree_to_J(x, reverse)
    real(real64), intent(in) :: x
    logical, intent(in), optional :: reverse
    logical :: do_reverse
    if (present(reverse)) then
        do_reverse = reverse
    else
        do_reverse = .false.
    end if
    if (do_reverse) then
        d_conv_hartree_to_J = x / E_hartree
    else
        d_conv_hartree_to_J = x * E_hartree
    end if
end function d_conv_hartree_to_J

! ======================================================================

elemental real(real32) function s_conv_dE_au_to_cgs(x, derord)
    !! Return energy derivatives from au to cm^-(derorder/2).
    real(real32), intent(in) :: x
    integer, intent(in), optional :: derord
    !! Derivative order
    integer :: dn
    if (present(derord)) then
        dn = max(0, derord)
    else
        dn = 2
    end if
    s_conv_dE_au_to_cgs = x &
        * real((E_hartree/hc)/(sqrt(factG)*bohr_radius)**dn, kind=real32)
end function s_conv_dE_au_to_cgs

! ======================================================================

elemental real(real64) function d_conv_dE_au_to_cgs(x, derord)
    !! Return energy derivatives from au to cm^-(derorder/2).
    real(real64), intent(in) :: x
    integer, intent(in), optional :: derord
    !! Derivative order
    integer :: dn
    if (present(derord)) then
        dn = max(0, derord)
    else
        dn = 2
    end if
    d_conv_dE_au_to_cgs = x * E_hartree/hc/(sqrt(factG)*bohr_radius)**dn
end function d_conv_dE_au_to_cgs

! ======================================================================

end module physics
