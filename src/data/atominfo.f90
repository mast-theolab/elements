module atominfo
    !! Module containing atomic data and properties
    !!
    !! Sources:
    !! masses: Wikipedia, Gaussian
    !! covalent radii: https://en.wikipedia.org/wiki/Covalent_radius
    !! vdW radii: multiple sources are used, see function `get_atom_rvdw`
    !! viewing radii (rvis): Molden internal database.
    !! color: inspired by JMol CPK coloring
    !!        https://en.wikipedia.org/wiki/CPK_coloring

    use iso_fortran_env, only: real64
    use physics, only: bohr => bohr_radius
    use string, only: locase
    use exception, only: BaseException, InitError, RaiseArgError

    implicit none
    
    private
    type :: AtomDB
        character(len=2) :: symbol           ! atomic symbol
        character(len=40) :: name            ! atom name
        integer :: number                    ! atomic number
        real(real64) :: mass                 ! atomic mass (in amu)
        real(real64), dimension(3) :: rcov   ! covalent radii
        ! Covalent radii are provided for single, double and triple bonds
        real(real64), dimension(4) :: rvdw   ! Van der Waals radii
        real(real64) :: rvis                 ! Visible radius, for display
        real(real64) :: rcovG                ! Covalent rad. used by Gaussian
        integer, dimension(3) :: color       ! Color as RGB vector.
    contains
        procedure :: get_rcov => get_atom_rcov
        procedure :: get_rvdw => get_atom_rvdw
    end type AtomDB

    type(AtomDB), dimension(118), public :: atdata

    data atdata(  1) / AtomDB( &
        symbol='H', name='Hydrogen', &
        number=1, &
        mass=1.00790_real64, &
        rcov=[32, 0, 0]/(100.0_real64*bohr), &
        rvdw=[120, 110, 120, 154]/(100.0_real64*bohr), &
        rvis=0.200_real64/bohr, &
        rcovG=0.643_real64, &
        color=[255, 255, 255] &
    ) /
    data atdata(  2) / AtomDB( &
        symbol='He', name='Helium', &
        number=2, &
        mass=4.00260_real64, &
        rcov=[46, 0, 0]/(100.0_real64*bohr), &
        rvdw=[140, 140, 143, 134]/(100.0_real64*bohr), &
        rvis=0.286_real64/bohr, &
        rcovG=0.643_real64, &
        color=[217, 255, 255] &
    ) /
    data atdata(  3) / AtomDB( &
        symbol='Li', name='Lithium', &
        number=3, &
        mass=6.94000_real64, &
        rcov=[133, 124, 0]/(100.0_real64*bohr), &
        rvdw=[181, 181, 212, 220]/(100.0_real64*bohr), &
        rvis=0.340_real64/bohr, &
        rcovG=0.457_real64, &
        color=[204, 128, 255] &
    ) /
    data atdata(  4) / AtomDB( &
        symbol='Be', name='Beryllium', &
        number=4, &
        mass=9.01218_real64, &
        rcov=[102, 90, 84]/(100.0_real64*bohr), &
        rvdw=[0, 153, 198, 219]/(100.0_real64*bohr), &
        rvis=0.589_real64/bohr, &
        rcovG=0.909_real64, &
        color=[194, 255, 0] &
    ) /
    data atdata(  5) / AtomDB( &
        symbol='B', name='Boron', &
        number=5, &
        mass=10.81000_real64, &
        rcov=[85, 78, 73]/(100.0_real64*bohr), &
        rvdw=[0, 192, 191, 205]/(100.0_real64*bohr), &
        rvis=0.415_real64/bohr, &
        rcovG=0.587_real64, &
        color=[255, 181, 181] &
    ) /
    data atdata(  6) / AtomDB( &
        symbol='C', name='Carbon', &
        number=6, &
        mass=12.01100_real64, &
        rcov=[75, 67, 60]/(100.0_real64*bohr), &
        rvdw=[170, 170, 177, 190]/(100.0_real64*bohr), &
        rvis=0.400_real64/bohr, &
        rcovG=0.436_real64, &
        color=[144, 144, 144] &
    ) /
    data atdata(  7) / AtomDB( &
        symbol='N', name='Nitrogen', &
        number=7, &
        mass=14.00670_real64, &
        rcov=[71, 60, 54]/(100.0_real64*bohr), &
        rvdw=[155, 155, 166, 179]/(100.0_real64*bohr), &
        rvis=0.400_real64/bohr, &
        rcovG=0.209_real64, &
        color=[48, 80, 248] &
    ) /
    data atdata(  8) / AtomDB( &
        symbol='O', name='Oxygen', &
        number=8, &
        mass=15.99940_real64, &
        rcov=[63, 57, 53]/(100.0_real64*bohr), &
        rvdw=[152, 152, 150, 171]/(100.0_real64*bohr), &
        rvis=0.400_real64/bohr, &
        rcovG=1.096_real64, &
        color=[255, 13, 13] &
    ) /
    data atdata(  9) / AtomDB( &
        symbol='F', name='Fluorine', &
        number=9, &
        mass=18.99840_real64, &
        rcov=[64, 59, 53]/(100.0_real64*bohr), &
        rvdw=[147, 147, 146, 163]/(100.0_real64*bohr), &
        rvis=0.320_real64/bohr, &
        rcovG=1.020_real64, &
        color=[144, 224,  80] &
    ) /
    data atdata( 10) / AtomDB( &
        symbol='Ne', name='Neon', &
        number=10, &
        mass=20.17900_real64, &
        rcov=[67, 96, 0]/(100.0_real64*bohr), &
        rvdw=[154, 154, 158, 156]/(100.0_real64*bohr), &
        rvis=0.423_real64/bohr, &
        rcovG=0.945_real64, &
        color=[179, 227, 245] &
    ) /
    data atdata( 11) / AtomDB( &
        symbol='Na', name='Sodium', &
        number=11, &
        mass=22.98977_real64, &
        rcov=[155, 160, 0]/(100.0_real64*bohr), &
        rvdw=[227, 227, 250, 225]/(100.0_real64*bohr), &
        rvis=0.485_real64/bohr, &
        rcovG=2.986_real64, &
        color=[171, 92, 242] &
    ) /
    data atdata( 12) / AtomDB( &
        symbol='Mg', name='Magnesium', &
        number=12, &
        mass=24.30500_real64, &
        rcov=[139, 132, 127]/(100.0_real64*bohr), &
        rvdw=[173, 173, 251, 240]/(100.0_real64*bohr), &
        rvis=0.550_real64/bohr, &
        rcovG=2.646_real64, &
        color=[138, 255, 0] &
    ) /
    data atdata( 13) / AtomDB( &
        symbol='Al', name='Aluminium', &
        number=13, &
        mass=26.98154_real64, &
        rcov=[126, 113, 111]/(100.0_real64*bohr), &
        rvdw=[0, 184, 225, 239]/(100.0_real64*bohr), &
        rvis=0.675_real64/bohr, &
        rcovG=2.400_real64, &
        color=[191, 166, 166] &
    ) /
    data atdata( 14) / AtomDB( &
        symbol='Si', name='Silicon', &
        number=14, &
        mass=28.08550_real64, &
        rcov=[116, 107, 102]/(100.0_real64*bohr), &
        rvdw=[210, 210, 219, 232]/(100.0_real64*bohr), &
        rvis=0.600_real64/bohr, &
        rcovG=2.192_real64, &
        color=[240, 200, 160] &
    ) /
    data atdata( 15) / AtomDB( &
        symbol='P', name='Phosphorus', &
        number=15, &
        mass=30.97376_real64, &
        rcov=[111, 102, 94]/(100.0_real64*bohr), &
        rvdw=[180, 180, 190, 223]/(100.0_real64*bohr), &
        rvis=0.525_real64/bohr, &
        rcovG=2.060_real64, &
        color=[255, 128, 0] &
    ) /
    data atdata( 16) / AtomDB( &
        symbol='S', name='Sulfur', &
        number=16, &
        mass=32.06000_real64, &
        rcov=[103, 94, 95]/(100.0_real64*bohr), &
        rvdw=[180, 180, 189, 214]/(100.0_real64*bohr), &
        rvis=0.510_real64/bohr, &
        rcovG=1.890_real64, &
        color=[255, 255, 48] &
    ) /
    data atdata( 17) / AtomDB( &
        symbol='Cl', name='Chlorine', &
        number=17, &
        mass=35.45300_real64, &
        rcov=[99, 95, 93]/(100.0_real64*bohr), &
        rvdw=[175, 175, 182, 206]/(100.0_real64*bohr), &
        rvis=0.495_real64/bohr, &
        rcovG=1.795_real64, &
        color=[31, 240, 31] &
    ) /
    data atdata( 18) / AtomDB( &
        symbol='Ar', name='Argon', &
        number=18, &
        mass=39.94800_real64, &
        rcov=[96, 107, 96]/(100.0_real64*bohr), &
        rvdw=[188, 188, 183, 197]/(100.0_real64*bohr), &
        rvis=0.508_real64/bohr, &
        rcovG=1.701_real64, &
        color=[128, 209, 227] &
    ) /
    data atdata( 19) / AtomDB( &
        symbol='K', name='Potassium', &
        number=19, &
        mass=39.09830_real64, &
        rcov=[196, 193, 0]/(100.0_real64*bohr), &
        rvdw=[275, 275, 273, 234]/(100.0_real64*bohr), &
        rvis=0.665_real64/bohr, &
        rcovG=3.836_real64, &
        color=[143, 64, 212] &
    ) /
    data atdata( 20) / AtomDB( &
        symbol='Ca', name='Calcium', &
        number=20, &
        mass=40.08000_real64, &
        rcov=[171, 147, 133]/(100.0_real64*bohr), &
        rvdw=[0, 231, 262, 270]/(100.0_real64*bohr), &
        rvis=0.495_real64/bohr, &
        rcovG=3.288_real64, &
        color=[61, 255, 0] &
    ) /
    data atdata( 21) / AtomDB( &
        symbol='Sc', name='Scandium', &
        number=21, &
        mass=44.95590_real64, &
        rcov=[148, 116, 114]/(100.0_real64*bohr), &
        rvdw=[0, 0, 258, 263]/(100.0_real64*bohr), &
        rvis=0.735_real64/bohr, &
        rcovG=2.721_real64, &
        color=[230, 230, 230] &
    ) /
    data atdata( 22) / AtomDB( &
        symbol='Ti', name='Titanium', &
        number=22, &
        mass=47.90000_real64, &
        rcov=[136, 117, 108]/(100.0_real64*bohr), &
        rvdw=[0, 0, 246, 257]/(100.0_real64*bohr), &
        rvis=0.720_real64/bohr, &
        rcovG=2.494_real64, &
        color=[191, 194, 199] &
    ) /
    data atdata( 23) / AtomDB( &
        symbol='V', name='Vanadium', &
        number=23, &
        mass=50.94150_real64, &
        rcov=[134, 112, 106]/(100.0_real64*bohr), &
        rvdw=[0, 0, 242, 252]/(100.0_real64*bohr), &
        rvis=0.665_real64/bohr, &
        rcovG=2.305_real64, &
        color=[166, 166, 171] &
    ) /
    data atdata( 24) / AtomDB( &
        symbol='Cr', name='Chromium', &
        number=24, &
        mass=51.99600_real64, &
        rcov=[122, 111, 103]/(100.0_real64*bohr), &
        rvdw=[0, 0, 245, 233]/(100.0_real64*bohr), &
        rvis=0.675_real64/bohr, &
        rcovG=2.230_real64, &
        color=[138, 153, 199] &
    ) /
    data atdata( 25) / AtomDB( &
        symbol='Mn', name='Manganese', &
        number=25, &
        mass=54.93800_real64, &
        rcov=[119, 105, 103]/(100.0_real64*bohr), &
        rvdw=[0, 0, 245, 242]/(100.0_real64*bohr), &
        rvis=0.675_real64/bohr, &
        rcovG=2.211_real64, &
        color=[156, 122, 199] &
    ) /
    data atdata( 26) / AtomDB( &
        symbol='Fe', name='Iron', &
        number=26, &
        mass=55.84700_real64, &
        rcov=[116, 109, 102]/(100.0_real64*bohr), &
        rvdw=[0, 0, 244, 237]/(100.0_real64*bohr), &
        rvis=0.670_real64/bohr, &
        rcovG=2.211_real64, &
        color=[224, 102, 51] &
    ) /
    data atdata( 27) / AtomDB( &
        symbol='Co', name='Cobalt', &
        number=27, &
        mass=58.93320_real64, &
        rcov=[111, 103, 96]/(100.0_real64*bohr), &
        rvdw=[0, 0, 240, 233]/(100.0_real64*bohr), &
        rvis=0.615_real64/bohr, &
        rcovG=2.192_real64, &
        color=[240, 144, 160] &
    ) /
    data atdata( 28) / AtomDB( &
        symbol='Ni', name='Nickel', &
        number=28, &
        mass=58.71000_real64, &
        rcov=[110, 101, 101]/(100.0_real64*bohr), &
        rvdw=[163, 0, 240, 229]/(100.0_real64*bohr), &
        rvis=0.750_real64/bohr, &
        rcovG=2.173_real64, &
        color=[80, 208, 80] &
    ) /
    data atdata( 29) / AtomDB( &
        symbol='Cu', name='Copper', &
        number=29, &
        mass=63.54600_real64, &
        rcov=[112, 115, 120]/(100.0_real64*bohr), &
        rvdw=[140, 0, 238, 217]/(100.0_real64*bohr), &
        rvis=0.760_real64/bohr, &
        rcovG=2.211_real64, &
        color=[200, 128, 51] &
    ) /
    data atdata( 30) / AtomDB( &
        symbol='Zn', name='Zinc', &
        number=30, &
        mass=65.38000_real64, &
        rcov=[118, 120, 0]/(100.0_real64*bohr), &
        rvdw=[139, 0, 239, 222]/(100.0_real64*bohr), &
        rvis=0.725_real64/bohr, &
        rcovG=2.362_real64, &
        color=[125, 128, 176] &
    ) /
    data atdata( 31) / AtomDB( &
        symbol='Ga', name='Gallium', &
        number=31, &
        mass=69.73500_real64, &
        rcov=[124, 117, 121]/(100.0_real64*bohr), &
        rvdw=[187, 187, 232, 233]/(100.0_real64*bohr), &
        rvis=0.610_real64/bohr, &
        rcovG=2.381_real64, &
        color=[194, 143, 143] &
    ) /
    data atdata( 32) / AtomDB( &
        symbol='Ge', name='Germanium', &
        number=32, &
        mass=72.59000_real64, &
        rcov=[121, 111, 114]/(100.0_real64*bohr), &
        rvdw=[0, 211, 229, 234]/(100.0_real64*bohr), &
        rvis=0.585_real64/bohr, &
        rcovG=2.305_real64, &
        color=[102, 143, 143] &
    ) /
    data atdata( 33) / AtomDB( &
        symbol='As', name='Arsenic', &
        number=33, &
        mass=74.92160_real64, &
        rcov=[121, 114, 106]/(100.0_real64*bohr), &
        rvdw=[185, 185, 188, 231]/(100.0_real64*bohr), &
        rvis=0.605_real64/bohr, &
        rcovG=2.268_real64, &
        color=[189, 128, 227] &
    ) /
    data atdata( 34) / AtomDB( &
        symbol='Se', name='Selenium', &
        number=34, &
        mass=78.96000_real64, &
        rcov=[116, 107, 107]/(100.0_real64*bohr), &
        rvdw=[190, 190, 182, 224]/(100.0_real64*bohr), &
        rvis=0.610_real64/bohr, &
        rcovG=2.192_real64, &
        color=[255, 161, 0] &
    ) /
    data atdata( 35) / AtomDB( &
        symbol='Br', name='Bromine', &
        number=35, &
        mass=79.90400_real64, &
        rcov=[114, 109, 110]/(100.0_real64*bohr), &
        rvdw=[183, 183, 186, 219]/(100.0_real64*bohr), &
        rvis=0.605_real64/bohr, &
        rcovG=2.154_real64, &
        color=[166,  41,  41] &
    ) /
    data atdata( 36) / AtomDB( &
        symbol='Kr', name='Krypton', &
        number=36, &
        mass=83.80000_real64, &
        rcov=[117, 121, 108]/(100.0_real64*bohr), &
        rvdw=[202, 202, 225, 212]/(100.0_real64*bohr), &
        rvis=0.524_real64/bohr, &
        rcovG=2.116_real64, &
        color=[92, 184, 209] &
    ) /
    data atdata( 37) / AtomDB( &
        symbol='Rb', name='Rubidium', &
        number=37, &
        mass=85.46780_real64, &
        rcov=[210, 202, 0]/(100.0_real64*bohr), &
        rvdw=[0, 303, 321, 240]/(100.0_real64*bohr), &
        rvis=0.735_real64/bohr, &
        rcovG=4.082_real64, &
        color=[112, 46, 176] &
    ) /
    data atdata( 38) / AtomDB( &
        symbol='Sr', name='Strontium', &
        number=38, &
        mass=87.62000_real64, &
        rcov=[185, 157, 139]/(100.0_real64*bohr), &
        rvdw=[0, 249, 284, 279]/(100.0_real64*bohr), &
        rvis=0.560_real64/bohr, &
        rcovG=3.609_real64, &
        color=[0, 255, 0] &
    ) /
    data atdata( 39) / AtomDB( &
        symbol='Y', name='Yttrium', &
        number=39, &
        mass=88.90590_real64, &
        rcov=[163, 130, 124]/(100.0_real64*bohr), &
        rvdw=[0, 0, 275, 274]/(100.0_real64*bohr), &
        rvis=0.890_real64/bohr, &
        rcovG=3.061_real64, &
        color=[148, 255, 255] &
    ) /
    data atdata( 40) / AtomDB( &
        symbol='Zr', name='Zirconium', &
        number=40, &
        mass=91.22000_real64, &
        rcov=[154, 127, 121]/(100.0_real64*bohr), &
        rvdw=[0, 0, 252, 269]/(100.0_real64*bohr), &
        rvis=0.780_real64/bohr, &
        rcovG=2.740_real64, &
        color=[148, 224, 224] &
    ) /
    data atdata( 41) / AtomDB( &
        symbol='Nb', name='Niobium', &
        number=41, &
        mass=92.90640_real64, &
        rcov=[147, 125, 116]/(100.0_real64*bohr), &
        rvdw=[0, 0, 256, 251]/(100.0_real64*bohr), &
        rvis=0.740_real64/bohr, &
        rcovG=2.532_real64, &
        color=[115, 194, 201] &
    ) /
    data atdata( 42) / AtomDB( &
        symbol='Mo', name='Molybdenum', &
        number=42, &
        mass=95.94000_real64, &
        rcov=[138, 121, 113]/(100.0_real64*bohr), &
        rvdw=[0, 0, 245, 244]/(100.0_real64*bohr), &
        rvis=0.735_real64/bohr, &
        rcovG=2.457_real64, &
        color=[84, 181, 181] &
    ) /
    data atdata( 43) / AtomDB( &
        symbol='Tc', name='Technetium', &
        number=43, &
        mass=98.90620_real64, &
        rcov=[128, 120, 110]/(100.0_real64*bohr), &
        rvdw=[0, 0, 244, 252]/(100.0_real64*bohr), &
        rvis=0.675_real64/bohr, &
        rcovG=2.400_real64, &
        color=[59, 158, 158] &
    ) /
    data atdata( 44) / AtomDB( &
        symbol='Ru', name='Ruthenium', &
        number=44, &
        mass=101.0700_real64, &
        rcov=[125, 114, 103]/(100.0_real64*bohr), &
        rvdw=[0, 0, 246, 237]/(100.0_real64*bohr), &
        rvis=0.700_real64/bohr, &
        rcovG=2.362_real64, &
        color=[36, 143, 143] &
    ) /
    data atdata( 45) / AtomDB( &
        symbol='Rh', name='Rhodium', &
        number=45, &
        mass=102.9055_real64, &
        rcov=[125, 110, 106]/(100.0_real64*bohr), &
        rvdw=[0, 0, 244, 233]/(100.0_real64*bohr), &
        rvis=0.725_real64/bohr, &
        rcovG=2.362_real64, &
        color=[10, 125, 140] &
    ) /
    data atdata( 46) / AtomDB( &
        symbol='Pd', name='Palladium', &
        number=46, &
        mass=106.4000_real64, &
        rcov=[120, 117, 112]/(100.0_real64*bohr), &
        rvdw=[163, 0, 215, 215]/(100.0_real64*bohr), &
        rvis=0.750_real64/bohr, &
        rcovG=2.419_real64, &
        color=[0, 105, 133] &
    ) /
    data atdata( 47) / AtomDB( &
        symbol='Ag', name='Silver', &
        number=47, &
        mass=107.8680_real64, &
        rcov=[128, 139, 137]/(100.0_real64*bohr), &
        rvdw=[172, 0, 253, 225]/(100.0_real64*bohr), &
        rvis=0.795_real64/bohr, &
        rcovG=2.532_real64, &
        color=[192, 192, 192] &
    ) /
    data atdata( 48) / AtomDB( &
        symbol='Cd', name='Cadmium', &
        number=48, &
        mass=112.4100_real64, &
        rcov=[136, 144, 0]/(100.0_real64*bohr), &
        rvdw=[158, 0, 249, 238]/(100.0_real64*bohr), &
        rvis=0.845_real64/bohr, &
        rcovG=2.797_real64, &
        color=[255, 217, 143] &
    ) /
    data atdata( 49) / AtomDB( &
        symbol='In', name='Indium', &
        number=49, &
        mass=114.8200_real64, &
        rcov=[142, 136, 146]/(100.0_real64*bohr), &
        rvdw=[193, 193, 243, 246]/(100.0_real64*bohr), &
        rvis=0.815_real64/bohr, &
        rcovG=2.721_real64, &
        color=[166, 117, 115] &
    ) /
    data atdata( 50) / AtomDB( &
        symbol='Sn', name='Tin', &
        number=50, &
        mass=118.6900_real64, &
        rcov=[140, 130, 132]/(100.0_real64*bohr), &
        rvdw=[217, 217, 242, 248]/(100.0_real64*bohr), &
        rvis=0.730_real64/bohr, &
        rcovG=2.665_real64, &
        color=[102, 128, 128] &
    ) /
    data atdata( 51) / AtomDB( &
        symbol='Sb', name='Antimony', &
        number=51, &
        mass=121.7500_real64, &
        rcov=[140, 133, 127]/(100.0_real64*bohr), &
        rvdw=[0, 206, 247, 246]/(100.0_real64*bohr), &
        rvis=0.730_real64/bohr, &
        rcovG=2.646_real64, &
        color=[158, 99, 181] &
    ) /
    data atdata( 52) / AtomDB( &
        symbol='Te', name='Tellurium', &
        number=52, &
        mass=127.6000_real64, &
        rcov=[136, 128, 121]/(100.0_real64*bohr), &
        rvdw=[206, 206, 199, 242]/(100.0_real64*bohr), &
        rvis=0.735_real64/bohr, &
        rcovG=2.570_real64, &
        color=[212, 122, 0] &
    ) /
    data atdata( 53) / AtomDB( &
        symbol='I', name='Iodine', &
        number=53, &
        mass=126.9045_real64, &
        rcov=[133, 129, 125]/(100.0_real64*bohr), &
        rvdw=[198, 198, 204, 238]/(100.0_real64*bohr), &
        rvis=0.700_real64/bohr, &
        rcovG=2.513_real64, &
        color=[148, 0, 148] &
    ) /
    data atdata( 54) / AtomDB( &
        symbol='Xe', name='Xenon', &
        number=54, &
        mass=131.3000_real64, &
        rcov=[131, 135, 122]/(100.0_real64*bohr), &
        rvdw=[216, 216, 206, 232]/(100.0_real64*bohr), &
        rvis=0.577_real64/bohr, &
        rcovG=2.476_real64, &
        color=[66, 158, 176] &
    ) /
    data atdata( 55) / AtomDB( &
        symbol='Cs', name='Caesium', &
        number=55, &
        mass=132.9054_real64, &
        rcov=[232, 209, 0]/(100.0_real64*bohr), &
        rvdw=[0, 343, 348, 249]/(100.0_real64*bohr), &
        rvis=0.835_real64/bohr, &
        rcovG=4.441_real64, &
        color=[87, 23, 143] &
    ) /
    data atdata( 56) / AtomDB( &
        symbol='Ba', name='Barium', &
        number=56, &
        mass=137.3300_real64, &
        rcov=[196, 161, 149]/(100.0_real64*bohr), &
        rvdw=[0, 268, 303, 293]/(100.0_real64*bohr), &
        rvis=0.670_real64/bohr, &
        rcovG=3.742_real64, &
        color=[0, 201, 0] &
    ) /
    data atdata( 57) / AtomDB( &
        symbol='La', name='Lanthanum', &
        number=57, &
        mass=138.9055_real64, &
        rcov=[180, 139, 139]/(100.0_real64*bohr), &
        rvdw=[0, 0, 298, 284]/(100.0_real64*bohr), &
        rvis=0.935_real64/bohr, &
        rcovG=3.194_real64, &
        color=[112, 212, 255] &
    ) /
    data atdata( 58) / AtomDB( &
        symbol='Ce', name='Cerium', &
        number=58, &
        mass=140.1200_real64, &
        rcov=[163, 137, 131]/(100.0_real64*bohr), &
        rvdw=[0, 0, 288, 282]/(100.0_real64*bohr), &
        rvis=0.915_real64/bohr, &
        rcovG=3.118_real64, &
        color=[255, 255, 199] &
    ) /
    data atdata( 59) / AtomDB( &
        symbol='Pr', name='Praseodymium', &
        number=59, &
        mass=140.9077_real64, &
        rcov=[176, 138, 128]/(100.0_real64*bohr), &
        rvdw=[0, 0, 292, 286]/(100.0_real64*bohr), &
        rvis=0.910_real64/bohr, &
        rcovG=3.118_real64, &
        color=[217, 255, 199] &
    ) /
    data atdata( 60) / AtomDB( &
        symbol='Nd', name='Neodymium', &
        number=60, &
        mass=144.2400_real64, &
        rcov=[174, 132, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 295, 284]/(100.0_real64*bohr), &
        rvis=0.905_real64/bohr, &
        rcovG=3.099_real64, &
        color=[199, 255, 199] &
    ) /
    data atdata( 61) / AtomDB( &
        symbol='Pm', name='Promethium', &
        number=61, &
        mass=0.0_real64, &
        rcov=[173, 135, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 290, 283]/(100.0_real64*bohr), &
        rvis=0.900_real64/bohr, &
        rcovG=3.080_real64, &
        color=[163, 255, 199] &
    ) /
    data atdata( 62) / AtomDB( &
        symbol='Sm', name='Samarium', &
        number=62, &
        mass=150.3600_real64, &
        rcov=[172, 134, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 290, 280]/(100.0_real64*bohr), &
        rvis=0.900_real64/bohr, &
        rcovG=3.061_real64, &
        color=[143, 255, 199] &
    ) /
    data atdata( 63) / AtomDB( &
        symbol='Eu', name='Europium', &
        number=63, &
        mass=151.9600_real64, &
        rcov=[168, 134, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 287, 280]/(100.0_real64*bohr), &
        rvis=0.995_real64/bohr, &
        rcovG=3.496_real64, &
        color=[97, 255, 199] &
    ) /
    data atdata( 64) / AtomDB( &
        symbol='Gd', name='Gadolinium', &
        number=64, &
        mass=157.2500_real64, &
        rcov=[169, 135, 132]/(100.0_real64*bohr), &
        rvdw=[0, 0, 283, 277]/(100.0_real64*bohr), &
        rvis=0.895_real64/bohr, &
        rcovG=3.042_real64, &
        color=[69, 255, 199] &
    ) /
    data atdata( 65) / AtomDB( &
        symbol='Tb', name='Terbium', &
        number=65, &
        mass=158.9253_real64, &
        rcov=[168, 135, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 279, 276]/(100.0_real64*bohr), &
        rvis=0.880_real64/bohr, &
        rcovG=3.005_real64, &
        color=[48, 255, 199] &
    ) /
    data atdata( 66) / AtomDB( &
        symbol='Dy', name='Dysprosium', &
        number=66, &
        mass=162.5000_real64, &
        rcov=[167, 133, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 287, 275]/(100.0_real64*bohr), &
        rvis=0.875_real64/bohr, &
        rcovG=3.005_real64, &
        color=[31, 255, 199] &
    ) /
    data atdata( 67) / AtomDB( &
        symbol='Ho', name='Holmium', &
        number=67, &
        mass=164.9303_real64, &
        rcov=[166, 133, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 281, 273]/(100.0_real64*bohr), &
        rvis=0.870_real64/bohr, &
        rcovG=2.986_real64, &
        color=[0, 255, 156] &
    ) /
    data atdata( 68) / AtomDB( &
        symbol='Er', name='Erbium', &
        number=68, &
        mass=167.2600_real64, &
        rcov=[165, 133, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 283, 272]/(100.0_real64*bohr), &
        rvis=0.865_real64/bohr, &
        rcovG=2.967_real64, &
        color=[0, 230, 117] &
    ) /
    data atdata( 69) / AtomDB( &
        symbol='Tm', name='Thulium', &
        number=69, &
        mass=168.9342_real64, &
        rcov=[164, 131, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 279, 271]/(100.0_real64*bohr), &
        rvis=0.860_real64/bohr, &
        rcovG=2.948_real64, &
        color=[0, 212, 82] &
    ) /
    data atdata( 70) / AtomDB( &
        symbol='Yb', name='Ytterbium', &
        number=70, &
        mass=173.0500_real64, &
        rcov=[170, 129, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 280, 277]/(100.0_real64*bohr), &
        rvis=0.970_real64/bohr, &
        rcovG=2.948_real64, &
        color=[0, 191, 56] &
    ) /
    data atdata( 71) / AtomDB( &
        symbol='Lu', name='Lutetium', &
        number=71, &
        mass=174.9670_real64, &
        rcov=[162, 131, 131]/(100.0_real64*bohr), &
        rvdw=[0, 0, 274, 270]/(100.0_real64*bohr), &
        rvis=0.860_real64/bohr, &
        rcovG=2.948_real64, &
        color=[0, 171, 36] &
    ) /
    data atdata( 72) / AtomDB( &
        symbol='Hf', name='Hafnium', &
        number=72, &
        mass=178.4900_real64, &
        rcov=[152, 128, 122]/(100.0_real64*bohr), &
        rvdw=[0, 0, 263, 264]/(100.0_real64*bohr), &
        rvis=0.785_real64/bohr, &
        rcovG=2.721_real64, &
        color=[77, 194, 255] &
    ) /
    data atdata( 73) / AtomDB( &
        symbol='Ta', name='Tantalum', &
        number=73, &
        mass=180.9479_real64, &
        rcov=[146, 126, 119]/(100.0_real64*bohr), &
        rvdw=[0, 0, 253, 258]/(100.0_real64*bohr), &
        rvis=0.715_real64/bohr, &
        rcovG=2.532_real64, &
        color=[77, 166, 255] &
    ) /
    data atdata( 74) / AtomDB( &
        symbol='W', name='Tungsten', &
        number=74, &
        mass=183.8400_real64, &
        rcov=[137, 120, 115]/(100.0_real64*bohr), &
        rvdw=[0, 0, 257, 253]/(100.0_real64*bohr), &
        rvis=0.685_real64/bohr, &
        rcovG=2.457_real64, &
        color=[33, 148, 214] &
    ) /
    data atdata( 75) / AtomDB( &
        symbol='Re', name='Rhenium', &
        number=75, &
        mass=186.2070_real64, &
        rcov=[131, 119, 110]/(100.0_real64*bohr), &
        rvdw=[0, 0, 249, 249]/(100.0_real64*bohr), &
        rvis=0.675_real64/bohr, &
        rcovG=2.419_real64, &
        color=[38, 125, 171] &
    ) /
    data atdata( 76) / AtomDB( &
        symbol='Os', name='Osmium', &
        number=76, &
        mass=190.2000_real64, &
        rcov=[129, 116, 109]/(100.0_real64*bohr), &
        rvdw=[0, 0, 248, 244]/(100.0_real64*bohr), &
        rvis=0.685_real64/bohr, &
        rcovG=2.381_real64, &
        color=[38, 102, 150] &
    ) /
    data atdata( 77) / AtomDB( &
        symbol='Ir', name='Iridium', &
        number=77, &
        mass=192.2200_real64, &
        rcov=[122, 115, 107]/(100.0_real64*bohr), &
        rvdw=[0, 0, 241, 240]/(100.0_real64*bohr), &
        rvis=0.660_real64/bohr, &
        rcovG=2.400_real64, &
        color=[23, 84, 135] &
    ) /
    data atdata( 78) / AtomDB( &
        symbol='Pt', name='Platinum', &
        number=78, &
        mass=195.0900_real64, &
        rcov=[123, 112, 110]/(100.0_real64*bohr), &
        rvdw=[172, 0, 229, 230]/(100.0_real64*bohr), &
        rvis=0.750_real64/bohr, &
        rcovG=2.457_real64, &
        color=[208, 208, 224] &
    ) /
    data atdata( 79) / AtomDB( &
        symbol='Au', name='Gold', &
        number=79, &
        mass=196.9665_real64, &
        rcov=[124, 121, 123]/(100.0_real64*bohr), &
        rvdw=[166, 0, 232, 226]/(100.0_real64*bohr), &
        rvis=0.750_real64/bohr, &
        rcovG=2.532_real64, &
        color=[255, 209, 35] &
    ) /
    data atdata( 80) / AtomDB( &
        symbol='Hg', name='Mercury', &
        number=80, &
        mass=200.5900_real64, &
        rcov=[133, 142, 0]/(100.0_real64*bohr), &
        rvdw=[170, 0, 245, 229]/(100.0_real64*bohr), &
        rvis=0.850_real64/bohr, &
        rcovG=2.816_real64, &
        color=[184, 184, 208] &
    ) /
    data atdata( 81) / AtomDB( &
        symbol='Tl', name='Thallium', &
        number=81, &
        mass=204.3700_real64, &
        rcov=[144, 142, 150]/(100.0_real64*bohr), &
        rvdw=[196, 196, 247, 242]/(100.0_real64*bohr), &
        rvis=0.775_real64/bohr, &
        rcovG=2.797_real64, &
        color=[166, 84, 77] &
    ) /
    data atdata( 82) / AtomDB( &
        symbol='Pb', name='Lead', &
        number=82, &
        mass=207.2000_real64, &
        rcov=[144, 135, 137]/(100.0_real64*bohr), &
        rvdw=[202, 202, 260, 249]/(100.0_real64*bohr), &
        rvis=0.770_real64/bohr, &
        rcovG=2.778_real64, &
        color=[87, 89, 97] &
    ) /
    data atdata( 83) / AtomDB( &
        symbol='Bi', name='Bismuth', &
        number=83, &
        mass=208.9804_real64, &
        rcov=[151, 141, 135]/(100.0_real64*bohr), &
        rvdw=[0, 207, 254, 250]/(100.0_real64*bohr), &
        rvis=0.770_real64/bohr, &
        rcovG=2.759_real64, &
        color=[158, 79, 181] &
    ) /
    data atdata( 84) / AtomDB( &
        symbol='Po', name='Polonium', &
        number=84, &
        mass=0.0_real64, &
        rcov=[145, 135, 129]/(100.0_real64*bohr), &
        rvdw=[0, 197, 0, 250]/(100.0_real64*bohr), &
        rvis=0.840_real64/bohr, &
        rcovG=2.759_real64, &
        color=[171, 92, 0] &
    ) /
    data atdata( 85) / AtomDB( &
        symbol='At', name='Astatine', &
        number=85, &
        mass=0.0_real64, &
        rcov=[147, 138, 138]/(100.0_real64*bohr), &
        rvdw=[0, 202, 0, 247]/(100.0_real64*bohr), &
        rvis=1.000_real64/bohr, &
        rcovG=2.740_real64, &
        color=[117, 79, 69] &
    ) /
    data atdata( 86) / AtomDB( &
        symbol='Rn', name='Radon', &
        number=86, &
        mass=0.0_real64, &
        rcov=[142, 145, 133]/(100.0_real64*bohr), &
        rvdw=[0, 220, 0, 243]/(100.0_real64*bohr), &
        rvis=1.000_real64/bohr, &
        rcovG=0.0_real64, &
        color=[66, 130, 150] &
    ) /
    data atdata( 87) / AtomDB( &
        symbol='Fr', name='Francium', &
        number=87, &
        mass=0.0_real64, &
        rcov=[223, 218, 0]/(100.0_real64*bohr), &
        rvdw=[0, 348, 0, 258]/(100.0_real64*bohr), &
        rvis=1.000_real64/bohr, &
        rcovG=0.0_real64, &
        color=[66, 0, 102] &
    ) /
    data atdata( 88) / AtomDB( &
        symbol='Ra', name='Radium', &
        number=88, &
        mass=0.0_real64, &
        rcov=[201, 173, 159]/(100.0_real64*bohr), &
        rvdw=[0, 283, 0, 292]/(100.0_real64*bohr), &
        rvis=0.950_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 125, 0] &
    ) /
    data atdata( 89) / AtomDB( &
        symbol='Ac', name='Actinium', &
        number=89, &
        mass=0.0_real64, &
        rcov=[186, 153, 140]/(100.0_real64*bohr), &
        rvdw=[0, 0, 280, 293]/(100.0_real64*bohr), &
        rvis=0.940_real64/bohr, &
        rcovG=0.0_real64, &
        color=[112, 171, 250] &
    ) /
    data atdata( 90) / AtomDB( &
        symbol='Th', name='Thorium', &
        number=90, &
        mass=232.038_real64, &
        rcov=[175, 143, 136]/(100.0_real64*bohr), &
        rvdw=[0, 0, 293, 288]/(100.0_real64*bohr), &
        rvis=0.895_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 186, 255] &
    ) /
    data atdata( 91) / AtomDB( &
        symbol='Pa', name='Protactinium', &
        number=91, &
        mass=231.0359_real64, &
        rcov=[169, 138, 129]/(100.0_real64*bohr), &
        rvdw=[0, 0, 288, 285]/(100.0_real64*bohr), &
        rvis=0.805_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 161, 255] &
    ) /
    data atdata( 92) / AtomDB( &
        symbol='U', name='Uranium', &
        number=92, &
        mass=238.0289_real64, &
        rcov=[170, 134, 118]/(100.0_real64*bohr), &
        rvdw=[186, 0, 271, 283]/(100.0_real64*bohr), &
        rvis=0.790_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 143, 255] &
    ) /
    data atdata( 93) / AtomDB( &
        symbol='Np', name='Neptunium', &
        number=93, &
        mass=0.0_real64, &
        rcov=[171, 133, 116]/(100.0_real64*bohr), &
        rvdw=[0, 0, 282, 281]/(100.0_real64*bohr), &
        rvis=0.775_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 128, 255] &
    ) /
    data atdata( 94) / AtomDB( &
        symbol='Pu', name='Plutonium', &
        number=94, &
        mass=0.0_real64, &
        rcov=[172, 135, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 281, 278]/(100.0_real64*bohr), &
        rvis=1.000_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 107, 255] &
    ) /
    data atdata( 95) / AtomDB( &
        symbol='Am', name='Americium', &
        number=95, &
        mass=0.0_real64, &
        rcov=[166, 135, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 283, 276]/(100.0_real64*bohr), &
        rvis=0.755_real64/bohr, &
        rcovG=0.0_real64, &
        color=[84, 92, 242] &
    ) /
    data atdata( 96) / AtomDB( &
        symbol='Cm', name='Curium', &
        number=96, &
        mass=0.0_real64, &
        rcov=[166, 136, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 305, 264]/(100.0_real64*bohr), &
        rvis=1.000_real64/bohr, &
        rcovG=0.0_real64, &
        color=[120, 92, 227] &
    ) /
    data atdata( 97) / AtomDB( &
        symbol='Bk', name='Berkelium', &
        number=97, &
        mass=0.0_real64, &
        rcov=[168, 139, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 340, 0]/(100.0_real64*bohr), &
        rvis=1.000_real64/bohr, &
        rcovG=0.0_real64, &
        color=[138, 79, 227] &
    ) /
    data atdata( 98) / AtomDB( &
        symbol='Cf', name='Californium', &
        number=98, &
        mass=0.0_real64, &
        rcov=[168, 140, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 305, 0]/(100.0_real64*bohr), &
        rvis=0.765_real64/bohr, &
        rcovG=0.0_real64, &
        color=[161, 54, 212] &
    ) /
    data atdata( 99) / AtomDB( &
        symbol='Es', name='Einsteinium', &
        number=99, &
        mass=0.0_real64, &
        rcov=[165, 140, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 270, 0]/(100.0_real64*bohr), &
        rvis=0.100_real64/bohr, &
        rcovG=0.0_real64, &
        color=[179, 31, 212] &
    ) /
    data atdata(100) / AtomDB( &
        symbol='Fm', name='Fermium', &
        number=100, &
        mass=0.0_real64, &
        rcov=[167, 0, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.900_real64/bohr, &
        rcovG=0.0_real64, &
        color=[179, 31, 186] &
    ) /
    data atdata(101) / AtomDB( &
        symbol='Md', name='Mendelevium', &
        number=101, &
        mass=0.0_real64, &
        rcov=[173, 139, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[179, 13, 166] &
    ) /
    data atdata(102) / AtomDB( &
        symbol='No', name='Nobelium', &
        number=102, &
        mass=0.0_real64, &
        rcov=[176, 0, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[189, 13, 135] &
    ) /
    data atdata(103) / AtomDB( &
        symbol='Lr', name='Lawrencium', &
        number=103, &
        mass=0.0_real64, &
        rcov=[161, 141, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[199, 0, 102] &
    ) /
    data atdata(104) / AtomDB( &
        symbol='Rf', name='Rutherfordium', &
        number=104, &
        mass=0.0_real64, &
        rcov=[157, 140, 131]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[204, 0, 89] &
    ) /
    data atdata(105) / AtomDB( &
        symbol='Db', name='Dubnium', &
        number=105, &
        mass=0.0_real64, &
        rcov=[149, 136, 126]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[209, 0, 79] &
    ) /
    data atdata(106) / AtomDB( &
        symbol='Sg', name='Seaborgium', &
        number=106, &
        mass=0.0_real64, &
        rcov=[143, 128, 121]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[217, 0, 69] &
    ) /
    data atdata(107) / AtomDB( &
        symbol='Bh', name='Bohrium', &
        number=107, &
        mass=0.0_real64, &
        rcov=[141, 128, 119]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[224, 0, 56] &
    ) /
    data atdata(108) / AtomDB( &
        symbol='Hs', name='Hassium', &
        number=108, &
        mass=0.0_real64, &
        rcov=[134, 125, 118]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[230, 0, 46] &
    ) /
    data atdata(109) / AtomDB( &
        symbol='Mt', name='Meitnerium', &
        number=109, &
        mass=0.0_real64, &
        rcov=[129, 125, 118]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[235, 0, 38] &
    ) /
    data atdata(110) / AtomDB( &
        symbol='Ds', name='Darmstadtium', &
        number=110, &
        mass=0.0_real64, &
        rcov=[128, 116, 112]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 0, 0] &
    ) /
    data atdata(111) / AtomDB( &
        symbol='Rg', name='Roentgenium', &
        number=111, &
        mass=0.0_real64, &
        rcov=[121, 116, 118]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 0, 0] &
    ) /
    data atdata(112) / AtomDB( &
        symbol='Cn', name='Copernicium', &
        number=112, &
        mass=0.0_real64, &
        rcov=[122, 116, 118]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 0, 0] &
    ) /
    data atdata(113) / AtomDB( &
        symbol='Nh', name='Nihonium', &
        number=113, &
        mass=0.0_real64, &
        rcov=[136, 0, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 0, 0] &
    ) /
    data atdata(114) / AtomDB( &
        symbol='Fl', name='Flerovium', &
        number=114, &
        mass=0.0_real64, &
        rcov=[143, 0, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 0, 0] &
    ) /
    data atdata(115) / AtomDB( &
        symbol='Mc', name='Moscovium', &
        number=115, &
        mass=0.0_real64, &
        rcov=[162, 0, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 0, 0] &
    ) /
    data atdata(116) / AtomDB( &
        symbol='Lv', name='Livermorium', &
        number=116, &
        mass=0.0_real64, &
        rcov=[175, 0, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 0, 0] &
    ) /
    data atdata(117) / AtomDB( &
        symbol='Ts', name='Tennessine', &
        number=117, &
        mass=0.0_real64, &
        rcov=[165, 0, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 0, 0] &
    ) /
    data atdata(118) / AtomDB( &
        symbol='Og', name='Oganesson', &
        number=118, &
        mass=0.0_real64, &
        rcov=[157, 0, 0]/(100.0_real64*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_real64*bohr), &
        rvis=0.0_real64/bohr, &
        rcovG=0.0_real64, &
        color=[0, 0, 0] &
    ) /

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
    real(real64) :: rcov
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
            rcov = 0.0_real64
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
    real(real64) :: rvdw
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
            rvdw = 0.0_real64
        end if
    end select

end function get_atom_rvdw

! ======================================================================

end module atominfo
