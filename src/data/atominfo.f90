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

    use numeric, only: realwp
    use datatypes, only: AtomDB
    use physics, only: bohr => bohr_radius

    implicit none

    type(AtomDB), dimension(118), public :: atdata

    data atdata(  1) / AtomDB( &
        symbol='H', name='Hydrogen', &
        number=1, &
        mass=1.00790_realwp, &
        rcov=[32, 0, 0]/(100.0_realwp*bohr), &
        rvdw=[120, 110, 120, 154]/(100.0_realwp*bohr), &
        rvis=0.200_realwp/bohr, &
        rcovG=0.643_realwp, &
        color=[255, 255, 255] &
    ) /
    data atdata(  2) / AtomDB( &
        symbol='He', name='Helium', &
        number=2, &
        mass=4.00260_realwp, &
        rcov=[46, 0, 0]/(100.0_realwp*bohr), &
        rvdw=[140, 140, 143, 134]/(100.0_realwp*bohr), &
        rvis=0.286_realwp/bohr, &
        rcovG=0.643_realwp, &
        color=[217, 255, 255] &
    ) /
    data atdata(  3) / AtomDB( &
        symbol='Li', name='Lithium', &
        number=3, &
        mass=6.94000_realwp, &
        rcov=[133, 124, 0]/(100.0_realwp*bohr), &
        rvdw=[181, 181, 212, 220]/(100.0_realwp*bohr), &
        rvis=0.340_realwp/bohr, &
        rcovG=0.457_realwp, &
        color=[204, 128, 255] &
    ) /
    data atdata(  4) / AtomDB( &
        symbol='Be', name='Beryllium', &
        number=4, &
        mass=9.01218_realwp, &
        rcov=[102, 90, 84]/(100.0_realwp*bohr), &
        rvdw=[0, 153, 198, 219]/(100.0_realwp*bohr), &
        rvis=0.589_realwp/bohr, &
        rcovG=0.909_realwp, &
        color=[194, 255, 0] &
    ) /
    data atdata(  5) / AtomDB( &
        symbol='B', name='Boron', &
        number=5, &
        mass=10.81000_realwp, &
        rcov=[85, 78, 73]/(100.0_realwp*bohr), &
        rvdw=[0, 192, 191, 205]/(100.0_realwp*bohr), &
        rvis=0.415_realwp/bohr, &
        rcovG=0.587_realwp, &
        color=[255, 181, 181] &
    ) /
    data atdata(  6) / AtomDB( &
        symbol='C', name='Carbon', &
        number=6, &
        mass=12.01100_realwp, &
        rcov=[75, 67, 60]/(100.0_realwp*bohr), &
        rvdw=[170, 170, 177, 190]/(100.0_realwp*bohr), &
        rvis=0.400_realwp/bohr, &
        rcovG=0.436_realwp, &
        color=[144, 144, 144] &
    ) /
    data atdata(  7) / AtomDB( &
        symbol='N', name='Nitrogen', &
        number=7, &
        mass=14.00670_realwp, &
        rcov=[71, 60, 54]/(100.0_realwp*bohr), &
        rvdw=[155, 155, 166, 179]/(100.0_realwp*bohr), &
        rvis=0.400_realwp/bohr, &
        rcovG=0.209_realwp, &
        color=[48, 80, 248] &
    ) /
    data atdata(  8) / AtomDB( &
        symbol='O', name='Oxygen', &
        number=8, &
        mass=15.99940_realwp, &
        rcov=[63, 57, 53]/(100.0_realwp*bohr), &
        rvdw=[152, 152, 150, 171]/(100.0_realwp*bohr), &
        rvis=0.400_realwp/bohr, &
        rcovG=1.096_realwp, &
        color=[255, 13, 13] &
    ) /
    data atdata(  9) / AtomDB( &
        symbol='F', name='Fluorine', &
        number=9, &
        mass=18.99840_realwp, &
        rcov=[64, 59, 53]/(100.0_realwp*bohr), &
        rvdw=[147, 147, 146, 163]/(100.0_realwp*bohr), &
        rvis=0.320_realwp/bohr, &
        rcovG=1.020_realwp, &
        color=[144, 224,  80] &
    ) /
    data atdata( 10) / AtomDB( &
        symbol='Ne', name='Neon', &
        number=10, &
        mass=20.17900_realwp, &
        rcov=[67, 96, 0]/(100.0_realwp*bohr), &
        rvdw=[154, 154, 158, 156]/(100.0_realwp*bohr), &
        rvis=0.423_realwp/bohr, &
        rcovG=0.945_realwp, &
        color=[179, 227, 245] &
    ) /
    data atdata( 11) / AtomDB( &
        symbol='Na', name='Sodium', &
        number=11, &
        mass=22.98977_realwp, &
        rcov=[155, 160, 0]/(100.0_realwp*bohr), &
        rvdw=[227, 227, 250, 225]/(100.0_realwp*bohr), &
        rvis=0.485_realwp/bohr, &
        rcovG=2.986_realwp, &
        color=[171, 92, 242] &
    ) /
    data atdata( 12) / AtomDB( &
        symbol='Mg', name='Magnesium', &
        number=12, &
        mass=24.30500_realwp, &
        rcov=[139, 132, 127]/(100.0_realwp*bohr), &
        rvdw=[173, 173, 251, 240]/(100.0_realwp*bohr), &
        rvis=0.550_realwp/bohr, &
        rcovG=2.646_realwp, &
        color=[138, 255, 0] &
    ) /
    data atdata( 13) / AtomDB( &
        symbol='Al', name='Aluminium', &
        number=13, &
        mass=26.98154_realwp, &
        rcov=[126, 113, 111]/(100.0_realwp*bohr), &
        rvdw=[0, 184, 225, 239]/(100.0_realwp*bohr), &
        rvis=0.675_realwp/bohr, &
        rcovG=2.400_realwp, &
        color=[191, 166, 166] &
    ) /
    data atdata( 14) / AtomDB( &
        symbol='Si', name='Silicon', &
        number=14, &
        mass=28.08550_realwp, &
        rcov=[116, 107, 102]/(100.0_realwp*bohr), &
        rvdw=[210, 210, 219, 232]/(100.0_realwp*bohr), &
        rvis=0.600_realwp/bohr, &
        rcovG=2.192_realwp, &
        color=[240, 200, 160] &
    ) /
    data atdata( 15) / AtomDB( &
        symbol='P', name='Phosphorus', &
        number=15, &
        mass=30.97376_realwp, &
        rcov=[111, 102, 94]/(100.0_realwp*bohr), &
        rvdw=[180, 180, 190, 223]/(100.0_realwp*bohr), &
        rvis=0.525_realwp/bohr, &
        rcovG=2.060_realwp, &
        color=[255, 128, 0] &
    ) /
    data atdata( 16) / AtomDB( &
        symbol='S', name='Sulfur', &
        number=16, &
        mass=32.06000_realwp, &
        rcov=[103, 94, 95]/(100.0_realwp*bohr), &
        rvdw=[180, 180, 189, 214]/(100.0_realwp*bohr), &
        rvis=0.510_realwp/bohr, &
        rcovG=1.890_realwp, &
        color=[255, 255, 48] &
    ) /
    data atdata( 17) / AtomDB( &
        symbol='Cl', name='Chlorine', &
        number=17, &
        mass=35.45300_realwp, &
        rcov=[99, 95, 93]/(100.0_realwp*bohr), &
        rvdw=[175, 175, 182, 206]/(100.0_realwp*bohr), &
        rvis=0.495_realwp/bohr, &
        rcovG=1.795_realwp, &
        color=[31, 240, 31] &
    ) /
    data atdata( 18) / AtomDB( &
        symbol='Ar', name='Argon', &
        number=18, &
        mass=39.94800_realwp, &
        rcov=[96, 107, 96]/(100.0_realwp*bohr), &
        rvdw=[188, 188, 183, 197]/(100.0_realwp*bohr), &
        rvis=0.508_realwp/bohr, &
        rcovG=1.701_realwp, &
        color=[128, 209, 227] &
    ) /
    data atdata( 19) / AtomDB( &
        symbol='K', name='Potassium', &
        number=19, &
        mass=39.09830_realwp, &
        rcov=[196, 193, 0]/(100.0_realwp*bohr), &
        rvdw=[275, 275, 273, 234]/(100.0_realwp*bohr), &
        rvis=0.665_realwp/bohr, &
        rcovG=3.836_realwp, &
        color=[143, 64, 212] &
    ) /
    data atdata( 20) / AtomDB( &
        symbol='Ca', name='Calcium', &
        number=20, &
        mass=40.08000_realwp, &
        rcov=[171, 147, 133]/(100.0_realwp*bohr), &
        rvdw=[0, 231, 262, 270]/(100.0_realwp*bohr), &
        rvis=0.495_realwp/bohr, &
        rcovG=3.288_realwp, &
        color=[61, 255, 0] &
    ) /
    data atdata( 21) / AtomDB( &
        symbol='Sc', name='Scandium', &
        number=21, &
        mass=44.95590_realwp, &
        rcov=[148, 116, 114]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 258, 263]/(100.0_realwp*bohr), &
        rvis=0.735_realwp/bohr, &
        rcovG=2.721_realwp, &
        color=[230, 230, 230] &
    ) /
    data atdata( 22) / AtomDB( &
        symbol='Ti', name='Titanium', &
        number=22, &
        mass=47.90000_realwp, &
        rcov=[136, 117, 108]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 246, 257]/(100.0_realwp*bohr), &
        rvis=0.720_realwp/bohr, &
        rcovG=2.494_realwp, &
        color=[191, 194, 199] &
    ) /
    data atdata( 23) / AtomDB( &
        symbol='V', name='Vanadium', &
        number=23, &
        mass=50.94150_realwp, &
        rcov=[134, 112, 106]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 242, 252]/(100.0_realwp*bohr), &
        rvis=0.665_realwp/bohr, &
        rcovG=2.305_realwp, &
        color=[166, 166, 171] &
    ) /
    data atdata( 24) / AtomDB( &
        symbol='Cr', name='Chromium', &
        number=24, &
        mass=51.99600_realwp, &
        rcov=[122, 111, 103]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 245, 233]/(100.0_realwp*bohr), &
        rvis=0.675_realwp/bohr, &
        rcovG=2.230_realwp, &
        color=[138, 153, 199] &
    ) /
    data atdata( 25) / AtomDB( &
        symbol='Mn', name='Manganese', &
        number=25, &
        mass=54.93800_realwp, &
        rcov=[119, 105, 103]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 245, 242]/(100.0_realwp*bohr), &
        rvis=0.675_realwp/bohr, &
        rcovG=2.211_realwp, &
        color=[156, 122, 199] &
    ) /
    data atdata( 26) / AtomDB( &
        symbol='Fe', name='Iron', &
        number=26, &
        mass=55.84700_realwp, &
        rcov=[116, 109, 102]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 244, 237]/(100.0_realwp*bohr), &
        rvis=0.670_realwp/bohr, &
        rcovG=2.211_realwp, &
        color=[224, 102, 51] &
    ) /
    data atdata( 27) / AtomDB( &
        symbol='Co', name='Cobalt', &
        number=27, &
        mass=58.93320_realwp, &
        rcov=[111, 103, 96]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 240, 233]/(100.0_realwp*bohr), &
        rvis=0.615_realwp/bohr, &
        rcovG=2.192_realwp, &
        color=[240, 144, 160] &
    ) /
    data atdata( 28) / AtomDB( &
        symbol='Ni', name='Nickel', &
        number=28, &
        mass=58.71000_realwp, &
        rcov=[110, 101, 101]/(100.0_realwp*bohr), &
        rvdw=[163, 0, 240, 229]/(100.0_realwp*bohr), &
        rvis=0.750_realwp/bohr, &
        rcovG=2.173_realwp, &
        color=[80, 208, 80] &
    ) /
    data atdata( 29) / AtomDB( &
        symbol='Cu', name='Copper', &
        number=29, &
        mass=63.54600_realwp, &
        rcov=[112, 115, 120]/(100.0_realwp*bohr), &
        rvdw=[140, 0, 238, 217]/(100.0_realwp*bohr), &
        rvis=0.760_realwp/bohr, &
        rcovG=2.211_realwp, &
        color=[200, 128, 51] &
    ) /
    data atdata( 30) / AtomDB( &
        symbol='Zn', name='Zinc', &
        number=30, &
        mass=65.38000_realwp, &
        rcov=[118, 120, 0]/(100.0_realwp*bohr), &
        rvdw=[139, 0, 239, 222]/(100.0_realwp*bohr), &
        rvis=0.725_realwp/bohr, &
        rcovG=2.362_realwp, &
        color=[125, 128, 176] &
    ) /
    data atdata( 31) / AtomDB( &
        symbol='Ga', name='Gallium', &
        number=31, &
        mass=69.73500_realwp, &
        rcov=[124, 117, 121]/(100.0_realwp*bohr), &
        rvdw=[187, 187, 232, 233]/(100.0_realwp*bohr), &
        rvis=0.610_realwp/bohr, &
        rcovG=2.381_realwp, &
        color=[194, 143, 143] &
    ) /
    data atdata( 32) / AtomDB( &
        symbol='Ge', name='Germanium', &
        number=32, &
        mass=72.59000_realwp, &
        rcov=[121, 111, 114]/(100.0_realwp*bohr), &
        rvdw=[0, 211, 229, 234]/(100.0_realwp*bohr), &
        rvis=0.585_realwp/bohr, &
        rcovG=2.305_realwp, &
        color=[102, 143, 143] &
    ) /
    data atdata( 33) / AtomDB( &
        symbol='As', name='Arsenic', &
        number=33, &
        mass=74.92160_realwp, &
        rcov=[121, 114, 106]/(100.0_realwp*bohr), &
        rvdw=[185, 185, 188, 231]/(100.0_realwp*bohr), &
        rvis=0.605_realwp/bohr, &
        rcovG=2.268_realwp, &
        color=[189, 128, 227] &
    ) /
    data atdata( 34) / AtomDB( &
        symbol='Se', name='Selenium', &
        number=34, &
        mass=78.96000_realwp, &
        rcov=[116, 107, 107]/(100.0_realwp*bohr), &
        rvdw=[190, 190, 182, 224]/(100.0_realwp*bohr), &
        rvis=0.610_realwp/bohr, &
        rcovG=2.192_realwp, &
        color=[255, 161, 0] &
    ) /
    data atdata( 35) / AtomDB( &
        symbol='Br', name='Bromine', &
        number=35, &
        mass=79.90400_realwp, &
        rcov=[114, 109, 110]/(100.0_realwp*bohr), &
        rvdw=[183, 183, 186, 219]/(100.0_realwp*bohr), &
        rvis=0.605_realwp/bohr, &
        rcovG=2.154_realwp, &
        color=[166,  41,  41] &
    ) /
    data atdata( 36) / AtomDB( &
        symbol='Kr', name='Krypton', &
        number=36, &
        mass=83.80000_realwp, &
        rcov=[117, 121, 108]/(100.0_realwp*bohr), &
        rvdw=[202, 202, 225, 212]/(100.0_realwp*bohr), &
        rvis=0.524_realwp/bohr, &
        rcovG=2.116_realwp, &
        color=[92, 184, 209] &
    ) /
    data atdata( 37) / AtomDB( &
        symbol='Rb', name='Rubidium', &
        number=37, &
        mass=85.46780_realwp, &
        rcov=[210, 202, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 303, 321, 240]/(100.0_realwp*bohr), &
        rvis=0.735_realwp/bohr, &
        rcovG=4.082_realwp, &
        color=[112, 46, 176] &
    ) /
    data atdata( 38) / AtomDB( &
        symbol='Sr', name='Strontium', &
        number=38, &
        mass=87.62000_realwp, &
        rcov=[185, 157, 139]/(100.0_realwp*bohr), &
        rvdw=[0, 249, 284, 279]/(100.0_realwp*bohr), &
        rvis=0.560_realwp/bohr, &
        rcovG=3.609_realwp, &
        color=[0, 255, 0] &
    ) /
    data atdata( 39) / AtomDB( &
        symbol='Y', name='Yttrium', &
        number=39, &
        mass=88.90590_realwp, &
        rcov=[163, 130, 124]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 275, 274]/(100.0_realwp*bohr), &
        rvis=0.890_realwp/bohr, &
        rcovG=3.061_realwp, &
        color=[148, 255, 255] &
    ) /
    data atdata( 40) / AtomDB( &
        symbol='Zr', name='Zirconium', &
        number=40, &
        mass=91.22000_realwp, &
        rcov=[154, 127, 121]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 252, 269]/(100.0_realwp*bohr), &
        rvis=0.780_realwp/bohr, &
        rcovG=2.740_realwp, &
        color=[148, 224, 224] &
    ) /
    data atdata( 41) / AtomDB( &
        symbol='Nb', name='Niobium', &
        number=41, &
        mass=92.90640_realwp, &
        rcov=[147, 125, 116]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 256, 251]/(100.0_realwp*bohr), &
        rvis=0.740_realwp/bohr, &
        rcovG=2.532_realwp, &
        color=[115, 194, 201] &
    ) /
    data atdata( 42) / AtomDB( &
        symbol='Mo', name='Molybdenum', &
        number=42, &
        mass=95.94000_realwp, &
        rcov=[138, 121, 113]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 245, 244]/(100.0_realwp*bohr), &
        rvis=0.735_realwp/bohr, &
        rcovG=2.457_realwp, &
        color=[84, 181, 181] &
    ) /
    data atdata( 43) / AtomDB( &
        symbol='Tc', name='Technetium', &
        number=43, &
        mass=98.90620_realwp, &
        rcov=[128, 120, 110]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 244, 252]/(100.0_realwp*bohr), &
        rvis=0.675_realwp/bohr, &
        rcovG=2.400_realwp, &
        color=[59, 158, 158] &
    ) /
    data atdata( 44) / AtomDB( &
        symbol='Ru', name='Ruthenium', &
        number=44, &
        mass=101.0700_realwp, &
        rcov=[125, 114, 103]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 246, 237]/(100.0_realwp*bohr), &
        rvis=0.700_realwp/bohr, &
        rcovG=2.362_realwp, &
        color=[36, 143, 143] &
    ) /
    data atdata( 45) / AtomDB( &
        symbol='Rh', name='Rhodium', &
        number=45, &
        mass=102.9055_realwp, &
        rcov=[125, 110, 106]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 244, 233]/(100.0_realwp*bohr), &
        rvis=0.725_realwp/bohr, &
        rcovG=2.362_realwp, &
        color=[10, 125, 140] &
    ) /
    data atdata( 46) / AtomDB( &
        symbol='Pd', name='Palladium', &
        number=46, &
        mass=106.4000_realwp, &
        rcov=[120, 117, 112]/(100.0_realwp*bohr), &
        rvdw=[163, 0, 215, 215]/(100.0_realwp*bohr), &
        rvis=0.750_realwp/bohr, &
        rcovG=2.419_realwp, &
        color=[0, 105, 133] &
    ) /
    data atdata( 47) / AtomDB( &
        symbol='Ag', name='Silver', &
        number=47, &
        mass=107.8680_realwp, &
        rcov=[128, 139, 137]/(100.0_realwp*bohr), &
        rvdw=[172, 0, 253, 225]/(100.0_realwp*bohr), &
        rvis=0.795_realwp/bohr, &
        rcovG=2.532_realwp, &
        color=[192, 192, 192] &
    ) /
    data atdata( 48) / AtomDB( &
        symbol='Cd', name='Cadmium', &
        number=48, &
        mass=112.4100_realwp, &
        rcov=[136, 144, 0]/(100.0_realwp*bohr), &
        rvdw=[158, 0, 249, 238]/(100.0_realwp*bohr), &
        rvis=0.845_realwp/bohr, &
        rcovG=2.797_realwp, &
        color=[255, 217, 143] &
    ) /
    data atdata( 49) / AtomDB( &
        symbol='In', name='Indium', &
        number=49, &
        mass=114.8200_realwp, &
        rcov=[142, 136, 146]/(100.0_realwp*bohr), &
        rvdw=[193, 193, 243, 246]/(100.0_realwp*bohr), &
        rvis=0.815_realwp/bohr, &
        rcovG=2.721_realwp, &
        color=[166, 117, 115] &
    ) /
    data atdata( 50) / AtomDB( &
        symbol='Sn', name='Tin', &
        number=50, &
        mass=118.6900_realwp, &
        rcov=[140, 130, 132]/(100.0_realwp*bohr), &
        rvdw=[217, 217, 242, 248]/(100.0_realwp*bohr), &
        rvis=0.730_realwp/bohr, &
        rcovG=2.665_realwp, &
        color=[102, 128, 128] &
    ) /
    data atdata( 51) / AtomDB( &
        symbol='Sb', name='Antimony', &
        number=51, &
        mass=121.7500_realwp, &
        rcov=[140, 133, 127]/(100.0_realwp*bohr), &
        rvdw=[0, 206, 247, 246]/(100.0_realwp*bohr), &
        rvis=0.730_realwp/bohr, &
        rcovG=2.646_realwp, &
        color=[158, 99, 181] &
    ) /
    data atdata( 52) / AtomDB( &
        symbol='Te', name='Tellurium', &
        number=52, &
        mass=127.6000_realwp, &
        rcov=[136, 128, 121]/(100.0_realwp*bohr), &
        rvdw=[206, 206, 199, 242]/(100.0_realwp*bohr), &
        rvis=0.735_realwp/bohr, &
        rcovG=2.570_realwp, &
        color=[212, 122, 0] &
    ) /
    data atdata( 53) / AtomDB( &
        symbol='I', name='Iodine', &
        number=53, &
        mass=126.9045_realwp, &
        rcov=[133, 129, 125]/(100.0_realwp*bohr), &
        rvdw=[198, 198, 204, 238]/(100.0_realwp*bohr), &
        rvis=0.700_realwp/bohr, &
        rcovG=2.513_realwp, &
        color=[148, 0, 148] &
    ) /
    data atdata( 54) / AtomDB( &
        symbol='Xe', name='Xenon', &
        number=54, &
        mass=131.3000_realwp, &
        rcov=[131, 135, 122]/(100.0_realwp*bohr), &
        rvdw=[216, 216, 206, 232]/(100.0_realwp*bohr), &
        rvis=0.577_realwp/bohr, &
        rcovG=2.476_realwp, &
        color=[66, 158, 176] &
    ) /
    data atdata( 55) / AtomDB( &
        symbol='Cs', name='Caesium', &
        number=55, &
        mass=132.9054_realwp, &
        rcov=[232, 209, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 343, 348, 249]/(100.0_realwp*bohr), &
        rvis=0.835_realwp/bohr, &
        rcovG=4.441_realwp, &
        color=[87, 23, 143] &
    ) /
    data atdata( 56) / AtomDB( &
        symbol='Ba', name='Barium', &
        number=56, &
        mass=137.3300_realwp, &
        rcov=[196, 161, 149]/(100.0_realwp*bohr), &
        rvdw=[0, 268, 303, 293]/(100.0_realwp*bohr), &
        rvis=0.670_realwp/bohr, &
        rcovG=3.742_realwp, &
        color=[0, 201, 0] &
    ) /
    data atdata( 57) / AtomDB( &
        symbol='La', name='Lanthanum', &
        number=57, &
        mass=138.9055_realwp, &
        rcov=[180, 139, 139]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 298, 284]/(100.0_realwp*bohr), &
        rvis=0.935_realwp/bohr, &
        rcovG=3.194_realwp, &
        color=[112, 212, 255] &
    ) /
    data atdata( 58) / AtomDB( &
        symbol='Ce', name='Cerium', &
        number=58, &
        mass=140.1200_realwp, &
        rcov=[163, 137, 131]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 288, 282]/(100.0_realwp*bohr), &
        rvis=0.915_realwp/bohr, &
        rcovG=3.118_realwp, &
        color=[255, 255, 199] &
    ) /
    data atdata( 59) / AtomDB( &
        symbol='Pr', name='Praseodymium', &
        number=59, &
        mass=140.9077_realwp, &
        rcov=[176, 138, 128]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 292, 286]/(100.0_realwp*bohr), &
        rvis=0.910_realwp/bohr, &
        rcovG=3.118_realwp, &
        color=[217, 255, 199] &
    ) /
    data atdata( 60) / AtomDB( &
        symbol='Nd', name='Neodymium', &
        number=60, &
        mass=144.2400_realwp, &
        rcov=[174, 132, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 295, 284]/(100.0_realwp*bohr), &
        rvis=0.905_realwp/bohr, &
        rcovG=3.099_realwp, &
        color=[199, 255, 199] &
    ) /
    data atdata( 61) / AtomDB( &
        symbol='Pm', name='Promethium', &
        number=61, &
        mass=0.0_realwp, &
        rcov=[173, 135, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 290, 283]/(100.0_realwp*bohr), &
        rvis=0.900_realwp/bohr, &
        rcovG=3.080_realwp, &
        color=[163, 255, 199] &
    ) /
    data atdata( 62) / AtomDB( &
        symbol='Sm', name='Samarium', &
        number=62, &
        mass=150.3600_realwp, &
        rcov=[172, 134, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 290, 280]/(100.0_realwp*bohr), &
        rvis=0.900_realwp/bohr, &
        rcovG=3.061_realwp, &
        color=[143, 255, 199] &
    ) /
    data atdata( 63) / AtomDB( &
        symbol='Eu', name='Europium', &
        number=63, &
        mass=151.9600_realwp, &
        rcov=[168, 134, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 287, 280]/(100.0_realwp*bohr), &
        rvis=0.995_realwp/bohr, &
        rcovG=3.496_realwp, &
        color=[97, 255, 199] &
    ) /
    data atdata( 64) / AtomDB( &
        symbol='Gd', name='Gadolinium', &
        number=64, &
        mass=157.2500_realwp, &
        rcov=[169, 135, 132]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 283, 277]/(100.0_realwp*bohr), &
        rvis=0.895_realwp/bohr, &
        rcovG=3.042_realwp, &
        color=[69, 255, 199] &
    ) /
    data atdata( 65) / AtomDB( &
        symbol='Tb', name='Terbium', &
        number=65, &
        mass=158.9253_realwp, &
        rcov=[168, 135, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 279, 276]/(100.0_realwp*bohr), &
        rvis=0.880_realwp/bohr, &
        rcovG=3.005_realwp, &
        color=[48, 255, 199] &
    ) /
    data atdata( 66) / AtomDB( &
        symbol='Dy', name='Dysprosium', &
        number=66, &
        mass=162.5000_realwp, &
        rcov=[167, 133, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 287, 275]/(100.0_realwp*bohr), &
        rvis=0.875_realwp/bohr, &
        rcovG=3.005_realwp, &
        color=[31, 255, 199] &
    ) /
    data atdata( 67) / AtomDB( &
        symbol='Ho', name='Holmium', &
        number=67, &
        mass=164.9303_realwp, &
        rcov=[166, 133, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 281, 273]/(100.0_realwp*bohr), &
        rvis=0.870_realwp/bohr, &
        rcovG=2.986_realwp, &
        color=[0, 255, 156] &
    ) /
    data atdata( 68) / AtomDB( &
        symbol='Er', name='Erbium', &
        number=68, &
        mass=167.2600_realwp, &
        rcov=[165, 133, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 283, 272]/(100.0_realwp*bohr), &
        rvis=0.865_realwp/bohr, &
        rcovG=2.967_realwp, &
        color=[0, 230, 117] &
    ) /
    data atdata( 69) / AtomDB( &
        symbol='Tm', name='Thulium', &
        number=69, &
        mass=168.9342_realwp, &
        rcov=[164, 131, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 279, 271]/(100.0_realwp*bohr), &
        rvis=0.860_realwp/bohr, &
        rcovG=2.948_realwp, &
        color=[0, 212, 82] &
    ) /
    data atdata( 70) / AtomDB( &
        symbol='Yb', name='Ytterbium', &
        number=70, &
        mass=173.0500_realwp, &
        rcov=[170, 129, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 280, 277]/(100.0_realwp*bohr), &
        rvis=0.970_realwp/bohr, &
        rcovG=2.948_realwp, &
        color=[0, 191, 56] &
    ) /
    data atdata( 71) / AtomDB( &
        symbol='Lu', name='Lutetium', &
        number=71, &
        mass=174.9670_realwp, &
        rcov=[162, 131, 131]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 274, 270]/(100.0_realwp*bohr), &
        rvis=0.860_realwp/bohr, &
        rcovG=2.948_realwp, &
        color=[0, 171, 36] &
    ) /
    data atdata( 72) / AtomDB( &
        symbol='Hf', name='Hafnium', &
        number=72, &
        mass=178.4900_realwp, &
        rcov=[152, 128, 122]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 263, 264]/(100.0_realwp*bohr), &
        rvis=0.785_realwp/bohr, &
        rcovG=2.721_realwp, &
        color=[77, 194, 255] &
    ) /
    data atdata( 73) / AtomDB( &
        symbol='Ta', name='Tantalum', &
        number=73, &
        mass=180.9479_realwp, &
        rcov=[146, 126, 119]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 253, 258]/(100.0_realwp*bohr), &
        rvis=0.715_realwp/bohr, &
        rcovG=2.532_realwp, &
        color=[77, 166, 255] &
    ) /
    data atdata( 74) / AtomDB( &
        symbol='W', name='Tungsten', &
        number=74, &
        mass=183.8400_realwp, &
        rcov=[137, 120, 115]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 257, 253]/(100.0_realwp*bohr), &
        rvis=0.685_realwp/bohr, &
        rcovG=2.457_realwp, &
        color=[33, 148, 214] &
    ) /
    data atdata( 75) / AtomDB( &
        symbol='Re', name='Rhenium', &
        number=75, &
        mass=186.2070_realwp, &
        rcov=[131, 119, 110]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 249, 249]/(100.0_realwp*bohr), &
        rvis=0.675_realwp/bohr, &
        rcovG=2.419_realwp, &
        color=[38, 125, 171] &
    ) /
    data atdata( 76) / AtomDB( &
        symbol='Os', name='Osmium', &
        number=76, &
        mass=190.2000_realwp, &
        rcov=[129, 116, 109]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 248, 244]/(100.0_realwp*bohr), &
        rvis=0.685_realwp/bohr, &
        rcovG=2.381_realwp, &
        color=[38, 102, 150] &
    ) /
    data atdata( 77) / AtomDB( &
        symbol='Ir', name='Iridium', &
        number=77, &
        mass=192.2200_realwp, &
        rcov=[122, 115, 107]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 241, 240]/(100.0_realwp*bohr), &
        rvis=0.660_realwp/bohr, &
        rcovG=2.400_realwp, &
        color=[23, 84, 135] &
    ) /
    data atdata( 78) / AtomDB( &
        symbol='Pt', name='Platinum', &
        number=78, &
        mass=195.0900_realwp, &
        rcov=[123, 112, 110]/(100.0_realwp*bohr), &
        rvdw=[172, 0, 229, 230]/(100.0_realwp*bohr), &
        rvis=0.750_realwp/bohr, &
        rcovG=2.457_realwp, &
        color=[208, 208, 224] &
    ) /
    data atdata( 79) / AtomDB( &
        symbol='Au', name='Gold', &
        number=79, &
        mass=196.9665_realwp, &
        rcov=[124, 121, 123]/(100.0_realwp*bohr), &
        rvdw=[166, 0, 232, 226]/(100.0_realwp*bohr), &
        rvis=0.750_realwp/bohr, &
        rcovG=2.532_realwp, &
        color=[255, 209, 35] &
    ) /
    data atdata( 80) / AtomDB( &
        symbol='Hg', name='Mercury', &
        number=80, &
        mass=200.5900_realwp, &
        rcov=[133, 142, 0]/(100.0_realwp*bohr), &
        rvdw=[170, 0, 245, 229]/(100.0_realwp*bohr), &
        rvis=0.850_realwp/bohr, &
        rcovG=2.816_realwp, &
        color=[184, 184, 208] &
    ) /
    data atdata( 81) / AtomDB( &
        symbol='Tl', name='Thallium', &
        number=81, &
        mass=204.3700_realwp, &
        rcov=[144, 142, 150]/(100.0_realwp*bohr), &
        rvdw=[196, 196, 247, 242]/(100.0_realwp*bohr), &
        rvis=0.775_realwp/bohr, &
        rcovG=2.797_realwp, &
        color=[166, 84, 77] &
    ) /
    data atdata( 82) / AtomDB( &
        symbol='Pb', name='Lead', &
        number=82, &
        mass=207.2000_realwp, &
        rcov=[144, 135, 137]/(100.0_realwp*bohr), &
        rvdw=[202, 202, 260, 249]/(100.0_realwp*bohr), &
        rvis=0.770_realwp/bohr, &
        rcovG=2.778_realwp, &
        color=[87, 89, 97] &
    ) /
    data atdata( 83) / AtomDB( &
        symbol='Bi', name='Bismuth', &
        number=83, &
        mass=208.9804_realwp, &
        rcov=[151, 141, 135]/(100.0_realwp*bohr), &
        rvdw=[0, 207, 254, 250]/(100.0_realwp*bohr), &
        rvis=0.770_realwp/bohr, &
        rcovG=2.759_realwp, &
        color=[158, 79, 181] &
    ) /
    data atdata( 84) / AtomDB( &
        symbol='Po', name='Polonium', &
        number=84, &
        mass=0.0_realwp, &
        rcov=[145, 135, 129]/(100.0_realwp*bohr), &
        rvdw=[0, 197, 0, 250]/(100.0_realwp*bohr), &
        rvis=0.840_realwp/bohr, &
        rcovG=2.759_realwp, &
        color=[171, 92, 0] &
    ) /
    data atdata( 85) / AtomDB( &
        symbol='At', name='Astatine', &
        number=85, &
        mass=0.0_realwp, &
        rcov=[147, 138, 138]/(100.0_realwp*bohr), &
        rvdw=[0, 202, 0, 247]/(100.0_realwp*bohr), &
        rvis=1.000_realwp/bohr, &
        rcovG=2.740_realwp, &
        color=[117, 79, 69] &
    ) /
    data atdata( 86) / AtomDB( &
        symbol='Rn', name='Radon', &
        number=86, &
        mass=0.0_realwp, &
        rcov=[142, 145, 133]/(100.0_realwp*bohr), &
        rvdw=[0, 220, 0, 243]/(100.0_realwp*bohr), &
        rvis=1.000_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[66, 130, 150] &
    ) /
    data atdata( 87) / AtomDB( &
        symbol='Fr', name='Francium', &
        number=87, &
        mass=0.0_realwp, &
        rcov=[223, 218, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 348, 0, 258]/(100.0_realwp*bohr), &
        rvis=1.000_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[66, 0, 102] &
    ) /
    data atdata( 88) / AtomDB( &
        symbol='Ra', name='Radium', &
        number=88, &
        mass=0.0_realwp, &
        rcov=[201, 173, 159]/(100.0_realwp*bohr), &
        rvdw=[0, 283, 0, 292]/(100.0_realwp*bohr), &
        rvis=0.950_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 125, 0] &
    ) /
    data atdata( 89) / AtomDB( &
        symbol='Ac', name='Actinium', &
        number=89, &
        mass=0.0_realwp, &
        rcov=[186, 153, 140]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 280, 293]/(100.0_realwp*bohr), &
        rvis=0.940_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[112, 171, 250] &
    ) /
    data atdata( 90) / AtomDB( &
        symbol='Th', name='Thorium', &
        number=90, &
        mass=232.038_realwp, &
        rcov=[175, 143, 136]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 293, 288]/(100.0_realwp*bohr), &
        rvis=0.895_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 186, 255] &
    ) /
    data atdata( 91) / AtomDB( &
        symbol='Pa', name='Protactinium', &
        number=91, &
        mass=231.0359_realwp, &
        rcov=[169, 138, 129]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 288, 285]/(100.0_realwp*bohr), &
        rvis=0.805_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 161, 255] &
    ) /
    data atdata( 92) / AtomDB( &
        symbol='U', name='Uranium', &
        number=92, &
        mass=238.0289_realwp, &
        rcov=[170, 134, 118]/(100.0_realwp*bohr), &
        rvdw=[186, 0, 271, 283]/(100.0_realwp*bohr), &
        rvis=0.790_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 143, 255] &
    ) /
    data atdata( 93) / AtomDB( &
        symbol='Np', name='Neptunium', &
        number=93, &
        mass=0.0_realwp, &
        rcov=[171, 133, 116]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 282, 281]/(100.0_realwp*bohr), &
        rvis=0.775_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 128, 255] &
    ) /
    data atdata( 94) / AtomDB( &
        symbol='Pu', name='Plutonium', &
        number=94, &
        mass=0.0_realwp, &
        rcov=[172, 135, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 281, 278]/(100.0_realwp*bohr), &
        rvis=1.000_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 107, 255] &
    ) /
    data atdata( 95) / AtomDB( &
        symbol='Am', name='Americium', &
        number=95, &
        mass=0.0_realwp, &
        rcov=[166, 135, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 283, 276]/(100.0_realwp*bohr), &
        rvis=0.755_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[84, 92, 242] &
    ) /
    data atdata( 96) / AtomDB( &
        symbol='Cm', name='Curium', &
        number=96, &
        mass=0.0_realwp, &
        rcov=[166, 136, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 305, 264]/(100.0_realwp*bohr), &
        rvis=1.000_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[120, 92, 227] &
    ) /
    data atdata( 97) / AtomDB( &
        symbol='Bk', name='Berkelium', &
        number=97, &
        mass=0.0_realwp, &
        rcov=[168, 139, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 340, 0]/(100.0_realwp*bohr), &
        rvis=1.000_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[138, 79, 227] &
    ) /
    data atdata( 98) / AtomDB( &
        symbol='Cf', name='Californium', &
        number=98, &
        mass=0.0_realwp, &
        rcov=[168, 140, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 305, 0]/(100.0_realwp*bohr), &
        rvis=0.765_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[161, 54, 212] &
    ) /
    data atdata( 99) / AtomDB( &
        symbol='Es', name='Einsteinium', &
        number=99, &
        mass=0.0_realwp, &
        rcov=[165, 140, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 270, 0]/(100.0_realwp*bohr), &
        rvis=0.100_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[179, 31, 212] &
    ) /
    data atdata(100) / AtomDB( &
        symbol='Fm', name='Fermium', &
        number=100, &
        mass=0.0_realwp, &
        rcov=[167, 0, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.900_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[179, 31, 186] &
    ) /
    data atdata(101) / AtomDB( &
        symbol='Md', name='Mendelevium', &
        number=101, &
        mass=0.0_realwp, &
        rcov=[173, 139, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[179, 13, 166] &
    ) /
    data atdata(102) / AtomDB( &
        symbol='No', name='Nobelium', &
        number=102, &
        mass=0.0_realwp, &
        rcov=[176, 0, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[189, 13, 135] &
    ) /
    data atdata(103) / AtomDB( &
        symbol='Lr', name='Lawrencium', &
        number=103, &
        mass=0.0_realwp, &
        rcov=[161, 141, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[199, 0, 102] &
    ) /
    data atdata(104) / AtomDB( &
        symbol='Rf', name='Rutherfordium', &
        number=104, &
        mass=0.0_realwp, &
        rcov=[157, 140, 131]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[204, 0, 89] &
    ) /
    data atdata(105) / AtomDB( &
        symbol='Db', name='Dubnium', &
        number=105, &
        mass=0.0_realwp, &
        rcov=[149, 136, 126]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[209, 0, 79] &
    ) /
    data atdata(106) / AtomDB( &
        symbol='Sg', name='Seaborgium', &
        number=106, &
        mass=0.0_realwp, &
        rcov=[143, 128, 121]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[217, 0, 69] &
    ) /
    data atdata(107) / AtomDB( &
        symbol='Bh', name='Bohrium', &
        number=107, &
        mass=0.0_realwp, &
        rcov=[141, 128, 119]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[224, 0, 56] &
    ) /
    data atdata(108) / AtomDB( &
        symbol='Hs', name='Hassium', &
        number=108, &
        mass=0.0_realwp, &
        rcov=[134, 125, 118]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[230, 0, 46] &
    ) /
    data atdata(109) / AtomDB( &
        symbol='Mt', name='Meitnerium', &
        number=109, &
        mass=0.0_realwp, &
        rcov=[129, 125, 118]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[235, 0, 38] &
    ) /
    data atdata(110) / AtomDB( &
        symbol='Ds', name='Darmstadtium', &
        number=110, &
        mass=0.0_realwp, &
        rcov=[128, 116, 112]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 0, 0] &
    ) /
    data atdata(111) / AtomDB( &
        symbol='Rg', name='Roentgenium', &
        number=111, &
        mass=0.0_realwp, &
        rcov=[121, 116, 118]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 0, 0] &
    ) /
    data atdata(112) / AtomDB( &
        symbol='Cn', name='Copernicium', &
        number=112, &
        mass=0.0_realwp, &
        rcov=[122, 116, 118]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 0, 0] &
    ) /
    data atdata(113) / AtomDB( &
        symbol='Nh', name='Nihonium', &
        number=113, &
        mass=0.0_realwp, &
        rcov=[136, 0, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 0, 0] &
    ) /
    data atdata(114) / AtomDB( &
        symbol='Fl', name='Flerovium', &
        number=114, &
        mass=0.0_realwp, &
        rcov=[143, 0, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 0, 0] &
    ) /
    data atdata(115) / AtomDB( &
        symbol='Mc', name='Moscovium', &
        number=115, &
        mass=0.0_realwp, &
        rcov=[162, 0, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 0, 0] &
    ) /
    data atdata(116) / AtomDB( &
        symbol='Lv', name='Livermorium', &
        number=116, &
        mass=0.0_realwp, &
        rcov=[175, 0, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 0, 0] &
    ) /
    data atdata(117) / AtomDB( &
        symbol='Ts', name='Tennessine', &
        number=117, &
        mass=0.0_realwp, &
        rcov=[165, 0, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 0, 0] &
    ) /
    data atdata(118) / AtomDB( &
        symbol='Og', name='Oganesson', &
        number=118, &
        mass=0.0_realwp, &
        rcov=[157, 0, 0]/(100.0_realwp*bohr), &
        rvdw=[0, 0, 0, 0]/(100.0_realwp*bohr), &
        rvis=0.0_realwp/bohr, &
        rcovG=0.0_realwp, &
        color=[0, 0, 0] &
    ) /

end module atominfo
