module atomic
    use iso_fortran_env, only: real64

    type, private :: atom
        character(len=2) :: symbol
        real(real64) :: mass
    end type atom

    type(atom), dimension(110) :: atdata

    data atdata(  1)/atom('H',    1.00790_real64 )/
    data atdata(  2)/atom('He',   4.00260_real64 )/
    data atdata(  3)/atom('Li',   6.94000_real64 )/
    data atdata(  4)/atom('Be',   9.01218_real64 )/
    data atdata(  5)/atom('B',   10.81000_real64 )/
    data atdata(  6)/atom('C',   12.01100_real64 )/
    data atdata(  7)/atom('N',   14.00670_real64 )/
    data atdata(  8)/atom('O',   15.99940_real64 )/
    data atdata(  9)/atom('F',   18.99840_real64 )/
    data atdata( 10)/atom('Ne',  20.17900_real64 )/
    data atdata( 11)/atom('Na',  22.98977_real64 )/
    data atdata( 12)/atom('Mg',  24.30500_real64 )/
    data atdata( 13)/atom('Al',  26.98154_real64 )/
    data atdata( 14)/atom('Si',  28.08550_real64 )/
    data atdata( 15)/atom('P',   30.97376_real64 )/
    data atdata( 16)/atom('S',   32.06000_real64 )/
    data atdata( 17)/atom('Cl',  35.45300_real64 )/
    data atdata( 18)/atom('Ar',  39.94800_real64 )/
    data atdata( 19)/atom('K',   39.09830_real64 )/
    data atdata( 20)/atom('Ca',  40.08000_real64 )/
    data atdata( 21)/atom('Sc',  44.95590_real64 )/
    data atdata( 22)/atom('Ti',  47.90000_real64 )/
    data atdata( 23)/atom('V',   50.94150_real64 )/
    data atdata( 24)/atom('Cr',  51.99600_real64 )/
    data atdata( 25)/atom('Mn',  54.93800_real64 )/
    data atdata( 26)/atom('Fe',  55.84700_real64 )/
    data atdata( 27)/atom('Co',  58.93320_real64 )/
    data atdata( 28)/atom('Ni',  58.71000_real64 )/
    data atdata( 29)/atom('Cu',  63.54600_real64 )/
    data atdata( 30)/atom('Zn',  65.38000_real64 )/
    data atdata( 31)/atom('Ga',  69.73500_real64 )/
    data atdata( 32)/atom('Ge',  72.59000_real64 )/
    data atdata( 33)/atom('As',  74.92160_real64 )/
    data atdata( 34)/atom('Se',  78.96000_real64 )/
    data atdata( 35)/atom('Br',  79.90400_real64 )/
    data atdata( 36)/atom('Kr',  83.80000_real64 )/
    data atdata( 37)/atom('Rb',  85.46780_real64 )/
    data atdata( 38)/atom('Sr',  87.62000_real64 )/
    data atdata( 39)/atom('Y',   88.90590_real64 )/
    data atdata( 40)/atom('Zr',  91.22000_real64 )/
    data atdata( 41)/atom('Nb',  92.90640_real64 )/
    data atdata( 42)/atom('Mo',  95.94000_real64 )/
    data atdata( 43)/atom('Tc',  98.90620_real64 )/
    data atdata( 44)/atom('Ru', 101.07000_real64 )/
    data atdata( 45)/atom('Rh', 102.90550_real64 )/
    data atdata( 46)/atom('Pd', 106.40000_real64 )/
    data atdata( 47)/atom('Ag', 107.86800_real64 )/
    data atdata( 48)/atom('Cd', 112.41000_real64 )/
    data atdata( 49)/atom('In', 114.82000_real64 )/
    data atdata( 50)/atom('Sn', 118.69000_real64 )/
    data atdata( 51)/atom('Sb', 121.75000_real64 )/
    data atdata( 52)/atom('Te', 127.60000_real64 )/
    data atdata( 53)/atom('I',  126.90450_real64 )/
    data atdata( 54)/atom('Xe', 131.30000_real64 )/
    data atdata( 55)/atom('Cs', 132.90540_real64 )/
    data atdata( 56)/atom('Ba', 137.33000_real64 )/
    data atdata( 57)/atom('La',   0.00000_real64 )/
    data atdata( 58)/atom('Ce',   0.00000_real64 )/
    data atdata( 59)/atom('Pr',   0.00000_real64 )/
    data atdata( 60)/atom('Nd',   0.00000_real64 )/
    data atdata( 61)/atom('Pm',   0.00000_real64 )/
    data atdata( 62)/atom('Sm',   0.00000_real64 )/
    data atdata( 63)/atom('Eu',   0.00000_real64 )/
    data atdata( 64)/atom('Gd',   0.00000_real64 )/
    data atdata( 65)/atom('Tb',   0.00000_real64 )/
    data atdata( 66)/atom('Dy',   0.00000_real64 )/
    data atdata( 67)/atom('Ho',   0.00000_real64 )/
    data atdata( 68)/atom('Er',   0.00000_real64 )/
    data atdata( 69)/atom('Tm',   0.00000_real64 )/
    data atdata( 70)/atom('Yb',   0.00000_real64 )/
    data atdata( 71)/atom('Lu',   0.00000_real64 )/
    data atdata( 72)/atom('Hf', 178.49000_real64 )/
    data atdata( 73)/atom('Ta', 180.94790_real64 )/
    data atdata( 74)/atom('W',  183.85000_real64 )/
    data atdata( 75)/atom('Re', 186.20700_real64 )/
    data atdata( 76)/atom('Os', 190.20000_real64 )/
    data atdata( 77)/atom('Ir', 192.22000_real64 )/
    data atdata( 78)/atom('Pt', 195.09000_real64 )/
    data atdata( 79)/atom('Au', 196.96650_real64 )/
    data atdata( 80)/atom('Hg', 200.59000_real64 )/
    data atdata( 81)/atom('Tl', 204.37000_real64 )/
    data atdata( 82)/atom('Pb', 207.20000_real64 )/
    data atdata( 83)/atom('Bi', 208.98040_real64 )/
    data atdata( 84)/atom('Po',   0.00000_real64 )/
    data atdata( 85)/atom('At',   0.00000_real64 )/
    data atdata( 86)/atom('Rn',   0.00000_real64 )/
    data atdata( 87)/atom('Fr',   0.00000_real64 )/
    data atdata( 88)/atom('Ra',   0.00000_real64 )/
    data atdata( 89)/atom('Ac',   0.00000_real64 )/
    data atdata( 90)/atom('Th',   0.00000_real64 )/
    data atdata( 91)/atom('Pa',   0.00000_real64 )/
    data atdata( 92)/atom('U',    0.00000_real64 )/
    data atdata( 93)/atom('Np',   0.00000_real64 )/
    data atdata( 94)/atom('Pu',   0.00000_real64 )/
    data atdata( 95)/atom('Am',   0.00000_real64 )/
    data atdata( 96)/atom('Cm',   0.00000_real64 )/
    data atdata( 97)/atom('Bk',   0.00000_real64 )/
    data atdata( 98)/atom('Cf',   0.00000_real64 )/
    data atdata( 99)/atom('Es',   0.00000_real64 )/
    data atdata(100)/atom('Fm',   0.00000_real64 )/
    data atdata(101)/atom('Md',   0.00000_real64 )/
    data atdata(102)/atom('No',   0.00000_real64 )/
    data atdata(103)/atom('Lr',   0.00000_real64 )/
    data atdata(104)/atom('Rf',   0.00000_real64 )/
    data atdata(105)/atom('Db',   0.00000_real64 )/
    data atdata(106)/atom('Sg',   0.00000_real64 )/
    data atdata(107)/atom('Bh',   0.00000_real64 )/
    data atdata(108)/atom('Hs',   0.00000_real64 )/
    data atdata(109)/atom('Mt',   0.00000_real64 )/
    data atdata(110)/atom('Ds',   0.00000_real64 )/


end module atomic
