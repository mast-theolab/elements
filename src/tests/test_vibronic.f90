program test_vibronic

    use numeric, only: f0, f1, realwp, is_close_to, small
    use blas_drv, only: xgemm
    use input, only: DataFile
    use datatypes, only: MoleculeDB, PropertyDB, VibrationsDB
    use geometry, only: Eckart_orient, superpose
    use output, only: prt_coord, prt_mat, prt_vec, sec_header
    use vibrational, only: build_modes, set_orientation
    use vibronic, only: Duschinsky_matrix, Duschinsky_shift
    use physics, only: phys_conv

    implicit none

    type(DataFile) :: dfile1, dfile2, dfile3
    type(MoleculeDB) :: mol1, mol2, mol3
    type(VibrationsDB) :: vib1, vib2, vib3
    class(PropertyDB), allocatable :: prop

    integer :: i, n_at3, n_vib, ia
    real(realwp), dimension(3,3) :: rot1, rot2, rot3
    real(realwp), dimension(:), allocatable :: kvec_ref, kvec
    real(realwp), dimension(:,:), allocatable :: jmat_ref, jmat
    real(realwp), dimension(:, :), allocatable :: grad
    real(realwp) :: x

    character(len=*), parameter :: S0_file = 'H2CO_S0_frq.fchk'
    character(len=*), parameter :: S2_opt_file = 'H2CO_S2_frq.fchk'
    character(len=*), parameter :: S2_file = 'H2CO.vac.S2.B3LYP.6-31+Gd.fchk'


    1000 format(' * Test on "',a,'" -> ',a)
    1010 format(' * Test on "',a,'" -----------> ',a)
    1011 format(' * Test on "',a,'" (',a7,') -> ',a)

    call sec_header(-1, 'Small Program to Test Vibronic Routines')

    call sec_header(1, 'Reading input data')
    dfile1 = DataFile(S0_file)
    dfile2 = DataFile(S2_opt_file)
    dfile3 = DataFile(S2_file)
    mol1 = dfile1%get_mol_data()
    mol2 = dfile2%get_mol_data()
    mol3 = dfile3%get_mol_data()
    n_at3 = 3*mol1%n_at
    
    call sec_header(1, 'Construction of normal modes')

    prop = dfile1%get_data(1, derorder=2)
    call build_modes(prop%data, mol1, vib1)
    
    prop = dfile2%get_data(1, derorder=2)
    call build_modes(prop%data, mol2, vib2)

    prop = dfile3%get_data(1, derorder=2)
    call build_modes(prop%data, mol3, vib3)

    prop = dfile3%get_data(1, derorder=1)

    call sec_header(1, 'Superposition procedure')
    ! Now we move mol1 to Eckart
    call sec_header(2, 'Initial state (Eckart)')
    call Eckart_orient(mol1, .true., rot_mat=rot1)
    call prt_coord(mol1%n_at, mol1%at_lab, mol1%at_crd)

    ! Now we superpose mol2 on it
    call sec_header(2, 'Final state (opt) (superposed)')
    call superpose(mol2, mol1%at_crd, rot_mat=rot2)
    call prt_coord(mol2%n_at, mol2%at_lab, mol2%at_crd)

    ! Now we superpose mol3 on it
    call sec_header(2, 'Final state (superposed)')
    call superpose(mol3, mol1%at_crd, rot_mat=rot3)
    call prt_coord(mol3%n_at, mol3%at_lab, mol3%at_crd)

    n_vib = vib1%n_vib
    call set_orientation(vib1)
    call set_orientation(vib2)
    call set_orientation(vib3)
    do i = 1, n_vib
        do ia = 1, n_at3, 3
            vib1%L_mat(ia:ia+2,i) = matmul(vib1%L_mat(ia:ia+2,i), rot1)
            vib2%L_mat(ia:ia+2,i) = matmul(vib2%L_mat(ia:ia+2,i), rot2)
            vib3%L_mat(ia:ia+2,i) = matmul(vib3%L_mat(ia:ia+2,i), rot3)
        end do
    end do

    call sec_header(1, 'Test on computation of the Duschinsky matrix')

    allocate(jmat_ref(n_vib,n_vib))

    call sec_header(2, 'Identity Duschinsky matrix')

    jmat_ref = 0.00_realwp
    do i = 1, n_vib
        jmat_ref(i,i) = 1.00_realwp
    end do

    jmat = Duschinsky_matrix(vib1%n_vib)

    if (all(is_close_to(jmat, jmat_ref, small, small))) then
        print 1000, '_dim', 'PASSED'
    else
        print 1000, '_dim', 'FAILED'
        stop 1
    end if

    jmat = Duschinsky_matrix(vib1%L_mat, vib2%L_mat, is_identity=.true.)

    if (all(is_close_to(jmat, jmat_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    jmat = Duschinsky_matrix(vib1, vib2, is_identity=.true.)

    if (all(is_close_to(jmat, jmat_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if

    call sec_header(2, 'Adiabatic Hessian Orthogonal Duschinsky matrix')

    jmat_ref = reshape([ &
         1.0217344924819350E-008_realwp,  0.9996450682024693E-000_realwp, -6.1931466839793724E-007_realwp, &    
        -2.3999036781250733E-007_realwp, -3.3569636310913962E-007_realwp,  2.5063029960540818E-002_realwp, &
         0.9999895363577169E-000_realwp, -1.0242417997169467E-008_realwp, -6.2625422185520692E-005_realwp, &     
         2.1707738316915286E-005_realwp, -4.4960651258961674E-007_realwp, -4.8291722454212944E-010_realwp, &
         6.6254903338511042E-005_realwp,  5.2661007593392972E-007_realwp,  0.9545631669302159E-000_realwp, &      
        -0.2978125091006597E-000_realwp,  1.0810428456230446E-002_realwp, -7.5655835554896773E-008_realwp, &
        -2.0396735536209910E-006_realwp,  3.8982783785573090E-007_realwp,  0.2978075036275402E-000_realwp, &
         0.9519584212714282E-000_realwp, -7.1315180186145896E-002_realwp, -1.2614175196741850E-007_realwp, &
        -4.2110744249167931E-007_realwp,  3.8628728919196828E-007_realwp,  1.0947474395588429E-002_realwp, &     
         7.1294271106128737E-002_realwp,  0.9973952474850436E-000_realwp, -1.2127969620050255E-006_realwp, &
        -2.4924684893224996E-010_realwp,  2.4963201275030343E-002_realwp, -1.3808764285922695E-007_realwp, &    
        -1.9350929923299694E-007_realwp, -1.2134422443411243E-006_realwp, -0.9996278421912512E-000_realwp], &
         [n_vib,n_vib])

    jmat = Duschinsky_matrix(vib1%L_mat, vib2%L_mat)

    if (all(is_close_to(jmat, jmat_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    jmat = Duschinsky_matrix(vib1, vib2)

    if (all(is_close_to(jmat, jmat_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if

    jmat = Duschinsky_matrix(vib1%L_mat, vib2%L_mat, at_mass1=mol1%at_mas, at_mass2=mol2%at_mas) 

    if (all(is_close_to(jmat, jmat_ref, small, small))) then
        print 1000, '_arr_mwg', 'PASSED'
    else
        print 1000, '_arr_mwg', 'FAILED'
        stop 1
    end if

    jmat = Duschinsky_matrix(vib1, vib2, mol1, mol2)

    if (all(is_close_to(jmat, jmat_ref, small, small))) then
        print 1000, '_arr_mwg', 'PASSED'
    else
        print 1000, '_arr_mwg', 'FAILED'
        stop 1
    end if

    call sec_header(1, 'Test on computation of the Duschinsky shift')

    call sec_header(2, 'Adiabatic Duschinsky shift')

    kvec_ref = [-1.5436864398590016E-011_realwp, &
                -1.5274448372792904E-012_realwp, &
                 4.6669595600321134E-000_realwp, &
                 1.5117772584274694E-000_realwp, &
                 1.4639506894651646E-000_realwp, &
                 6.1661786787681194E-013_realwp]

    kvec = Duschinsky_shift(vib1%L_mat, &
                            mol1%at_mas, &
                            coord1=mol1%at_crd, &
                            coord2=mol2%at_crd)

    if (all(is_close_to(kvec, kvec_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    kvec = Duschinsky_shift(vib1, mol1, mol2=mol2)

    if (all(is_close_to(kvec, kvec_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if

    call sec_header(2, 'Vertical Gradient Duschinsky shift')

    kvec_ref = [ 3.3140190430981111E-012_realwp, &
                -2.1917745396393684E-012_realwp, &
                 2.9851243022072689E-000_realwp, &
                 2.0984462594019946E-000_realwp, &
                 1.3071141884221709E-000_realwp, &
                -2.2534925314676713E-013_realwp]   

    ! Get the gradient data and rotate to the intermediate state orientation
    grad = matmul(rot3, reshape(prop%data, [3, mol3%n_at]))

    ! Important the frequecy of the initial state
    kvec = Duschinsky_shift(vib1%L_mat, &
                            mol1%at_mas, &
                            mode='VG', &
                            red_freq2=vib1%red_freq, &
                            grad2=grad)

    if (all(is_close_to(kvec, kvec_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    kvec = Duschinsky_shift(vib1, mol1, mode='VG', grad2=grad)

    if (all(is_close_to(kvec, kvec_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if
    call sec_header(2, 'Vertical Hessian Duschinsky shift')

    kvec_ref = [  4.1114568116253733E-012_realwp, &
                 -5.5780380314729427E-012_realwp, &
                  4.6718294300157908E-000_realwp, &
                  1.6649934616623288E-000_realwp, &
                  1.3923875612603047E-000_realwp, &
                 -4.1423158306241881E-013_realwp]

    jmat = Duschinsky_matrix(vib1, vib3)

    kvec = Duschinsky_shift(vib1%L_mat, &
                            mol1%at_mas, &
                            mode='VH', &
                            red_freq2=vib3%red_freq, &
                            Jmat=jmat, &
                            grad2=grad) 

    if (all(is_close_to(kvec, kvec_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    kvec = Duschinsky_shift(vib1, mol1, mode='VH', vib2=vib3, Jmat=jmat, grad2=grad)

    if (all(is_close_to(kvec, kvec_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if

end program test_vibronic
