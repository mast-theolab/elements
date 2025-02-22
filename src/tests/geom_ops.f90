program test_geom_ops
    !! Program to test some standard structural operations
    use numeric, only: small, realwp, is_close_to
    use output, only: sec_header
    use datatypes, only: MoleculeDB
    use geometry

    implicit none

    integer, parameter :: n_at_H2CO = 4, n_at_naph = 18, n_at_meox = 10
    real(realwp), parameter :: crd_H2CO(3,n_at_H2CO,2) = reshape([ &
        [[-0.571255_realwp,  -0.287251_realwp,   0.208084_realwp], &
         [ 0.446748_realwp,   0.224644_realwp,  -0.162731_realwp], &
         [ 0.944775_realwp,   1.026406_realwp,   0.416946_realwp], &
         [ 0.944777_realwp,  -0.076260_realwp,  -1.105229_realwp]], &
        [[ 0.375270_realwp,   0.284320_realwp,   0.460849_realwp], &
         [-0.297423_realwp,  -0.225320_realwp,  -0.365275_realwp], &
         [-0.633992_realwp,   0.341542_realwp,  -1.277301_realwp], &
         [-0.627284_realwp,  -1.297252_realwp,  -0.271499_realwp]]], &
         [3, n_at_H2CO, 2])
    real(realwp), parameter :: crd_meox(3,n_at_meox) = reshape([ &  ! meox.anh.log
        [-0.006714_realwp,   0.104984_realwp,  -0.400089_realwp], &
        [-0.205008_realwp,   0.045326_realwp,   1.008884_realwp], &
        [ 1.140780_realwp,  -0.001355_realwp,   0.440905_realwp], &
        [-0.572666_realwp,   0.960676_realwp,   1.463977_realwp], &
        [-0.659099_realwp,  -0.872850_realwp,   1.371991_realwp], &
        [ 1.932216_realwp,  -1.270083_realwp,   0.350377_realwp], &
        [ 1.726204_realwp,   0.915477_realwp,   0.483638_realwp], &
        [ 2.547932_realwp,  -1.276691_realwp,  -0.551417_realwp], &
        [ 2.596846_realwp,  -1.367924_realwp,   1.211914_realwp], &
        [ 1.270751_realwp,  -2.136540_realwp,   0.320091_realwp]], &
        [3, n_at_meox])
    real(realwp), parameter :: crd_naph(3,n_at_naph) = reshape([ &
        [ 0.02781033_realwp,  0.01441865_realwp,  -9.37094265_realwp], &
        [-0.79403778_realwp, -0.41167996_realwp,  -8.27852964_realwp], &
        [-0.16878179_realwp, -0.87520457_realwp,  -7.09016523_realwp], &
        [-2.20638673_realwp, -0.35623467_realwp,  -8.42067811_realwp], &
        [-0.59744567_realwp,  0.47794326_realwp, -10.55930706_realwp], &
        [ 1.44015927_realwp, -0.04102664_realwp,  -9.22879418_realwp], &
        [ 1.20256312_realwp, -0.91673788_realwp,  -6.98368369_realwp], &
        [-2.78166665_realwp,  0.09802964_realwp,  -9.58530135_realwp], &
        [ 2.01543919_realwp, -0.49529095_realwp,  -8.06417094_realwp], &
        [-1.96879057_realwp,  0.51947658_realwp, -10.66578860_realwp], &
        [-2.82896745_realwp, -0.68013906_realwp,  -7.59026615_realwp], &
        [-0.79425646_realwp, -1.19837188_realwp,  -6.26164297_realwp], &
        [ 2.06273999_realwp,  0.28287776_realwp, -10.05920614_realwp], &
        [ 0.02802900_realwp,  0.80111058_realwp, -11.38782932_realwp], &
        [ 1.66755644_realwp, -1.27349386_realwp,  -6.06904811_realwp], &
        [-3.86315338_realwp,  0.13515601_realwp,  -9.68048450_realwp], &
        [ 3.09692592_realwp, -0.53241732_realwp,  -7.96898779_realwp], &
        [-2.43378389_realwp,  0.87623255_realwp, -11.58042418_realwp]], &
        [3, n_at_naph])
    real(realwp) :: wgt_H2CO(n_at_H2CO), wgt_meox(n_at_meox), &
        wgt_naph(n_at_naph)

    type(MoleculeDB) :: mol_H2CO(2), mol_meox, mol_naph

    real(realwp), dimension(3) :: com, com_ref
    real(realwp), dimension(3,3) :: tensor, tensor_ref
    real(realwp), dimension(:,:), allocatable :: crd_new, crd_ref, crd_tmp

    data wgt_meox/15.99940_realwp, 2*12.01100_realwp, 2*1.00790_realwp, &
                  12.01100_realwp, 4*1.00790_realwp/, &
        wgt_naph/10*12.01100_realwp, 8*1.00790_realwp/, &
        wgt_H2CO/15.99940_realwp, 12.01100_realwp, 2*1.00790_realwp/

    1000 format(' * Test on "',a,'" -> ',a)
    1010 format(' * Test on "',a,'" -----------> ',a)
    1011 format(' * Test on "',a,'" (',a7,') -> ',a)

    ! Initialization of some components of MoleculeDB
    mol_H2CO(1)%n_at = n_at_H2CO
    mol_H2CO(2)%n_at = n_at_H2CO
    mol_H2CO(1)%at_mas = wgt_H2CO
    mol_H2CO(2)%at_mas = wgt_H2CO
    mol_H2CO(1)%at_crd = crd_H2CO(:,:,1)
    mol_H2CO(2)%at_crd = crd_H2CO(:,:,2)
    mol_meox%n_at = n_at_meox
    mol_meox%at_mas = wgt_meox
    mol_meox%at_crd = crd_meox
    mol_naph%n_at = n_at_naph
    mol_naph%at_mas = wgt_naph
    mol_naph%at_crd = crd_naph

    ! Start of the test
    call sec_header(-1, 'Test Program on Structural Operations')

    ! =====================================

    call sec_header(1, 'Test on computation of center of mass')

    call sec_header(2, 'Methyloxirane')

    com_ref = [0.711168782324E+00_realwp, -0.290201537423E+00_realwp, &
               0.336688063871E+00_realwp]

    com = center_of_mass(crd_meox, wgt_meox)
    ! print '(3e20.12)', com
    if (all(is_close_to(com, com_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    com = center_of_mass(n_at_meox, crd_meox, wgt_meox)
    if (all(is_close_to(com, com_ref, small, small))) then
        print 1000, '_dim', 'PASSED'
    else
        print 1000, '_dim', 'FAILED'
        stop 1
    end if

    com = center_of_mass(mol_meox)
    if (all(is_close_to(com, com_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if

    ! -------------------------------------

    call sec_header(2, 'Naphthalene')

    com_ref = [-0.383113728047E+00_realwp, -0.198630653906E+00_realwp, &
               -0.882473614500E+01_realwp]

    com = center_of_mass(crd_naph, wgt_naph)
    ! print '(3e20.12)', com
    if (all(is_close_to(com, com_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    com = center_of_mass(n_at_naph, crd_naph, wgt_naph)
    if (all(is_close_to(com, com_ref, small, small))) then
        print 1000, '_dim', 'PASSED'
    else
        print 1000, '_dim', 'FAILED'
        stop 1
    end if

    com = center_of_mass(mol_naph)
    if (all(is_close_to(com, com_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if

    ! =====================================

    call sec_header(1, 'Test on computation of inertia tensor')

    call sec_header(2, 'Formaldehyde (structure for S0)')

    tensor_ref = reshape([ &
        [ 0.541121209577E+01_realwp, -0.473557936782E+01_realwp, &
          0.343044172056E+01_realwp], &
        [-0.473557936782E+01_realwp,  0.118348514614E+02_realwp, &
          0.879114824411E+00_realwp], &
        [ 0.343044172056E+01_realwp,  0.879114824411E+00_realwp, &
          0.124116162657E+02_realwp]], &
        [3,3])

    tensor = inertia_moments(crd_H2CO(:,:,1), wgt_H2CO)
    ! print '(3e20.12)', tensor
    if (all(is_close_to(tensor, tensor_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    tensor = inertia_moments(n_at_H2CO, crd_H2CO(:,:,1), wgt_H2CO)
    if (all(is_close_to(tensor, tensor_ref, small, small))) then
        print 1000, '_dim', 'PASSED'
    else
        print 1000, '_dim', 'FAILED'
        stop 1
    end if

    tensor = inertia_moments(mol_H2CO(1))
    if (all(is_close_to(tensor, tensor_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if

    ! -------------------------------------

    call sec_header(2, 'Formaldehyde (structure for S2)')

    tensor_ref = reshape([ &
        [ 0.104361145144E+02_realwp, -0.311393410859E+01_realwp, &
         -0.505971881945E+01_realwp], &
        [-0.311393410859E+01_realwp,  0.108366096864E+02_realwp, &
         -0.300021597612E+01_realwp], &
        [-0.505971881945E+01_realwp, -0.300021597612E+01_realwp, &
          0.783424483607E+01_realwp]], &
        [3,3])

    tensor = inertia_moments(crd_H2CO(:,:,2), wgt_H2CO)
    ! print '(3e20.12)', tensor
    if (all(is_close_to(tensor, tensor_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    tensor = inertia_moments(n_at_H2CO, crd_H2CO(:,:,2), wgt_H2CO)
    if (all(is_close_to(tensor, tensor_ref, small, small))) then
        print 1000, '_dim', 'PASSED'
    else
        print 1000, '_dim', 'FAILED'
        stop 1
    end if

    tensor = inertia_moments(mol_H2CO(2))
    if (all(is_close_to(tensor, tensor_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if

    ! -------------------------------------

    call sec_header(2, 'Methyloxirane')

    tensor_ref = reshape([ &
        [ 0.550276255667E+02_realwp, 0.375946418606E+02_realwp, &
         -0.129824386512E+02_realwp], &
        [ 0.375946418606E+02_realwp, 0.104497391220E+03_realwp, &
          0.646880968529E+01_realwp], &
        [-0.129824386512E+02_realwp, 0.646880968529E+01_realwp, &
          0.109966938189E+03_realwp]], &
        [3,3])

    tensor = inertia_moments(crd_meox, wgt_meox)
    ! print '(3e20.12)', tensor
    if (all(is_close_to(tensor, tensor_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    tensor = inertia_moments(n_at_meox, crd_meox, wgt_meox)
    if (all(is_close_to(tensor, tensor_ref, small, small))) then
        print 1000, '_dim', 'PASSED'
    else
        print 1000, '_dim', 'FAILED'
        stop 1
    end if

    tensor = inertia_moments(mol_meox)
    if (all(is_close_to(tensor, tensor_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if

    ! -------------------------------------

    call sec_header(12, 'Naphthalene')

    tensor_ref = reshape([ &
        [ 0.102306701297E+05_realwp,  0.327524889239E+02_realwp, &
         -0.542313580512E+03_realwp], &
        [ 0.327524889239E+02_realwp,  0.105411116246E+05_realwp, &
         -0.142064494936E+03_realwp], &
        [-0.542313580512E+03_realwp, -0.142064494936E+03_realwp, &
          0.384996460382E+03_realwp]], &
        [3, 3])

    tensor = inertia_moments(crd_naph, wgt_naph)
    ! print '(3e20.12)', tensor
    if (all(is_close_to(tensor, tensor_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    tensor = inertia_moments(n_at_naph, crd_naph, wgt_naph)
    if (all(is_close_to(tensor, tensor_ref, small, small))) then
        print 1000, '_dim', 'PASSED'
    else
        print 1000, '_dim', 'FAILED'
        stop 1
    end if

    tensor = inertia_moments(mol_naph)
    if (all(is_close_to(tensor, tensor_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if

    ! =====================================

    call sec_header(1, 'Test on Eckart orientation')

    call sec_header(2, 'Formaldehyde (structure for S0)')

    allocate(crd_ref(3,n_at_H2CO), crd_tmp(3,n_at_H2CO))

    crd_ref = reshape([ &
        [-0.599133873225E+00_realwp, -0.498904758596E-07_realwp, &
         -0.130429458247E-07_realwp], &
        [ 0.599144209383E+00_realwp, -0.134664658737E-07_realwp, &
          0.528876946591E-07_realwp], &
        [ 0.118536710434E+01_realwp, -0.939798576528E+00_realwp, &
         -0.211605819259E-06_realwp], &
        [ 0.118536550045E+01_realwp,  0.939799528967E+00_realwp, &
         -0.211605603509E-06_realwp]], &
        [3, n_at_H2CO])
    
    crd_tmp = crd_H2CO(:,:,1)
    call Eckart_orient(crd_tmp, wgt_H2CO, new_crd=crd_tmp)
    ! print '(3e20.12)', crd_tmp
    if (all(is_close_to(crd_tmp, crd_ref, small, small))) then
        print 1010, '_arr', 'PASSED'
    else
        print 1010, '_arr', 'FAILED'
        stop 1
    end if

    crd_tmp = crd_H2CO(:,:,1)
    call Eckart_orient(n_at_H2CO, crd_tmp, wgt_H2CO, new_crd=crd_tmp)
    if (all(is_close_to(crd_tmp, crd_ref, small, small))) then
        print 1010, '_dim', 'PASSED'
    else
        print 1010, '_dim', 'FAILED'
        stop 1
    end if

    call Eckart_orient(mol_H2CO(1), new_crd=crd_tmp)
    if (all(is_close_to(crd_tmp, crd_ref, small, small))) then
        print 1011, '_mol', 'new_crd', 'PASSED'
    else
        print 1011, '_mol', 'new_crd', 'FAILED'
        stop 1
    end if

    call Eckart_orient(mol_H2CO(1), new_mol=mol_H2CO(1))
    if (all(is_close_to(mol_H2CO(1)%at_crd, crd_ref, small, small))) then
        print 1011, '_mol', 'new_mol', 'PASSED'
    else
        print 1011, '_mol', 'new_mol', 'FAILED'
        stop 1
    end if
    mol_H2CO(1)%at_crd = crd_H2CO(:,:,1)

    ! -------------------------------------

    call sec_header(2, 'Formaldehyde (structure for S2)')

    crd_ref = reshape([ &
        [ 0.590972097803E+00_realwp, -0.284050401051E-05_realwp, &
         -0.825099777851E-06_realwp], &
        [-0.590014635862E+00_realwp, -0.135679736780E-04_realwp, &
          0.331800980789E-05_realwp], &
        [-0.117507481974E+01_realwp, -0.961318527102E+00_realwp, &
         -0.132229903102E-04_realwp], &
        [-0.117489361984E+01_realwp,  0.961525304849E+00_realwp, &
         -0.132196274258E-04_realwp]], &
        [3, n_at_H2CO])

    crd_tmp = crd_H2CO(:,:,2)
    call Eckart_orient(crd_tmp, wgt_H2CO, new_crd=crd_tmp)
    ! print '(3e20.12)', crd_tmp
    if (all(is_close_to(crd_tmp, crd_ref, small, small))) then
        print 1010, '_arr', 'PASSED'
    else
        print 1010, '_arr', 'FAILED'
        stop 1
    end if

    crd_tmp = crd_H2CO(:,:,2)
    call Eckart_orient(n_at_H2CO, crd_tmp, wgt_H2CO, new_crd=crd_tmp)
    if (all(is_close_to(crd_tmp, crd_ref, small, small))) then
        print 1010, '_dim', 'PASSED'
    else
        print 1010, '_dim', 'FAILED'
        stop 1
    end if

    call Eckart_orient(mol_H2CO(2), new_crd=crd_tmp)
    if (all(is_close_to(crd_tmp, crd_ref, small, small))) then
        print 1011, '_mol', 'new_crd', 'PASSED'
    else
        print 1011, '_mol', 'new_crd', 'FAILED'
        stop 1
    end if

    call Eckart_orient(mol_H2CO(2), new_mol=mol_H2CO(2))
    if (all(is_close_to(mol_H2CO(2)%at_crd, crd_ref, small, small))) then
        print 1011, '_mol', 'new_mol', 'PASSED'
    else
        print 1011, '_mol', 'new_mol', 'FAILED'
        stop 1
    end if
    mol_H2CO(2)%at_crd = crd_H2CO(:,:,2)

    ! -------------------------------------
    
    call sec_header(2, 'Methyloxirane')

    deallocate(crd_ref, crd_tmp)
    allocate(crd_ref(3,n_at_meox), crd_tmp(3,n_at_meox))

    crd_ref = reshape([ &
        [-0.812642056488E+00_realwp,  0.730992060505E+00_realwp, &
         -0.140111510931E+00_realwp], &
        [-0.959851381212E+00_realwp, -0.685241704373E+00_realwp, &
         -0.113732730384E+00_realwp], &
        [ 0.210571542913E+00_realwp, -0.442685416136E-01_realwp, &
          0.482243933968E+00_realwp], &
        [-0.175999502771E+01_realwp, -0.105593332794E+01_realwp, &
          0.520700505902E+00_realwp], &
        [-0.861338164259E+00_realwp, -0.117694275745E+01_realwp, &
         -0.107789958171E+01_realwp], &
        [ 0.155457896863E+01_realwp, -0.488629679165E-01_realwp, &
         -0.179465919169E+00_realwp], &
        [ 0.221095880372E+00_realwp,  0.529077955474E-01_realwp, &
          0.156648296809E+01_realwp], &
        [ 0.208709400482E+01_realwp,  0.882661746051E+00_realwp, &
          0.231956713015E-01_realwp], &
        [ 0.216360245814E+01_realwp, -0.873309027546E+00_realwp, &
          0.198596224572E+00_realwp], &
        [ 0.145278270222E+01_realwp, -0.157386909723E+00_realwp, &
         -0.125977189206E+01_realwp]], &
        [3, n_at_meox])
    
    crd_tmp = crd_meox
    call Eckart_orient(crd_tmp, wgt_meox, new_crd=crd_tmp)
    ! print '(3e20.12)', crd_tmp
    if (all(is_close_to(crd_tmp, crd_ref, small, small))) then
        print 1010, '_arr', 'PASSED'
    else
        print 1010, '_arr', 'FAILED'
        stop 1
    end if

    crd_tmp = crd_meox
    call Eckart_orient(n_at_meox, crd_tmp, wgt_meox, new_crd=crd_tmp)
    if (all(is_close_to(crd_tmp, crd_ref, small, small))) then
        print 1010, '_dim', 'PASSED'
    else
        print 1010, '_dim', 'FAILED'
        stop 1
    end if

    call Eckart_orient(mol_meox, new_crd=crd_tmp)
    if (all(is_close_to(crd_tmp, crd_ref, small, small))) then
        print 1011, '_mol', 'new_crd', 'PASSED'
    else
        print 1011, '_mol', 'new_crd', 'FAILED'
        stop 1
    end if

    call Eckart_orient(mol_meox, new_mol=mol_meox)
    if (all(is_close_to(mol_meox%at_crd, crd_ref, small, small))) then
        print 1011, '_mol', 'new_mol', 'PASSED'
    else
        print 1011, '_mol', 'new_mol', 'FAILED'
        stop 1
    end if
    mol_meox%at_crd = crd_meox

    ! -------------------------------------

    call sec_header(2, 'Naphthalene')

    deallocate(crd_ref, crd_tmp)
    allocate(crd_ref(3,n_at_naph), crd_tmp(3,n_at_naph))

    crd_ref = reshape([ &
        [-0.516258640199E-07_realwp,  0.715954002350E+00_realwp, &
         -0.823327918881E-09_realwp], &
        [ 0.570729596093E-07_realwp, -0.715953999503E+00_realwp, &
         -0.121577448944E-08_realwp], &
        [ 0.124412914036E+01_realwp, -0.140163097314E+01_realwp, &
         -0.738289827425E-09_realwp], &
        [-0.124412891950E+01_realwp, -0.140163104733E+01_realwp, &
         -0.248367655359E-08_realwp], &
        [-0.124412914310E+01_realwp,  0.140163097025E+01_realwp, &
         -0.130081199803E-08_realwp], &
        [ 0.124412891676E+01_realwp,  0.140163104444E+01_realwp, &
          0.444574755890E-09_realwp], &
        [ 0.243270710882E+01_realwp, -0.708137923960E+00_realwp, &
          0.349660267852E-08_realwp], &
        [-0.243270699777E+01_realwp, -0.708138082189E+00_realwp, &
         -0.814448799847E-09_realwp], &
        [ 0.243270699503E+01_realwp,  0.708138079296E+00_realwp, &
         -0.122465302561E-08_realwp], &
        [-0.243270710546E+01_realwp,  0.708137929782E+00_realwp, &
          0.378067417486E-08_realwp], &
        [-0.124236183293E+01_realwp, -0.248887605537E+01_realwp, &
          0.753562593182E-08_realwp], &
        [ 0.124236222076E+01_realwp, -0.248887597147E+01_realwp, &
          0.718739853979E-08_realwp], &
        [ 0.124236182810E+01_realwp,  0.248887605545E+01_realwp, &
         -0.258348200577E-09_realwp], &
        [-0.124236222559E+01_realwp,  0.248887597155E+01_realwp, &
          0.898795142027E-10_realwp], &
        [ 0.337696317237E+01_realwp, -0.124519690355E+01_realwp, &
         -0.234612538392E-08_realwp], &
        [-0.337696297235E+01_realwp, -0.124519711157E+01_realwp, &
          0.220858712669E-08_realwp], &
        [ 0.337696296960E+01_realwp,  0.124519710868E+01_realwp, &
         -0.424768889664E-08_realwp], &
        [-0.337696316693E+01_realwp,  0.124519690640E+01_realwp, &
          0.307023104059E-09_realwp]], &
        [3, n_at_naph])
    
    crd_tmp = crd_naph
    call Eckart_orient(crd_tmp, wgt_naph, new_crd=crd_tmp)
    ! print '(3e20.12)', crd_tmp
    if (all(is_close_to(crd_tmp, crd_ref, small, small))) then
        print 1010, '_arr', 'PASSED'
    else
        print 1010, '_arr', 'FAILED'
        stop 1
    end if

    crd_tmp = crd_naph
    call Eckart_orient(n_at_naph, crd_tmp, wgt_naph, new_crd=crd_tmp)
    if (all(is_close_to(crd_tmp, crd_ref, small, small))) then
        print 1010, '_dim', 'PASSED'
    else
        print 1010, '_dim', 'FAILED'
        stop 1
    end if

    call Eckart_orient(mol_naph, new_crd=crd_tmp)
    if (all(is_close_to(crd_tmp, crd_ref, small, small))) then
        print 1011, '_mol', 'new_crd', 'PASSED'
    else
        print 1011, '_mol', 'new_crd', 'FAILED'
        stop 1
    end if

    call Eckart_orient(mol_naph, new_mol=mol_naph)
    if (all(is_close_to(mol_naph%at_crd, crd_ref, small, small))) then
        print 1011, '_mol', 'new_mol', 'PASSED'
    else
        print 1011, '_mol', 'new_mol', 'FAILED'
        stop 1
    end if
    mol_naph%at_crd = crd_naph

    ! =====================================

    call sec_header(1, 'Test on superposition')

    call sec_header(2, 'Formaldehyde - S0 in its original geometry')

    deallocate(crd_ref, crd_tmp)
    allocate(crd_new(3,n_at_H2CO), crd_ref(3,n_at_H2CO), crd_tmp(3,n_at_H2CO))

    crd_ref = reshape([ &
        [-0.564321560889E+00_realwp, -0.283770377403E+00_realwp, &
          0.205548769218E+00_realwp], &
        [ 0.438993698933E+00_realwp,  0.220757774463E+00_realwp, &
         -0.159881796838E+00_realwp], &
        [ 0.936004614661E+00_realwp,  0.103464975607E+01_realwp, &
          0.437581704480E+00_realwp], &
        [ 0.935892900180E+00_realwp, -0.934435499567E-01_realwp, &
         -0.111957400878E+01_realwp]], &
        [3, n_at_H2CO])

    crd_new = crd_H2CO(:,:,2)
    call superpose(crd_new, crd_H2CO(:,:,1), wgt_H2CO)
    ! print '(3e20.12)', crd_new
    if (all(is_close_to(crd_new, crd_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    crd_new = crd_H2CO(:,:,2)
    call superpose(n_at_H2CO, crd_new, crd_H2CO(:,:,1), wgt_H2CO)
    if (all(is_close_to(crd_new, crd_ref, small, small))) then
        print 1000, '_dim', 'PASSED'
    else
        print 1000, '_dim', 'FAILED'
        stop 1
    end if

    call superpose(mol_H2CO(2), crd_H2CO(:,:,1))
    if (all(is_close_to(mol_H2CO(2)%at_crd, crd_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if
    mol_H2CO(2)%at_crd = crd_H2CO(:,:,2)

    ! -------------------------------------

    call sec_header(2, 'Formaldehyde - S0 in Eckart orientation')

    crd_ref = reshape([ &
        [-0.590972097704E+00_realwp,  0.111906383143E-04_realwp, &
          0.816187199813E-06_realwp], &
        [ 0.590014635373E+00_realwp, -0.275763834744E-04_realwp, &
         -0.330911169769E-05_realwp], &
        [ 0.117505199533E+01_realwp, -0.961346426019E+00_realwp, &
          0.132390632404E-04_realwp], &
        [ 0.117491644850E+01_realwp,  0.961497409692E+00_realwp, &
          0.132389952132E-04_realwp]], &
        [3, n_at_H2CO])

    crd_tmp = crd_H2CO(:,:,1)
    call Eckart_orient(crd_tmp, wgt_H2CO, new_crd=crd_tmp)

    crd_new = crd_H2CO(:,:,2)
    call superpose(crd_new, crd_tmp, wgt_H2CO)
    ! print '(3e20.12)', crd_new
    if (all(is_close_to(crd_new, crd_ref, small, small))) then
        print 1000, '_arr', 'PASSED'
    else
        print 1000, '_arr', 'FAILED'
        stop 1
    end if

    crd_new = crd_H2CO(:,:,2)
    call superpose(n_at_H2CO, crd_new, crd_tmp, wgt_H2CO)
    if (all(is_close_to(crd_new, crd_ref, small, small))) then
        print 1000, '_dim', 'PASSED'
    else
        print 1000, '_dim', 'FAILED'
        stop 1
    end if

    call superpose(mol_H2CO(2), crd_tmp)
    if (all(is_close_to(mol_H2CO(2)%at_crd, crd_ref, small, small))) then
        print 1000, '_mol', 'PASSED'
    else
        print 1000, '_mol', 'FAILED'
        stop 1
    end if
    mol_H2CO(2)%at_crd = crd_H2CO(:,:,2)

end program test_geom_ops