module electronic
    !! Module related to electronic structure calculations
    !!
    !! Module managing different levels of operations related to
    !!   electronic calculations:
    !! * orbitals-related constants and operations
    !! * integral calculations
    use iso_fortran_env, only: real64

    use output, only: iu_out
    use exception, only: ArgumentError, BaseException, InitError, &
        RaiseArgError, RaiseError
    use math, only: build_PascalTriangle, cross, itri_pa, pi, phii_xn_phij
    use basisset, only: max_nxyz, PrimitiveFunction, set_primC_comp, &
        transfo_cart2pure

    type, public :: ovij_1e
        integer :: n_ao = 0
        !! number of atomic orbitals
        real(real64), dimension(:,:), allocatable :: i_j
        !! Store the overlap integral < i | j >
        real(real64), dimension(:,:,:), allocatable :: i_r_j
        !! Store the integral < i | r | j > / r=x,y,z
        real(real64), dimension(:,:,:), allocatable :: i_rr_j_s
        !! Store the integral < i | r^2 | j >  (symmetric form)
        !! The elements are stored as: xx,xy,xz,yy,yz,zz
        real(real64), dimension(:,:,:), allocatable :: i_p_j
        !! Store the integral < i | p | j > / p=d/dx,d/dy,d/dz
        real(real64), dimension(:,:,:), allocatable :: i_rxp_j
        !! Store the integral < i | r x p | j >
        real(real64), dimension(:,:,:), allocatable :: i_rrjxp_j
        !! Store the integral < i | (r-Rj) x p | j >
        real(real64), dimension(:,:,:), allocatable :: i_rrijxp_j
        !! Store the integral < i | (r-(Ri+Rj)/2) x p | j >
        real(real64), dimension(:,:,:), allocatable :: i_rjxr_j
        !! Store the integral < i | Rj x r | j >
        real(real64), dimension(:,:,:), allocatable :: i_rixr_j
        !! Store the integral < i | Ri x r | j >
    end type ovij_1e

    integer, private :: MAXAO = 5000
    !! Maximum number of atomic orbitals, for storage

    interface convert_AO2MO
        module procedure convert_AO2MO_1, convert_AO2MO_N
    end interface convert_AO2MO

contains

! ======================================================================

subroutine overlap_ao_1e(iout, n_at, n_ao, qty_flag, ondisk, inmem, &
                         at_crd, nprim_per_at, bsetDB, ovij, err)
    !! Compute overlap one-electron integrals between atomic orbitals
    !!
    !! Computes different types of 1-electron overlap integrals between
    !!   atomic orbitals upon request:
    !! * < i | j >
    !! * < i | r | j > (e.g., electric dipole)
    !! * < i | r^2 | j > (e.g., electric quadrupole)
    !! * < i | p | j > (e.g., gradient, velocity gauge)
    !! * < i | r x p | j > (e.g., angular momentum)
    !! * < i | (r-Rj) x p | j >  (alternative form for angular momentum)
    !! * < i | (r-(Rj+Ri)/2) x p | j >  (alt. form for angular momentum)
    !! * < i | Ri x r | j >  (GIAO-like form)
    !! * < i | Rj x r | j >  (GIAO-like form)
    implicit none

    integer, intent(in) :: iout
    !! Unit for output
    integer, intent(in) :: n_at
    !! Number of atoms
    integer, intent(in) :: n_ao
    !! Number of atomic orbitals
    integer, intent(in) :: qty_flag
    !! Quantity flag, as bits array (see below)
    logical, intent(in) :: ondisk
    !! If True, the integrals are stored in dedicated binary files.
    logical, intent(in) :: inmem
    !! If True, the work arrays are stored in `ovij`
    real(real64), dimension(:,:), intent(in) :: at_crd
    !! Atomic coordinates (in au)
    integer, dimension(:), intent(in) :: nprim_per_at
    !! Number of basis primitives per atom
    class(PrimitiveFunction), dimension(:,:), intent(in) :: bsetDB
    !! Basis set DB
    class(ovij_1e), allocatable, intent(out) :: ovij
    !! Array used to store overlap integrals.
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance

    logical :: do_ij, do_p, do_r, do_r2, do_rixr, do_rjxr, do_rrijxp, &
        do_rrjxp, do_rxp
    
    integer, parameter :: &
        dimPa = 24  ! maximum used dimension for Pascal's triangle
    real(real64), parameter :: &
        tol_expm = 100.0_real64, &  ! Max. value of x for e^-x to be relevant
        sqpi3 = sqrt(pi**3)
    integer :: iu_ij, iu_p, iu_r, iu_rixr, iu_rjxr, iu_rrijxp, iu_rrjxp, &
        iu_rxp, iu_r2, lrec, ndi, ndj
    integer :: iprim, jprim 
    integer :: i, idi, idj, ix, j, ia, ia0, ii, ja, ja0, jj
    integer, dimension(:,:), allocatable :: ldi, ldj
    real(real64), parameter :: f2 = 2.0_real64, pt5 = 0.5_real64
    real(real64) :: a, ai, aj, cier, ebase, r2ij, x
    real(real64), dimension(3) :: fij_p, fij_r, fij_r2, fij, fijm1, &
        fijp1, ri, rj
    ! Temporary arrays for one couple of primitives, store data wrt Cart. basis
    real(real64), dimension(max_nxyz,max_nxyz) :: oij
    real(real64), dimension(3,max_nxyz,max_nxyz) :: oij_p, oij_r, oij_rixr, &
        oij_rjxr, oij_rrijxp, oij_rrjxp, oij_rxp
    real(real64), dimension(6,max_nxyz,max_nxyz) :: oij_r2
    ! Longer-term temporary arrays, for one initial-set primitive.
    real(real64), dimension(MAXAO,max_nxyz) :: tij
    real(real64), dimension(3,MAXAO,max_nxyz) :: tij_p, tij_r, tij_rixr, &
        tij_rjxr, tij_rrijxp, tij_rrjxp, tij_rxp
    real(real64), dimension(6,MAXAO,max_nxyz) :: tij_r2
    real(real64), dimension(max_nxyz,max_nxyz), target :: ident_mat
    real(real64), dimension(:), allocatable :: ci, cj
    real(real64), dimension(:,:), allocatable, target :: c2p_D, c2p_F, c2p_G, &
        c2p_H, c2p_I
    real(real64), dimension(:,:), pointer :: c2pi, c2pj
    class(BaseException), allocatable :: suberr

    err = InitError()

    ident_mat = 0.0_real64
    do i = 1, max_nxyz
        ident_mat(i,i) = 1.0_real64
    end do

    if (n_ao > MAXAO) then
        call RaiseError(err, 'Too many atomic orbitals')
        return
    end if

    write(iout, *)  'Entering overlap_ao_1e'

    allocate(ovij)
    tij        = 0.0_real64
    tij_r      = 0.0_real64
    tij_r2     = 0.0_real64
    tij_p      = 0.0_real64
    tij_rxp    = 0.0_real64
    tij_rrjxp  = 0.0_real64
    tij_rrijxp = 0.0_real64
    tij_rjxr   = 0.0_real64
    tij_rixr   = 0.0_real64

    ! analyse qty_flag to find which quantity to retrieve
    ! if null, nothing to do
    if (qty_flag == 0) then
        call RaiseError(err, 'Empty qty_flag, nothing to do.')
        return
    end if
    do_ij     = btest(qty_flag, 0)
    do_r      = btest(qty_flag, 1)
    do_r2     = btest(qty_flag, 2)
    do_p      = btest(qty_flag, 3)
    do_rxp    = btest(qty_flag, 4)
    do_rrjxp  = btest(qty_flag, 5)
    do_rrijxp = btest(qty_flag, 6)
    do_rjxr   = btest(qty_flag, 7)
    do_rixr   = btest(qty_flag, 8)

    ! Now, initialize the ovij structure
    ovij%n_ao = n_ao
    if (inmem) then
        if (do_ij)     allocate(ovij%i_j(n_ao,n_ao))
        if (do_r)      allocate(ovij%i_r_j(3,n_ao,n_ao))
        if (do_r2)     allocate(ovij%i_rr_j_s(6,n_ao,n_ao))
        if (do_p)      allocate(ovij%i_p_j(3,n_ao,n_ao))
        if (do_rxp)    allocate(ovij%i_rxp_j(3,n_ao,n_ao))
        if (do_rrjxp)  allocate(ovij%i_rrjxp_j(3,n_ao,n_ao))
        if (do_rrijxp) allocate(ovij%i_rrijxp_j(3,n_ao,n_ao))
        if (do_rjxr)   allocate(ovij%i_rjxr_j(3,n_ao,n_ao))
        if (do_rixr)   allocate(ovij%i_rixr_j(3,n_ao,n_ao))
    end if
    if (ondisk) then
        inquire(iolength=lrec) tij(:n_ao,1)
        if (do_ij) &
            open(newunit=iu_ij, file='SAO.SCR', access='direct', &
                 action='write', recl=lrec)
        if (do_r) &
            open(newunit=iu_r, file='P.SCR', access='direct', &
                 action='write', recl=lrec)
        if (do_r2) &
            open(newunit=iu_r2, file='Q.SCR', access='direct', &
                 action='write', recl=lrec)
        if (do_p) &
            open(newunit=iu_p, file='V.SCR', access='direct', &
                 action='write', recl=lrec)
        if (do_rxp) &
            open(newunit=iu_rxp, file='L.SCR', access='direct', &
                 action='write', recl=lrec)
        if (do_rrjxp) &
            open(newunit=iu_rrjxp, file='W.SCR', access='direct', &
                 action='write', recl=lrec)
        if (do_rrijxp) &
            open(newunit=iu_rrijxp, file='T.SCR', access='direct', &
                 action='write', recl=lrec)
        if (do_rjxr) &
            open(newunit=iu_rjxr, file='D.SCR', access='direct', &
                 action='write', recl=lrec)
        if (do_rixr) &
            open(newunit=iu_rixr, file='E.SCR', access='direct', &
                 action='write', recl=lrec)
    end if

    if (.not.allocated(itri_pa)) then
        call build_PascalTriangle(dimPa)
    else
        if (size(itri_pa, 1) < 24) then
            deallocate(itri_pa)
            call build_PascalTriangle(dimPa)
        end if
    end if

    ia0 = 0
    do ia = 1, n_at
        ri = at_crd(:,ia)
        do iprim = 1, nprim_per_at(ia)
            ai = bsetDB(ia,iprim)%alpha
            call set_primC_comp(bsetDB(ia,iprim), ndi, ldi, ci, suberr)
            if (suberr%raised()) then
                select type(suberr)
                    class is (ArgumentError)
                        call RaiseError(err, &
                                        'Unrecognized shell type for iprim')
                        return
                    class default
                        call RaiseError(err, 'Generic error')
                        return
                end select
            end if
            if (bsetDB(ia,iprim)%pure) then
                select case (bsetDB(ia,iprim)%shelltype)
                    case ('D')
                        if (.not.allocated(c2p_D)) &
                            c2p_D = transfo_cart2pure(bsetDB(ia,iprim)%L)
                        c2pi => c2p_D
                    case ('F')
                        if (.not.allocated(c2p_F)) &
                            c2p_F = transfo_cart2pure(bsetDB(ia,iprim)%L)
                        c2pi => c2p_F
                    case ('G')
                        if (.not.allocated(c2p_G)) &
                            c2p_G = transfo_cart2pure(bsetDB(ia,iprim)%L)
                        c2pi => c2p_G
                    case ('H')
                        if (.not.allocated(c2p_H)) &
                            c2p_H = transfo_cart2pure(bsetDB(ia,iprim)%L)
                        c2pi => c2p_H
                    case ('I')
                        if (.not.allocated(c2p_I)) &
                            c2p_I = transfo_cart2pure(bsetDB(ia,iprim)%L)
                        c2pi => c2p_I
                    case default
                        c2pi => ident_mat
                end select
            else
                c2pi => ident_mat
            end if
            ja0 = 0
            do ja = 1, n_at
                rj = at_crd(:,ja)
                r2ij = sum((rj-ri)**2)
                do jprim = 1, nprim_per_at(ja)
                    aj = bsetDB(ja,jprim)%alpha
                    call set_primC_comp(bsetDB(ja,jprim), ndj, ldj, cj, suberr)
                    if (suberr%raised()) then
                        select type(suberr)
                            class is (ArgumentError)
                                call RaiseError(&
                                    err, 'Unrecognized shell type for jprim')
                                return
                            class default
                                call RaiseError(err, 'Generic error')
                                return
                        end select
                    end if
                    if (bsetDB(ja,jprim)%pure) then
                        select case (bsetDB(ja,jprim)%shelltype)
                            case ('D')
                                if (.not.allocated(c2p_D)) &
                                    c2p_D = transfo_cart2pure(bsetDB(ja,jprim)%L)
                                c2pj => c2p_D
                            case ('F')
                                if (.not.allocated(c2p_F)) &
                                    c2p_F = transfo_cart2pure(bsetDB(ja,jprim)%L)
                                c2pj => c2p_F
                            case ('G')
                                if (.not.allocated(c2p_G)) &
                                    c2p_G = transfo_cart2pure(bsetDB(ja,jprim)%L)
                                c2pj => c2p_G
                            case ('H')
                                if (.not.allocated(c2p_H)) &
                                    c2p_H = transfo_cart2pure(bsetDB(ja,jprim)%L)
                                c2pj => c2p_H
                            case ('I')
                                if (.not.allocated(c2p_I)) &
                                    c2p_I = transfo_cart2pure(bsetDB(ja,jprim)%L)
                                c2pj => c2p_I
                            case default
                                c2pj => ident_mat
                        end select
                    else
                        c2pj => ident_mat
                    end if
                    a = 1.0_real64/(ai+aj)
                    x = ai*aj*r2ij*a
                    if (x < tol_expm) then
                        ebase = sqpi3*exp(-x)
                    else
                        if (bsetDB(ja,jprim)%shell_last) &
                            ja0 = ja0 + bsetDB(ja,jprim)%ndim
                        cycle
                    end if

                    ! Now loop over orbitals to which primitives contribute
                    do idi = 1, ndi
                        cier = ebase*ci(idi)
                        do idj = 1, ndj
                            ! Note: for optimization, oij* are stored in a
                            !       transposed way.
                            a = cier*cj(idj)
                            do ix = 1, 3
                                fij(ix) = phii_xn_phij(.true., ldi(ix,idi), &
                                    ri(ix), ai, ldj(ix,idj), rj(ix), aj, 0)
                                fijp1(ix) = phii_xn_phij(.true., ldi(ix,idi), &
                                    ri(ix), ai, ldj(ix,idj)+1, rj(ix), aj, 0)
                                if (ldj(ix,idj) > 0) then
                                    fijm1(ix) = phii_xn_phij(.true., &
                                        ldi(ix,idi), ri(ix), ai, &
                                        ldj(ix,idj)-1, rj(ix), aj, 0)
                                else
                                    fijm1(ix) = 0.0_real64
                                end if
                            end do
                            fij_r = phii_xn_phij(.true., ldi(:,idi), ri, ai, &
                                                 ldj(:,idj), rj, aj, 1)
                            fij_r2 = phii_xn_phij(.true., ldi(:,idi), ri, &
                                                  ai, ldj(:,idj), rj, aj, 2)
                            fij_p = ldj(:,idj)*fijm1 - f2*aj*fijp1
                            if (do_ij) oij(idj,idi) = a*product(fij)
                            if (do_r) then
                                oij_r(1,idj,idi) = a*fij_r(1)*fij(2)*fij(3)
                                oij_r(2,idj,idi) = a*fij(1)*fij_r(2)*fij(3)
                                oij_r(3,idj,idi) = a*fij(1)*fij(2)*fij_r(3)
                            end if
                            if (do_r2) then
                                oij_r2(1,idj,idi) = a*fij_r2(1)*fij(2)*fij(3)
                                oij_r2(2,idj,idi) = a*fij_r(1)*fij_r(2)*fij(3)
                                oij_r2(3,idj,idi) = a*fij_r(1)*fij(2)*fij_r(3)
                                oij_r2(4,idj,idi) = a*fij(1)*fij_r2(2)*fij(3)
                                oij_r2(5,idj,idi) = a*fij(1)*fij_r(2)*fij_r(3)
                                oij_r2(6,idj,idi) = a*fij(1)*fij(2)*fij_r2(3)
                            end if
                            if (do_p) then
                                oij_p(1,idj,idi) = a*fij_p(1)*fij(2)*fij(3)
                                oij_p(2,idj,idi) = a*fij(1)*fij_p(2)*fij(3)
                                oij_p(3,idj,idi) = a*fij(1)*fij(2)*fij_p(3)
                            end if
                            if (do_rxp) &
                                oij_rxp(:,idj,idi) = &
                                    a*cross(fij_r, fij_p)*fij*pt5
                            if (do_rrjxp) &
                                oij_rrjxp(:,idj,idi) = &
                                    a*cross(fij_r-rj, fij_p)*fij*pt5
                            if (do_rrijxp) &
                                oij_rrijxp(:,idj,idi) = &
                                    a*fij*pt5*cross(fij_r-(ri+rj)*pt5, fij_p)
                            if (do_rjxr) &
                                oij_rjxr(:,idj,idi) = &
                                    a*cross(fij_r, rj)*fij*pt5
                            if (do_rixr) &
                                oij_rixr(:,idj,idi) = &
                                    a*cross(fij_r, ri)*fij*pt5
                        end do
                    end do

                    ! Check now if we need to convert to pure basis
                    if (bsetDB(ia,iprim)%pure .or. bsetDB(ja,jprim)%pure) then
                        do i = 1, bsetDB(ia,iprim)%ndim
                            do j = 1, bsetDB(ja,jprim)%ndim
                                jj = ja0 + j
                                do idi = 1, ndi
                                    x = c2pi(i,idi)
                                    do idj = 1, ndj
                                        a = x*c2pj(j,idj)
                                        if (do_ij) &
                                            tij(jj,i) = tij(jj,i) &
                                                + a*oij(idj,idi)
                                        if (do_r) &
                                            tij_r(:,jj,i) = tij_r(:,jj,i) &
                                                + a*oij_r(:,idj,idi)
                                        if (do_r2) &
                                            tij_r2(:,jj,i) = tij_r2(:,jj,i) &
                                                + a*oij_r2(:,idj,idi)
                                        if (do_p) &
                                            tij_p(:,jj,i) = tij_p(:,jj,i) &
                                                + a*oij_p(:,idj,idi)
                                        if (do_rxp) &
                                            tij_rxp(:,jj,i) = tij_rxp(:,jj,i) &
                                                + a*oij_rxp(:,idj,idi)
                                        if (do_rrjxp) &
                                            tij_rrjxp(:,jj,i) = &
                                                tij_rrjxp(:,jj,i) &
                                                + a*oij_rrjxp(:,idj,idi)
                                        if (do_rrijxp) &
                                            tij_rrijxp(:,jj,i) = &
                                                tij_rrijxp(:,jj,i) &
                                                + a*oij_rrijxp(:,idj,idi)
                                        if (do_rjxr) &
                                            tij_rjxr(:,jj,i) = &
                                                tij_rjxr(:,jj,i) &
                                                + a*oij_rjxr(:,idj,idi)
                                        if (do_rixr) &
                                            tij_rixr(:,jj,i) = &
                                                tij_rixr(:,jj,i) &
                                                + a*oij_rixr(:,idj,idi)
                                    end do
                                end do
                            end do
                        end do
                    else
                        do idi = 1, ndi
                            do idj = 1, ndj
                                jj = ja0 + idj
                                if (do_ij) &
                                    tij(jj,idi) = tij(jj,idi) &
                                        + oij(idj,idi)
                                if (do_r) &
                                    tij_r(:,jj,idi) = tij_r(:,jj,idi) &
                                        + oij_r(:,idj,idi)
                                if (do_r2) &
                                    tij_r2(:,jj,idi) = tij_r2(:,jj,idi) &
                                        + oij_r2(:,idj,idi)
                                if (do_p) &
                                    tij_p(:,jj,idi) = tij_p(:,jj,idi) &
                                        + oij_p(:,idj,idi)
                                if (do_rxp) &
                                    tij_rxp(:,jj,idi) = tij_rxp(:,jj,idi) &
                                        + oij_rxp(:,idj,idi)
                                if (do_rrjxp) &
                                    tij_rrjxp(:,jj,idi) = tij_rrjxp(:,jj,idi) &
                                        + oij_rrjxp(:,idj,idi)
                                if (do_rrijxp) &
                                    tij_rrijxp(:,jj,idi) = &
                                        tij_rrijxp(:,jj,idi) &
                                        + oij_rrijxp(:,idj,idi)
                                if (do_rjxr) &
                                    tij_rjxr(:,jj,idi) = tij_rjxr(:,jj,idi) &
                                        + oij_rjxr(:,idj,idi)
                                if (do_rixr) &
                                    tij_rixr(:,jj,idi) = tij_rixr(:,jj,idi) &
                                        + oij_rixr(:,idj,idi)
                            end do
                        end do
                    end if
                    if (bsetDB(ja,jprim)%shell_last) &
                        ja0 = ja0 + bsetDB(ja,jprim)%ndim
                end do
            end do
            if (bsetDB(ia,iprim)%shell_last) then
                ndi = bsetDB(ia,iprim)%ndim
                if (ondisk) then
                    if (do_ij) then
                        do i = 1, ndi
                            write(iu_ij, rec=ia0+i) (tij(j,i),j=1,n_ao)
                        end do
                    end if
                    if (do_r) then
                        do ix = 1, 3
                            ii = (ix-1)*n_ao + ia0
                            do i = 1, ndi
                                write(iu_r, rec=ii+i)  (tij_r(ix,j,i),j=1,n_ao)
                            end do
                        end do
                    end if
                    if (do_r2) then
                        do ix = 1, 6
                            ii = (ix-1)*n_ao + ia0
                            do i = 1, ndi
                                write(iu_r2, rec=ii+i) (tij_r2(ix,j,i),j=1,n_ao)
                            end do
                        end do
                    end if
                    if (do_p) then
                        do ix = 1, 3
                            ii = (ix-1)*n_ao + ia0
                            do i = 1, ndi
                                write(iu_p, rec=ii+i) &
                                    (tij_p(ix,j,i),j=1,n_ao)
                            end do
                        end do
                    end if
                    if (do_rxp) then
                        do ix = 1, 3
                            ii = (ix-1)*n_ao + ia0
                            do i = 1, ndi
                                write(iu_rxp, rec=ii+i) &
                                    (tij_rxp(ix,j,i),j=1,n_ao)
                            end do
                        end do
                    end if
                    if (do_rrjxp) then
                        do ix = 1, 3
                            ii = (ix-1)*n_ao + ia0
                            do i = 1, ndi
                                write(iu_rrjxp, rec=ii+i) &
                                    (tij_rrjxp(ix,j,i),j=1,n_ao)
                            end do
                        end do
                    end if
                    if (do_rrijxp) then
                        do ix = 1, 3
                            ii = (ix-1)*n_ao + ia0
                            do i = 1, ndi
                                write(iu_rrijxp, rec=ii+i) &
                                    (tij_rrijxp(ix,j,i),j=1,n_ao)
                            end do
                        end do
                    end if
                    if (do_rjxr) then
                        do ix = 1, 3
                            ii = (ix-1)*n_ao + ia0
                            do i = 1, ndi
                                write(iu_rjxr, rec=ii+i) &
                                    (tij_rjxr(ix,j,i),j=1,n_ao)
                            end do
                        end do
                    end if
                    if (do_rixr) then
                        do ix = 1, 3
                            ii = (ix-1)*n_ao + ia0
                            do i = 1, ndi
                                write(iu_rixr, rec=ii+i) &
                                    (tij_rixr(ix,j,i),j=1,n_ao)
                            end do
                        end do
                    end if
                end if
                if (inmem) then
                    if (do_ij) &
                        ovij%i_j(:,ia0+1:ia0+ndi) &
                            = tij(:n_ao,:ndi)
                    if (do_r) &
                        ovij%i_r_j(:,:,ia0+1:ia0+ndi) &
                            = tij_r(:,:n_ao,:ndi)
                    if (do_r2) &
                        ovij%i_rr_j_s(:,:,ia0+1:ia0+ndi) &
                            = tij_r2(:,:n_ao,:ndi)
                    if (do_p) &
                        ovij%i_p_j(:,:,ia0+1:ia0+ndi) &
                            = tij_p(:,:n_ao,:ndi)
                    if (do_rxp) &
                        ovij%i_rxp_j(:,:,ia0+1:ia0+ndi) &
                            = tij_rxp(:,:n_ao,:ndi)
                    if (do_rrjxp) &
                        ovij%i_rrjxp_j(:,:,ia0+1:ia0+ndi) &
                            = tij_rrjxp(:,:n_ao,:ndi)
                    if (do_rrijxp) &
                        ovij%i_rrijxp_j(:,:,ia0+1:ia0+ndi) &
                            = tij_rrijxp(:,:n_ao,:ndi)
                    if (do_rjxr) &
                        ovij%i_rjxr_j(:,:,ia0+1:ia0+ndi) &
                            = tij_rjxr(:,:n_ao,:ndi)
                    if (do_rixr) &
                        ovij%i_rixr_j(:,:,ia0+1:ia0+ndi) &
                            = tij_rixr(:,:n_ao,:ndi)
                end if
            
                ia0 = ia0 + bsetDB(ia,iprim)%ndim
                if (do_ij)     tij        = 0.0_real64
                if (do_r)      tij_r      = 0.0_real64
                if (do_r2)     tij_r2     = 0.0_real64
                if (do_p)      tij_p      = 0.0_real64
                if (do_rxp)    tij_rxp    = 0.0_real64
                if (do_rrjxp)  tij_rrjxp  = 0.0_real64
                if (do_rrijxp) tij_rrijxp = 0.0_real64
                if (do_rjxr)   tij_rjxr   = 0.0_real64
                if (do_rixr)   tij_rixr   = 0.0_real64
            end if
        end do
    end do

    if (ondisk) then
        if (do_ij)     close(iu_ij)
        if (do_r)      close(iu_r)
        if (do_r2)     close(iu_r2)
        if (do_p)      close(iu_p)
        if (do_rxp)    close(iu_rxp)
        if (do_rrjxp)  close(iu_rrjxp)
        if (do_rrijxp) close(iu_rrijxp)
        if (do_rjxr)   close(iu_rjxr)
        if (do_rixr)   close(iu_rixr)
    end if

    return

end subroutine overlap_ao_1e

! ======================================================================

subroutine convert_AO2MO_1(n_ao, n_mo, c_ia, q_ao, q_mo, tmp_arr)
    !! Convert a scalar quantity from atomic to molecular orbitals
    !!
    !! Takes a quantity in atomic orbitals (`q_ao`) to molecular orbitals
    !!   (`q_mo`).
    !! Note: The quantity must be scalar (see convert_AO2MO_N otherwise)
    implicit none
    integer, intent(in) :: n_ao
    !! Number of atomic orbitals
    integer, intent(in) :: n_mo
    !! Number of molecular orbitals
    real(real64), dimension(:,:), intent(in) :: c_ia
    !! Coefficients of MOs in AOs basis (geometry: n_mo, n_ao)
    real(real64), dimension(:,:), intent(in) :: q_ao
    !! Quantity in AO basis
    real(real64), dimension(:,:), intent(out) :: q_mo
    !! Quantity in MO basis
    real(real64), dimension(:,:) :: tmp_arr
    !! Temporary array

    integer :: a, b, i, j

    !$omp parallel do collapse(2)
    do j = 1, n_mo
        do a = 1, n_ao
            tmp_arr(a,j) = 0.0_real64
            do b = 1, n_ao
                tmp_arr(a,j) = tmp_arr(a,j) + c_ia(j,b)*q_ao(a,b)
            end do
        end do
    end do
    !$omp end parallel do
    
    !$omp parallel do collapse(2)
    do i = 1, n_mo
        do j = 1, n_mo
            q_mo(i,j) = 0.0_real64
            do a = 1, n_ao
                q_mo(i,j) = q_mo(i,j) + c_ia(i,a)*tmp_arr(a,j)
            end do
        end do
    end do
    !$omp end parallel do

end subroutine convert_AO2MO_1

! ======================================================================

subroutine convert_AO2MO_N(n_ao, n_mo, c_ia, q_ao, q_mo, tmp_arr)
    !! Convert a vector from atomic to molecular orbitals
    !!
    !! Takes a quantity in atomic orbitals (`q_ao`) to molecular orbitals
    !!   (`q_mo`).
    !! Note: The quantity is expected to be a vector.
    implicit none
    integer, intent(in) :: n_ao
    !! Number of atomic orbitals
    integer, intent(in) :: n_mo
    !! Number of molecular orbitals
    real(real64), dimension(:,:), intent(in) :: c_ia
    !! Coefficients of MOs in AOs basis (geometry: n_mo, n_ao)
    real(real64), dimension(:,:,:), intent(in) :: q_ao
    !! Quantity in AO basis
    real(real64), dimension(:,:,:), intent(out) :: q_mo
    !! Quantity in MO basis
    real(real64), dimension(:,:,:) :: tmp_arr
    !! Temporary array

    integer :: a, b, i, j

    !$omp parallel do collapse(2)
    do j = 1, n_mo
        do a = 1, n_ao
            tmp_arr(:,a,j) = 0.0_real64
            do b = 1, n_ao
                tmp_arr(:,a,j) = tmp_arr(:,a,j) + c_ia(j,b)*q_ao(:,a,b)
            end do
        end do
    end do
    !$omp end parallel do

    !$omp parallel do collapse(2)
    do i = 1, n_mo
        do j = 1, n_mo
            q_mo(:,i,j) = 0.0_real64
            do a = 1, n_ao
                q_mo(:,i,j) = q_mo(:,i,j) + c_ia(i,a)*tmp_arr(:,a,j)
            end do
        end do
    end do
    !$omp end parallel do

end subroutine convert_AO2MO_N

! ======================================================================

function eltrans_amp(n_ab, n_ao, n_mos, ovlp_ao, trans_el_dens, tmp_arr, &
                     to_MO_, c_ia_)
    !! Build and return transition amplitudes in the AO or MO basis.
    !!
    !! Builds the transition amplitudes in the atomic orbital basis
    !!   from the electronic transition density matrix of a given
    !!   transition.  The amplitudes are returned as a array.
    !! Default is to return in AO basis.
    implicit none
    integer, intent(in) :: n_ab
    !! Number of unique MO sets (1 for closed-shell, 2 for open-shell)
    integer, intent(in) :: n_ao
    !! Number of atomic orbitals
    integer, dimension(:), intent(in) :: n_mos
    !! Number of molecular orbitals.
    real(real64), dimension(n_ao, n_ao), intent(in) :: ovlp_ao
    !! Overlap integrals between atomic orbitals.
    real(real64), dimension(n_ao, n_ao, n_ab), intent(in) :: trans_el_dens
    !! Transition electronic density matrix
    real(real64), dimension(n_ao, n_ao):: tmp_arr
    !! Temporary array
    logical, intent(in), optional :: to_MO_
    !! Convert transition amplitudes to the MO basis.
    real(real64), dimension(:,:,:), intent(in), optional :: c_ia_
    !! Coefficients of MOs in AOs basis
    real(real64), dimension(:,:,:), allocatable, target :: eltrans_amp
    !! Electronic transition amplitudes

    integer :: i, iab, j, k, l, max_nmo
    real(real64), dimension(:,:,:), allocatable, target :: tmp
    real(real64), dimension(:,:,:), pointer :: amp_ao => null()
    logical :: to_MO


    if (present(to_MO_)) then
        to_MO = to_MO_
    else
        to_MO = .False.
    end if

    if (to_MO) then
        max_nmo = maxval(n_mos(:n_ab))
        allocate(eltrans_amp(max_nmo,max_nmo,n_ab), tmp(n_ao,n_ao,n_ab))
        amp_ao => tmp
    else
        allocate(eltrans_amp(n_ao,n_ao,n_ab))
        amp_ao => eltrans_amp
    end if

    do iab = 1, n_ab
        !$omp parallel do collapse(2)
        do i = 1, n_ao
            do k = 1, n_ao
                tmp_arr(k,i) = 0.0_real64
                do l = 1, n_ao
                    tmp_arr(k,i) = tmp_arr(k,i) + &
                        trans_el_dens(l,k,iab)*ovlp_ao(l,i)
                end do
            end do
        end do
        !$omp end parallel do

        !$omp parallel do collapse(2)
        do i = 1, n_ao
            do j = 1, n_ao
                amp_ao(j,i,iab) = 0.0_real64
                do k = 1, n_ao
                    amp_ao(j,i,iab) = amp_ao(j,i,iab) + &
                        ovlp_ao(j,k)*tmp_arr(k,i)
                end do
            end do
        end do
        !$omp end parallel do

        if (to_MO) then
            if (.not.present(c_ia_)) then
                write(iu_out, '(a)') &
                    'DEVERR: Missing c_ia for the conversion from AO to MO'
                stop
            end if
            call convert_AO2MO(n_ao, n_mos(iab), c_ia_(:,:,iab), tmp(:,:,iab), &
                               eltrans_amp(:,:,iab), tmp_arr)
        end if
    end do
    if (to_MO) deallocate(tmp)

end function eltrans_amp

! ======================================================================

end module electronic
