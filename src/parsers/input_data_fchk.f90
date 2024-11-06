submodule (input:input_data) input_data_fchk
    !! Submodule containing the definition of procedures related to the
    !! data extraction.
    use numeric, only: f2
    use arrays, only: symm_tri_array
    use parsefchk, only: fchkdata, fchkparser
    use basisset, only: build_bset_DB
    use atominfo, only: atdata
    use propinfo, only: load_property_info
    use exception, only: BaseException, Error, InitError, &
        RaiseError, RaiseFileError, RaiseQuantityError

    implicit none

contains

! ======================================================================
! MODULE COMPONENTS
! ======================================================================

module procedure build_mol_data_fchk
    !! Builds molecular specs. from Gaussian formatted checkpoint file.
    !!
    !! Parses a Gaussian formatted checkpoint file and extracts basic
    !! molecular specification data.
    !! An instance of MoleculeDB is returned.

    character(len=42), dimension(9) :: fchk_keys = [ &
        'Number of atoms                           ', &  ! 1
        'Number of electrons                       ', &  ! 2
        'Charge                                    ', &  ! 3
        'Multiplicity                              ', &  ! 4
        'Current cartesian coordinates             ', &  ! 5
        'Atomic numbers                            ', &  ! 6
        'Nuclear charges                           ', &  ! 7
        'Real atomic weights                       ', &  ! 8
        'Total Energy                              '  &  ! 9
    ]

    integer :: ia
    logical :: ok
    type(fchkparser) :: dfchk
    type(fchkdata), dimension(:), allocatable :: dbase

    dfile%error = InitError()

    dfchk = fchkparser(dfile%name)
    dbase = dfchk%read(fchk_keys)

    ok = dfchk%close()
    if (.not. ok) then
        call RaiseFileError(dfile%error, dfile%name, 'closing')
        return
    end if

    ! Basic molecular information
    ! Start with standard dimensions that may be used for different applications
    mol%n_at = dbase(1)%idata(1)
    mol%n_el = dbase(2)%idata(1)

    ! Basic molecular specifications
    mol%charge = dbase(3)%idata(1)
    mol%multip = dbase(4)%idata(1)
    allocate(mol%at_chg(mol%n_at), mol%at_lab(mol%n_at), mol%at_num(mol%n_at))
    do ia = 1, mol%n_at
        mol%at_num(ia) = dbase(6)%idata(ia)
        mol%at_chg(ia) = dbase(7)%rdata(ia)
        mol%at_lab(ia) = atdata(mol%at_num(ia))%symbol
    end do
    mol%at_mas = dbase(8)%rdata
    mol%at_crd = reshape(dbase(5)%rdata, [3, mol%n_at])
    deallocate(dbase(5)%rdata, dbase(6)%idata, dbase(7)%rdata, dbase(8)%rdata)

    ! Other basic quantities
    ! -- energy
    mol%energy = dbase(9)%rdata(1)

    mol%loaded = .true.

    deallocate(dbase)

end procedure build_mol_data_fchk

! ======================================================================

module procedure build_bset_data_fchk
    !! Builds basis set information from Gaussian formatted checkpoint file.
    !!
    !! Parses a Gaussian formatted checkpoint file and extracts basic
    !! set information data.
    !! An instance of BasisSetDB is returned.

    character(len=42), dimension(12) :: fchk_keys = [ &
        'Number of atoms                           ', &  !  1
        'Number of basis functions                 ', &  !  2
        'Number of independent functions           ', &  !  3
        'Number of contracted shells               ', &  !  4
        'Shell types                               ', &  !  5
        'Number of primitives per shell            ', &  !  6
        'Shell to atom map                         ', &  !  7
        'Primitive exponents                       ', &  !  8
        'Contraction coefficients                  ', &  !  9
        'P(S=P) Contraction coefficients           ', &  ! 10
        'Pure/Cartesian d shells                   ', &  ! 11
        'Pure/Cartesian f shells                   '  &  ! 12
    ]

    integer :: n_at
    integer, dimension(:), allocatable :: &
        shell_types, &   ! shell types
        prim_per_sh, &   ! number of primitives per shell
        shell_to_at      ! shell to atom mapping
    real(realwp), dimension(:), allocatable :: &
        coef_contr, &    ! contraction coefficients
        coef_contrSP, &  ! contraction coefficients (P(S=P))
        prim_exp         ! primitives exponents
    logical :: ok
    character(len=512) :: errmsg
    type(fchkparser) :: dfchk
    type(fchkdata), dimension(:), allocatable :: dbase
    class(BaseException), allocatable :: suberr

    dfile%error = InitError()

    dfchk = fchkparser(dfile%name)
    dbase = dfchk%read(fchk_keys)

    ok = dfchk%close()
    if (.not. ok) then
        call RaiseFileError(dfile%error, dfile%name, 'closing')
        return
    end if

    ! Dimensions
    n_at = dbase(1)%idata(1)
    bset%n_basis = dbase(2)%idata(1)
    bset%n_basok = dbase(3)%idata(1)
    bset%n_shells = dbase(4)%idata(1)

    ! Basis function and shells
    bset%pureD = dbase(11)%idata(1) == 0
    bset%pureF = dbase(12)%idata(1) == 0
    shell_types = dbase(5)%idata
    prim_per_sh = dbase(6)%idata
    shell_to_at = dbase(7)%idata
    deallocate(dbase(5)%idata, dbase(6)%idata, dbase(7)%idata)
    prim_exp = dbase(8)%rdata
    deallocate(dbase(8)%rdata)
    coef_contr = dbase(9)%rdata
    deallocate(dbase(9)%rdata)
    if (dbase(10)%dtype /= '0') then
        coef_contrSP = dbase(10)%rdata
        deallocate(dbase(10)%rdata)
    end if

    ! Build DB of basis functions
    call build_bset_DB(n_at, bset%n_shells, bset%pureD, bset%pureF, &
                       shell_types, prim_per_sh, shell_to_at, coef_contr, &
                       coef_contrSP, prim_exp, bset%nprim_per_at, bset%info, &
                       suberr)
    if (suberr%raised()) then
        select type(suberr)
            class is (Error)
                write(errmsg, '(a,a,"Original error:", a)') &
                    'Error when parsing the basis set', new_line(' '), &
                    trim(suberr%msg())
                call RaiseError(dfile%error, errmsg)
                return
            class default
                call RaiseError(dfile%error, 'Generic error')
                return
        end select
    end if
    bset%loaded = .true.

    deallocate(dbase)

end procedure build_bset_data_fchk

! ======================================================================

module procedure build_orb_data_fchk
    !! Builds molecular specs. from Gaussian formatted checkpoint file.
    !!
    !! Parses a Gaussian formatted checkpoint file and extracts basic
    !! molecular specification data.
    !! An instance of OrbitalsDB is returned.

    character(len=42), dimension(7) :: fchk_keys = [ &
        'Number of alpha electrons                 ', &  ! 1
        'Number of beta electrons                  ', &  ! 2
        'Number of basis functions                 ', &  ! 3
        'Alpha MO coefficients                     ', &  ! 4
        'Beta MO coefficients                      ', &  ! 5
        'Alpha Orbital Energies                    ', &  ! 6
        'Beta Orbital Energies                     '  &  ! 7
    ]

    integer :: i, ij, j
    logical :: ok
    type(fchkparser) :: dfchk
    type(fchkdata), dimension(:), allocatable :: dbase

    dfile%error = InitError()

    dfchk = fchkparser(dfile%name)
    dbase = dfchk%read(fchk_keys)

    ok = dfchk%close()
    if (.not. ok) then
        call RaiseFileError(dfile%error, dfile%name, 'closing')
        return
    end if

    ! Dimensions
    ! -- Number of electrons
    orb%n_els(1) = dbase(1)%idata(1)
    orb%n_els(2) = dbase(2)%idata(1)
    ! -- Number of atomic orbitals
    orb%n_ao = dbase(3)%idata(1)

    orb%n_mos(1) = dbase(6)%len
    if (dbase(7)%dtype /= '0') then
        orb%n_ab = 2
        orb%openshell = .true.
        orb%n_mos(2) = dbase(7)%len
    else
        orb%n_ab = 1
        orb%openshell = .false.
    end if
    orb%n_mo = maxval(orb%n_mos(:orb%n_ab))
    ! -- Save orbital energies and clear redundant memory
    allocate(orb%en_mos(orb%n_mo,orb%n_ab))
    orb%en_mos(:,1) = dbase(6)%rdata
    deallocate(dbase(6)%rdata)
    if (orb%openshell) then
        orb%en_mos(:,2) = dbase(7)%rdata
        deallocate(dbase(7)%rdata)
    end if
    ! -- Save orbital coefficients and clear redundant memory
    allocate(orb%coef_mos(orb%n_mo,orb%n_ao,orb%n_ab))
    ij = 0
    do i = 1, orb%n_mos(1)
        do j = 1, orb%n_ao
            ij = ij + 1
            orb%coef_mos(i,j,1) = dbase(4)%rdata(ij)
        end do
    end do
    deallocate(dbase(4)%rdata)
    if (orb%openshell) then
        ij = 0
        do i = 1, orb%n_mos(2)
            do j = 1, orb%n_ao
                ij = ij + 1
                orb%coef_mos(i,j,2) = dbase(5)%rdata(ij)
            end do
        end do
        deallocate(dbase(5)%rdata)
    end if
    
    orb%loaded = .true.

    deallocate(dbase)

end procedure build_orb_data_fchk

! ======================================================================

module procedure build_exc_data_fchk
    !! Builds excited-state data from Gaussian formatted chk file.
    !!
    !! Parses a Gaussian formatted checkpoint file and extracts all
    !! relevant information on electronic excited-state data and
    !! transition data.
    !! An instance of ExcitationDB is returned.

    character(len=42), dimension(12), parameter :: fchk_keys = [ &
        'Number of basis functions            ', &  !  1.
        'Beta Orbital Energies                ', &  !  2.
        'Total CI Density                     ', &  !  3.
        'ETran scalars                        ', &  !  4.
        'ETran spin                           ', &  !  5.
        'ETran sym                            ', &  !  6.
        'ETran state values                   ', &  !  7.
        'Excited state NLR                    ', &  !  8.
        'Number of g2e trans dens             ', &  !  9.
        'G to E trans densities               ', &  ! 10.
        'Excited state densities              ', &  ! 11.
        'SCF Energy                           '  &  ! 12.
    ]

    integer :: &
        i,  i_ab, ioff, n_ab, nbas_lt, n_basis, &
        lblock_ETran, &  ! Length of a data block for a given state
        NLR              ! indicate if left and right trans. matrix data stored
    logical :: ok
    type(fchkparser) :: dfchk
    type(fchkdata), dimension(:), allocatable :: dbase

    dfile%error = InitError()

    dfchk = fchkparser(dfile%name)
    dbase = dfchk%read(fchk_keys)

    ok = dfchk%close()
    if (.not. ok) then
        call RaiseFileError(dfile%error, dfile%name, 'closing')
        return
    end if

    ! Check if excited-states properties available
    ! First check that this exists, since a fchk may be missing it.
    if (dbase(4)%dtype == '0') then
        call RaiseQuantityError(dfile%error, &
            msg='Missing excited-states data in file.')
        return
    end if

    ! -- Basic information
    ! Number of basis functions
    n_basis = dbase(1)%idata(1)
    if (dbase(2)%dtype /= '0') then
        n_ab = 2
    else
        n_ab = 1
    end if

    ! First check that NLR == 1. We do not support the other case for now
    NLR = dbase(8)%idata(1)
    if (NLR > 1) then
        call RaiseError(dfile%error, 'NLR /= 1 not yet supported.  Sorry.')
        return
    end if

    ! Now let us process the information of interest in the file.
    ! The "scalars" array provides information to parse "ETran state values"
    exc%n_states = dbase(4)%idata(1)
    exc%id_state = dbase(4)%idata(5)
    lblock_ETran = dbase(4)%idata(2)

    ! Extract excited-states data
    exc%ispin_exc = dbase(5)%idata

    ! Extract transition data values
    allocate(exc%g2e_energy(exc%n_states), exc%exc_energy(exc%n_states), &
             exc%g2e_eldip(3,exc%n_states), exc%g2e_magdip(3,exc%n_states))
    exc%gs_energy = dbase(12)%rdata(1)
    do i = 1, exc%n_states
        ioff = (i-1)*lblock_ETran
        exc%exc_energy(i) = dbase(7)%rdata(1+ioff)
        exc%g2e_energy(i) = exc%exc_energy(i) - exc%gs_energy
        exc%g2e_eldip(:,i) = dbase(7)%rdata(2+ioff:4+ioff)
        exc%g2e_magdip(:,i) = dbase(7)%rdata(8+ioff:10+ioff)
    end do

    ! Now extract the ground to excited transition moments
    ! structure: n_basis, n_basis, alpha/beta, n_states
    allocate(exc%g2e_dens(n_basis,n_basis,2,exc%n_states), &
             exc%exc_dens(n_basis,n_basis,n_ab,exc%n_states))
    ! WARNING: We need to be careful here!
    ! Gaussian up to G16 included saved the transition densities with embedded
    ! coefficients, sqrt(2) for TD and it seems 1/sqrt(2) for CI.
    ! We need to do the inverse operation to get the correct coefficients.
    if (dfile%check_version(major='G16')) then
        exc%g2e_dens = reshape(dbase(10)%rdata, &
                               [n_basis,n_basis,2,exc%n_states]) / sqrt(f2)
    else
        exc%g2e_dens = reshape(dbase(10)%rdata, &
                               [n_basis,n_basis,2,exc%n_states])
    end if

    ! Excited-state densities are stored in lower-triangular forms.
    ! We need to unpack them.
    ! They are always stored for both a and b, but in case of closed-shell.
    nbas_lt = n_basis*(n_basis+1)/2
    do i = 1, exc%n_states
        ioff = (i-1)*2*nbas_lt
        do i_ab = 1, n_ab
            call symm_tri_array(n_basis, &
                dbase(11)%rdata(ioff+(i_ab-1)*nbas_lt+1:), linear=.true., &
                lower=.true., anti_symm=.false., arr_new=exc%exc_dens(:,:,i_ab,i))
        end do
    end do

    exc%dens_loaded = .true.
    exc%prop_loaded = .true.

    deallocate(dbase)

end procedure build_exc_data_fchk

! ======================================================================

module procedure build_vib_data_fchk
    !! Builds vibrational data from Gaussian formatted chk file.
    !!
    !! Parses a Gaussian formatted checkpoint file and extracts all
    !! relevant information related to vibrations.
    !! An instance of VibrationsDB is returned.
    !!
    !! See [[build_vibdata(proc)]] for details.

    character(len=42), dimension(4), parameter :: fchk_keys = [ &
        'Number of Normal Modes               ', &  !  1.
        'Vib-AtMass                           ', &  !  2.
        'Vib-E2                               ', &  !  3.
        'Vib-Modes                            '  &  !  4.
    ]

    integer :: i, ia, n_at3
    real(realwp), dimension(:,:), allocatable :: evec
    logical :: build_Lmat, build_Lmweig, ok
    type(fchkparser) :: dfchk
    type(fchkdata), dimension(:), allocatable :: dbase

    dfile%error = InitError()

    if (present(get_Lmat)) then
        build_Lmat = get_Lmat
    else
        build_Lmat = .true.
    end if

    if (present(get_Lmweig)) then
        build_Lmweig = get_Lmweig
    else
        build_Lmweig = .true.
    end if

    dfchk = fchkparser(dfile%name)
    dbase = dfchk%read(fchk_keys)

    ok = dfchk%close()
    if (.not. ok) then
        call RaiseFileError(dfile%error, dfile%name, 'closing')
        return
    end if

    ! First check that this exists, since a fchk may be missing it.
    if (dbase(1)%dtype == '0') then
        call RaiseQuantityError(dfile%error, &
            msg='Missing vibrational data in file.')
        return
    end if

    ! Get number of normal modes, necessary for the rest of the parsing
    vib%n_vib = dbase(1)%idata(1)
    ! Get number of Cartesian coordinates
    n_at3 = 3*size(dbase(2)%rdata)

    ! Extract frequencies and reduced mass
    ! They are stored in Vib-E2, starting with freq, then red. mass...
    vib%freq = dbase(3)%rdata(1:vib%n_vib)
    vib%red_mass = dbase(3)%rdata(vib%n_vib+1:2*vib%n_vib)

    ! Now extract eigenvectors and do necessary operations/conversions
    evec = reshape(dbase(4)%rdata, [n_at3,vib%n_vib])
    if (build_Lmweig) then
        allocate(vib%L_mwg(n_at3,vib%n_vib))
        do i = 1, vib%n_vib
            vib%L_mwg(:,i) = evec(:,i)/sqrt(vib%red_mass(i))
        end do
    end if

    if (build_Lmat) then
        allocate(vib%L_mat(n_at3,vib%n_vib))
        do i = 1, vib%n_vib
            do ia = 1, n_at3
                vib%L_mat(ia,i) = evec(ia,i)*sqrt(dbase(2)%rdata((ia+2)/3)) &
                    /sqrt(vib%red_mass(i))
            end do
            vib%L_mat(:,i) = vib%L_mat(:,i) / sqrt(sum(vib%L_mat(:,i)**2))
        end do
    end if

    deallocate(dbase, evec)

end procedure build_vib_data_fchk

! ======================================================================

module procedure get_data_from_id_fchk
    !! Extract data for a specific quantity from Gaussian fchk file.
    !!
    !! Extracts data for a specific quantity from the Gaussian fchk file.
    !! The procedure also takes care of doing necessary adjustments if
    !! needed.
    integer :: der_ord, i, ioff, LP
    real(realwp), dimension(:), allocatable :: tmpvec
    logical :: ok
    character(len=42), dimension(:), allocatable :: fchk_keys
    type(fchkparser) :: dfchk
    type(fchkdata), dimension(:), allocatable :: dbase

    ! Check if reference/starting state makes sense and set it in DB
    if (present(start_state)) then
        if (start_state < 0) then
            prop%istat = 99
            return
        else
            prop%states(1) = start_state
        end if
    else
        prop%states(2) = 0
    end if

    ! Check if end state provided and consistent with start state
    if (present(end_state)) then
        if (end_state == -1) then
            prop%states(2) = -1
        else if (end_state < start_state) then
            prop%istat = 99
            return
        else
            prop%states(2) = end_state
        end if
    else
        prop%states(2) = prop%states(1)
    end if

    ! Check derivative order.
    ! Since does not necessarily makes sense for all properties, we
    ! store it in temporary variables and will use later.
    if (present(derorder)) then
        if (derorder < 0) then
            prop%istat = 99
            return
        else
            der_ord = derorder
        end if
    else
        der_ord = 0
    end if

    ! load basic information about property
    call load_property_info(prop, identifier)

    ! Load parsing capabilities
    dfchk = fchkparser(dfile%name)

    ! Now extract information
    select case(identifier)
    case(101)
        prop%order = der_ord
        if (prop%states(1) /= prop%states(2)) then
            ! Electronic transition moment
            if (prop%states(1) == 0) then
                fchk_keys = [ &
                    'ETran scalars                        ', &  !  1.
                    'ETran state values                   '  &  !  2.
                ]
                dbase = dfchk%read(fchk_keys)
                if (dbase(1)%dtype == '0') then
                    prop%istat = 2
                    return
                end if
                LP = 3  ! we could use prop%pdim but not much sense.
                if (prop%order == 0) then
                    block
                    integer :: n_states, lblock
                    n_states = dbase(1)%idata(1)
                    lblock = dbase(1)%idata(2)
                    if (n_states < prop%states(2)) then
                        prop%istat = 2
                        return
                    end if
                    if (prop%states(2) == -1) then
                        allocate(prop%data(LP*n_states))
                        do i = 1, n_states
                            ioff = (i-1)*lblock
                            prop%data((i-1)*LP+1:i*LP) = &
                                dbase(2)%rdata(2+ioff:4+ioff)
                        end do
                        prop%loaded = .true.
                        prop%shape = [LP, n_states]
                        prop%dim_shape = [1, 1]
                    else
                        allocate(prop%data(LP))
                        ioff = (prop%states(2)-1)*lblock
                        prop%data = dbase(2)%rdata(2+ioff:4+ioff)
                        prop%loaded = .true.
                        prop%shape = [LP]
                        prop%dim_shape = [1]
                    end if
                    end block
                else if (prop%order == 1) then
                    block
                    integer :: exc_state, lblock, n_at3, n_LR, n_states
                    n_states = dbase(1)%idata(1)
                    exc_state = dbase(1)%idata(5)
                    lblock = dbase(1)%idata(2)
                    n_LR = dbase(1)%idata(3)
                    n_at3 = dbase(1)%idata(6) - 3  ! 3*natoms + 3 (elec. field)
                    if (exc_state /= prop%states(2) .and. &
                            prop%states(2) /= -1) then
                        prop%istat = 2
                        return
                    end if
                    allocate(prop%data(LP*n_at3))
                    ioff = lblock*n_states*n_LR + 3*lblock + 2
                    do i = 0, n_at3 - 1
                        prop%data(i*LP+1:i*LP+3) = &
                            dbase(2)%rdata(ioff+i*lblock:ioff+i*lblock+2)
                    end do
                    prop%loaded = .true.
                    prop%shape = [LP, n_at3]
                    prop%dim_shape = [1, 1]
                    end block
                else
                    prop%istat = 1
                end if
            else
                prop%istat = 1
            end if
        else
            ! Reference state
            stop 'El. dipole NYI'
        end if
    case(111)
        prop%order = der_ord
        if (prop%states(1) /= prop%states(2)) then
            ! Electronic transition moment
            if (prop%states(1) == 0) then
                fchk_keys = [ &
                    'ETran scalars                        ', &  !  1.
                    'ETran state values                   '  &  !  2.
                ]
                dbase = dfchk%read(fchk_keys)
                if (dbase(1)%dtype == '0') then
                    prop%istat = 2
                    return
                end if
                LP = 3  ! we could use prop%pdim but not much sense.
                if (prop%order == 0) then
                    block
                    integer :: n_states, lblock
                    n_states = dbase(1)%idata(1)
                    lblock = dbase(1)%idata(2)
                    if (n_states < prop%states(2)) then
                        prop%istat = 2
                        return
                    end if
                    if (prop%states(2) == -1) then
                        allocate(prop%data(LP*n_states))
                        do i = 1, n_states
                            ioff = (i-1)*lblock
                            prop%data((i-1)*LP+1:i*LP) = &
                                dbase(2)%rdata(5+ioff:7+ioff)
                        end do
                        prop%loaded = .true.
                        prop%shape = [LP, n_states]
                        prop%dim_shape = [1, 1]
                    else
                        allocate(prop%data(LP))
                        ioff = (prop%states(2)-1)*lblock
                        prop%data = dbase(2)%rdata(5+ioff:7+ioff)
                        prop%loaded = .true.
                        prop%shape = [LP]
                        prop%dim_shape = [1]
                    end if
                    end block
                else if (prop%order == 1) then
                    block
                    integer :: exc_state, lblock, n_at3, n_LR, n_states
                    n_states = dbase(1)%idata(1)
                    exc_state = dbase(1)%idata(5)
                    lblock = dbase(1)%idata(2)
                    n_LR = dbase(1)%idata(3)
                    n_at3 = dbase(1)%idata(6) - 3  ! 3*natoms + 3 (elec. field)
                    if (exc_state /= prop%states(2) .and. &
                            prop%states(2) /= -1) then
                        prop%istat = 2
                        return
                    end if
                    allocate(prop%data(LP*n_at3))
                    ioff = lblock*n_states*n_LR + 3*lblock + 5
                    do i = 0, n_at3 - 1
                        prop%data(i*LP+1:i*LP+3) = &
                            dbase(2)%rdata(ioff+i*lblock:ioff+i*lblock+2)
                    end do
                    prop%loaded = .true.
                    prop%shape = [LP, n_at3]
                    prop%dim_shape = [1, 1]
                    end block
                else
                    prop%istat = 1
                end if
            else
                prop%istat = 1
            end if
        else
            ! Reference state
            stop 'El. dipole (vel.) NYI'
        end if
    case(102)
        prop%order = der_ord
        if (prop%states(1) /= prop%states(2)) then
            ! Electronic transition moment
            if (prop%states(1) == 0) then
                fchk_keys = [ &
                    'ETran scalars                        ', &  !  1.
                    'ETran state values                   '  &  !  2.
                ]
                dbase = dfchk%read(fchk_keys)
                if (dbase(1)%dtype == '0') then
                    prop%istat = 2
                    return
                end if
                LP = 3  ! we could use prop%pdim but not much sense.
                if (prop%order == 0) then
                    block
                    integer :: n_states, lblock
                    n_states = dbase(1)%idata(1)
                    lblock = dbase(1)%idata(2)
                    if (n_states < prop%states(2)) then
                        prop%istat = 2
                        return
                    end if
                    if (prop%states(2) == -1) then
                        allocate(prop%data(LP*n_states))
                        do i = 1, n_states
                            ioff = (i-1)*lblock
                            prop%data((i-1)*LP+1:i*LP) = &
                                dbase(2)%rdata(8+ioff:10+ioff)
                        end do
                        prop%loaded = .true.
                        prop%shape = [LP, n_states]
                        prop%dim_shape = [1, 1]
                    else
                        allocate(prop%data(LP))
                        ioff = (prop%states(2)-1)*lblock
                        prop%data = dbase(2)%rdata(8+ioff:10+ioff)
                        prop%loaded = .true.
                        prop%shape = [LP]
                        prop%dim_shape = [1, 1]
                    end if
                    end block
                else if (prop%order == 1) then
                    block
                    integer :: exc_state, lblock, n_at3, n_LR, n_states
                    n_states = dbase(1)%idata(1)
                    exc_state = dbase(1)%idata(5)
                    lblock = dbase(1)%idata(2)
                    n_LR = dbase(1)%idata(3)
                    n_at3 = dbase(1)%idata(6) - 3  ! 3*natoms + 3 (elec. field)
                    if (exc_state /= prop%states(2) .and. &
                            prop%states(2) /= -1) then
                        prop%istat = 2
                        return
                    end if
                    allocate(prop%data(LP*n_at3))
                    ioff = lblock*n_states*n_LR + 3*lblock + 8
                    do i = 0, n_at3 - 1
                        prop%data(i*LP+1:i*LP+3) = &
                            dbase(2)%rdata(ioff+i*lblock:ioff+i*lblock+2)
                    end do
                    prop%loaded = .true.
                    prop%shape = [LP, n_at3]
                    prop%dim_shape = [1, 1]
                    end block
                else
                    prop%istat = 1
                end if
            else
                prop%istat = 1
            end if
        else
            ! Reference state
            stop 'Mag. dipole NYI'
        end if
    case(107)
        prop%order = der_ord
        if (prop%states(1) /= prop%states(2)) then
            ! Electronic transition moment
            if (prop%states(1) == 0) then
                fchk_keys = [ &
                    'ETran scalars                        ', &  !  1.
                    'ETran state values                   '  &  !  2.
                ]
                dbase = dfchk%read(fchk_keys)
                if (dbase(1)%dtype == '0') then
                    prop%istat = 2
                    return
                end if
                LP = 9  ! length stored in memeory.
                if (prop%order == 0) then
                    block
                    integer :: n_states, lblock
                    allocate(tmpvec(LP))
                    n_states = dbase(1)%idata(1)
                    lblock = dbase(1)%idata(2)
                    if (n_states < prop%states(2)) then
                        prop%istat = 2
                        return
                    end if
                    if (prop%states(2) == -1) then
                        allocate(prop%data(LP*n_states))
                        do i = 1, n_states
                            ioff = (i-1)*lblock
                            tmpvec = dbase(2)%rdata(11+ioff:16+ioff)
                            prop%data((i-1)*LP+1:i*LP) = &
                                reshape(elquad_LT_to_2D(tmpvec), [LP])
                        end do
                        prop%loaded = .true.
                        prop%shape = [prop%pdim, n_states]
                        prop%dim_shape = [1, 1, 1]
                    else
                        allocate(prop%data(LP))
                        ioff = (prop%states(2)-1)*lblock
                        tmpvec = dbase(2)%rdata(11+ioff:16+ioff)
                        prop%data = reshape(elquad_LT_to_2D(tmpvec), [LP])
                        prop%loaded = .true.
                        prop%shape = prop%pdim
                        prop%dim_shape = [1, 1]
                    end if
                    end block
                else if (prop%order == 1) then
                    block
                    integer :: exc_state, lblock, n_at3, n_LR, n_states
                    n_states = dbase(1)%idata(1)
                    exc_state = dbase(1)%idata(5)
                    lblock = dbase(1)%idata(2)
                    n_LR = dbase(1)%idata(3)
                    n_at3 = dbase(1)%idata(6) - 3  ! 3*natoms + 3 (elec. field)
                    if (exc_state /= prop%states(2) .and. &
                            prop%states(2) /= -1) then
                        prop%istat = 2
                        return
                    end if
                    allocate(prop%data(LP*n_at3))
                    ioff = lblock*n_states*n_LR + 3*lblock + 11
                    do i = 0, n_at3 - 1
                        tmpvec = dbase(2)%rdata(ioff+i*lblock:ioff+i*lblock+5)
                        prop%data(i*LP+1:i*LP+9) = &
                            reshape(elquad_LT_to_2D(tmpvec), [LP])
                    end do
                    prop%loaded = .true.
                    prop%shape = [prop%pdim, n_at3]
                    prop%dim_shape = [1, 1, 1]
                    end block
                else
                    prop%istat = 1
                end if
            else
                prop%istat = 1
            end if
        else
            ! Reference state
            stop 'El. quadrupole NYI'
        end if
    case default
        stop 'Unsupported property'
    end select
    ok = dfchk%close()
    if (.not. ok) then
        call RaiseFileError(dfile%error, dfile%name, 'closing')
        return
    end if

end procedure get_data_from_id_fchk

! ======================================================================

module procedure get_data_from_tag_fchk
    !! Extract data for a specific quantity from Gaussian fchk file.
    !!
    !! Extracts data for a specific quantity from the Gaussian fchk file.
    !! The procedure also takes care of doing necessary adjustments if
    !! needed.
    integer :: der_ord, i, ioff, LP
    real(realwp), dimension(:), allocatable :: tmpvec
    logical :: ok
    character(len=42), dimension(:), allocatable :: fchk_keys
    type(fchkparser) :: dfchk
    type(fchkdata), dimension(:), allocatable :: dbase

    ! Check if reference/starting state makes sense and set it in DB
    if (present(start_state)) then
        if (start_state < 0) then
            prop%istat = 99
            return
        else
            prop%states(1) = start_state
        end if
    else
        prop%states(2) = 0
    end if

    ! Check if end state provided and consistent with start state
    if (present(end_state)) then
        if (end_state == -1) then
            prop%states(2) = -1
        else if (end_state < start_state) then
            prop%istat = 99
            return
        else
            prop%states(2) = end_state
        end if
    else
        prop%states(2) = prop%states(1)
    end if

    ! Check derivative order.
    ! Since does not necessarily makes sense for all properties, we
    ! store it in temporary variables and will use later.
    if (present(derorder)) then
        if (derorder < 0) then
            prop%istat = 99
            return
        else
            der_ord = derorder
        end if
    else
        der_ord = 0
    end if

    ! load basic information about property
    call load_property_info(prop, name, tag)

    ! Load parsing capabilities
    dfchk = fchkparser(dfile%name)

    ! Now extract information
    select case(name)
    case default
        stop 'Unsupported property'
    end select
    ok = dfchk%close()
    if (.not. ok) then
        call RaiseFileError(dfile%error, dfile%name, 'closing')
        return
    end if

end procedure get_data_from_tag_fchk

! ======================================================================

function elquad_LT_to_2D(equad_LT) result(equad_2D)
    !! Sort elements of electric quadrupole and returns full 2D tensor.
    !!
    !! Gaussian stores the electric quadrupole in a particular way,
    !! starting from the diagonal elements (trace), in the following
    !! sequence:
    !! XX, YY, ZZ, XY, XZ, ZZ
    !! The function builds a complete, symmetric 2D tensor, reordering
    !! the elements.
    real(realwp), dimension(:), intent(in) :: equad_LT
    !! Electric quadrupole, in lower triangular form as stored by Gaussian.
    real(realwp), dimension(3,3) :: equad_2D
    !! Electric quadrupole, as symmetric tensor
    equad_2D(1,1) = equad_LT(1)
    equad_2D(2,2) = equad_LT(2)
    equad_2D(3,3) = equad_LT(3)
    equad_2D(2,1) = equad_LT(4)
    equad_2D(1,2) = equad_LT(4)
    equad_2D(3,1) = equad_LT(5)
    equad_2D(1,3) = equad_LT(5)
    equad_2D(3,2) = equad_LT(6)
    equad_2D(2,3) = equad_LT(6)
end function elquad_LT_to_2D

end
