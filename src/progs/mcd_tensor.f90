program mcd_tensor
    use iso_fortran_env, only: output_unit
    use numeric, only: realwp, f0, f1, f2
    use string, only: timestamp, labXYZ_1D, locase
    use input, only: DataFile
    use parse_cmdline, only: CmdArgDB
    use output, only: iu_out, sec_header, len_int, prt_coord, prt_mat, &
      write_err
    use exception, only: BaseException, Error, AllocateError, ArgumentError, &
        FileError, ValueError
    use gmcd_legacy, only: write_control, build_MOs
    use basisset, only: fix_norm_AOs, chk_bset_redundancy
    use electronic, only: convert_AO2MO, eltrans_amp, overlap_ao_1e, ovij_1e
    use exc_sos, only: sos_eiOg, sos_ejOei, sos_prefac_ejOei, &
        sos_MCD_tensor_LORG_corr
    use workdata, only: MoleculeDB, BasisSetDB, OrbitalsDB, ExcitationDB

    implicit none

    type DebugType
        logical :: user_set ! user has activated some debugging features.
        logical :: print    ! enable debug printing
        logical :: timer    ! add timer
        logical :: add_ijaa ! add <ia||ja> terms in pfac for exc-exc integrals
    end type DebugType

    integer :: fstate, i, iab, istate, ispin_gs, ix, jstate, jx, max_state, &
        model_sos, qty_flag
    real(realwp) :: de, ef, ei
    real(realwp) :: e_gamma = 1.0e-4
    real(realwp), dimension(3) :: &
        r_gg, &      ! < g | r | g >
        p_gg, &      ! < g | p | g >
        rxp_gg, &    ! < g | r x p | g >
        r_fg, &      ! < f | r | g >
        p_fg, &      ! < f | p | g >
        rxp_fg, &    ! < f | r x p | g >
        r_fk, &      ! < f | r | k >
        rxp_fk, &    ! < f | r x p | k >
        r_kg, &      ! < f | r | g >
        rxp_kg, &    ! < f | r x p | g >
        edip_nuc     ! Nuclear contribution to electric dipoles
    real(realwp), dimension(3,2) :: &
        r_gg_ab, &   ! < g | r | g >
        p_gg_ab, &   ! < g | p | g >
        rxp_gg_ab    ! < g | r x p | g >
    real(realwp), dimension(3,3) :: G_if
    real(realwp), dimension(:), allocatable :: &
        ov_eieg, &   ! Inverse of (E_i-E_g)
        ov_eief      ! Inverse of (E_i-E_f)
    real(realwp), dimension(:,:), allocatable :: &
        ov_eiej      ! Inverse of (E_i-E_j)
    real(realwp), dimension(:,:,:), allocatable :: &
        r_lk, &      ! Stores all combinations < l | r | k >
        p_lk         ! Stores all combinations < l | p | k >
    real(realwp), dimension(:,:,:,:), allocatable :: &
        t_mo, &      ! Transition amplitude in MO basis
        pfac_r, &    ! Stores the prefactor for < e_j | r | e_i >
        pfac_p, &    ! Stores the prefactor for < e_j | p | e_i >
        pfac_rxp, &  ! Stores the prefactor for < e_j | r x p | e_i >
        Smo_ipj, &   ! MO-basis integral < i | p | j >
        Smo_irj, &   ! MO-basis integral < i | r | j >
        Smo_irxpj    ! MO-basis integral < i | r x p | j >
    real(realwp), allocatable :: tmp_ao_arr1(:,:), tmp_ao_arrN(:,:,:)
    ! logical, parameter :: DEBUG = .true., TIMEIT = .false.
    logical :: do_giao, do_guvcde, exists, forbid, in_mem, show_timestamp, &
        use_gamma
    character(len=512) :: fmt_elstate
    character(len=:), allocatable :: infile, outfile, datetime
    character(len=*), parameter :: PROGTITLE = 'MCD Tensor Calculator'
    type(DataFile) :: dfile
    class(ovij_1e), allocatable :: ao_int
    class(BaseException), allocatable :: err
    type(DebugType) :: debug
    type(MoleculeDB) :: moldb
    type(BasisSetDB) :: bsetdb
    type(OrbitalsDB) :: orbdb
    type(ExcitationDB) :: excdb

    interface write_param
        procedure write_param_bool, write_param_int, write_param_real, &
            write_param_char
    end interface write_param

    1200 format(/, &
        ' > Number of excited states : ',i0,/, &
        ' > Reference excited state  : ',i0)

    call parse_argopts(infile, outfile, e_gamma, fstate, max_state, &
        use_gamma, do_giao, do_guvcde, in_mem, model_sos, show_timestamp, &
        debug)
    dfile = DataFile(infile)
    if (dfile%has_error()) then
        call write_err('std', 'Error found while initializing data file', &
                       dfile%get_error())
        stop 1
    end if

    if (show_timestamp) then
        datetime = timestamp()
        write(iu_out, '(/," Calculations started on ", a)') datetime
    end if

    call sec_header(1, 'Data on Molecular System')
    if (debug%timer) call write_time('Molecular data')
    moldb = dfile%build_mol_data()
    if (dfile%has_error()) then
        call write_err('std', &
            'Error found while parsing molecular specifications in file', &
            dfile%get_error())
            stop 1
    end if
    bsetdb = dfile%build_bset_data()
    if (dfile%has_error()) then
        call write_err('std', &
            'Error found while parsing basis set data in file', &
            dfile%get_error())
            stop 1
    end if
    orbdb = dfile%build_orb_data()
    if (dfile%has_error()) then
        call write_err('std', &
            'Error found while parsing molecular orbitals data in file', &
            dfile%get_error())
            stop 1
    end if
    excdb = dfile%build_exc_data()
    if (dfile%has_error()) then
        call write_err('std', &
            'Error found while parsing electronic excitations data in file', &
            dfile%get_error())
            stop 1
    end if
    ! if (openshell) then
    !     call write_err('gen', &
    !         'Sorry, open-shell systems not yet supported.  Working on it.' &
    !         )
    !     stop
    ! end if
    call write_moldata(moldb, bsetdb, orbdb)
    ! Post-processing to build additional data
    ispin_gs = moldb%multip - 1

    if (do_guvcde) then
        if (debug%timer) call write_time('GUVCDE write control')
        call write_control(iu_out, moldb, bsetdb, orbdb, err)
        if (err%raised()) then
            select type(err)
                class is (Error)
                    call write_err('std', &
                        'Error while printing ROAAI-like control output', &
                        err%msg() &
                        )
                    stop
                class default
                    call write_err('gen', &
                        'Something went wrong inside write_control.' &
                    )
                    stop
            end select
        end if
    end if

    ! normalization coeffs
    if (debug%timer) call write_time('Check AO normalization')
    call fix_norm_AOs(iu_out, moldb%n_at, orbdb%n_ao, moldb%at_crd, &
                      bsetdb%nprim_per_at, bsetdb%info, err, debug%print)
    if (err%raised()) then
        select type(err)
            class is (Error)
                call write_err('std', &
                    'Error while normalizing the coefficients for the AOs', &
                    err%msg() &
                    )
                stop
            class default
                call write_err('gen', &
                    'Something went wrong inside fix_norm_AOs.' &
                )
                stop
        end select
    end if

    call sec_header(2, 'Computation of one-electron integrals')
    if (debug%timer) call write_time('1-e AO integral')
    ! overlap6
    qty_flag = 2**10 - 1
    ! qty_flag = ibset(qty_flag, 0)
    call overlap_ao_1e(iu_out, moldb%n_at, orbdb%n_ao, qty_flag, do_guvcde, &
                       in_mem, moldb%at_crd, bsetdb%nprim_per_at, &
                       bsetdb%info, ao_int, err)
    if (err%raised()) then
        select type(err)
            class is (Error)
                call write_err('std', &
                    'Error while computing the AO 1-electron integrals', &
                    err%msg() &
                )
                stop
            class default
                call write_err('gen', &
                    'Something went wrong inside overlap_ao_1e.' &
                )
                stop
        end select
    end if
    ! call prt_mat(ao_int%i_j, n_ao, n_ao)

    ! call chk_bset_redundancy(n_ao, ao_int%i_j)

    call sec_header(2, 'Definition of Molecular Orbitals')
    ! convert AO -> MO quantities of interest
    if (do_guvcde) then
        if (debug%timer) call write_time('GUVCDE build MOs')
        call build_MOs(orbdb%n_ab, orbdb%n_ao, orbdb%n_mos, orbdb%coef_mos, &
                       .true.)
    endif
    ! call convert_AO2MO(n_ao, n_mo, coef_mos(1), ao_int%i_j, Smo_ij)
    allocate(Smo_irj(3,orbdb%n_mo,orbdb%n_mo,orbdb%n_ab), &
             Smo_ipj(3,orbdb%n_mo,orbdb%n_mo,orbdb%n_ab), &
             Smo_irxpj(3,orbdb%n_mo,orbdb%n_mo,orbdb%n_ab), &
             tmp_ao_arr1(orbdb%n_ao,orbdb%n_ao), &
             tmp_ao_arrN(3,orbdb%n_ao,orbdb%n_ao))
    if (debug%timer) call write_time('AO to MO conversion')
    do iab = 1, orbdb%n_ab
        if (debug%timer) call write_time('<i|r|j>')
        call convert_AO2MO(orbdb%n_ao, orbdb%n_mos(iab), &
                           orbdb%coef_mos(:,:,iab), ao_int%i_r_j, &
                           Smo_irj(:,:,:,iab), tmp_ao_arrN)
        if (debug%timer) call write_time('<i|p|j>')
        call convert_AO2MO(orbdb%n_ao, orbdb%n_mos(iab), &
                           orbdb%coef_mos(:,:,iab), ao_int%i_p_j, &
                           Smo_ipj(:,:,:,iab), tmp_ao_arrN)
        if (debug%timer) call write_time('<i|rxp|j>')
        call convert_AO2MO(orbdb%n_ao, orbdb%n_mos(iab), &
                           orbdb%coef_mos(:,:,iab), ao_int%i_rxp_j, &
                           Smo_irxpj(:,:,:,iab), tmp_ao_arrN)
    end do
    ! call convert_AO2MO(n_ao, n_mo, coef_mos(1), ao_int%i_j)

    call sec_header(1, 'Transition Data')
    if (debug%timer) call write_time('Transition data')
    excdb = dfile%build_exc_data()
    if (dfile%has_error()) then
        call write_err('std', &
            'Error found while parsing excited-states data in file', &
            dfile%get_error())
            stop 1
    end if
    write(iu_out, 1200) excdb%n_states, excdb%id_state

    call sec_header(2, 'User overrides')
    exists = .false.
    if (fstate > 0 .and. fstate /= excdb%id_state) then
        exists = .true.
        write(iu_out, '(" > New final state, num. ",i0,".")') &
            fstate
        excdb%id_state = fstate
    end if
    if (max_state > excdb%n_states) then
        exists = .true.
        write(iu_out, '(1x,"NOTE: ",a,/,1x,a)') &
            'Requested upper bound exceeds available states.', &
            'Reverting to highest available state.'
    else if (max_state > 0) then
        exists = .true.
        write(iu_out, '(" > Setting highest electronic state to ",i0,".")') &
            max_state
        excdb%n_states = max_state
    end if
    if (.not.exists) &
        write(iu_out, '(1x,a)') &
            'No modifications requested. Proceeding with data found.'

    ! call prt_mat(g2e_dens(:,:,1,1), n_basis, n_basis)

    call sec_header(2, 'Definition of Transition Amplitudes')
    if (debug%timer) call write_time('Build trans. amplitudes')
    allocate(t_mo(orbdb%n_mo,orbdb%n_mo,orbdb%n_ab,excdb%n_states))
    ! "Now doing state num. ",i3," out of 123
    write(fmt_elstate, '("('' Now doing state num. '',i",i0,",'' out of '',&
        &i0,''.'')")') len_int(excdb%n_states)
    do istate = 1, excdb%n_states
        write(iu_out, fmt_elstate) istate, excdb%n_states
        t_mo(:,:,:,istate) = &
            eltrans_amp(orbdb%n_ab, orbdb%n_ao, orbdb%n_mos, ao_int%i_j, &
                        excdb%g2e_dens(:,:,:,istate), tmp_ao_arr1, .true., &
                        orbdb%coef_mos)
        if (.not.orbdb%openshell) &
            t_mo(:,:,:,istate) = t_mo(:,:,:,istate)*sqrt(f2)
    end do

    call sec_header(1, 'Computation of Properties')
    call sec_header(2, 'Ground-state properties')

    if (debug%timer) call write_time('Ground-state moments')
    p_gg_ab = f0
    r_gg_ab = f0
    rxp_gg_ab = f0
    if (orbdb%openshell) then
        do iab = 1, orbdb%n_ab
            do i = 1, orbdb%n_els(iab)
                p_gg_ab(:,iab) = p_gg_ab(:,iab) + Smo_ipj(:,i,i,iab)
                r_gg_ab(:,iab) = r_gg_ab(:,iab) + Smo_irj(:,i,i,iab)
                rxp_gg_ab(:,iab) = rxp_gg_ab(:,iab) + Smo_irxpj(:,i,i,iab)
            end do
        end do
    else
        do i = 1, orbdb%n_els(1)
            p_gg_ab(:,1)   = p_gg_ab(:,1) + Smo_ipj(:,i,i,1)
            r_gg_ab(:,1)   = r_gg_ab(:,1) + Smo_irj(:,i,i,1)
            rxp_gg_ab(:,1) = rxp_gg_ab(:,1) + Smo_irxpj(:,i,i,1)
        end do
        p_gg_ab(:,1)   = f2*p_gg_ab(:,1)
        r_gg_ab(:,1)   = f2*r_gg_ab(:,1)
        rxp_gg_ab(:,1) = f2*rxp_gg_ab(:,1)
    end if
    p_gg   = p_gg_ab(:,1)
    r_gg   = r_gg_ab(:,1)
    rxp_gg = rxp_gg_ab(:,1)
    if (orbdb%n_ab == 2) then
        p_gg   = p_gg   + p_gg_ab(:,2)
        r_gg   = r_gg   + r_gg_ab(:,2)
        rxp_gg = rxp_gg + rxp_gg_ab(:,2)
    end if

    edip_nuc(1) = sum(moldb%at_crd(1,:)*moldb%at_chg)
    edip_nuc(2) = sum(moldb%at_crd(2,:)*moldb%at_chg)
    edip_nuc(3) = sum(moldb%at_crd(3,:)*moldb%at_chg)

    if (debug%print) then
        write(iu_out, '(a)') 'Electric dipole'
        write(iu_out, '("Nucl.: X = ",f12.6,", Y = ",f12.6,", Z = ",f12.6,/&
            &"Elec.: X = ",f12.6,", Y = ",f12.6,", Z = ",f12.6,/,&
            &"Total: X = ",f12.6,", Y = ",f12.6,", Z = ",f12.6)') &
            edip_nuc, -r_gg, edip_nuc-r_gg
        write(iu_out, '("rxp",3f12.6)') rxp_gg
    end if

    call sec_header(2, 'Transition Energies Factors')

    if (debug%timer) call write_time('Transition energies')
    allocate(ov_eieg(excdb%n_states), ov_eief(excdb%n_states))
    ef = excdb%g2e_energy(excdb%id_state)
    do istate = 1, excdb%n_states
        ei = excdb%g2e_energy(istate)
        if (use_gamma) then
            ov_eieg(istate) = ei/(ei**2 + e_gamma**2)
        else
            ov_eieg(istate) = f1/ei
        end if
        if (istate /= excdb%id_state) then
            de = ei - ef
            if (use_gamma) then
                ov_eief(istate) = de/(de**2 + e_gamma**2)
            else
                ov_eief(istate) = f1/de
            end if
        end if
    end do

    call sec_header(1, 'MCD Tensor')
    write(iu_out, '(/,4(a,/),6(/,a))') &
        ' Computing tensor G_if:', &
        '           --      <f|u|k><k|m|i>    --       <f|m|k><k|u|i>', &
        '    G_if = \       --------------  + \        --------------', &
        '           /_ k!=i    E_k - E_i      /_ k!=f     E_k - E_f', &
        '    i:   initial state', &
        '    f:   final state', &
        '    k:   intermediate state', &
        '    u:   electric dipole', &
        '    m:   magnetic dipole', &
        '    E_i: energy of state i'

    ! For GIAO, we do some pre-processing by building the combinations
    !   < l | O | k > and < l | O | g > to be used later.
    if (do_giao) then
        if (debug%timer) call write_time('GIAO terms')
        allocate(r_lk(3,excdb%n_states,0:excdb%n_states), &
                 p_lk(3,excdb%n_states,0:excdb%n_states), &
                 ov_eiej(0:excdb%n_states,0:excdb%n_states))
        ov_eiej(0,0) = f0
        ov_eiej(excdb%id_state,0) = ov_eieg(excdb%id_state)
        ov_eiej(0,excdb%id_state) = -ov_eieg(excdb%id_state)
        do istate = 1, excdb%n_states
            ! This should never be used, since NaN, filled only for display
            ov_eiej(istate,istate) = f0
            ! We exclude id_states so we can treat last and preserve the
            !   prefactors
            if (istate /= excdb%id_state) then
                ov_eiej(istate,0) = ov_eieg(istate)
                ov_eiej(0,istate) = -ov_eieg(istate)
                pfac_r = &
                    sos_prefac_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                                     t_mo(:,:,:,istate), r_gg_ab, Smo_irj, &
                                     model_sos, debug%add_ijaa)
                pfac_p = &
                    sos_prefac_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                                     t_mo(:,:,:,istate), p_gg_ab, Smo_ipj, &
                                     model_sos, debug%add_ijaa)
                r_lk(:,istate,istate) = &
                    sos_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                              t_mo(:,:,:,istate), pfac_r, model=model_sos)
                p_lk(:,istate,istate) = &
                    sos_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                              t_mo(:,:,:,istate), pfac_p, model=model_sos)
                forbid = .not.orbdb%openshell .and. &
                          excdb%ispin_exc(istate) /= ispin_gs
                r_lk(:,istate,0) = &
                    sos_eiOg(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                             t_mo(:,:,:,istate), Smo_irj, forbid, &
                             model=model_sos)
                p_lk(:,istate,0) = &
                    sos_eiOg(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                             t_mo(:,:,:,istate), Smo_ipj, forbid, &
                             model=model_sos)
                ei = excdb%g2e_energy(istate)
                do jstate = 1, excdb%n_states
                    if (jstate /= istate) then
                        forbid = .not.orbdb%openshell .and. &
                            excdb%ispin_exc(jstate) /= excdb%ispin_exc(istate)
                        r_lk(:,istate,jstate) = &
                            sos_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els,&
                                      t_mo(:,:,:,jstate), pfac_r, forbid, &
                                      model=model_sos)
                        p_lk(:,istate,jstate) = &
                            sos_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els,&
                                      t_mo(:,:,:,jstate), pfac_p, forbid, &
                                      model=model_sos)
                        de = ei - excdb%g2e_energy(jstate)
                        if (use_gamma) then
                            ov_eiej(istate,jstate) = de/(de**2 + e_gamma**2)
                        else
                            ov_eiej(istate,jstate) = f1/de
                        end if
                    end if
                end do
            end if
        end do
    end if

    ! The prefactor are used to compute < f | O | e_i > and < f | O | g >
    !   (the latter for GIAO)
    ! For this reason, they are pre-computed once to speed up later
    if (debug%timer) call write_time('Prefactors')
    pfac_r = sos_prefac_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                              t_mo(:,:,:,excdb%id_state), r_gg_ab, Smo_irj, &
                              model_sos, debug%add_ijaa)
    pfac_p = sos_prefac_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                              t_mo(:,:,:,excdb%id_state), p_gg_ab, Smo_ipj, &
                              model_sos, debug%add_ijaa)
    pfac_rxp = sos_prefac_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                                t_mo(:,:,:,excdb%id_state), rxp_gg_ab, &
                                Smo_irxpj, model_sos, debug%add_ijaa)
    if (do_giao) then
        r_lk(:,excdb%id_state,excdb%id_state) = &
            sos_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                      t_mo(:,:,:,excdb%id_state), pfac_r, model=model_sos)
        p_lk(:,excdb%id_state,excdb%id_state) = &
            sos_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                      t_mo(:,:,:,excdb%id_state), pfac_p, model=model_sos)
        r_lk(:,excdb%id_state,0) = &
            sos_eiOg(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                     t_mo(:,:,:,excdb%id_state), Smo_irj, model=model_sos)
        p_lk(:,excdb%id_state,0) = &
            sos_eiOg(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                     t_mo(:,:,:,excdb%id_state), Smo_ipj, model=model_sos)
        do istate = 1, excdb%n_states
            if (istate /= excdb%id_state) then
                forbid = .not.orbdb%openshell .and. &
                    excdb%ispin_exc(istate) /= excdb%ispin_exc(excdb%id_state)
                r_lk(:,excdb%id_state,istate) = &
                    sos_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                              t_mo(:,:,:,istate), pfac_r, forbid, &
                              model=model_sos)
                p_lk(:,excdb%id_state,istate) = &
                    sos_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                              t_mo(:,:,:,istate), pfac_p, forbid, &
                              model=model_sos)
                ov_eiej(istate,excdb%id_state) = ov_eief(istate)
                ov_eiej(excdb%id_state,istate) = - ov_eief(istate)
            end if
        end do
    end if

    ! Compute the transition moment < f | O | g >
    ! we need to correct the sign of p and divided by the energy
    if (debug%timer) call write_time('Transition moments')
    forbid = .not.orbdb%openshell .and. &
        excdb%ispin_exc(excdb%id_state) /= ispin_gs
    if (do_giao) then
        r_fg = r_lk(:,excdb%id_state,0)
        p_fg = -p_lk(:,excdb%id_state,0) / excdb%g2e_energy(excdb%id_state)
    else
        r_fg = sos_eiOg(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                        t_mo(:,:,:,excdb%id_state), Smo_irj, forbid, &
                        model=model_sos)
        p_fg = -sos_eiOg(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                         t_mo(:,:,:,excdb%id_state), Smo_ipj, forbid, &
                         model=model_sos) / excdb%g2e_energy(excdb%id_state)
    endif
    rxp_fg = sos_eiOg(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                      t_mo(:,:,:,excdb%id_state), Smo_irxpj, forbid, &
                      model=model_sos)

    G_if = f0
    if (debug%print) write(iu_out, '(a)') 'NOW ON G_IF'
    do istate = 1, excdb%n_states
        if (istate /= excdb%id_state) then
            forbid = .not.orbdb%openshell .and. &
                excdb%ispin_exc(istate) /= ispin_gs
            r_kg = sos_eiOg(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                            t_mo(:,:,:,istate), Smo_irj, forbid, &
                            model=model_sos)
            rxp_kg = sos_eiOg(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                              t_mo(:,:,:,istate), Smo_irxpj, forbid, &
                              model=model_sos)
            forbid = .not.orbdb%openshell .and. &
                excdb%ispin_exc(istate) /= excdb%ispin_exc(excdb%id_state)
            r_fk = sos_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                             t_mo(:,:,:,istate), pfac_r, forbid, &
                             model=model_sos)
            rxp_fk = sos_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                               t_mo(:,:,:,istate), pfac_rxp, forbid, &
                               model=model_sos)
            if (debug%print) then
                write(iu_out, '(i4,"r_kg   ",3f12.6)') istate, r_kg
                write(iu_out, '(4x,"rxp_kg ",3f12.6)') rxp_kg
                write(iu_out, '(4x,"r_fk   ",3f12.6)') r_fk
                write(iu_out, '(4x,"rxp_fk ",3f12.6)') rxp_fk
                write(iu_out, '(4x,"ov_eiej",2f12.6)') &
                    ov_eieg(istate), ov_eief(istate)
            end if
            do jx = 1, 3
                do ix = 1, 3
                    G_if(ix,jx) = G_if(ix,jx) &
                        + r_fk(ix)*rxp_kg(jx)*ov_eieg(istate) &
                        + r_kg(ix)*rxp_fk(jx)*ov_eief(istate)
                end do
            end do
            if (do_giao) then
                call sos_MCD_tensor_LORG_corr( &
                    excdb%n_states, orbdb%n_els, excdb%id_state, istate, &
                    r_gg, r_lk, p_lk, ov_eieg, ov_eiej, G_if)
            end if
        end if
    end do
    ! Add special cases for k=g or k=f
    forbid = .not.orbdb%openshell .and. &
        excdb%ispin_exc(excdb%id_state) /= ispin_gs
    rxp_kg = sos_eiOg(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                      t_mo(:,:,:,excdb%id_state), Smo_irxpj, forbid, &
                      model=model_sos)
    r_fk = sos_ejOei(3, orbdb%n_ab, orbdb%n_mos, orbdb%n_els, &
                     t_mo(:,:,:,excdb%id_state), pfac_r, forbid, &
                     model=model_sos)
    ! Add case k=f for first term of G
    do jx = 1, 3
        do ix = 1, 3
            G_if(ix,jx) = G_if(ix,jx) &
                + r_fk(ix)*rxp_kg(jx)*ov_eieg(excdb%id_state)
        end do
    end do
    if (debug%timer) call write_time('LORG correction')
    if (do_giao) then
        call sos_MCD_tensor_LORG_corr( &
            excdb%n_states, orbdb%n_els, excdb%id_state, excdb%id_state, &
            r_gg, r_lk, p_lk, ov_eieg, ov_eiej, G_if)
    end if
    ! Add case k=g for second term of G
    do jx = 1, 3
        do ix = 1, 3
            G_if(ix,jx) = G_if(ix,jx) &
                - r_gg(ix)*rxp_kg(jx)*ov_eieg(excdb%id_state)
        end do
    end do
    if (do_giao) then
        call sos_MCD_tensor_LORG_corr( &
            excdb%n_states, orbdb%n_els, excdb%id_state, 0, r_gg, r_lk, p_lk, &
            ov_eieg, ov_eiej, G_if)
    end if

    ! Correct using a factor of 1/2
    G_if = G_if / f2

    call sec_header(2, 'Final Value')
    write(iu_out, '(10X,"X              Y              Z")')
    do i = 1, 3
        write(iu_out, '(1x,a,3es15.6)') labXYZ_1D(i), G_if(i,:)
    end do

    if (show_timestamp) then
        datetime = timestamp()
        write(iu_out, '(/," Calculations ended on ", a)') datetime
    end if

contains
    function fmt_param(label, sub) result(res)
        !! Formats the parameter label for output
        implicit none

        integer, parameter :: maxchar = 44

        character(len=*), intent(in) :: label
        !! Label of the parameter
        logical, intent(in), optional :: sub
        !! Build parameter label as suboption.
        character(len=maxchar) :: res
        !! Result of the function

        logical :: main

        if (present(sub)) then
            main = .not.sub
        else
            main = .true.
        end if

        if (main) then
            res = ' - ' // trim(label) // ' ' // repeat('-', maxchar)
            res(maxchar:maxchar) = '>'
        else
            res = '   > ' // trim(label)
            res(maxchar:maxchar) = ':'
        end if
    end function fmt_param

    subroutine parse_argopts(infile, outfile, e_gamma, f_state, max_state, &
                             use_gamma, do_giao, do_guvcde, in_mem, &
                             model_sos, show_timestamp, debug)
        !! Parses commandline options and updates information.
        !!
        !! Builds options parser and parse user-options, setting
        !! default values where necessary.
        implicit none

        character(len=:), allocatable, intent(out) :: infile
        !! Input filename.
        character(len=:), allocatable, intent(out) :: outfile
        !! Optional output filename.
        real(realwp), intent(out) :: e_gamma
        !! Energy shift to avoid accidental degeneracies in summation.
        integer, intent(out) :: f_state
        !! Final state of interest.
        integer, intent(out) :: max_state
        !! Highest electronic state to consider in the summations.
        logical, intent(out) :: use_gamma
        !! Use energy shift correction gamma.
        logical, intent(out) :: do_giao
        !! Include GIAO correction.
        logical, intent(out) :: do_guvcde
        !! Add special routines to generate data compatible with GUVCDE.
        !! This is intended for debugging purpose.
        logical, intent(out) :: in_mem
        !! Store integrals in memory (default, alternative NYI)
        integer, intent(out) :: model_sos
        !! Model for the SOS excitation:
        !! 0: auto
        !! 1: use pure TD-DFT amplitudes
        !! 2: use pure Slater determinants
        !! 3: hybrid, all TD-DFT amplitudes + Slater permutations
        logical, intent(out) :: show_timestamp
        !! Show timestamp in output.
        type(DebugType), intent(out) :: debug
        !! Debug flags

        integer :: i
        character(len=1024) :: string
        character(len=:), dimension(:), allocatable :: svalues
        class(CmdArgDB), allocatable :: opts

        ! Basic initialization
        e_gamma = 1.0e-4_realwp
        do_guvcde = .false.  ! No testing.  Code does not try to match SOS/GUVCDE.
        do_giao = .true.  ! Include GIAO correction
        in_mem = .true.  ! Store everything in memory
        use_gamma = .true.  ! Use energy correction gamma
        f_state = -1  ! Final state is chosen automatically
        max_state = -1  ! Highest state in summation chosen automatically
        show_timestamp = .true.  ! By default, show timestamp.
        ! Debug mode
        debug%user_set = .false.
        debug%print = .false.
        debug%timer = .false.
        debug%add_ijaa = .true.

        ! Build parser and check arguments
        opts = CmdArgDB(progname='mcd_tensor')
        if (opts%has_error()) then
            call sec_header(-1, PROGTITLE)
            write(*, '(a)') 'Failed to initialize the command-line parser'
            select type (err => opts%exception())
                class is (AllocateError)
                    write(*, '(a)') trim(err%msg())
                class is (ArgumentError)
                    write(*, '(a)') trim(err%msg())
                class is (Error)
                    write(*, '(a)') trim(err%msg())
                class default
                    write(*, '(a)') 'Unknown error.'
            end select
            stop
        end if
        call opts%add_arg_char( &
            'string', label='filename', &
            help='Gaussian formatted checkpoint file')
        call opts%add_arg_real( &
            'real', longname='--gamma', &
            help='Include a shift term, gamma, to avoid a divergence of the &
                &denominator in the summation with close states (unit: Hartree).',&
            def_value=real(e_gamma), min_value=0.0)
        call opts%add_arg_int( &
            'scalar', shortname='-f', longname='--final', &
            help='Final excited electronic state (starting from 1)', &
            min_value=1)
        call opts%add_arg_char( &
            'string', longname='--exc-model', &
            help='Model to represent excited wave functions.  &
            &Possible choices are: &
            &amplitudes (use TD-DFT amplitudes coefficients), &
            &slater (use Slater determinants with amplitudes, only &
            &considering excitations occupied->virtual), &
            &hybrid (default, use TD-DFT amplitudes and add Slater-like &
            &permutations)')
        call opts%add_arg_bool( &
            'store_true', longname='--giao', &
            help='Include GIAO corrections (default)', &
            def_value=.True.)
        call opts%add_arg_bool( &
            'store_true', longname='--guvcde', &
            help='Add printing similar to SOS/GUVCDE for control.')
        call opts%add_arg_int( &
            'scalar', shortname='-m', longname='--max-state', &
            help='Highest electronic state to include (upper bound of summation)', &
            min_value=1)
        call opts%add_arg_bool( &
            'store_false', longname='--no-giao', &
            help='Include GIAO corrections')
        call opts%add_arg_char( &
            'string', label='output', shortname='-o', longname='--output', &
            help='Name of the file to store the output. Existing content will &
            &be overwritten.')
        call opts%add_arg_bool( &
            'store_false', longname='--no-timestamp', &
            help='Do not show date and time in output (to facilitate regression tests).')
        call opts%add_arg_char( &
            'list', label='debug', longname='--debug', &
            help='Debug flags. Supported: "print", "ijaa", "ijaa", "timer".  &
            &Use "no" before keyword to deactivate.')

        call opts%parse_args()
        if (opts%has_error()) then
            call sec_header(-1, PROGTITLE)
            select type (err => opts%exception())
                class is (ValueError)
                    write(*, '(a)') trim(opts%get_error())
                class is (Error)
                    write(*, '(a)') trim(opts%get_error())
                class default
                    write(*, '(a)') 'Unknown error while reading "gamma"'
            end select
            stop
        end if

        ! Check debug flags
        if (opts%is_user_set('debug')) then
            debug%user_set = .true.
            call opts%get_value('debug', svalues)
            do i = 1, size(svalues)
                select case(locase(svalues(i)))
                    case ('print')
                        debug%print = .true.
                    case ('noprint')
                        debug%print = .false.
                    case ('timer')
                        debug%timer = .true.
                    case ('notimer')
                        debug%timer = .false.
                    case ('ijaa')
                        debug%add_ijaa = .true.
                    case ('noijaa')
                        debug%add_ijaa = .false.
                end select
            end do
        end if
        
        ! Output file
        if (opts%is_user_set('output')) then
            call opts%get_value('output', string)
            outfile = trim(string)
            open(newunit=iu_out, file=outfile, action='write')
        end if
        call sec_header(-1, PROGTITLE)
        if (debug%timer) call write_time('Option parser')
        call sec_header(1, 'Simulation Parameters')
        call sec_header(2, 'Standard parameters')

        ! Check requirement to use/not use GIAO correction
        if (opts%is_user_set('giao') &
            .and. opts%is_user_set('no-giao')) then
            write(*, '(" Error: conflicting option for the definition of GIAO")')
            stop
        end if

        if (opts%is_user_set('no-giao')) then
            do_giao = .false.
        else if (opts%is_user_set('giao')) then
            do_giao = .true.
        end if
        call write_param('GIAO correction', do_giao)

        ! Check requirement for GUVCDE compatibility
        if (opts%is_user_set('guvcde')) do_guvcde = .true.
        call write_param('Compatibility mode with GUVCDE', do_guvcde)

        ! Check for value of e_gamma
        call opts%get_value('gamma', e_gamma)
        use_gamma = abs(e_gamma) > tiny(e_gamma)
        call write_param('Correction term against degeneracies', use_gamma)
        if (use_gamma) &
            call write_param('Value of the term (in Hartrees)', e_gamma, .true.)

        ! Get input file name
        call opts%get_value('filename', string)
        infile = trim(string)
        inquire(file=infile, exist=exists)
        if (.not.exists) then
            write(*, '("Error: File ",a," does not exist.")') trim(infile)
            stop
        end if
        call write_param('Input filename', infile)

        ! Set default final state.
        if (opts%is_user_set('final')) then
            call opts%get_value('final', f_state)
            call write_param('Final electronic state', f_state)
        else
            call write_param('Final electronic state', 'automatic')
        end if

        if (opts%is_user_set('max-state')) then
            call opts%get_value('max-state', max_state)
            call write_param('Highest excited state', max_state)
        else
            call write_param('Highest excited state', 'include all')
        end if

        if (opts%is_user_set('exc-model')) then
            call opts%get_value('exc-model', string)
            select case (trim(locase(string)))
            case ('ampl', 'amplitude', 'amplitudes')
                model_sos = 1
            case ('slater')
                model_sos = 2
            case ('hybrid')
                model_sos = 3
            case default
                write(*, '(" Error: Unrecognized model for SOS excitation.")')
                stop 1
            end select
        else
            model_sos = 3
        end if

        select case (model_sos)
        case (1)
            call write_param('Model for SOS excitations', 'TD amplitudes')
        case (2)
            call write_param('Model for SOS excitations', 'Slater determinant')
        case (3)
            call write_param('Model for SOS excitations', 'Hybrid model')
        case default
            write(*, '(" Error: List of supported SOS models needs updating")')
            stop 2
        end select

        if (opts%is_user_set('no-timestamp')) show_timestamp = .false.

        ! debugging mode - print information
        if (debug%user_set) then
            call sec_header(2, 'Debugging parameters')
            call write_param("Extra printing", debug%print)
            call write_param("Add timer", debug%timer)
            call write_param("include <ia||ja> terms", debug%add_ijaa)
        end if

    end subroutine parse_argopts

    subroutine write_moldata(moldb, bsetdb, orbdb)
        !! Write molecular data.
        !!
        !! Writes some of the molecular data stored in the moldata module
        !! for mcd_tensor.
        implicit none

        class(MoleculeDB), intent(in) :: moldb
        !! Molecule specifications database.
        class(BasisSetDB), intent(in) :: bsetdb
        !! Basis set specifications database.
        class(OrbitalsDB), intent(in) :: orbdb
        !! Molecular orbitals database.

        character(len=:), allocatable :: tag_bfunD, tag_bfunF, tag_open

        1100 format(/, &
            ' > Number of atoms           : ',i0,/, &
            ' > Charge                    : ',i0,/, &
            ' > Multiplicity              : ',i0,/, &
            ' > Open-/Closed-shell        : ',a,//, &
            ' > Number of basis functions : ',i0,/, &
            ' > Type of D functions       : ',a,/, &
            ' > Type of F functions       : ',a,/, &
            ' > Number of atomic orbitals : ',i0)

        if (orbdb%openshell) then
            tag_open = 'open'
        else
            tag_open = 'closed'
        end if

        if (bsetdb%pureD) then
            tag_bfunD = 'pure'
        else
            tag_bfunD = 'Cartesian'
        end if
        if (bsetdb%pureF) then
            tag_bfunF = 'pure'
        else
            tag_bfunF = 'Cartesian'
        end if

        write(iu_out, 1100) moldb%n_at, moldb%charge, moldb%multip, tag_open, &
            bsetdb%n_basis, tag_bfunD, tag_bfunF, orbdb%n_ao

        call sec_header(2, 'Atomic coordinates')
        call prt_coord(moldb%n_at, moldb%at_lab, moldb%at_crd, moldb%at_mas)
    end subroutine write_moldata

    subroutine write_param_bool(param, status, subparam)
        !! Writes status for logical-type parameter.
        !!
        !! Formats and writes one parameter of simulation with the
        !! specified status.
        implicit none

        character(len=*), intent(in) :: param
        !! Parameter name (with info like unit).
        logical, intent(in) :: status
        !! Status
        logical, intent(in), optional :: subparam
        !! Parameter is a sub-parameter type

        1000 format(a,1x,"enabled")
        1010 format(a,1x,"disabled")

        if (status) then
            write(iu_out, 1000) fmt_param(param, subparam)
        else
            write(iu_out, 1010) fmt_param(param, subparam)
        end if
    end subroutine write_param_bool

    subroutine write_param_int(param, value, subparam)
        !! Writes value of integer-type parameter.
        !!
        !! Formats and writes one parameter of simulation with the
        !! indicated value.
        implicit none

        character(len=*), intent(in) :: param
        !! Parameter name (with info like unit).
        integer, intent(in) :: value
        !! Value associated to parameter
        logical, intent(in), optional :: subparam
        !! Parameter is a sub-parameter type

        1000 format(a,1x,i0)

        write(iu_out, 1000) fmt_param(param, subparam), value
    end subroutine write_param_int

    subroutine write_param_real(param, value, subparam)
        !! Writes value of real-type parameter.
        !!
        !! Formats and writes one parameter of simulation with the
        !! indicated value.
        implicit none

        character(len=*), intent(in) :: param
        !! Parameter name (with info like unit).
        real(realwp), intent(in) :: value
        !! Value associated to parameter
        logical, intent(in), optional :: subparam
        !! Parameter is a sub-parameter type

        1000 format(a,1x,f0.4)
        1010 format(a,es13.6)

        if (abs(value) < 1.0e-2_realwp .or. abs(value) > 1.0e4_realwp) then
            write(iu_out, 1010) fmt_param(param, subparam), value
        else
            write(iu_out, 1000) fmt_param(param, subparam), value
        end if
    end subroutine write_param_real

    subroutine write_param_char(param, value, subparam)
        !! Writes value of character-type parameter.
        !!
        !! Formats and writes one parameter of simulation with the
        !! indicated value.
        implicit none

        character(len=*), intent(in) :: param
        !! Parameter name (with info like unit).
        character(len=*), intent(in) :: value
        !! Value associated to parameter
        logical, intent(in), optional :: subparam
        !! Parameter is a sub-parameter type

        1000 format(a,1x,a)

        write(iu_out, 1000) fmt_param(param, subparam), trim(value)
    end subroutine write_param_char

    subroutine write_time(label)
        !! Writes time information.
        !!
        !! Writes time data at call time, displaying the label as well.
        implicit none

        character(len=*), intent(in) :: label
        !! Label to display in printing information
        
        character(len=:), allocatable :: dtime

        1000 format('Entering ',a,' on ',a)
        dtime = timestamp()
        write(iu_out, 1000) trim(label), dtime
    end subroutine write_time

end program mcd_tensor
