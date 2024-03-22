program mcd_tensor
    use iso_fortran_env, only: real64, output_unit
    use string, only: labXYZ_1D
    use input, only: build_moldata, build_transdata
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
    use moldata
    use transdata

    implicit none

    integer :: fstate, i, iab, istate, ix, jstate, jx, max_states, qty_flag
    real(real64) :: de, ef, ei
    real(real64) :: e_gamma = 1.0e-4
    real(real64), dimension(3) :: &
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
    real(real64), dimension(3,2) :: &
        r_gg_ab, &   ! < g | r | g >
        p_gg_ab, &   ! < g | p | g >
        rxp_gg_ab    ! < g | r x p | g >
    real(real64), dimension(3,3) :: G_if
    real(real64), dimension(:), allocatable :: &
        ov_eieg, &   ! Inverse of (E_i-E_g)
        ov_eief      ! Inverse of (E_i-E_f)
    real(real64), dimension(:,:), allocatable :: &
        ov_eiej      ! Inverse of (E_i-E_j)
    real(real64), dimension(:,:,:), allocatable :: &
        r_lk, &      ! Stores all combinations < l | r | k >
        p_lk         ! Stores all combinations < l | p | k >
    real(real64), dimension(:,:,:,:), allocatable :: &
        t_mo, &      ! Transition amplitude in MO basis
        pfac_r, &    ! Stores the prefactor for < e_j | r | e_i >
        pfac_p, &    ! Stores the prefactor for < e_j | p | e_i >
        pfac_rxp, &  ! Stores the prefactor for < e_j | r x p | e_i >
        Smo_ipj, &   ! MO-basis integral < i | p | j >
        Smo_irj, &   ! MO-basis integral < i | r | j >
        Smo_irxpj    ! MO-basis integral < i | r x p | j >
    real(real64), allocatable :: tmp_ao_arr1(:,:), tmp_ao_arrN(:,:,:)
    ! for_guvcde = true to print information to compare with SOS/GUVCDE
    logical, parameter :: DEBUG = .False., TIMEIT = .False.
    logical :: &
        for_guvcde = .False.,  & ! Code behavior matches SOS/GUVCDE for testing
        use_gamma = .True., &    ! Use gamma to "fuzzy" singularity condition
        use_giao = .True., &     ! Use GIAO
        in_mem = .True., &       ! Store integrals in memory
        exists
    character(len=512) :: fmt_elstate, fname
    class(ovij_1e), allocatable :: ao_int
    class(BaseException), allocatable :: err

    interface write_param
        procedure write_param_bool, write_param_int, write_param_real, &
            write_param_char
    end interface write_param

    1200 format(/, &
        ' > Number of excited states : ',i0,/, &
        ' > Reference excited state  : ',i0)

    ! Write title
    call sec_header(-1, 'MCD Tensor Calculator')

    call parse_opts

    call sec_header(1, 'Data on Molecular System')
    if (TIMEIT) call write_time('Molecular data')
    call build_moldata(fname, err)
    if (err%raised()) then
        select type(err)
            class is (FileError)
                call write_err('std', 'Error found with input file', err%msg())
                stop
            class is (Error)
                call write_err('std', &
                    'Error found while parsing basic molecular data', &
                    err%msg() &
                )
                stop
            class default
                call write_err('gen', &
                    'Something went wrong while parsing molecular data.' &
                )
                stop
        end select
    end if
    if (openshell) then
        call write_err('gen', &
            'Sorry, open-shell systems not yet supported.  Working on it.' &
            )
        stop
    end if
    call write_moldata

    if (for_guvcde) then
        if (TIMEIT) call write_time('GUVCDE write control')
        call write_control(iu_out, err)
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
    if (TIMEIT) call write_time('Check AO normalization')
    call fix_norm_AOs(iu_out, n_at, n_ao, at_crd, nprim_per_at, bset_info, &
                      err, DEBUG)
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
    if (TIMEIT) call write_time('1-e AO integral')
    ! overlap6
    qty_flag = 2**10 - 1
    ! qty_flag = ibset(qty_flag, 0)
    call overlap_ao_1e(iu_out, n_at, n_ao, qty_flag, for_guvcde, in_mem, &
                       at_crd, nprim_per_at, bset_info, ao_int, err)
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
    if (for_guvcde) then
        if (TIMEIT) call write_time('GUVCDE build MOs')
        call build_MOs(n_ab, n_ao, n_mos, coef_mos, .true.)
    endif
    ! call convert_AO2MO(n_ao, n_mo, coef_mos(1), ao_int%i_j, Smo_ij)
    allocate(Smo_irj(3,n_mo,n_mo,n_ab), Smo_ipj(3,n_mo,n_mo,n_ab), &
             Smo_irxpj(3,n_mo,n_mo,n_ab))
    allocate(tmp_ao_arr1(n_ao,n_ao), tmp_ao_arrN(3,n_ao,n_ao))
    if (TIMEIT) call write_time('AO to MO conversion')
    do iab = 1, n_ab
        call convert_AO2MO(n_ao, n_mos(iab), coef_mos(:,:,iab), ao_int%i_r_j, &
                           Smo_irj(:,:,:,iab), tmp_ao_arrN)
        if (TIMEIT) call write_time('<i|p|j>')
        call convert_AO2MO(n_ao, n_mos(iab), coef_mos(:,:,iab), ao_int%i_p_j, &
                           Smo_ipj(:,:,:,iab), tmp_ao_arrN)
        if (TIMEIT) call write_time('<i|rxp|j>')
        call convert_AO2MO(n_ao, n_mos(iab), coef_mos(:,:,iab), &
                           ao_int%i_rxp_j, Smo_irxpj(:,:,:,iab), tmp_ao_arrN)
    end do
    ! call convert_AO2MO(n_ao, n_mo, coef_mos(1), ao_int%i_j)

    call sec_header(1, 'Transition Data')
    if (TIMEIT) call write_time('Transition data')
    call build_transdata(fname, n_ab, n_basis, err)
    if (err%raised()) then
        select type(err)
            class is (FileError)
                call write_err('std', 'Error found with input file', err%msg())
                stop
            class is (Error)
                call write_err('std', &
                    'Error while parsing transition data', &
                    err%msg() &
                )
                stop
            class default
                call write_err('gen', &
                    'Something went wrong while parsing transition data.' &
                )
                stop
        end select
    end if
    write(iu_out, 1200) n_states, id_state

    call sec_header(2, 'User overrides')
    exists = .false.
    if (fstate > 0 .and. fstate /= id_state) then
        exists = .true.
        write(iu_out, '(" > New final state, num. ",i0,".")') &
            fstate
        id_state = fstate
    end if
    if (max_states > n_states) then
        exists = .true.
        write(iu_out, '(1x,"NOTE: ",a,/,1x,a)') &
            'Requested upper bound exceeds available states.', &
            'Reverting to highest available state.'
    else if (max_states > 0) then
        exists = .true.
        write(iu_out, '(" > Setting highest electronic state to ",i0,".")') &
            max_states
        n_states = max_states
    end if
    if (.not.exists) &
        write(iu_out, '(1x,a)') &
            'No modifications requested. Proceeding with data found.'

    ! call prt_mat(g2e_dens(:,:,1,1), n_basis, n_basis)

    call sec_header(2, 'Definition of Transition Amplitudes')
    if (TIMEIT) call write_time('Build trans. amplitudes')
    allocate(t_mo(n_mo, n_mo, n_ab, n_states))
    ! "Now doing state num. ",i3," out of 123
    write(fmt_elstate, '("('' Now doing state num. '',i",i0,",'' out of '',&
        &i0,''.'')")') len_int(n_states)
    do istate = 1, n_states
        write(iu_out, fmt_elstate) istate, n_states
        t_mo(:,:,:,istate) = eltrans_amp(n_ab, n_ao, n_mos, ao_int%i_j, &
                                         g2e_dens(:,:,:,istate), tmp_ao_arr1, &
                                         .True., coef_mos)
        ! ! call prt_mat(transpose(g2e_dens(:,:,1,istate)), n_ao, n_ao)
        ! do i = 1, n_basis
        !     do j = 1, n_basis
        !         tmp(j,i) = 0.0_real64
        !         do k = 1, n_basis
        !             do l = 1, n_basis
        !                 tmp(j,i) = tmp(j,i) + &
        !                     ao_int%i_j(j,k)*g2e_dens(l,k,1,istate)*ao_int%i_j(l,i)
        !             end do
        !         end do
        !     end do
        ! end do
        ! ! call prt_mat(tmp, n_ao, n_ao, thresh_=1.0e-6_real64)
        ! call convert_AO2MO(n_ao, n_mo, coef_a_mo, tmp, g2e_MOs)
        ! call prt_mat(t_mo(:,:,1,istate), n_mos(1), n_mos(1), &
        !              thresh_=1.0e-6_real64)

    end do

    call sec_header(1, 'Computation of Properties')
    call sec_header(2, 'Ground-state properties')

    if (TIMEIT) call write_time('Ground-state moments')
    p_gg_ab = 0.0_real64
    r_gg_ab = 0.0_real64
    rxp_gg_ab = 0.0_real64
    if (openshell) then
        do iab = 1, n_ab
            do i = 1, n_els(iab)
                p_gg_ab(:,iab) = p_gg_ab(:,iab) + Smo_ipj(:,i,i,iab)
                r_gg_ab(:,iab) = r_gg_ab(:,iab) + Smo_irj(:,i,i,iab)
                rxp_gg_ab(:,iab) = rxp_gg_ab(:,iab) + Smo_irxpj(:,i,i,iab)
            end do
        end do
    else
        do i = 1, n_els(1)
            p_gg_ab(:,1)   = p_gg_ab(:,1) + Smo_ipj(:,i,i,1)
            r_gg_ab(:,1)   = r_gg_ab(:,1) + Smo_irj(:,i,i,1)
            rxp_gg_ab(:,1) = rxp_gg_ab(:,1) + Smo_irxpj(:,i,i,1)
        end do
        p_gg_ab(:,1)   = 2.0_real64*p_gg_ab(:,1)
        r_gg_ab(:,1)   = 2.0_real64*r_gg_ab(:,1)
        rxp_gg_ab(:,1) = 2.0_real64*rxp_gg_ab(:,1)
    end if
    p_gg   = p_gg_ab(:,1)
    r_gg   = r_gg_ab(:,1)
    rxp_gg = rxp_gg_ab(:,1)
    if (n_ab == 2) then
        p_gg   = p_gg   + p_gg_ab(:,2)
        r_gg   = r_gg   + r_gg_ab(:,2)
        rxp_gg = rxp_gg + rxp_gg_ab(:,2)
    end if

    edip_nuc(1) = sum(at_crd(1,:)*at_chg)
    edip_nuc(2) = sum(at_crd(2,:)*at_chg)
    edip_nuc(3) = sum(at_crd(3,:)*at_chg)

    if (DEBUG) then
        write(iu_out, '(a)') 'Electric dipole'
        write(iu_out, '("Nucl.: X = ",f12.6,", Y = ",f12.6,", Z = ",f12.6,/&
            &"Elec.: X = ",f12.6,", Y = ",f12.6,", Z = ",f12.6,/,&
            &"Total: X = ",f12.6,", Y = ",f12.6,", Z = ",f12.6)') &
            edip_nuc, -r_gg, edip_nuc-r_gg
        write(iu_out, '("rxp",3f12.6)') rxp_gg
    end if

    call sec_header(2, 'Transition Energies Factors')

    if (TIMEIT) call write_time('Transition energies')
    allocate(ov_eieg(n_states), ov_eief(n_states))
    ef = g2e_energy(id_state)
    do istate = 1, n_states
        ei = g2e_energy(istate)
        if (use_gamma) then
            ov_eieg(istate) = ei/(ei**2 + e_gamma**2)
        else
            ov_eieg(istate) = 1.0_real64/ei
        end if
        if (istate /= id_state) then
            de = ei - ef
            if (use_gamma) then
                ov_eief(istate) = de/(de**2 + e_gamma**2)
            else
                ov_eief(istate) = 1.0_real64/de
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
    if (use_giao) then
        if (TIMEIT) call write_time('GIAO terms')
        allocate(r_lk(3,n_states,0:n_states), p_lk(3,n_states,0:n_states), &
                 ov_eiej(0:n_states,0:n_states))
        do istate = 1, n_states
            ! We exclude id_states so we can treat last and preserve the
            !   prefactors
            if (istate /= id_state) then
                ei = g2e_energy(istate)
                if (use_gamma) then
                    ov_eiej(istate,0) = ei/(ei**2 + e_gamma**2)
                else
                    ov_eiej(istate,0) = 1.0_real64/ei
                end if
                ov_eiej(0,istate) = - ov_eiej(istate,0)
                pfac_r = sos_prefac_ejOei(3, n_ab, n_mos, n_els, &
                                          t_mo(:,:,:,istate), r_gg_ab, Smo_irj)
                pfac_p = sos_prefac_ejOei(3, n_ab, n_mos, n_els, &
                                          t_mo(:,:,:,istate), p_gg_ab, Smo_ipj)
                r_lk(:,istate,istate) = &
                    sos_ejOei(3, n_ab, n_mos, t_mo(:,:,:,istate), pfac_r)
                p_lk(:,istate,istate) = &
                    sos_ejOei(3, n_ab, n_mos, t_mo(:,:,:,istate), pfac_p)
                r_lk(:,istate,0) = &
                    sos_eiOg(3, n_ab, n_mos, t_mo(:,:,:,istate), Smo_irj)
                p_lk(:,istate,0) = &
                    sos_eiOg(3, n_ab, n_mos, t_mo(:,:,:,istate), Smo_ipj)
                do jstate = 1, n_states
                    if (jstate /= istate) then
                        r_lk(:,istate,jstate) = &
                            sos_ejOei(3, n_ab, n_mos, t_mo(:,:,:,jstate), &
                                      pfac_r)
                        p_lk(:,istate,jstate) = &
                            sos_ejOei(3, n_ab, n_mos, t_mo(:,:,:,jstate), &
                                      pfac_p)
                        de = ei - g2e_energy(jstate)
                        if (use_gamma) then
                            ov_eiej(istate,jstate) = de/(de**2 + e_gamma**2)
                        else
                            ov_eiej(istate,jstate) = 1.0_real64/de
                        end if
                    end if
                end do
            end if
        end do
    end if

    ! The prefactor are used to compute < f | O | e_i > and < f | O | g >
    !   (the latter for GIAO)
    ! For this reason, they are pre-computed once to speed up later
    if (TIMEIT) call write_time('Prefactors')
    pfac_r = sos_prefac_ejOei(3, n_ab, n_mos, n_els, t_mo(:,:,:,id_state), &
                              r_gg_ab, Smo_irj)
    pfac_p = sos_prefac_ejOei(3, n_ab, n_mos, n_els, t_mo(:,:,:,id_state), &
                              p_gg_ab, Smo_ipj)
    pfac_rxp = sos_prefac_ejOei(3, n_ab, n_mos, n_els, t_mo(:,:,:,id_state), &
                                rxp_gg_ab, Smo_irxpj)
    if (use_giao) then
        r_lk(:,id_state,id_state) = &
            sos_ejOei(3, n_ab, n_mos, t_mo(:,:,:,id_state), pfac_r)
        p_lk(:,id_state,id_state) = &
            sos_ejOei(3, n_ab, n_mos, t_mo(:,:,:,id_state), pfac_p)
        r_lk(:,id_state,0) = &
            sos_eiOg(3, n_ab, n_mos, t_mo(:,:,:,id_state), Smo_irj)
        p_lk(:,id_state,0) = &
            sos_eiOg(3, n_ab, n_mos, t_mo(:,:,:,id_state), Smo_ipj)
        do istate = 1, n_states
            if (istate /= id_state) then
                r_lk(:,id_state,istate) = &
                    sos_ejOei(3, n_ab, n_mos, t_mo(:,:,:,istate), &
                                pfac_r)
                p_lk(:,id_state,istate) = &
                    sos_ejOei(3, n_ab, n_mos, t_mo(:,:,:,istate), &
                                pfac_p)
                ov_eiej(istate,id_state) = ov_eief(istate)
                ov_eiej(id_state,istate) = - ov_eief(istate)
            end if
        end do
    end if

    ! Compute the transition moment < f | O | g >
    ! we need to correct the sign of p and divided by the energy
    if (TIMEIT) call write_time('Transition moments')
    if (use_giao) then
        r_fg = r_lk(:,id_state,0)
        p_fg = -p_lk(:,id_state,0) / g2e_energy(id_state)
    else
        r_fg = sos_eiOg(3, n_ab, n_mos, t_mo(:,:,:,id_state), Smo_irj)
        p_fg = -sos_eiOg(3, n_ab, n_mos, t_mo(:,:,:,id_state), Smo_ipj) &
            / g2e_energy(id_state)
    endif
    rxp_fg = sos_eiOg(3, n_ab, n_mos, t_mo(:,:,:,id_state), Smo_irxpj)

    G_if = 0.0_real64
    if (DEBUG) write(iu_out, '(a)') 'NOW ON G_IF'
    do istate = 1, n_states
        if (istate /= id_state) then
            r_kg = sos_eiOg(3, n_ab, n_mos, t_mo(:,:,:,istate), Smo_irj)
            rxp_kg = sos_eiOg(3, n_ab, n_mos, t_mo(:,:,:,istate), Smo_irxpj)
            r_fk = sos_ejOei(3, n_ab, n_mos, t_mo(:,:,:,istate), pfac_r)
            rxp_fk = sos_ejOei(3, n_ab, n_mos, t_mo(:,:,:,istate), pfac_rxp)
            if (DEBUG) then
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
            if (use_giao) then
                call sos_MCD_tensor_LORG_corr( &
                    n_states, n_els, id_state, istate, r_gg, r_lk, p_lk, &
                    ov_eieg, ov_eiej, G_if)
            end if
        end if
    end do
    ! Add special cases for k=g or k=f
    rxp_kg = sos_eiOg(3, n_ab, n_mos, t_mo(:,:,:,id_state), Smo_irxpj)
    r_fk = sos_ejOei(3, n_ab, n_mos, t_mo(:,:,:,id_state), pfac_r)
    ! Add case k=f for first term of G
    do jx = 1, 3
        do ix = 1, 3
            G_if(ix,jx) = G_if(ix,jx) &
                + r_fk(ix)*rxp_kg(jx)*ov_eieg(id_state)
        end do
    end do
    if (TIMEIT) call write_time('LORG correction')
    if (use_giao) then
        call sos_MCD_tensor_LORG_corr( &
            n_states, n_els, id_state, id_state, r_gg, r_lk, p_lk, ov_eieg, &
            ov_eiej, G_if)
    end if
    ! Add case k=g for second term of G
    do jx = 1, 3
        do ix = 1, 3
            G_if(ix,jx) = G_if(ix,jx) &
                - r_gg(ix)*rxp_kg(jx)*ov_eieg(id_state)
        end do
    end do
    if (use_giao) then
        call sos_MCD_tensor_LORG_corr( &
            n_states, n_els, id_state, 0, r_gg, r_lk, p_lk, ov_eieg, ov_eiej, &
            G_if)
    end if

    ! Correct using a factor of 1/2
    G_if = G_if / 2.0_real64

    call sec_header(2, 'Final Value')
    write(iu_out, '(10X,"X              Y              Z")')
    do i = 1, 3
        write(iu_out, '(1x,a,3es15.6)') labXYZ_1D(i), G_if(i,:)
    end do

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

    subroutine parse_opts
        !! Parses commandline options and updates information.
        !!
        !! Builds options parser and parse user-options, setting
        !! default values where necessary
        implicit none

        character(len=512) :: outfile
        class(CmdArgDB), allocatable :: opts

        ! Build parser and check arguments
        opts = CmdArgDB(progname='mcd_tensor')
        if (opts%error%raised()) then
            write(*, '(a)') 'Failed to initialize the command-line parser'
            select type (err => opts%error)
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
            'string', err, label='filename', &
            help='Gaussian formatted checkpoint file')
        call opts%add_arg_real( &
            'real', err, longname='--gamma', &
            help='Include a shift term, gamma, to avoid a divergence of the &
                &denominator in the summation with close states (unit: Hartree).',&
            def_value=real(e_gamma), min_value=0.0)
        call opts%add_arg_int( &
            'scalar', err, shortname='-f', longname='--final', &
            help='Final excited electronic state (starting from 1)', &
            min_value=1)
        call opts%add_arg_bool( &
            'store_true', err, longname='--giao', &
            help='Include GIAO corrections (default)', &
            def_value=.True.)
        call opts%add_arg_bool( &
            'store_true', err, longname='--guvcde', &
            help='Add printing similar to SOS/GUVCDE for control.')
        call opts%add_arg_int( &
            'scalar', err, shortname='-m', longname='--max-state', &
            help='Highest electronic state to include (upper bound of summation)', &
            min_value=1)
        call opts%add_arg_bool( &
            'store_false', err, longname='--no-giao', &
            help='Include GIAO corrections')
        call opts%add_arg_char( &
            'string', err, label='output', shortname='-o', longname='--output', &
            help='Name of the file to store the output. Existing content will &
            &be overwritten.')
        
        if (TIMEIT) call write_time('Option parser')
        call opts%parse_args(err)
        if (err%raised()) then
            select type (err)
                class is (ValueError)
                    write(*, '(a)') trim(err%msg())
                class is (Error)
                    write(*, '(a)') trim(err%msg())
                class default
                    write(*, '(a)') 'Unknown error while reading "gamma"'
            end select
            stop
        end if
        if (opts%is_user_set('output', err)) then
            call opts%get_value('output', outfile, err)
            open(newunit=iu_out, file=outfile, action='write')
            call sec_header(-1, 'MCD Tensor Calculator')
        end if
        call sec_header(1, 'Simulation Parameters')
        write(iu_out, '(1x)')
        if(DEBUG) call write_param("Debugging mode", .true.)
        if(TIMEIT) call write_param("Timer", .true.)


        if (opts%is_user_set('giao', err) &
            .and. opts%is_user_set('no-giao', err)) then
            write(*, '(" Error: conflicting option for the definition of GIAO")')
            stop
        end if
        
        if (opts%is_user_set('no-giao', err)) use_giao = .False.
        call write_param('GIAO correction', use_giao)

        if (opts%is_user_set('guvcde', err)) for_guvcde = .True.
        call write_param('Compatibility mode with GUVCDE', for_guvcde)

        call opts%get_value('gamma', e_gamma, err)
        use_gamma = abs(e_gamma) > tiny(e_gamma)
        call write_param('Correction term against degeneracies', use_gamma)
        if (use_gamma) &
            call write_param('Value of the term (in Hartrees)', e_gamma, .true.)

        call opts%get_value('filename', fname, err)
        inquire(file=fname, exist=exists)
        if (.not.exists) then
            write(*, '("Error: File ",a," does not exist.")') trim(fname)
            stop
        end if
        call write_param('Input filename', fname)

        if (opts%is_user_set('final', err)) then
            call opts%get_value('final', fstate, err)
            call write_param('Final electronic state', fstate)
        else
            fstate = -1
            call write_param('Final electronic state', 'automatic')
        end if

        if (opts%is_user_set('max-state', err)) then
            call opts%get_value('max-state', max_states, err)
            call write_param('Highest excited state', max_states)
        else
            max_states = -1
            call write_param('Highest excited state', 'include all')
        end if

    end subroutine parse_opts

    subroutine write_moldata
        !! Writes molecular data.
        !!
        !! Write some of the molecular data stored in the moldata module
        !! for mcd_tensor.
        implicit none

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

        if (openshell) then
            tag_open = 'open'
        else
            tag_open = 'closed'
        end if

        if (pureD) then
            tag_bfunD = 'pure'
        else
            tag_bfunD = 'Cartesian'
        end if
        if (pureF) then
            tag_bfunF = 'pure'
        else
            tag_bfunF = 'Cartesian'
        end if

        write(iu_out, 1100) n_at, charge, multip, tag_open, n_basis, &
            tag_bfunD, tag_bfunF, n_ao

        call sec_header(2, 'Atomic coordinates')
        call prt_coord(n_at, at_lab, at_crd, at_mas)
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
        real(real64), intent(in) :: value
        !! Value associated to parameter
        logical, intent(in), optional :: subparam
        !! Parameter is a sub-parameter type

        1000 format(a,1x,f0.4)
        1010 format(a,es13.6)

        if (abs(value) < 1.0e-2_real64 .or. abs(value) > 1.0e4_real64) then
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

        integer, dimension(8) :: ia_dtime

        1000 format('Entering: ',a,' - Date: ',i4,'/',i2.2,'/',i2.2,' at ', &
            i2.2,':',i2.2,':',i2.2)
        call date_and_time(values=ia_dtime)
        write(iu_out, 1000) trim(label), ia_dtime(1:3), ia_dtime(5:7)
    end subroutine write_time

end program mcd_tensor
