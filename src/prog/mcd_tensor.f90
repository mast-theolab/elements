program mcd_tensor
    use iso_fortran_env, only: real64, output_unit
    use input, only: build_moldata, build_transdata
    use parse_cmdline, only: CmdArgDB
    use output, only: prt_mat, iu_out, write_err, num_digits_int
    use exception, only: BaseException, Error, ValueError
    use gmcd_legacy, only: write_control, build_MOs
    use basisset, only: fix_norm_AOs, chk_bset_redundancy
    use electronic, only: convert_AO2MO, eltrans_amp, overlap_ao_1e, ovij_1e
    use exc_sos, only: sos_eiOg, sos_ejOei, sos_prefac_ejOei, &
        sos_MCD_tensor_LORG_corr
    use moldata
    use transdata

    implicit none

    integer :: i, iab, istate, ix, jstate, jx, qty_flag
    integer, dimension(8) :: ia_dtime
    real(real64) :: de, ef, ei
    real(real64) :: e_gamma = .02_real64
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
    logical :: for_guvcde = .False., use_gamma = .True., use_giao = .True., &
        exists
    character(len=512) :: fmt_elstate, fname
    character(len=*), parameter :: fmt_dtime = &
        '("Entering: ",a," - Date: ",i4,"/",i2.2,"/",i2.2," at &
            &",i2.2,":",i2.2,":",i2.2)'
    class(ovij_1e), allocatable :: ao_int
    class(BaseException), allocatable :: err
    class(CmdArgDB), allocatable :: opts

    ! Build parser and check arguments
    opts = CmdArgDB(progname='mcd_tensor')
    call opts%add_arg_char( &
        'string', err, label='filename', &
        help='Gaussian formatted checkpoint file')
    call opts%add_arg_bool( &
        'store_true', err, longname='--giao', &
        help='Include GIAO corrections (default)', &
        def_value=.True.)
    call opts%add_arg_bool( &
        'store_false', err, longname='--no-giao', &
        help='Include GIAO corrections')
    call opts%add_arg_real( &
        'real', err, longname='--gamma', &
        help='Include a shift term, gamma, to avoid a divergence of the &
            &denominator in the summation with close states (unit: Hartree).',&
        def_value=0.02, min_value=0.0)
    call opts%add_arg_bool( &
        'store_true', err, longname='--guvcde', &
        help='Add printing similar to SOS/GUVCDE for control.')
    
    if (TIMEIT) then
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "Option parser", ia_dtime(1:3), &
            ia_dtime(5:7)
    end if
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

    if (opts%is_user_set('giao', err) &
        .and. opts%is_user_set('no-giao', err)) then
        write(*, '(" Error: conflicting option for the definition of GIAO")')
        stop
    end if
    
    if (opts%is_user_set('no-giao', err)) use_giao = .False.

    if (opts%is_user_set('guvcde', err)) for_guvcde = .True.

    call opts%get_value('gamma', e_gamma, err)
    use_gamma = abs(e_gamma) > tiny(e_gamma)

    call opts%get_value('filename', fname, err)
    inquire(file=fname, exist=exists)
    if (.not.exists) then
        write(*, '("Error: File ",a," does not exist.")') trim(fname)
        stop
    end if

    if (TIMEIT) then
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "Molecular data", ia_dtime(1:3), &
            ia_dtime(5:7)
    end if
    call build_moldata(fname, err)
    if (err%raised()) then
        select type(err)
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

    if (for_guvcde) then
        if (TIMEIT) then
            call date_and_time(values=ia_dtime)
            write(output_unit, fmt_dtime) "GUVCDE write control", &
                ia_dtime(1:3), ia_dtime(5:7)
        end if
        call write_control(output_unit, err)
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
    if (TIMEIT) then
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "Check AO normalization", ia_dtime(1:3), &
            ia_dtime(5:7)
    end if
    call fix_norm_AOs(output_unit, n_at, n_ao, at_crd, nprim_per_at, bset_info, &
                      err)
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

    if (TIMEIT) then
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "1-e AO integral", ia_dtime(1:3), &
            ia_dtime(5:7)
    end if
    ! overlap6
    qty_flag = 2**10 - 1
    ! qty_flag = ibset(qty_flag, 0)
    call overlap_ao_1e(output_unit, n_at, n_ao, qty_flag, .True., .True., &
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

    ! convert AO -> MO quantities of interest
    if (for_guvcde) then
        if (TIMEIT) then
            call date_and_time(values=ia_dtime)
            write(output_unit, fmt_dtime) "GUVCDE build MOs", ia_dtime(1:3), &
                ia_dtime(5:7)
        end if
        call build_MOs(n_ab, n_ao, n_mos, coef_mos, .true.)
    endif
    ! call convert_AO2MO(n_ao, n_mo, coef_mos(1), ao_int%i_j, Smo_ij)
    allocate(Smo_irj(3,n_mo,n_mo,n_ab), Smo_ipj(3,n_mo,n_mo,n_ab), &
             Smo_irxpj(3,n_mo,n_mo,n_ab))
    allocate(tmp_ao_arr1(n_ao,n_ao), tmp_ao_arrN(3,n_ao,n_ao))
    if (TIMEIT) then
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "AO 2 MO conversion", ia_dtime(1:3), &
            ia_dtime(5:7)
    end if
    do iab = 1, n_ab
        call convert_AO2MO(n_ao, n_mos(iab), coef_mos(:,:,iab), ao_int%i_r_j, &
                           Smo_irj(:,:,:,iab), tmp_ao_arrN)
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "<i|p|j>", ia_dtime(1:3), &
            ia_dtime(5:7)
        call convert_AO2MO(n_ao, n_mos(iab), coef_mos(:,:,iab), ao_int%i_p_j, &
                           Smo_ipj(:,:,:,iab), tmp_ao_arrN)
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "<i|rxp|j>", ia_dtime(1:3), &
            ia_dtime(5:7)
        call convert_AO2MO(n_ao, n_mos(iab), coef_mos(:,:,iab), ao_int%i_rxp_j, &
                           Smo_irxpj(:,:,:,iab), tmp_ao_arrN)
    end do
    ! call convert_AO2MO(n_ao, n_mo, coef_mos(1), ao_int%i_j)

    if (TIMEIT) then
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "Transition data", ia_dtime(1:3), &
            ia_dtime(5:7)
    end if
    call build_transdata(fname, n_ab, n_basis, err)
    if (err%raised()) then
        select type(err)
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

    ! call prt_mat(g2e_dens(:,:,1,1), n_basis, n_basis)

    if (TIMEIT) then
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "Build trans. amplitudes", ia_dtime(1:3), &
            ia_dtime(5:7)
    end if
    write(iu_out, '(a)') 'Building transition amplitudes'
    allocate(t_mo(n_mo, n_mo, n_ab, n_states))
    ! "Now doing state num. ",i3," out of 123
    write(fmt_elstate, '("(''Now doing state num. '',i",i0,",'' out of '',&
        &i0,''.'')")') num_digits_int(n_states)
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

    write(iu_out, '(a)') 'Computing the ground-states moments'

    if (TIMEIT) then
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "Ground-state moments", ia_dtime(1:3), &
            ia_dtime(5:7)
    end if
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

    write(iu_out, '(a)') 'Computing the transition energies factors'

    if (TIMEIT) then
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "Transition energies", ia_dtime(1:3), &
            ia_dtime(5:7)
    end if
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

    write(iu_out, '(4(a,/),6(/,a))') &
        'Computing tensor G_if:', &
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
        if (TIMEIT) then
            call date_and_time(values=ia_dtime)
            write(output_unit, fmt_dtime) "GIAO terms", ia_dtime(1:3), &
                ia_dtime(5:7)
        end if
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
    if (TIMEIT) then
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "Prefactors", ia_dtime(1:3), &
            ia_dtime(5:7)
    end if
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
    if (TIMEIT) then
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "Transition moments", ia_dtime(1:3), &
            ia_dtime(5:7)
    end if
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
    if (TIMEIT) then
        call date_and_time(values=ia_dtime)
        write(output_unit, fmt_dtime) "LORG correction", ia_dtime(1:3), &
            ia_dtime(5:7)
    end if
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

    write(*, '("MCD G tensor",/)')
    write(*, '(3es15.6)') G_if

end program mcd_tensor
