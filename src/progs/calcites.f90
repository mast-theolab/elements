program calcites
    !! Calculate Contribution to Intensity of Thermally Excited States
    !!
    !! Evaluates the contributions to intensity of thermally excited
    !! vibrational states, either using exact analytic forms or by
    !! estimating the total intensity with approximate models.
    !! The intensity can depend on the spectroscopy of interest and the
    !! description model
    use iso_fortran_env, only: int64

    use blas_drv, only: xgemm
    use datatypes, only: ExcitationDB, MoleculeDB, PropertyDB, VibrationsDB
    use input, only: DataFile
    use numeric, only: f0, f1, f2, f10m1, f10p2, realwp
    use output, only: iu_out, sec_header, write_err
    use parse_cmdline, only: CmdLineArgsDB
    use physics, only: spec_conv
    use run_env, only: run
    use string, only: locase, num_chars_int, upcase
    use vibrational, only: boltz_pop_max_quanta, build_modes, full_boltz_pop

    implicit none

    type Params
        character(len=:), allocatable :: file_prp
        character(len=:), allocatable :: file_vib
        character(len=:), allocatable :: file_out
        character(len=6) :: spec = 'OPA'
        character(len=12) :: level = 'FCHT'
        integer :: max_modes = 3
        real(realwp) :: temp = 298.15_realwp
        real(realwp) :: rho_min = 0.1_realwp
        real(realwp) :: thresh_int = 0.1_realwp
    end type Params

    type Property
        integer :: id
        real(realwp), dimension(:), allocatable :: p_ref
        real(realwp), dimension(:,:), allocatable :: p_d1
    end type Property

    integer :: len_mode, n, n_vib
    integer(int64), dimension(2) :: n_states
    real(realwp) :: gs_int, prop_cst, rho_tot, thresh_int
    real(realwp), dimension(2) :: int_sum, rho_sum
    real(realwp), dimension(:), allocatable :: freq, prop_nq1
    character(len=*), parameter :: PROGNAME = 'CalCITES'
    character(len=80) :: fmt_state, fmt_line
    character(len=1024) :: line

    type(Params) :: opts
    type(Property), dimension(:), allocatable :: prop

    call parse_options(opts)

    ! Start program
    call sec_header(-1, PROGNAME)
1000 format(/,"This program builds a list of all vibrational states with a &
        &population above",/, &
        "a set threshold, compared to the vibrational ground state (set to 1) &
        &and a",/, &
        "relative intensity with respect to that from the ground state, as an &
        &absolute",/, &
        "value",//, &
        "NOTE: The Boltzmann populations are given relative to ground state's &
        &and not",/,6x,"normalized.")
    write(iu_out, 1000)

    call extract_data(opts, n_vib, freq, prop)

    call get_spectro_coeffs(opts, freq, prop, prop_cst, prop_nq1)

    n_vib = size(freq)

    ! Construct format to display states and populations
    len_mode = num_chars_int(n_vib)
    write(fmt_state, '("(1x,",i0,"(i0,""^"",i0,:,"", ""))")') opts%max_modes
    ! We make a rough estimate of how many characters we will need, considering
    ! that it is unlikely we have more than 99 quanta per mode
    ! format: i^n, j^m... so len_mode + 1 + 2(quanta) + 2 (", ")
    ! 2 removed since n_modes-1 separators
    n = opts%max_modes*(len_mode + 5) - 2
    write(fmt_line, '("(a",i0,","" | "",f12.4,"" | "",f10.8,"" | "",&
        &es13.6)")') n
    line = 'State'
    write(iu_out,'(/,1x,a," |  Energy",6x,"| Population |  Intensity")') &
        line(:n-1)

    call sec_header(1, "Vibrational ground state")
    line = ' '
    gs_int = prop_cst
    write(line, fmt_state) 0, 0
    write(iu_out, fmt_line) line, f0, f1, gs_int
    n_states = 1
    rho_sum = f1
    int_sum = gs_int
    thresh_int = abs(gs_int) * opts%thresh_int

    call build_states(opts, n_vib, freq, thresh_int, prop_cst, prop_nq1, &
                      n_states, rho_sum, int_sum)


    call sec_header(1, 'Final result')

    rho_tot = full_boltz_pop(freq, opts%temp, .false.)
    write(iu_out, &
        '("                         | Pop. criterion |  All criteria  |&
        & Fraction 2/1")')
    write(iu_out, &
        '("> Number of states       | ",i14," | ",i14," | ",f10.6,"%")') &
        n_states, real(n_states(2), kind=realwp) &
                  / real(n_states(1), kind=realwp)*f10p2
    write(iu_out, &
        '("> Computed population    | ",es14.8," | ",es14.8," | ",f10.6,"%")')&
        rho_sum, rho_sum(2)/rho_sum(1)*f10p2
    write(iu_out, &
        '("> Computed intensity     | ",es14.8," | ",es14.8," | ",f10.6,"%")')&
        int_sum, int_sum(2)/int_sum(1)*f10p2
    write(iu_out, &
        '("> Recovered population   | ",f13.6,"% | ",f13.6,"% |")') &
        rho_sum(1)/rho_tot*f10p2, rho_sum(2)/rho_tot*f10p2

    write(iu_out, '(" Total analytical population: ",es14.8)') rho_tot

contains

! ======================================================================

subroutine build_states(opts_db, nvib, omega, int_min, coef_cst, coef_nq1, &
                        nstates, calc_rho, calc_int)
    !! Build vibrational states and compute intensity.
    !!
    !! Builds vibrational states and computes the relative intensity.
    !!
    !! @note
    !! Each of the counters (`nstates`, `calc_rho`, `calc_int`) stores 2
    !! internal counters:
    !! 1. counter considering only the thresholds on the Boltzmann pop.
    !! 2. counter considering both selection criteria.
    !! @endnote
    !!
    !! @warning
    !! The ground-state contributions are ignored in this routine,
    !! `nstates_bz`, `nstates_ok`, `calc_rho` and `calc_int` are assumed
    !! already initialized with sensible values and updated outside.
    !! @endwarning
    type(Params), intent(in) :: opts_db
        !! Options parameters.
    integer, intent(in) :: nvib
        !! Number of normal modes.
    real(realwp), dimension(:), intent(in) :: omega
        !! Harmonic wavenumbers, in cm<sup>-1</sup>.
    real(realwp), intent(in) :: int_min
        !! Minimum intensity, related to the ground-state value.
    real(realwp), intent(in) :: coef_cst
        !! Constant coefficient in intensity calculation, independent of |v>.
    real(realwp), dimension(:), intent(in) :: coef_nq1
        !! Coefficients dependent on |v>.
    integer(int64), dimension(2), intent(inout) :: nstates
        !! Number of states treated, updated.
    real(realwp), dimension(2), intent(inout) :: calc_rho
        !! Sum of the relative Boltzmann population, updated.
    real(realwp), dimension(2), intent(inout) :: calc_int
        !! Sum of the intensities, updated.

    integer :: i, id_spec, id_level, iq, j, max_modes, n_combis, n_modes
    integer, dimension(:), allocatable :: nq_modes, nq_max
    integer, dimension(:,:), allocatable :: nq_list
    real(realwp) :: val, contrib
    real(realwp), dimension(:), allocatable :: bz_list

    select case (upcase(trim(opts_db%spec)))
        case ('OPA')
            id_spec = 1
        case ('OPE')
            id_spec = 2
        case ('ECD')
            id_spec = 3
        case ('CPL')
            id_spec = 4
        case ('RR')
            id_spec = 5
        case ('RROA')
            id_spec = 6
        case default
            call run%error%raise_argerror('val', 'Unsupported spectroscopy')
    end select
    select case (upcase(trim(opts_db%level)))
        case ('FC')
            id_level = 0
        case ('FCHT')
            id_level = 1
        case ('HT')
            id_level = 2
        case default
            call run%error%raise_argerror('val', 'Unsupported level')
    end select
    max_modes = opts_db%max_modes

    allocate(nq_modes(max_modes), nq_max(max_modes))

    do n_modes = 1, opts_db%max_modes
        if (n_modes <= 1) then
            write(line, '(i0," Excited mode")') n_modes
        else
            write(line, '(i0," Excited modes")') n_modes
        end if
        call sec_header(1, trim(line))
        do i = 1, n_modes
            nq_modes(i) = i
        end do
        iq = n_modes
        mainloop: do
            call boltz_pop_max_quanta(nvib, n_modes, omega, nq_modes, nq_max, &
                                      opts_db%temp, opts_db%rho_min, .true., &
                                      nq_list, bz_list)
            if (allocated(nq_list)) then
                n_combis = size(bz_list)
                nstates(1) = nstates(1) + n_combis
                do i = 1, n_combis
                    calc_rho(1) = calc_rho(1) + bz_list(i)
                    contrib = get_spectro_contrib( &
                        nq_modes(:n_modes), nq_list(:,i), id_spec, id_level, &
                        coef_cst, coef_nq1)
                    val = bz_list(i)*contrib
                    calc_int(1) = calc_int(1) + val
                    if (abs(val) >= int_min) then
                        nstates(2) = n_states(2) + 1
                        line = ' '
                        write(line, fmt_state) (nq_modes(j), nq_list(j,i), &
                                                j=1, n_modes)
                        write(iu_out, fmt_line) line, &
                            sum(nq_list(:,i)*freq(nq_modes(:n_modes))), &
                            bz_list(i), val
                        calc_rho(2) = calc_rho(2) + bz_list(i)
                        calc_int(2) = calc_int(2) + val
                    end if
                end do
            end if
            do
                nq_modes(iq) = nq_modes(iq) + 1
                if (nq_modes(iq) > n_vib - n_modes + iq) then
                    iq = iq - 1
                    if (iq == 0) exit mainloop
                else
                    exit
                end if
            end do
            if (iq < n_modes) then
                do i = iq+1, n_modes
                    nq_modes(i) = nq_modes(iq) + i - iq
                end do
                iq = n_modes
            end if
        end do mainloop

    end do
end subroutine build_states

! ======================================================================

subroutine extract_data(opts_db, nvib, omega, props)
    !! Extract all relevant data.
    !!
    !! Extracts all relevant data
    type(Params), intent(in) :: opts_db
        !! Options parameters.
    integer, intent(out) :: nvib
        !! Number of normal modes.
    real(realwp), dimension(:), allocatable, intent(out) :: omega
        !! Harmonic wavenumbers, in cm<sup>-1</sup>.
    type(Property), dimension(:), allocatable, intent(out) :: props
        !! List of extracted properties of interest.

    integer :: istate_exc, n_at3
    real(realwp), dimension(:), allocatable :: tmp_vec
    real(realwp), dimension(:,:), allocatable :: tmp_mat

    type(DataFile) :: dfile_prp, dfile_vib
    type(ExcitationDB) :: exc_db
    type(VibrationsDB) :: vib_db
    type(MoleculeDB) :: mol_db0
    type(PropertyDB) :: ffx_db, p0_db, p_d1_db

    ! Extract data for the definition of the vibrational states
    dfile_vib = DataFile(opts_db%file_vib)

    vib_db = dfile_vib%get_vib_data()

    mol_db0 = dfile_vib%get_mol_data()

    n_at3 = 3*mol_db0%n_at
    nvib = vib_db%n_vib
    ffx_db = dfile_vib%get_data(1, derorder=2)
    call build_modes(ffx_db%data, mol_db0, vib_db, set_Lmweigh=.true.)
    omega = vib_db%freq
    call ffx_db%clear()

    ! Initialize file storing properties of interest.
    dfile_prp = DataFile(opts_db%file_prp)

    select case (upcase(trim(opts_db%spec)))
        case ('OPA', 'RR')
            allocate(props(1))
            props(1)%id = 101
            exc_db = dfile_prp%get_exc_data()
            istate_exc = exc_db%id_state
            if (upcase(opts_db%level(:2)) == 'FC') then
                p0_db = dfile_prp%get_data(props(1)%id, 0, istate_exc)
                if (.not.p0_db%loaded) then
                    call run%check(dfile_prp%error, &
                                   'Unable to extract property')
                end if
                ! Now convert property quantities
                props(1)%p_ref = p0_db%data
                ! Clear memory
                call p0_db%clear()
            else
                allocate(props(1)%p_ref(3))
                props(1)%p_ref = f0
            end if
            if (index(upcase(trim(opts_db%level)), 'HT') > 0) then
                p_d1_db = dfile_prp%get_data(props(1)%id, 0, istate_exc, 1)
                if (.not.p_d1_db%loaded) then
                    call run%check(dfile_prp%error, 'Unable to extract property')
                end if
                tmp_mat = reshape(p_d1_db%data, [p_d1_db%shape(1), p_d1_db%shape(2)])
                allocate(props(1)%p_d1(3,nvib))
                call xgemm('N', 'N', 3, nvib, n_at3, f1, tmp_mat, 3, &
                        vib_db%L_mwg, n_at3, f0, props(1)%p_d1, 3)
                ! Clear memory
                deallocate(tmp_mat)
                call p_d1_db%clear()
                ! Now correct to be function of q instead of Q
                tmp_vec = spec_conv%mwq2q(f1, omega, .true.)
                props(1)%p_d1(1,:) = props(1)%p_d1(1,:) * tmp_vec
                props(1)%p_d1(2,:) = props(1)%p_d1(2,:) * tmp_vec
                props(1)%p_d1(3,:) = props(1)%p_d1(3,:) * tmp_vec
                deallocate(tmp_vec)
            end if
        case default
            call run%error%raise_error('opt', 'value', &
                'unsupported spectroscopy.')
    end select

end subroutine extract_data

! ======================================================================

subroutine get_spectro_coeffs(opts_db, omega, props, coef_cst, coef_nq)
    !! Get spectroscopy-related coefficients.
    !!
    !! Gets the spectroscopy-related coefficients.
    type(Params), intent(in) :: opts_db
        !! Options parameters.
    real(realwp), dimension(:), allocatable, intent(out) :: omega
        !! Harmonic wavenumbers, in cm<sup>-1</sup>.
    type(Property), dimension(:), intent(in) :: props
        !! List of properties of interest.
    real(realwp), intent(out) :: coef_cst
        !! Constant coefficient in intensity calculation, independent of |v>
    real(realwp), dimension(:), allocatable, intent(out) :: coef_nq
        !! Coefficients dependent on |v>.

    integer :: i, nvib

    nvib = size(omega)

    select case(opts_db%spec)
        case ('OPA')
            if (props(1)%id /= 101) &
                call run%error%raise_error('data', 'missing', &
                    'unexpected property found. Aborting.')
            select case(opts_db%level)
                case ('FC')
                    coef_cst = sum(props(1)%p_ref**2)
                case ('FCHT')
                    coef_cst = sum(props(1)%p_ref**2) &
                        + sum(props(1)%p_d1**2)/f2
                    allocate(coef_nq(nvib))
                    do i = 1, nvib
                        coef_nq(i) = sum(props(1)%p_d1(:,i)**2)/f2
                    end do
                case ('HT')
                    coef_cst = sum(props(1)%p_d1**2)/f2
                    allocate(coef_nq(nvib))
                    do i = 1, nvib
                        coef_nq(i) = sum(props(1)%p_d1(:,i)**2)/f2
                    end do
                case default
                    call run%error%raise_deverror('nyi', &
                        'unsupported level for OPA.')
            end select
        case default
            call run%error%raise_error('opt', 'value', &
                'unsupported spectroscopy.')
    end select

end subroutine get_spectro_coeffs

! ======================================================================

function get_spectro_contrib(nq_modes, nq_quanta, id_spectro, id_level, &
                             coef_cst, coef_nq1) result(contrib)
    !! Get spectroscopic contribution.
    !!
    !! Gets spectroscopic contributions for a given vibrational state.
    !!
    !! @note
    !! The contribution ignore the Boltzmann population, which needs to
    !! be added separately.
    !! @endnote
    integer, dimension(:), intent(in) :: nq_modes
        !! List of excited modes.
    integer, dimension(:), intent(in) :: nq_quanta
        !! List of quanta for the excited modes.
    integer, intent(in) :: id_spectro
        !! Identifier of the spectroscopy.
    integer, intent(in) :: id_level
        !! Identifier of the level of approximation
    real(realwp), intent(in) :: coef_cst
        !! Constant coefficient
    real(realwp), dimension(:), intent(in) :: coef_nq1
        !! Coefficients dependent on |v>.
    real(realwp) :: contrib
        !! Computed contribution to intensity

    select case (id_spectro)
        case (1:4)
            if (id_level == 0) then
                contrib = coef_cst
            else
                contrib = coef_cst + sum(coef_nq1(nq_modes)*(2*nq_quanta))
            end if
        case default
            call run%error%raise_deverror('nyi', 'unsupported spectroscopy')
    end select

end function get_spectro_contrib

! ======================================================================

subroutine parse_options(opts_db)
    type(Params), intent(out) :: opts_db

    character(len=1024) :: argval
    type(CmdLineArgsDB) :: parser

    ! Build option parser for commandline
    parser = CmdLineArgsDB(progname=locase(PROGNAME))
    call run%check(parser%error, 'Unable to initialize the commandline parser')
    call parser%add_arg_char( &
        'string', label='file_vib', &
        help='Gaussian formatted checkpoint file containing the description &
            &of the vibrational states.')
    call parser%add_arg_char( &
        'string', label='file_prp', &
        help='Gaussian formatted checkpoint file containing data on the &
            &properties of interest.')
    call parser%add_arg_real( &
        'real', shortname='-m', longname='--min-intensity', &
        def_value=opts%thresh_int, &
        help='Minimum intensity with respect to the ground-state contribution&
            & (ex: 0.1 -> 10% of the intensity from |0>)')
    call parser%add_arg_int( &
        'int', shortname='-n', longname='--num-modes', &
        def_value=opts%max_modes, &
        help='Maximum number of excited modes')
    call parser%add_arg_char( &
        'string', shortname='-o', longname='--output', &
        help='Output filename.')
    call parser%add_arg_real( &
        'real', shortname='-p', longname='--population', &
        def_value=opts%rho_min, &
        help='Population with respect to ground vibrational state, as absolute&
        & value (0.1 -> 10%).')
    call parser%add_arg_real( &
        'real', shortname='-t', longname='--temperature', &
        def_value=opts%temp, help='Temperature.')

    call parser%parse_args()
    call run%check(parser%error, 'Failure to parse commandline options')

    ! Check commandline arguments and set information
    call parser%get_value('file_vib', argval)
    opts_db%file_vib = trim(argval)
    call parser%get_value('file_prp', argval)
    opts_db%file_prp = trim(argval)
    if (parser%is_user_set('output')) then
        call parser%get_value('output', argval)
        opts_db%file_out = trim(argval)
        open(newunit=iu_out, file=opts_db%file_out, action='write')
    end if

    call parser%get_value('temperature', opts_db%temp)
    call parser%get_value('population', opts_db%rho_min)
    call parser%get_value('min-intensity', opts_db%thresh_int)
    call parser%get_value('num-modes', opts_db%max_modes)

    if (opts_db%temp <= f0) then
        call run%error%raise_error('opt', 'value', &
            'temperature must be strictly positive')
    end if
    if (opts_db%rho_min <= f0) then
        call run%error%raise_error('opt', 'value', &
            'minimum population must be strictly positive')
    end if
    if (opts_db%thresh_int <= f0) then
        call run%error%raise_error('opt', 'value', &
            'minimum relative intensity must be strictly positive')
    end if
    if (opts_db%max_modes <= 0) then
        call run%error%raise_error('opt', 'value', &
            'maximum number of modes must be strictly positive')
    end if

end subroutine parse_options

! ======================================================================

end program calcites
