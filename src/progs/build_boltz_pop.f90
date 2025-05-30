program build_boltz_pop
    use numeric, only: f0, f1, f10p2, realwp
    use input, only: DataFile
    use string, only: num_chars_int
    use parse_cmdline, only: CmdArgDB
    use output, only: iu_out, sec_header, write_err
    use exception, only: BaseException, runstat
    use datatypes, only: VibrationsDB
    use vibrational, only: boltz_pop_max_quanta, full_boltz_pop

    integer :: i, iq, j, len_mode, max_modes, n, n_states
    integer, dimension(:), allocatable :: nq_modes, nq_max
    integer, dimension(:,:), allocatable :: nq_list
    real(realwp) :: rho_min, rho_sum, rho_tot, T_ref
    real(realwp), dimension(:), allocatable :: bz_list
    character(len=80) :: fmt_state, fmt_line
    character(len=1024) :: line
    character(len=:), allocatable :: infile, outfile
    type(CmdArgDB) :: opts
    type(DataFile) :: dfile
    type(VibrationsDB) :: vibDB
    class(BaseException), allocatable :: err

    ! Build option parser for commandline
    opts = CmdArgDB(progname='build_boltz_pop')
    if (opts%has_error()) then
        err = opts%exception()
        call runstat%raise_error( &
            'Unable to initialize the commandline parser', &
            err%msg())
    end if
    call opts%add_arg_char( &
        'string', label='filename', &
        help='Gaussian formatted checkpoint file.')
    call opts%add_arg_int( &
        'int', shortname='-n', longname='--num-modes', &
        def_value=3, &
        help='Maximum number of excited modes')
    call opts%add_arg_char( &
        'string', shortname='-o', longname='--output', &
        help='Output filename.')
    call opts%add_arg_real( &
        'real', shortname='-p', longname='--population', &
        def_value=0.1_realwp, &
        help='Population with respect to ground vibrational state, as absolute&
        & value (0.1 -> 10%).')
    call opts%add_arg_real( &
        'real', shortname='-t', longname='--temperature', &
        def_value=298.15_realwp, &
        help='Temperature.')

    call opts%parse_args()
    if (opts%has_error()) then
        err = opts%exception()
        call runstat%raise_error( &
            'Failure to parse commandline options', &
            err%msg())
    end if

    ! Check commandline arguments and set information
    call opts%get_value('filename', line)
    infile = trim(line)
    if (opts%is_user_set('output')) then
        call opts%get_value('output', line)
        outfile = trim(line)
        open(newunit=iu_out, file=outfile, action='write')
    end if

    call opts%get_value('temperature', T_ref)
    call opts%get_value('population', rho_min)
    call opts%get_value('num-modes', max_modes)

    if (T_ref <= f0) then
        call runstat%raise_error(&
            'Temperature must be strictly positive')
    end if
    if (rho_min <= f0) then
        call runstat%raise_error(&
            'Minimum population must be strictly positive')
    end if
    if (max_modes <= 0) then
        call runstat%raise_error(&
            'Maximum number of modes must be strictly positive')
    end if

    ! Start program
    call sec_header(-1, 'Test Program to Build Boltzmann population')
1000 format(/,"This program builds a list of all vibrational states with a &
            &population above",/,"a set threshold, compared to the &
            &vibrational ground state (set to 1).",//, &
            "NOTE: The Boltzmann populations are given relative to ground &
            &state's and not",/,6x,"normalized.")
    write(iu_out, 1000)

    ! Extract data from input file.
    dfile = DataFile(infile)
    if (dfile%has_error()) then
        call runstat%raise_error( &
            'Failed to initialize data file', dfile%get_error())
    end if

    vibDB = dfile%get_vib_data()
    if (dfile%has_error()) then
        call runstat%raise_error( &
            'Unable to parse vibrational data, check file.', dfile%get_error())
    end if

    ! Check consistency between max. num. excited modes and num. vibrations.
    if (max_modes > vibDB%n_vib) then
        write(line, '("Reducing max. number of modes from ",i0," to ",i0)') &
            max_modes, vibDB%n_vib
        call runstat%raise_warning('Too many modes chosen', line)
        max_modes = vibDB%n_vib
    end if

    ! Construct format to display states and populations
    len_mode = num_chars_int(vibDB%n_vib)
    write(fmt_state, '("(1x,",i0,"(i0,""^"",i0,:,"", ""))")') max_modes
    ! We make a rough estimate of how many characters we will need, considering
    ! that it is unlikely we have more than 99 quanta per mode
    ! format: i^n, j^m... so len_mode + 1 + 2(quanta) + 2 (", ")
    ! 2 removed since n_modes-1 separators
    n = max_modes*(len_mode + 5) - 2
    write(fmt_line, '("(a",i0,","" | "",f12.4,"" | "",f10.8)")') n
    line = 'State'
    write(iu_out,'(/,1x,a," |  Energy",6x,"| Population")') line(:n-1)

    allocate(nq_modes(max_modes), nq_max(max_modes))
    call sec_header(1, "Vibrational ground state")
    line = ' '
    write(line, fmt_state) 0, 0
    write(iu_out, fmt_line) line, f0, f1
    n_states = 1
    rho_sum = f1
    do n_modes = 1, max_modes
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
            call boltz_pop_max_quanta(vibDB, n_modes, nq_modes, nq_max, &
                                      T_ref, rho_min, .true., nq_list, bz_list)
            if (allocated(nq_list)) then
                n = size(bz_list)
                n_states = n_states + n
                do i = 1, n
                    line = ' '
                    write(line, fmt_state) (nq_modes(j), nq_list(j,i), &
                                            j=1, n_modes)
                    write(iu_out, fmt_line) line, &
                        sum(nq_list(:,i)*vibDB%freq(nq_modes(:n_modes))), bz_list(i)
                    rho_sum = rho_sum + bz_list(i)
                end do
            end if
            do
                nq_modes(iq) = nq_modes(iq) + 1
                if (nq_modes(iq) > vibDB%n_vib - n_modes + iq) then
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

    call sec_header(1, 'Final result')
    write(iu_out, '("> Total number of states : ",i0)') n_states
    write(iu_out, '("> Computed population : ",es14.8)') rho_sum
    rho_tot = full_boltz_pop(vibDB, T_ref, .false.)
    write(iu_out, '("> Analytical population : ",es14.8)') rho_tot
    write(iu_out, '("> Recovered percentage : ",f0.6,"%")') rho_sum/rho_tot*f10p2

end program build_boltz_pop
