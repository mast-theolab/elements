program test_getdata
    use iso_fortran_env, only: output_unit
    use datatypes, only: PropertyDB
    use input, only: DataFile
    use output, only: iu_out, prt_mat, write_err, sec_header
    use parse_cmdline, only: CmdLineArgsDB
    use run_env, only: run

    integer :: i, iprp, ider, LP
    integer, parameter :: MAXFILES = 1
    logical :: auto, do_g2e
    character(len=:), dimension(:), allocatable :: files
    character(len=256) :: outfile, string, fmt
    character(len=3), dimension(4) :: &
        ordinal = ['1st', '2nd', '3rd', '4th']
    type(DataFile) :: dfile
    class(PropertyDB), allocatable :: prop
    class(CmdLineArgsDB), allocatable :: opts

    opts = CmdLineArgsDB(progname='test_read_vib')
    call run%check(opts%error, 'Failed to initialize the command-line parser')

    call opts%add_arg_int( &
        'scalar', label='derord', shortname='-d', &
        longname='--derorder', &
        help='Property derivative (0: no derivative)')
    call opts%add_arg_char( &
        'list', label='file', shortname='-f', longname='--file', &
        required=.true., max_nvals=MAXFILES, help='Filename')
    call opts%add_arg_bool( &
        'store_true', label='eltrans', longname='--g2e', &
        help='Return the ground to excited transition moment')
    call opts%add_arg_char( &
        'string', label='output', shortname='-o', &
        longname='--output', help='Optional output')
    call opts%add_arg_int( &
        'scalar', label='propid', longname='--propid', &
        help='Property id to read (as integer)')

    call opts%parse_args()
    call run%check(opts%error, 'Failed to parse command-line arguments')
    if (opts%is_user_set('output')) then
        call opts%get_value('output', string)
        outfile = trim(string)
        open(newunit=iu_out, file=outfile, action='write')
    else
        iu_out = output_unit
    end if
    if (opts%is_user_set('file')) then
        call opts%get_value('file', files)
    end if

    iprp = 0
    if (opts%is_user_set('propid')) then
        call opts%get_value('propid', iprp)
        if (opts%is_user_set('derord')) then
            call opts%get_value('derord', ider)
        else
            ider = 0
        end if
        if (opts%is_user_set('eltrans')) then
            do_g2e = .true.
        else
            do_g2e = .false.
        end if
        auto = .false.
    else
        auto = .true.
    end if

    dfile = DataFile(files(1))
    call run%check(dfile%error, 'Error found while initializing data file')

    call sec_header(0, 'Test Program for DataFile%GetData')

    write(iu_out, '(/," Filename: ",a)') trim(dfile%get_filename())

    if (.not.auto) then
        if (iprp > 0) then
            fmt = '("Reading property num. ",i0,:," (",a," derivative)")'
            if (ider > 0) then
                write(string, fmt) iprp, ordinal(ider)
            else
                write(string, fmt) iprp
            end if
            call sec_header(1, trim(string))
            if (do_g2e) then
                prop = dfile%get_data(iprp, 0, -1)
            else
                prop = dfile%get_data(iprp, 0)
            end if
            write(iu_out, '("Shape: ",*(i0,:,", "))') &
                (prop%shape(i), i=1, size(prop%shape))
            write(iu_out, '(5es16.8)') prop%data
        else
            write(iu_out, '(a)') 'Unsupported property'
            stop
        end if
    else
        call sec_header(1, 'Ground to excited states electric dipole (len.)')
        call sec_header(2, 'Reference quantity')
        prop = dfile%get_data(101, 0, -1)
        LP = product(prop%pdim)
        do i = 1, prop%shape(2)
            write(iu_out, '(i6, *(es16.8))') i, prop%data((i-1)*LP+1:i*LP)
        end do

        call sec_header(2, '1st derivative')
        prop = dfile%get_data(101, 0, -1, 1)
        LP = product(prop%pdim)
        do i = 1, prop%shape(2)
            write(iu_out, '(i6, *(es16.8))') i, prop%data((i-1)*LP+1:i*LP)
        end do

        call sec_header(1, 'Ground to excited states electric dipole (vel.)')
        call sec_header(2, 'Reference quantity')
        prop = dfile%get_data(111, 0, -1)
        LP = product(prop%pdim)
        do i = 1, prop%shape(2)
            write(iu_out, '(i6, *(es16.8))') i, prop%data((i-1)*LP+1:i*LP)
        end do

        call sec_header(2, '1st derivative')
        prop = dfile%get_data(111, 0, -1, 1)
        LP = product(prop%pdim)
        do i = 1, prop%shape(2)
            write(iu_out, '(i6, *(es16.8))') i, prop%data((i-1)*LP+1:i*LP)
        end do

        call sec_header(1, 'Ground to excited states magnetic dipole')
        call sec_header(2, 'Reference quantity')
        prop = dfile%get_data(102, 0, -1)
        LP = product(prop%pdim)
        do i = 1, prop%shape(2)
            write(iu_out, '(i6, *(es16.8))') i, prop%data((i-1)*LP+1:i*LP)
        end do

        call sec_header(2, '1st derivative')
        prop = dfile%get_data(102, 0, -1, 1)
        LP = product(prop%pdim)
        do i = 1, prop%shape(2)
            write(iu_out, '(i6, *(es16.8))') i, prop%data((i-1)*LP+1:i*LP)
        end do

        call sec_header(1, 'Ground to excited states electric quadrupole')
        call sec_header(2, 'Reference quantity')
        prop = dfile%get_data(107, 0, -1)
        LP = product(prop%pdim)
        do i = 1, prop%shape(3)
            write(iu_out, '(i6, *(3es16.8,:,/,6x))') &
                i, prop%data((i-1)*LP+1:i*LP)
        end do

        call sec_header(2, '1st derivative')
        prop = dfile%get_data(107, 0, -1, 1)
        LP = product(prop%pdim)
        do i = 1, prop%shape(3)
            write(iu_out, '(i6, *(3es16.8,:/,6x))') &
                i, prop%data((i-1)*LP+1:i*LP)
        end do
    end if

end program test_getdata
