program test_getdata
    use iso_fortran_env, only: output_unit
    use exception, only: AllocateError, ArgumentError, Error, ValueError
    use input, only: DataFile
    use parse_cmdline, only: CmdArgDB
    use output, only: iu_out, prt_mat, write_err, sec_header
    use datatypes, only: PropertyDB

    integer :: i, iprp, ider, LP
    integer, parameter :: MAXFILES = 1
    real :: rval
    logical :: auto, do_g2e
    character(len=:), dimension(:), allocatable :: files
    character(len=256) :: outfile, string, fmt
    character(len=3), dimension(4) :: &
        ordinal = ['1st', '2nd', '3rd', '4th']
    type(DataFile) :: dfile
    class(PropertyDB), allocatable :: prop
    class(CmdArgDB), allocatable :: opts

    opts = CmdArgDB(progname='test_read_vib')
    if (opts%has_error()) then
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
    if (opts%has_error()) then
        select type (err => opts%exception())
            class is (ValueError)
                write(*, '(a)') trim(opts%get_error())
            class is (Error)
                write(*, '(a)') trim(opts%get_error())
            class default
                write(*, '(a)') 'Unknown error while parsing options'
        end select
        stop
    end if
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
    if (dfile%has_error()) then
        call write_err('std', 'Error found while initializing data file', &
                       dfile%get_error())
        stop 1
    end if

    call sec_header(0, 'Test Program for DataFile%GetData')

    write(iu_out, '(/," Filename: ",a)') trim(dfile%get_name())
    
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
