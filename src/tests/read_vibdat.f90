program test_read_vibdata
    use iso_fortran_env, only: output_unit
    use exception, only: AllocateError, ArgumentError, Error, ValueError
    use input, only: DataFile
    use parse_cmdline, only: CmdArgDB
    use output, only: prt_mat, write_err
    use workdata, only: MoleculeDB, VibrationsDB

    integer :: ifile, iu_out, nfiles
    integer, parameter :: MAXFILES = 2
    character(len=:), dimension(:), allocatable :: files
    character(len=256) :: outfile, string
    type(DataFile) :: dfile
    type(MoleculeDB), dimension(:), allocatable :: mols
    type(VibrationsDB), dimension(:), allocatable :: vibs
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
    call opts%add_arg_char( &
        'list', label='file', shortname='-f', longname='--file', &
        required=.true., max_nvals=MAXFILES, help='Filename')
    call opts%add_arg_char( &
        'string', label='output', shortname='-o', &
        longname='--output', help='Optional output')

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

    nfiles = size(files)
    allocate(mols(nfiles), vibs(nfiles))
    do ifile = 1, size(files)
        dfile = DataFile(files(ifile))
        if (dfile%has_error()) then
            call write_err('std', 'Error found while initializing data file', &
                        dfile%get_error())
            stop 1
        end if
        mols(ifile) = dfile%get_mol_data()
        if (dfile%has_error()) then
            call write_err('std', &
                'Error found while parsing molecular data in file', &
                dfile%get_error())
                stop 1
        end if
        vibs(ifile) = dfile%get_vib_data()
        if (dfile%has_error()) then
            call write_err('std', &
                'Error found while parsing vibrational data in file', &
                dfile%get_error())
                stop 1
        end if

        write(iu_out, '(/,"L for file: ",a)') trim(files(ifile))
        call prt_mat(vibs(ifile)%L_mat, 3*mols(ifile)%n_at, &
                     vibs(ifile)%n_vib, iunit=iu_out)

        if (nfiles == 2) then
            if (vibs(1)%n_vib == vibs(2)%n_vib) then
                write(iu_out, '(/,"Overlap")')
                call prt_mat(matmul(transpose(vibs(1)%L_mat), vibs(2)%L_mat), &
                            vibs(1)%n_vib, vibs(1)%n_vib, iunit=iu_out)
            end if
        end if
    end do

end program test_read_vibdata
