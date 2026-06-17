program test_read_vibdata
    use iso_fortran_env, only: output_unit
    use datatypes, only: MoleculeDB, VibrationsDB
    use input, only: DataFile
    use output, only: prt_mat, write_err
    use parse_cmdline, only: CmdLineArgsDB
    use run_env, only: run

    integer :: ifile, iu_out, nfiles
    integer, parameter :: MAXFILES = 2
    character(len=:), dimension(:), allocatable :: files
    character(len=256) :: outfile, string
    type(DataFile) :: dfile
    type(MoleculeDB), dimension(:), allocatable :: mols
    type(VibrationsDB), dimension(:), allocatable :: vibs
    class(CmdLineArgsDB), allocatable :: opts

    opts = CmdLineArgsDB(progname='test_read_vib')
    call run%check(opts%error, 'Failed to initialize the command-line parser')
    call opts%add_arg_char( &
        'list', label='file', shortname='-f', longname='--file', &
        required=.true., max_nvals=MAXFILES, help='Filename')
    call opts%add_arg_char( &
        'string', label='output', shortname='-o', &
        longname='--output', help='Optional output')

    call opts%parse_args()
    call run%check(opts%error, 'Unable to parse options for read_vibdat')
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
        call run%check(dfile%error, 'Error found while initializing data file')
        mols(ifile) = dfile%get_mol_data()
        call run%check(dfile%error, &
            'Error found while parsing molecular data in file')
        vibs(ifile) = dfile%get_vib_data()
        call run%check(dfile%error, &
            'Error found while parsing vibrational data in file')

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
