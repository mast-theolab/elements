program test_read_vibdata
    use input, only : build_moldata, build_vibdata, DataFile
    use output, only: prt_mat, write_err
    use moldata
    use vibdata

    character(len=256) :: file1, file2 = ' ', progname
    type(DataFile) :: dfile1, dfile2

    call get_command_argument(0, progname)

    if (command_argument_count() < 1) then
        print '("Usage: ",a," file1 [file2]")', trim(progname)
        stop
    end if

    call get_command_argument(1, file1)
    ! if (command_argument_count() >= 2) call get_command_argument(2, file1)

    dfile1 = DataFile(file1)
    if (dfile1%has_error()) then
        call write_err('std', 'Error found while initializing data file', &
                       dfile1%get_error())
        stop 1
    end if
    call dfile1%build_moldata()
    if (dfile1%has_error()) then
        call write_err('std', &
            'Error found while parsing molecular data in file', &
            dfile1%get_error())
            stop 1
    end if
    call dfile1%build_vibdata()
    if (dfile1%has_error()) then
        call write_err('std', &
            'Error found while parsing vibrational data in file', &
            dfile1%get_error())
            stop 1
    end if

    call prt_mat(L_mat, 3*n_at, n_vib)

end program test_read_vibdata
