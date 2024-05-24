program test_read_file
    use input, only: DataFile

    character(len=256) :: infile, progname
    type(DataFile) :: dfile

    call get_command_argument(0, progname)

    if (command_argument_count() == 0) then
        print '("Usage: ",a," filename")', trim(progname)
        stop
    end if

    call get_command_argument(1, infile)
    dfile = DataFile(infile)

end program test_read_file