program autoclave
    !! Automated Computation from Logfiles of Anharmonic Vibrational Energies
    !!
    !! A simple tool to compute the anharmonic energies levels of one or
    !! more vibrational states taking data from logfiles (Gaussian).
    !! The script uses an ad hoc parser since some quantities like the
    !! anharmonic X matrix are not easily available otherwise.
    use arrays, only: ij2lin => ij2lin_lt
    use exception, only: BaseException, runstat
    use numeric, only: f0, f1, realwp
    use output, only: iu_out, sec_header
    use parse_cmdline, only: CmdArgDB
    use string, only: locase
    use vibrational_PT2, only: calc_en_vib

    implicit none

    type Params
        character(len=:), allocatable :: file_log
        character(len=:), allocatable :: file_out
    end type Params

    integer :: istat, iu_in, mode, n_mode, n_vib
    integer, dimension(:), allocatable :: nq_i, nq_n
    real(realwp), dimension(:), allocatable :: h_freq
    real(realwp), dimension(:), allocatable :: a_Xmat
    logical :: exists
    character(len=256) :: fmt, line
    character(len=10), dimension(:), allocatable :: state
    
    type(Params) :: opts
    
    ! Basic parameters for the program
    ! NCOLS_XXX: maximum number of data of interest in a given row
    integer, parameter :: NCOLS_FREQ = 3, NCOLS_H2A = 9, NCOLS_XMAT = 5
    character(len=*), parameter :: PROGNAME = 'AutoCLAVE'


    ! Start program
    call sec_header(-1, PROGNAME)
1000 format(/, &
        'This program computes the energies of vibrational states using the', &
        /, &
        'vibrational perturbation theory at the 2nd order.',//, &
        'A (Gaussian) logfile is expected to extract the harmonic energies &
        &and',/, &
        'anharmonic X matrix.',/, &
        'Note that the harmonic frequency block is used to set the &
        &problem''s',/, &
        'dimension but the block actually used is in the anharmonic &
        &calculations',/, &
        'output, as "QUADRATIC FORCE CONSTANTS".',/, &
        'The program is able to handle different numbering orders in the',/, &
        'anharmonic block')
    write(iu_out, 1000)

    call sec_header(1, 'Parsing user options')

    call parse_options(opts)

    ! Check if file exists
    inquire(file=opts%file_log, exist=exists)
    if (.not.exists) &
        call runstat%raise_error('Logfile does not exist.  Aborting.')

    call sec_header(1, 'Reading input data')

    open(newunit=iu_in, file=opts%file_log)
    read(iu_in, '(a)') line
    if (index(line, 'Entering Gaussian System') == 0) then
        call runstat%raise_error( &
            'Logfile does not seem to be a Gaussian file.  Aborting.')
    else
        call parse_glog(iu_in, n_vib, h_freq, a_Xmat)
    end if

    ! Now parse user state specifications
    call sec_header(1, 'Reading state specifications')
1001 format( &
    'States should be provide as:',/, &
    'mode1 num_quanta1 [mode2 num_quanta2 [...]]',/, &
    'One state specification per line.',/, &
    'Ex:',/, &
    '1 1',/, &
    '2 1 3 1"',/, &
    '"q" to stop')
    write(iu_out, 1001)

    allocate(nq_i(n_vib), nq_n(n_vib), state(2*n_vib))
    do
        write(*, '("> ")', advance='no')
        read '(a)', line
        if (trim(line) == 'q' .or. trim(line) == 'Q' &
            .or. len_trim(line) == 0) exit
        istat = 0
        n_mode = 0
        do while(istat == 0)
            n_mode = n_mode + 1  ! We increase by 1 to check if odd number.
            read(line, *, iostat=istat) state(:n_mode)
        end do
        ! We check if n_mode is even since we have overshot the number by 1
        ! So, the true number found is n_mode + 1.
        if (mod(n_mode, 2) == 0) then
            call runstat%raise_error('Wrong state specification')
        end if
        n_mode = (n_mode - 1)/2  ! correct number since we went overboard
        read(line, *, iostat=istat) (nq_i(mode), nq_n(mode), mode=1, n_mode)
        if (istat /= 0) then
            call runstat%raise_error( &
                'Integers expected as state specification')
        else if(any(nq_i(:n_mode) <= 0) .or. any(nq_n(:n_mode) <= 0)) then
            call runstat%raise_error( &
                'Only excited modes should be listed, with positive quanta')
        end if
1100 format('("state: "',i0,'(i0,"(",i0,") ")," Energy: ",f0.6," cm^-1")')
        write(fmt, 1100) n_mode
        write(iu_out, fmt) (nq_i(mode), nq_n(mode), mode=1, n_mode), &
            calc_en_vib(n_vib, n_mode, nq_i, nq_n, h_freq, a_Xmat)
    end do

    ! do i = 1, n_vib
    !     print '(i6,f15.6)', i, h_freq(i)
    ! end do

    ! do i = 1, n_vib
    !     do j = 1, i
    !         print '(2i6,e15.6)', i, j, a_Xmat(i*(i-1)/2+j)
    !     end do
    ! end do

contains

! ======================================================================

subroutine parse_glog(iu, nvib, freq, Xmat)
    !! Parse Gaussian logfile and extract relevant information.
    integer, intent(in) :: iu
        !! Unit identifier.
    integer, intent(out) :: nvib
        !! Number of vibrations.
    real(realwp), dimension(:), allocatable, intent(out) :: freq
        !! Harmonic wavenumbers, in cm^-1^.
    real(realwp), dimension(:), allocatable, intent(out) :: Xmat
        !! Anharmonic X matrix, in cm^-1^.

    integer :: i, ibloc, icol, ios, irow, j, n, nblocs
    integer, dimension(NCOLS_XMAT) :: idx_X_j
    integer, dimension(NCOLS_H2A,2) :: indexes_HA
    integer, dimension(:), allocatable :: id_A2H
    real(realwp) :: val
    real(realwp), dimension(NCOLS_XMAT) :: vals
    logical :: found
    character(len=60) :: blocks(max(NCOLS_FREQ, NCOLS_H2A, NCOLS_XMAT))

    ! Read up to the harmonic block
    do while (index(line, 'Frequencies -- ') == 0)
        read(iu, '(a)', iostat=ios) line
        if (ios /= 0) then
            call runstat%raise_error( &
                'end-of-file reached while looking for harmonic frequencies')
        end if
    end do

    ! We have reached harmonic block
    nvib = 0
    do while (index(line, 'Second-order Perturbative Anharmonic Analysis') &
              == 0)
        if (line(:16) == ' Frequencies -- ') then
            call parse_blocs(line, 2, NCOLS_FREQ, 'harmonic-frequency', n, &
                             blocks)
            nvib = nvib + n
        end if
        read(iu, '(a)', iostat=ios) line
        if (ios /= 0) then
            call runstat%raise_error( &
                'end-of-file reached while looking for anharmonic block')
        end if
    end do

    if (nvib == 0) then
        call runstat%raise_error('Failed to read the number of modes')
    end if

    ! allocate arrays
    allocate(id_A2H(nvib), freq(nvib), Xmat(nvib*(nvib+1)/2))

    ! Now read data from anharmonic block in this order:
    ! 1. equivalence table
    ! 2. quadratic force constants
    ! 3. anharmonic X matrix

    ! Equivalence table
    do while (index(line, 'Vibro-Rotational Analysis Based on Symmetry') == 0)
        read(iu, '(a)', iostat=ios) line
        if (ios /= 0) then
            call runstat%raise_error( &
                'end-of-file reached while looking for equivalency table')
        end if
    end do

    found = .false.
    do while (index(line, 'Coriolis Couplings') == 0)
        if (line(:4) == ' (H)') then
            found = .true.
            call parse_blocs(line, 2, NCOLS_H2A, 'harmonic indexes', n, &
                             blocks)
            call parse_id_HA_blocks(n, blocks, indexes_HA(:,1))
        else if (line(:4) == ' (A)') then
            call parse_blocs(line, 2, NCOLS_H2A, 'anharm. indexes', n, &
                             blocks)
            call parse_id_HA_blocks(n, blocks, indexes_HA(:,2))
            do icol = 1, n
                if (indexes_HA(icol,2) /= 0) &
                    id_A2H(indexes_HA(icol,2)) = indexes_HA(icol,1)
            end do
        end if
        read(iu, '(a)', iostat=ios) line
        if (ios /= 0) then
            call runstat%raise_error( &
                'end-of-file reached while parsing equivalency table')
        end if
    end do
    ! No equivalence table, the indexes are matched
    if (.not.found) then
        do i = 1, nvib
            id_A2H(i) = i
        end do
    end if

    ! Harmonic frequencies
    do while (index(line, 'QUADRATIC FORCE CONSTANTS IN NORMAL MODES') == 0)
        read(iu, '(a)', iostat=ios) line
        if (ios /= 0) then
            call runstat%raise_error( &
                'end-of-file reached while looking for harmonic freq.')
        end if
    end do

    do while (index(line, 'CUBIC FORCE CONSTANTS') == 0)
        read(line, *, iostat=ios) i, j, val, blocks(:2)
        if (ios == 0) then
            freq(id_A2H(i)) = val
        end if
        read(iu, '(a)', iostat=ios) line
        if (ios /= 0) then
            call runstat%raise_error( &
                'end-of-file reached while parsing harmonic freq.')
        end if
    end do

    ! Anharmonic X matrix
    do while (index(line, 'Total Anharmonic X Matrix') == 0)
        read(iu, '(a)', iostat=ios) line
        if (ios /= 0) then
            call runstat%raise_error( &
                'end-of-file reached while looking for anh. X matrix.')
        end if
    end do
    ! Read "underline" line
    read(iu, '(a)') line

    ! The structure is particular since lower-triangular printing.
    ! Count the number of blocks
    nblocs = ceiling(real(nvib)/NCOLS_XMAT)
    do ibloc = 1, nblocs
        read(iu, '(a)') line
        read(line, *) idx_X_j
        do irow = 1, nvib - (ibloc-1)*NCOLS_XMAT
            read(iu, '(a)') line
            n = min(irow, NCOLS_XMAT)
            read(line, *) i, (vals(i), i=1, n)
            do j = 1, n
                Xmat(ij2lin(id_A2H(i), id_A2H(idx_X_j(j)))) = vals(j)
            end do
        end do
    end do

end subroutine parse_glog

! ======================================================================

subroutine parse_options(opts_db)
    type(Params), intent(out) :: opts_db

    character(len=1024) :: argval
    class(BaseException), allocatable :: err
    type(CmdArgDB) :: parser

    ! Build option parser for commandline
    parser = CmdArgDB(progname=locase(PROGNAME))
    if (parser%has_error()) then
        err = parser%exception()
        call runstat%raise_error( &
            'Unable to initialize the commandline parser', &
            err%msg())
    end if
    call parser%add_arg_char( &
        'string', label='logfile', &
        help='Logfile with anharmonic data.')
    call parser%add_arg_char( &
        'string', shortname='-o', longname='--output', &
        help='Output filename.')

    call parser%parse_args()
    if (parser%has_error()) then
        err = parser%exception()
        call runstat%raise_error( &
            'Failure to parse commandline options', &
            err%msg())
    end if

    ! Check commandline arguments and set information
    call parser%get_value('logfile', argval)
    opts_db%file_log = trim(argval)
    if (parser%is_user_set('output')) then
        call parser%get_value('output', argval)
        opts_db%file_out = trim(argval)
        open(newunit=iu_out, file=opts_db%file_out, action='write')
    end if

end subroutine parse_options

! ======================================================================

subroutine parse_id_HA_blocks(n_data, blocks, indexes_HA)
    integer, intent(in) :: n_data
        !! Number of data blocks stored in `blocks`.
    character(len=*), dimension(:), intent(in) :: blocks
        !! Parsed blocks of indexes to process.
    integer, dimension(:), intent(out) :: indexes_HA
        !! Slide of indexes to fill.

    integer :: i, ivert

    do i = 1, n_data
        ivert = index(blocks(i), '|')
        read(blocks(i)(:ivert-1), *) indexes_HA(i)
    end do
    
end subroutine parse_id_HA_blocks

! ======================================================================

subroutine parse_blocs(cstring, n_offset, max_data, label, n_data, blocks)
    character(len=*), intent(in) :: cstring
        !! String to parse.
    integer, intent(in) :: n_offset
        !! Number of elements at the beginning to ignore.
    integer, intent(in) :: max_data
        !! Maximum number of data in a given row.
    character(len=*) :: label
        !! Label of the quantity to print in error message.
    integer, intent(out) :: n_data
        !! Number of elements actually extracted.
    character(len=*), dimension(:), intent(out) :: blocks
        !! Array to store extracted data blocks, as strings.

    integer :: i, ios, n
    character(len=1024) :: msg1, msg2

    n = n_offset + max_data
    do
        read(cstring, *, iostat=ios) blocks(:n)
        if (ios /= 0) then
            n = n - 1
            if (n == n_offset) then
                write(msg2, '("Failed to parse line: ",a)') trim(cstring)
                write(msg1, '("Failed to parse ",a," block")') trim(label)
                call runstat%raise_error(trim(msg1), details=trim(msg2))
            end if
        else
            exit
        end if
    end do
    n_data = n - n_offset
    ! Remove initial elements to be ignored.
    do i = 1, n
        blocks(i) = blocks(i+n_offset)
    end do

end subroutine parse_blocs

! ======================================================================

end program autoclave