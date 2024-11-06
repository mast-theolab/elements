program test_timestamp
    use iso_fortran_env, only: output_unit
    use string, only: timestamp
    use output, only: iu_out, sec_header

    character(len=*), parameter :: fname = 'test_timestamp.txt'
    character(len=:), allocatable :: datetime

    open(newunit=iu_out, file=fname, action='write')

    call sec_header(0, 'Test Program for TimeStamp')

    write(iu_out, '(/," Filename: ",a)') fname

    call sec_header(1, 'Date (ISO)')
    
    call sec_header(2, 'Full date')
    datetime = timestamp()
    write(iu_out, '(1x,a)') datetime
    
    call sec_header(2, 'Full date with ms')
    datetime = timestamp(add_ms=.true.)
    write(iu_out, '(1x,a)') datetime
    
    call sec_header(2, 'Only date')
    datetime = timestamp(time=.false.)
    write(iu_out, '(1x,a)') datetime
    
    call sec_header(2, 'Only time')
    datetime = timestamp(date=.false.)
    write(iu_out, '(1x,a)') datetime
    
    call sec_header(2, 'Only time with ms')
    datetime = timestamp(date=.false., add_ms=.true.)
    write(iu_out, '(1x,a)') datetime

    call sec_header(1, 'Date (US)')
    
    call sec_header(2, 'Full date')
    datetime = timestamp(US_date_order=.true.)
    write(iu_out, '(1x,a)') datetime
    
    call sec_header(2, 'Full date with ms')
    datetime = timestamp(add_ms=.true., US_date_order=.true.)
    write(iu_out, '(1x,a)') datetime
    
    call sec_header(2, 'Only date')
    datetime = timestamp(time=.false., US_date_order=.true.)
    write(iu_out, '(1x,a)') datetime
    
    call sec_header(1, 'Date (short format)')
    
    call sec_header(2, 'Full date')
    datetime = timestamp(shorten=.true.)
    write(iu_out, '(1x,a)') datetime
    
    call sec_header(2, 'Full date with ms')
    datetime = timestamp(add_ms=.true., shorten=.true.)
    write(iu_out, '(1x,a)') datetime
    
    call sec_header(2, 'Only date')
    datetime = timestamp(time=.false., shorten=.true.)
    write(iu_out, '(1x,a)') datetime
    
end program test_timestamp
