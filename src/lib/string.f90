module string
    !! Module with string operations
    !!
    !! The module contains string operations:
    !! - Conversion to lower case
    !! - Conversion to upper case
contains

! ======================================================================

function upcase(string)
    !! Converts a string to uppercase
    !!
    !! Given a string, converts to uppercase.  The chosen way makes
    !!   easier the extension to non-ASCII characters if needed.
    character(len=*), intent(in) :: string
    !! String to convert to uppercase
    character(len=len(string)) :: upcase

    integer :: i, j
    character(len=*), parameter :: &
        lcase = 'abcdefghijklmnopqrstuvwxyz', &
        ucase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len=1) :: chr
    
    do i = 1, len(string)
        chr = string(i:i)
        j = index(lcase, chr)
        if (j > 0) then
            upcase(i:i) = ucase(j:j)
        else
            upcase(i:i) = chr
        end if
    end do
    return
end function upcase

! ======================================================================

function locase(string)
    !! Converts a string to lowercase
    !!
    !! Given a string, converts to lowercase.  The chosen way makes
    !!   easier the extension to non-ASCII characters if needed.
    character(len=*), intent(in) :: string
    !! String to convert to lowercase
    character(len=len(string)) :: locase

    integer :: i, j
    character(len=*), parameter :: &
        lcase = 'abcdefghijklmnopqrstuvwxyz', &
        ucase = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len=1) :: chr
    
    do i = 1, len(string)
        chr = string(i:i)
        j = index(ucase, chr)
        if (j > 0) then
            locase(i:i) = lcase(j:j)
        else
            locase(i:i) = chr
        end if
    end do
    return
end function locase

! ======================================================================

end module string
