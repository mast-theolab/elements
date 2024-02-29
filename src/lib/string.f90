module string
    !! Module with string operations
    !!
    !! The module contains string operations:
    !! - Conversion to lower case
    !! - Conversion to upper case
    !! - find string in array
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

function findstr(array, string) result(pos)
    !! Find position of `string` in `array` of strings.
    !!
    !! Returns the position of `string` in `array`, 0 otherwise.
    !! Array is expected to be 1D.
    implicit none

    character(len=*), dimension(:), intent(in) :: array
    !! Array of strings to search.
    character(len=*), intent(in) :: string
    !! String to find.
    integer :: pos
    !! Position of `string` in `array`.

    integer :: i

    pos = 0
    do i = 1, size(array)
        if (array(i) == string) then
            pos = i
            exit
        end if
    end do

end function findstr

! ======================================================================

end module string
