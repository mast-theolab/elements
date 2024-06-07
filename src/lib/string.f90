module string
    !! Module with string operations
    !!
    !! The module contains string operations:
    !! - Conversion to lower case
    !! - Conversion to upper case
    !! - find string in array
    implicit none

    character(len=1), dimension(3), parameter :: labXYZ_1D = ['X', 'Y', 'Z']
    character(len=2), dimension(6), parameter :: &
        labXYZ_LT = ['XX', 'XY', 'YY', 'XZ', 'YZ', 'ZZ']
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

function str_equal(string1, string2, ignore_case, ignore_blanks) result(res)
    !! Check if string1 and string2 are equal
    !!
    !! Checks if two strings are equal.
    !! The test can be case sensitive or insensitive (default).
    !! Trailing and leading blanks can be ignored (default).

    character(len=*), intent(in) :: string1
    !! First string
    character(len=*), intent(in) :: string2
    !! Second string
    logical, intent(in), optional :: ignore_case
    !! Ignore case in comparison
    logical, intent(in), optional :: ignore_blanks
    !! Ignore trailing and leading blanks characters
    logical :: res
    !! Result of the check

    logical :: nocase, noblanks
    character(len=:), allocatable :: str1, str2

    if (present(ignore_case)) then
        nocase = ignore_case
    else
        nocase = .true.
    end if

    if (present(ignore_blanks)) then
        noblanks = ignore_blanks
    else
        noblanks = .true.
    end if

    if (nocase .and. noblanks) then
        str1 = locase(trim(string1))
        str2 = locase(trim(string2))
    else if (nocase) then
        str1 = locase(string1)
        str2 = locase(string2)
    else if (noblanks) then
        str1 = trim(string1)
        str2 = trim(string2)
    else
        str1 = string1
        str2 = string2
    end if

    res = len(str1) == len(str2) .and. str1 == str2
            
end function str_equal

! ======================================================================

end module string
