module parsefchk
    use iso_fortran_env, only: int64, real64
    use output, only: iu_out
    private

    integer, parameter :: LHEAD = 42, NCOLS_R = 5, NCOLS_I = 6
    type, public :: fchkdata
        character(len=:), allocatable :: key
        character(len=1) :: dtype
        integer :: len
        integer(int64), dimension(:), allocatable :: idata
        real(real64), dimension(:), allocatable :: rdata
        complex(real64), dimension(:), allocatable :: zdata
        logical, dimension(:), allocatable :: ldata
        character(len=:), dimension(:), allocatable :: cdata
    end type fchkdata
    
    ! TODO: We could add int_kind, real_kind for more versatile parser
    ! TODO: Add possibility to store keys for faster search
    type, public :: fchkparser
        private
        character(len=:), allocatable :: name
        integer :: unit
    contains
        procedure :: filename => get_fchkname
        procedure, private :: readitem, readitems
        generic :: read => readitem, readitems
    end type fchkparser

    ! The interface must have the same name as the class to allow
    !   creating the appearance of a constructor
    interface fchkparser
        module procedure init_fchk
    end interface fchkparser
contains

! ======================================================================

function init_fchk(fname, break) result(fchk)
    ! Pseudo-constructor for fchkparser
    implicit none
    type(fchkparser) :: fchk
    character(len=*), intent(in) :: fname
    logical, intent(in), optional :: break
    integer :: ios
    logical :: dobreak, exists

    if (.not.present(break)) then
        dobreak = .True.
    else
        dobreak = break
    end if

    inquire(file=fname, exist=exists)
    if (.not.exists) then
        write(iu_out, '("Error: File ",a," does not exist")') trim(fname)
        if (dobreak) stop
    else
        allocate(character(len=len_trim(fname)) :: fchk%name)
        fchk%name = trim(fname)
        open(newunit=fchk%unit, file=fchk%name, iostat=ios)
        if (ios /= 0) then
            deallocate(fchk%name)
            write(iu_out, '("Error: Could not open file ",a)') fchk%name
            if (dobreak) stop
        end if
    end if
end function init_fchk

! ======================================================================

function get_fchkname(this) result(fname)
    ! Method to retrieve the filename from the class
    character(len=:), pointer :: fname
    class(fchkparser), intent(in), target :: this
    fname => this%name
end function get_fchkname

! ======================================================================

function readitem(this, key) result(res)
    ! Read 1 item in the formatted checkpoint file
    type(fchkdata) :: res
    class(fchkparser) :: this
    character(len=*) :: key

    integer :: ios
    logical :: status
    character(len=256) :: line

    inquire(this%unit, opened=status)
    if (.not.status) then
        open(newunit=this%unit, file=this%name)
    else
        rewind(this%unit)
    end if

    do
        read(this%unit, '(a)', iostat=ios) line
        if (ios < 0) then
            res%dtype = '0'
            return
        end if
        if (line(1:LHEAD) == key) then
            res = readdata(line, this%unit, trim(key))
            exit
        end if
    end do

end function readitem

! ======================================================================

function readitems(this, keys) result(res)
    ! Read multiple items in the formatted checkpoint file
    type(fchkdata), dimension(:), allocatable :: res
    class(fchkparser) :: this
    character(len=*), dimension(:) :: keys

    integer :: idx, ios, nkeys
    integer, dimension(:), allocatable :: indexes ! store indexes still to find
    logical :: status
    character(len=256) :: line

    inquire(this%unit, opened=status)
    if (.not.status) then
        open(newunit=this%unit, file=this%name)
    else
        rewind(this%unit)
    end if

    nkeys = size(keys)
    allocate(res(nkeys), indexes(nkeys))
    indexes = [(i, i = 1, nkeys)]

    do
        read(this%unit, '(a)', iostat=ios) line
        if (ios < 0) exit
        do i = 1, nkeys
            if (line(1:LHEAD) == keys(indexes(i))) then
                idx = indexes(i)
                res(idx) = readdata(line, this%unit, trim(keys(idx)))
                indexes(i:) = eoshift(indexes(i:), 1)
                nkeys = nkeys - 1
                if (nkeys == 0) return
                exit
            end if
        end do
    end do

    if (nkeys > 0) then
        do i = 1, nkeys
            idx = indexes(i)
            res(idx)%dtype = '0'
            allocate(character(len=len_trim(keys(idx))) :: res(idx)%key)
            res(idx)%key = trim(keys(idx))
        end do
    end if

end function readitems

! ======================================================================

function readdata(line, unit, key) result(res)
    ! Read 1 dataset from an opened formatted checkpoint file
    type(fchkdata) :: res
    integer, intent(in) :: unit
    character(len=*), intent(in) :: key, line

    character(len=64) :: string

    allocate(character(len=len(key)) :: res%key)
    res%key = key

    if (index(line(LHEAD+1:), 'N=') > 0) then
        read(line(LHEAD+1:), *) res%dtype, string, res%len
        select case (res%dtype)
            case ('R')
                allocate(res%rdata(res%len))
                read(unit, *) (res%rdata(i), i=1, res%len)
            case ('I')
                allocate(res%idata(res%len))
                read(unit, *) (res%idata(i), i=1, res%len)
            case default
                write(iu_out, '(a)') 'DEVERR: Unsupported data type in FChk'
                stop
        end select
    else
        read(line(LHEAD+1:), *) res%dtype, string
        res%len = 1
        select case (res%dtype)
            case ('R')
                allocate(res%rdata(res%len))
                read(string, *) res%rdata(1)
            case ('I')
                allocate(res%idata(res%len))
                read(string, *) res%idata(1)
            case default
                write(iu_out, '(a)')  'DEVERR: Unsupported data type in FChk'
                stop
        end select
    end if

    return
end function readdata

! ======================================================================

end module parsefchk
