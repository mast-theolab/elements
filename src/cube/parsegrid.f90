module parsegrid
    use iso_fortran_env, only: int64, real64
    use output, only: iu_out
    use run_env, only: CoreExecObject, ErrorHandle
    private

    type, public :: griddata
        character(len=42) :: grid_type
        real(real64), dimension(3) :: grid_origin
        integer, dimension(3) :: n_points
        real(real64), dimension(3) :: step_size
    end type griddata

    type, public, extends(CoreExecObject) :: gridparser
        private
        character(len=:), allocatable :: name
        integer :: unit
    contains
        procedure :: filename => get_grid_fname
        procedure :: read => read_items
        procedure :: close => close_grid_fname
    end type gridparser

    interface gridparser
        module procedure init_gridparser
    end interface gridparser

contains

    function init_gridparser(fname) result(parser)
        implicit none
        type(gridparser) :: parser
        character(len=*), intent(in) :: fname

        integer :: ios
        character(len=1024) :: msg

        open(newunit=parser%unit, file=fname, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            write(msg, '("Could not open file ",a)') trim(fname)
            call parser%error%raise_error('file', 'missing', fname)
        end if
        parser%name = trim(fname)

    end function init_gridparser

    function get_grid_fname(this) result(fname)
        implicit none
        character(len=:), pointer :: fname
        class(gridparser), intent(in), target :: this
        fname => this%name
    end function get_grid_fname

    function close_grid_fname(this) result(status)
        implicit none
        logical :: status
        class(gridparser), intent(in), target :: this

        integer :: ios
        close(this%unit, iostat=ios)
        status = ios == 0
    end function close_grid_fname

    function read_items(this, keys) result(res)
        ! Read multiple items from a grid file
        ! Distances are in atomic units
        implicit none
        type(griddata) :: res
        class(gridparser) :: this
        character(len=*), dimension(:), intent(in) :: keys

        integer, dimension(:), allocatable :: indexes ! store indexes still to find
        character(len=256) :: line
        character(len=1024) :: msg
        integer :: i, ios, nkeys, idx, offset
        logical :: status

        inquire(this%unit, opened=status)
        if (.not.status) then
            open(newunit=this%unit, file=this%name)
        else
            rewind(this%unit)
        end if

        nkeys = size(keys)
        allocate(indexes(nkeys))

        indexes = [(i, i = 1, nkeys)]

        do
            read(this%unit, '(a)', iostat=ios) line
            if (ios < 0) exit
            do i = 1, nkeys
                if (index(trim(line), trim(keys(indexes(i)))) > 0) then
                    idx = indexes(i)
                    offset = len(trim(keys(idx))) + 2
                    line = trim(line(offset:))
                    if (keys(idx) == 'Grid type') then
                        read(line, '(a)') res%grid_type
                    else if (keys(idx) == 'Grid origin') then
                        read(line, *) res%grid_origin
                    else if (keys(idx) == 'N points') then
                        read(line, *) res%n_points
                    else if (keys(idx) == 'Step') then
                        read(line, *) res%step_size
                    end if
                    indexes(i:) = eoshift(indexes(i:), 1)
                    nkeys = nkeys - 1
                    if (nkeys == 0) return
                    exit
                end if
            end do
        end do

        if (nkeys > 0) then
            write(msg, '("Could not find all keys in input grid file ",a)') &
                trim(this%name)
            call this%error%raise_error('key', 'missing', msg)
        end if

   end function read_items

end module parsegrid
