module input
    !! Input-processing module
    !!
    !! Contains data and procedures to process input data/options.
    use iso_fortran_env, only: real64
    use string, only: locase
    use parsefchk, only: fchkdata, fchkparser
    use basisset, only: build_bset_DB
    use atomic, only: atdata
    use moldata
    use transdata
    use exception, only: BaseException, Error, InitError, RaiseFileError, &
        RaiseError

    public :: build_moldata, build_transdata
    private :: build_moldata_fchk, build_transdata_fchk

contains

! ======================================================================

subroutine build_moldata(fname, err, ftype)
    !! Main subroutine to extract molecule-related data
    !!
    !! This subroutine should be public and the wrapper to specialized
    !!   subroutines, which can remain hidden.
    !! Note: the subroutine only uses the section after the last "."
    !!   as extension.
    character(len=*), intent(in) :: fname
    !! File name containing data of interest
    character(len=*), optional, intent(in) :: ftype
    !! File type, superseeds the automatic search
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance

    integer :: pos
    character(len=4) :: ft

    err = InitError()

    if (present(ftype)) then
        select case (locase(ftype))
            case ('fchk')
                ft = 'fchk'
            case default
                call RaiseError(err, 'Unsupported type of file')
                return
        end select
    else
        pos = index(fname, '.', back=.True.)
        select case (locase(fname(pos+1:)))
            case ('fchk', 'fck', 'fch')
                ft = 'fchk'
            case default
                call RaiseError(err, 'ERROR: Unrecognized file extension')
                return
        end select
    end if

    select case (ft)
        case ('fchk')
            call build_moldata_fchk(fname, err)
    end select

end subroutine build_moldata

! ======================================================================

subroutine build_transdata(fname, n_ab, n_basis, err, ftype)
    !! Main subroutine to extract electronic transition data.
    !!
    !! This subroutine should be public and the wrapper to specialized
    !!   subroutines, which can remain hidden.
    !! Note: the subroutine only uses the section after the last "."
    !!   as extension.
    character(len=*), intent(in) :: fname
    !! File name containing data of interest
    integer, intent(in) :: n_ab
    !! Number of unique MO sets (1 for closed-shell, 2 for open-shell)
    integer, intent(in) :: n_basis
    !! Number of basis functions
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance
    character(len=*), optional, intent(in) :: ftype
    !! File type, superseeds the automatic search

    integer :: pos
    character(len=4) :: ft

    err = InitError()

    if (present(ftype)) then
        select case (locase(ftype))
            case ('fchk')
                ft = 'fchk'
            case default
                call RaiseError(err, 'Unsupported type of file')
                return
        end select
    else
        pos = index(fname, '.', back=.True.)
        select case (locase(fname(pos+1:)))
            case ('fchk', 'fck', 'fch')
                ft = 'fchk'
            case default
                call RaiseError(err, 'ERROR: Unrecognized file extension')
                return
        end select
    end if

    select case (ft)
        case ('fchk')
            call build_transdata_fchk(fname, n_ab, n_basis, err)
    end select

end subroutine build_transdata

! ======================================================================

subroutine build_moldata_fchk(fname, err)
    !! Build molecular data from Gaussian formatted checkpoint file
    !!
    !! Parses a Gaussian checkpoint file and extracts all relevant
    !!   information
    implicit none
    character(len=*), intent(in) :: fname
    !! Name of the formatted checkpoint file
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance

    character(len=42), dimension(*), parameter :: fchk_keys = [ &
        'Number of atoms                      ', &  !  1.
        'Charge                               ', &  !  2.
        'Multiplicity                         ', &  !  3.
        'Number of electrons                  ', &  !  4.
        'Number of alpha electrons            ', &  !  5.
        'Number of beta electrons             ', &  !  6.
        'Number of basis functions            ', &  !  7.
        'Number of independent functions      ', &  !  8.
        'Number of contracted shells          ', &  !  9.
        'Current cartesian coordinates        ', &  ! 10.
        'Nuclear charges                      ', &  ! 11.
        'Real atomic weights                  ', &  ! 12.
        'Alpha MO coefficients                ', &  ! 13.
        'Beta MO coefficients                 ', &  ! 14.
        'Alpha Orbital Energies               ', &  ! 15.
        'Beta Orbital Energies                ', &  ! 16.
        'Shell types                          ', &  ! 17.
        'Number of primitives per shell       ', &  ! 18.
        'Shell to atom map                    ', &  ! 19.
        'Primitive exponents                  ', &  ! 20.
        'Contraction coefficients             ', &  ! 21.
        'P(S=P) Contraction coefficients      ', &  ! 22.
        'Pure/Cartesian d shells              ', &  ! 23.
        'Pure/Cartesian f shells              '  &  ! 24.
    ]

    integer :: i, ia, ij, j
    integer, dimension(:), allocatable :: &
        shell_types, &   ! shell types
        prim_per_sh, &   ! number of primitives per shell
        shell_to_at      ! shell to atom mapping
    real(real64), dimension(:), allocatable :: &
        coef_contr, &    ! contraction coefficients
        coef_contrSP, &  ! contraction coefficients (P(S=P))
        prim_exp         ! primitives exponents
    logical :: ok
    character(len=512) :: errmsg
    type(fchkparser) :: datafile
    type(fchkdata), dimension(:), allocatable :: dbase
    class(BaseException), allocatable :: suberr
    
    err = InitError()

    datafile = fchkparser(fname)
    dbase = datafile%read(fchk_keys)

    ok = datafile%close()
    if (.not. ok) then
        call RaiseFileError(err, fname, 'closing')
        return
    end if

    ! Molecule information
    n_at = dbase(1)%idata(1)
    allocate(at_chg(n_at), at_lab(n_at))
    do ia = 1, n_at
        at_chg(ia) = dbase(11)%rdata(ia)
        at_lab(ia) = atdata(int(at_chg(ia)))%symbol
    end do
    at_mas = dbase(12)%rdata
    at_crd = reshape(dbase(10)%rdata, [3, n_at])
    deallocate(dbase(10)%rdata, dbase(11)%rdata, dbase(12)%rdata)
    charge = dbase(2)%idata(1)
    multip = dbase(3)%idata(1)

    ! Dimension
    ! -- Electrons and open shell
    n_el = dbase(4)%idata(1)
    n_els(1) = dbase(5)%idata(1)
    n_els(2) = dbase(6)%idata(1)

    ! -- Basis function, shells atomic orbitals
    n_basis = dbase(7)%idata(1)
    n_basok = dbase(8)%idata(1)
    n_ao = n_basis
    n_shells = dbase(9)%idata(1)

    ! Molecular orbitals
    n_mos(1) = dbase(15)%len
    if (dbase(16)%dtype /= '0') then
        n_ab = 2
        openshell = .True.
        n_mos(2) = dbase(16)%len
    else
        n_ab = 1
        openshell = .False.
    end if
    n_mo = maxval(n_mos(:n_ab))
    ! -- Save orbital energies and clear redundant memory
    allocate(en_mos(n_mo,n_ab))
    do i = 1, n_mos(1)
        en_mos(i,1) = dbase(15)%rdata(i)
    end do
    deallocate(dbase(15)%rdata)
    if (openshell) then
        do i = 1, n_mos(2)
            en_mos(i,2) = dbase(16)%rdata(i)
        end do
        deallocate(dbase(16)%rdata)
    end if
    ! -- Save orbital coefficients and clear redundant memory
    allocate(coef_mos(n_mo,n_ao,n_ab))
    ij = 0
    do i = 1, n_mos(1)
        do j = 1, n_ao
            ij = ij + 1
            coef_mos(i,j,1) = dbase(13)%rdata(ij)
        end do
    end do
    deallocate(dbase(13)%rdata)
    if (openshell) then
        ij = 0
        do i = 1, n_mos(2)
            do j = 1, n_ao
                ij = ij + 1
                coef_mos(i,j,2) = dbase(14)%rdata(ij)
            end do
        end do
        deallocate(dbase(14)%rdata)
    end if
    
    ! Basis function and shells
    pureD = dbase(23)%idata(1) == 0
    pureF = dbase(24)%idata(1) == 0
    shell_types = dbase(17)%idata
    prim_per_sh = dbase(18)%idata
    shell_to_at = dbase(19)%idata
    deallocate(dbase(17)%idata, dbase(18)%idata, dbase(19)%idata)
    prim_exp = dbase(20)%rdata
    deallocate(dbase(20)%rdata)
    coef_contr = dbase(21)%rdata
    deallocate(dbase(21)%rdata)
    if (dbase(22)%dtype /= '0') then
        coef_contrSP = dbase(22)%rdata
        deallocate(dbase(22)%rdata)
    end if

    ! Build DB of basis functions
    call build_bset_DB(n_at, n_shells, pureD, pureF, shell_types, &
                       prim_per_sh, shell_to_at, coef_contr, coef_contrSP, &
                       prim_exp, nprim_per_at, bset_info, suberr)
    if (suberr%raised()) then
        select type(suberr)
            class is (Error)
                write(errmsg, '(a,a,"Original error:", a)') &
                    'Error when parsing the basis set', new_line(' '), &
                    trim(suberr%msg())
                call RaiseError(err, errmsg)
                return
            class default
                call RaiseError(err, 'Generic error')
                return
        end select
    end if

end subroutine build_moldata_fchk

! ======================================================================

subroutine build_transdata_fchk(fname, n_ab, n_basis, err)
    !! Extract data on electronic transition from Gaussian fchk
    !!
    !! Parses a Gaussian checkpoint file and extracts all relevant
    !!   information connected to electronic transitions.
    implicit none
    character(len=*), intent(in) :: fname
    !! Name of the formatted checkpoint file
    integer, intent(in) :: n_ab
    !! Number of unique MO sets (1 for closed-shell, 2 for open-shell)
    integer, intent(in) :: n_basis
    !! Number of basis functions
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance

    character(len=42), dimension(*), parameter :: fchk_keys = [ &
        'Total CI Density                     ', &  !  1.
        'ETran scalars                        ', &  !  2.
        'ETran spin                           ', &  !  3.
        'ETran sym                            ', &  !  4.
        'ETran state values                   ', &  !  5.
        'Excited state NLR                    ', &  !  6.
        'Number of g2e trans dens             ', &  !  7.
        'G to E trans densities               ', &  !  8.
        'Excited state densities              ', &  !  9.
        'SCF Energy                           '  &  ! 10.
    ]

    integer :: &
        i, &
        ioff, &
        lblock_ETran, &  ! Length of a data block for a given state
        NLR              ! indicate if left and right trans. matrix data stored
    real(real64) :: rval
    logical :: ok
    type(fchkparser) :: datafile
    type(fchkdata), dimension(:), allocatable :: dbase
    
    err = InitError()

    datafile = fchkparser(fname)
    dbase = datafile%read(fchk_keys)

    ok = datafile%close()
    if (.not. ok) then
        call RaiseFileError(err, fname, 'closing')
        return
    end if

    ! First check that NLR == 1. We do not support the other case for now
    NLR = dbase(6)%idata(1)
    if (NLR > 1) then
        call RaiseError(err, 'NLR /= 1 not yet supported.  Sorry.')
        return
    end if

    ! Now let us process the information of interest in the file.
    ! The "scalars" array provides information to parse "ETran state values"
    n_states = dbase(2)%idata(1)
    id_state = dbase(2)%idata(5)
    lblock_ETran = dbase(2)%idata(2)

    ! Extract transition data values
    allocate(g2e_energy(n_states), g2e_eldip(3,n_states), g2e_magdip(3,n_states))
    rval = dbase(10)%rdata(1)
    do i = 1, n_states
        ioff = (i-1)*lblock_ETran
        g2e_energy(i) = dbase(5)%rdata(1+ioff) - rval
        g2e_eldip(:,i) = dbase(5)%rdata(2+ioff:4+ioff)
        g2e_magdip(:,i) = dbase(5)%rdata(5+ioff:7+ioff)
    end do

    ! Now extract the ground to excited transition moments
    ! structure: n_basis, n_basis, alpha/beta, n_states
    allocate(g2e_dens(n_basis,n_basis,2,n_states), &
             exc_dens(n_basis,n_basis,n_ab,n_states))
    g2e_dens = reshape(dbase(8)%rdata, [n_basis,n_basis,2,n_states])
    exc_dens = reshape(dbase(9)%rdata, [n_basis,n_basis,n_ab,n_states])

end subroutine build_transdata_fchk

! ======================================================================

end module input
