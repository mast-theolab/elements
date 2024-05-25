submodule (input) input_extract
    !! Submodule containing the definition of procedures related to the
    !! data extraction.
    use arrays, only: symm_tri_array
    use parsefchk, only: fchkdata, fchkparser
    use basisset, only: build_bset_DB
    use atomic, only: atdata
    use moldata
    use transdata
    use exception, only: BaseException, Error, InitError, RaiseArgError, &
        RaiseError, RaiseFileError

    implicit none

contains

! ======================================================================
! MODULE INTERFACE
! ======================================================================

module procedure build_moldata
    !! Main subroutine to load and store molecule-related data.
    !!
    !! This subroutine is intended to be bound to DataFile instances
    !! but can be used directly, providing either a DataFile instance
    !! or the filename.  In the latter case, an instance of DataFile
    !! is created internally to be transmitted to the internal routines.
    type(DataFile), target :: finf
    class(DataFile), pointer :: file

    if (present(err)) err = InitError()

    if (.not.present(dfile)) then
        finf = DataFile(fname, ftype)
        if (finf%error%raised()) then
            if (present(err)) then
                err = finf%error
                return
            else
                print '("Failed to build file information")'
                print '("Motive: ",a)', finf%error%msg()
                stop 1
            end if
        end if
        file => finf
    else
        file => dfile
    end if

    if (file%type == 'GFChk') then
        call build_moldata_fchk(file)
        if (file%error%raised() .and. .not.present(dfile)) then
            if (present(err)) then
                err = file%error
                return
            else
                print '("Error while parsing molecular data")'
                print '("Reason: ",a)', file%error%msg()
            end if
        end if
    else
        if (present(err)) then
            err = InitError()
            call RaiseArgError(err, 'Unsupported file type')
            return
        else
            print '("Unsupported file type: ",a)', file%type
            stop 1
        end if
    end if

end procedure build_moldata

! ----------------------------------------------------------------------

module procedure build_excdata
    !! Main subroutine to load and store excited-states data.
    !!
    !! This subroutine is intended to be bound to DataFile instances
    !! but can be used directly, providing either a DataFile instance
    !! or the filename.  In the latter case, an instance of DataFile
    !! is created internally to be transmitted to the internal routines.
    type(DataFile), target :: finf
    class(DataFile), pointer :: file

    if (present(err)) err = InitError()

    if (.not.present(dfile)) then
        finf = DataFile(fname, ftype)
        if (finf%error%raised()) then
            if (present(err)) then
                err = finf%error
                return
            else
                print '("Failed to build file information")'
                print '("Motive: ",a)', finf%error%msg()
                stop 1
            end if
        end if
        file => finf
    else
        file => dfile
    end if

    if (file%type == 'GFChk') then
        call build_excdata_fchk(file)
        if (file%error%raised() .and. .not.present(dfile)) then
            if (present(err)) then
                err = file%error
                return
            else
                print '("Error while parsing excited-state data")'
                print '("Reason: ",a)', file%error%msg()
            end if
        end if
    else
        if (present(err)) then
            err = InitError()
            call RaiseArgError(err, 'Unsupported file type')
            return
        else
            print '("Unsupported file type: ",a)', file%type
            stop 1
        end if
    end if

end procedure build_excdata

! ======================================================================
! SUB-MODULE COMPONENTS
! ======================================================================


subroutine build_moldata_fchk(dfile)
    !! Builds molecular data from Gaussian formatted checkpoint file.
    !!
    !! Parses a Gaussian formatted checkpoint file and extracts basic
    !! electronic structure and molecular specification data.  Data are
    !! stored in the moldata module.
    class(DataFile), intent(inout) :: dfile
    !! Name of the formatted checkpoint file

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

    dfile%error = InitError()

    datafile = fchkparser(dfile%name)
    dbase = datafile%read(fchk_keys)

    ok = datafile%close()
    if (.not. ok) then
        call RaiseFileError(dfile%error, dfile%name, 'closing')
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
    en_mos(:,1) = dbase(15)%rdata
    deallocate(dbase(15)%rdata)
    if (openshell) then
        en_mos(:,2) = dbase(16)%rdata
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
                call RaiseError(dfile%error, errmsg)
                return
            class default
                call RaiseError(dfile%error, 'Generic error')
                return
        end select
    end if

end subroutine build_moldata_fchk

! ======================================================================

subroutine build_excdata_fchk(dfile)
    !! Builds excited-state data from Gaussian formatted chk file.
    !!
    !! Parses a Gaussian formatted checkpoint file and extracts all
    !! relevant information on electronic excited-state data and
    !! transition data.  Data are stored in the transdata module.
    class(DataFile), intent(inout) :: dfile
    !! Name of the formatted checkpoint file

    character(len=42), dimension(*), parameter :: fchk_keys = [ &
        'Number of basis functions            ', &  !  1.
        'Beta Orbital Energies                ', &  !  2.
        'Total CI Density                     ', &  !  3.
        'ETran scalars                        ', &  !  4.
        'ETran spin                           ', &  !  5.
        'ETran sym                            ', &  !  6.
        'ETran state values                   ', &  !  7.
        'Excited state NLR                    ', &  !  8.
        'Number of g2e trans dens             ', &  !  9.
        'G to E trans densities               ', &  ! 10.
        'Excited state densities              ', &  ! 11.
        'SCF Energy                           '  &  ! 12.
    ]

    integer :: &
        i,  i_ab, ioff, n_ab, nbas_lt, n_basis, &
        lblock_ETran, &  ! Length of a data block for a given state
        NLR              ! indicate if left and right trans. matrix data stored
    real(real64) :: rval
    logical :: ok
    type(fchkparser) :: datafile
    type(fchkdata), dimension(:), allocatable :: dbase

    dfile%error = InitError()

    datafile = fchkparser(dfile%name)
    dbase = datafile%read(fchk_keys)

    ok = datafile%close()
    if (.not. ok) then
        call RaiseFileError(dfile%error, dfile%name, 'closing')
        return
    end if

    ! -- Basic information
    ! Number of basis functions
    n_basis = dbase(1)%idata(1)
    if (dbase(2)%dtype /= '0') then
        n_ab = 2
    else
        n_ab = 1
    end if

    ! First check that NLR == 1. We do not support the other case for now
    NLR = dbase(8)%idata(1)
    if (NLR > 1) then
        call RaiseError(dfile%error, 'NLR /= 1 not yet supported.  Sorry.')
        return
    end if

    ! Now let us process the information of interest in the file.
    ! The "scalars" array provides information to parse "ETran state values"
    n_states = dbase(4)%idata(1)
    id_state = dbase(4)%idata(5)
    lblock_ETran = dbase(4)%idata(2)

    ! Extract excited-states data
    ispin_exc = dbase(5)%idata

    ! Extract transition data values
    allocate(g2e_energy(n_states), g2e_eldip(3,n_states), &
             g2e_magdip(3,n_states))
    rval = dbase(12)%rdata(1)
    do i = 1, n_states
        ioff = (i-1)*lblock_ETran
        g2e_energy(i) = dbase(7)%rdata(1+ioff) - rval
        g2e_eldip(:,i) = dbase(7)%rdata(2+ioff:4+ioff)
        g2e_magdip(:,i) = dbase(7)%rdata(5+ioff:7+ioff)
    end do

    ! Now extract the ground to excited transition moments
    ! structure: n_basis, n_basis, alpha/beta, n_states
    allocate(g2e_dens(n_basis,n_basis,2,n_states), &
             exc_dens(n_basis,n_basis,n_ab,n_states))
    ! WARNING: We need to be careful here!
    ! Gaussian up to G16 included saved the transition densities with embedded
    ! coefficients, sqrt(2) for TD and it seems 1/sqrt(2) for CI.
    ! We need to do the inverse operation to get the correct coefficients.
    if (dfile%check_version(major='G16')) then
        g2e_dens = reshape(dbase(10)%rdata, [n_basis,n_basis,2,n_states]) &
        / sqrt(2.0_real64)
    else
        g2e_dens = reshape(dbase(10)%rdata, [n_basis,n_basis,2,n_states])
    end if

    ! Excited-state densities are stored in lower-triangular forms.
    ! We need to unpack them.
    ! They are always stored for both a and b, but in case of closed-shell.
    nbas_lt = n_basis*(n_basis+1)/2
    do i = 1, n_states
        ioff = (i-1)*2*nbas_lt
        do i_ab = 1, n_ab
            call symm_tri_array(n_basis, &
                dbase(11)%rdata(ioff+(i_ab-1)*nbas_lt+1:), linear=.True., &
                lower=.True., anti_symm=.False., arr_new=exc_dens(:,:,i_ab,i))
        end do
    end do

end subroutine build_excdata_fchk

! ======================================================================

end