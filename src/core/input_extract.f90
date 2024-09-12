submodule (input) input_extract
    !! Submodule containing the definition of procedures related to the
    !! data extraction.
    use arrays, only: symm_tri_array
    use parsefchk, only: fchkdata, fchkparser
    use basisset, only: build_bset_DB
    use atominfo, only: atdata
    use moldata
    use transdata
    use exception, only: BaseException, Error, InitError, RaiseArgError, &
        RaiseError, RaiseFileError, RaiseQuantityError

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
        call build_moldata_fchk(file, load_spec_data, load_bset_data, &
                                load_morb_data)
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


subroutine build_moldata_fchk(dfile, load_spec_data, load_bset_data, &
                              load_morb_data)
    !! Builds molecular data from Gaussian formatted checkpoint file.
    !!
    !! Parses a Gaussian formatted checkpoint file and extracts basic
    !! electronic structure and molecular specification data.  Data are
    !! stored in the moldata module.
    class(DataFile), intent(inout) :: dfile
    !! Name of the formatted checkpoint file
    logical, intent(in), optional :: load_spec_data
    !! Load basic molecular specifications data from file.
    logical, intent(in), optional :: load_bset_data
    !! Load basis set data from file.
    logical, intent(in), optional :: load_morb_data
    !! Load data related to molecular orbitals from file.

    character(len=42), dimension(:), allocatable :: fchk_keys

    integer :: i, ia, ij, ind_bset, ind_morb, ind_spec, j, n_fields
    integer, dimension(:), allocatable :: &
        shell_types, &   ! shell types
        prim_per_sh, &   ! number of primitives per shell
        shell_to_at      ! shell to atom mapping
    real(real64), dimension(:), allocatable :: &
        coef_contr, &    ! contraction coefficients
        coef_contrSP, &  ! contraction coefficients (P(S=P))
        prim_exp         ! primitives exponents
    logical :: ok, load_bset, load_spec, load_morb
    character(len=512) :: errmsg
    type(fchkparser) :: dfchk
    type(fchkdata), dimension(:), allocatable :: dbase
    class(BaseException), allocatable :: suberr

    dfile%error = InitError()

    ! check data to extract and set default
    ! By default, we assume all data should be loaded.
    if (present(load_spec_data)) then
        load_spec = load_spec_data
    else
        load_spec = .true.
    end if
    if (present(load_bset_data)) then
        load_bset = load_bset_data
    else
        load_bset = .true.
    end if
    if (present(load_morb_data)) then
        load_morb = load_morb_data
    else
        load_morb = .true.
    end if

    n_fields = 4  ! number of atoms and electrons always loaded

    ! Set indexes to blocks of data fields and compute total number of fields
    if (load_spec) then
        ind_spec = n_fields + 1
        n_fields = n_fields + 6
    else
        ind_spec = 0
    end if
    if (load_bset) then
        ind_bset = n_fields + 1
        n_fields = n_fields + 11
    else
        ind_bset = 0
    end if
    if (load_morb) then
        ind_morb = n_fields + 1
        n_fields = n_fields + 4
    else
        ind_morb = 0
    end if
    allocate(fchk_keys(n_fields))
    fchk_keys(1) = 'Number of atoms'
    fchk_keys(2) = 'Number of electrons'
    fchk_keys(3) = 'Number of alpha electrons'
    fchk_keys(4) = 'Number of beta electrons'
    if (load_spec) then
        fchk_keys(ind_spec  ) = 'Charge'
        fchk_keys(ind_spec+1) = 'Multiplicity'
        fchk_keys(ind_spec+2) = 'Current cartesian coordinates'
        fchk_keys(ind_spec+3) = 'Atomic numbers'
        fchk_keys(ind_spec+4) = 'Nuclear charges'
        fchk_keys(ind_spec+5) = 'Real atomic weights'
    end if
    if (load_bset) then
        fchk_keys(ind_bset   ) = 'Number of basis functions'
        fchk_keys(ind_bset+ 1) = 'Number of independent functions'
        fchk_keys(ind_bset+ 2) = 'Number of contracted shells'
        fchk_keys(ind_bset+ 3) = 'Shell types'
        fchk_keys(ind_bset+ 4) = 'Number of primitives per shell'
        fchk_keys(ind_bset+ 5) = 'Shell to atom map'
        fchk_keys(ind_bset+ 6) = 'Primitive exponents'
        fchk_keys(ind_bset+ 7) = 'Contraction coefficients'
        fchk_keys(ind_bset+ 8) = 'P(S=P) Contraction coefficients'
        fchk_keys(ind_bset+ 9) = 'Pure/Cartesian d shells'
        fchk_keys(ind_bset+10) = 'Pure/Cartesian f shells'
    end if
    if (load_morb) then
        fchk_keys(ind_morb) = 'Alpha MO coefficients'
        fchk_keys(ind_morb+1) = 'Beta MO coefficients'
        fchk_keys(ind_morb+2) = 'Alpha Orbital Energies'
        fchk_keys(ind_morb+3) = 'Beta Orbital Energies'
    end if

    dfchk = fchkparser(dfile%name)
    dbase = dfchk%read(fchk_keys)

    ok = dfchk%close()
    if (.not. ok) then
        call RaiseFileError(dfile%error, dfile%name, 'closing')
        return
    end if

    ! Basic molecular information
    ! Start with standard dimensions that may be used for different applications
    if (.not.dat_base_loaded) then
        n_at = dbase(1)%idata(1)
        n_el = dbase(2)%idata(1)
        n_els(1) = dbase(3)%idata(1)
        n_els(2) = dbase(4)%idata(1)
        dat_base_loaded = .true.
    else
        ! check if data are consistent
        if (n_at /= dbase(1)%idata(1) .or. n_el /= dbase(2)%idata(1)) then
            call RaiseQuantityError(dfile%error, &
                msg='Inconsistency found between molecular data sets.')
            return
        end if
    end if

    ! Basic molecular specifications
    if (load_spec) then
        charge = dbase(ind_spec)%idata(1)
        multip = dbase(ind_spec+1)%idata(1)
        allocate(at_chg(n_at), at_lab(n_at), at_num(n_at))
        do ia = 1, n_at
            at_num(ia) = dbase(ind_spec+3)%idata(ia)
            at_chg(ia) = dbase(ind_spec+4)%rdata(ia)
            at_lab(ia) = atdata(at_num(ia))%symbol
        end do
        at_mas = dbase(ind_spec+5)%rdata
        at_crd = reshape(dbase(ind_spec+2)%rdata, [3, n_at])
        deallocate(dbase(ind_spec+2)%rdata, dbase(ind_spec+3)%idata, &
                   dbase(ind_spec+4)%rdata, dbase(ind_spec+5)%rdata)
        dat_spec_loaded = .true.
    end if

    ! Basis functions
    if (load_bset) then
        ! -- Basis function, shells atomic orbitals
        n_basis = dbase(ind_bset)%idata(1)
        n_basok = dbase(ind_bset+1)%idata(1)
        n_ao = n_basis
        n_shells = dbase(ind_bset+2)%idata(1)

        ! Basis function and shells
        pureD = dbase(ind_bset+9)%idata(1) == 0
        pureF = dbase(ind_bset+10)%idata(1) == 0
        shell_types = dbase(ind_bset+3)%idata
        prim_per_sh = dbase(ind_bset+4)%idata
        shell_to_at = dbase(ind_bset+5)%idata
        deallocate(dbase(ind_bset+3)%idata, dbase(ind_bset+4)%idata, &
                   dbase(ind_bset+5)%idata)
        prim_exp = dbase(ind_bset+6)%rdata
        deallocate(dbase(ind_bset+6)%rdata)
        coef_contr = dbase(ind_bset+7)%rdata
        deallocate(dbase(ind_bset+7)%rdata)
        if (dbase(ind_bset+8)%dtype /= '0') then
            coef_contrSP = dbase(ind_bset+8)%rdata
            deallocate(dbase(ind_bset+8)%rdata)
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
        dat_bset_loaded = .true.
    end if

    ! Molecular orbitals
    if (load_morb) then
        n_mos(1) = dbase(ind_morb+2)%len
        if (dbase(ind_morb+3)%dtype /= '0') then
            n_ab = 2
            openshell = .True.
            n_mos(2) = dbase(ind_morb+3)%len
        else
            n_ab = 1
            openshell = .False.
        end if
        n_mo = maxval(n_mos(:n_ab))
        ! -- Save orbital energies and clear redundant memory
        allocate(en_mos(n_mo,n_ab))
        en_mos(:,1) = dbase(ind_morb+2)%rdata
        deallocate(dbase(ind_morb+2)%rdata)
        if (openshell) then
            en_mos(:,2) = dbase(ind_morb+3)%rdata
            deallocate(dbase(ind_morb+3)%rdata)
        end if
        ! -- Save orbital coefficients and clear redundant memory
        allocate(coef_mos(n_mo,n_ao,n_ab))
        ij = 0
        do i = 1, n_mos(1)
            do j = 1, n_ao
                ij = ij + 1
                coef_mos(i,j,1) = dbase(ind_morb)%rdata(ij)
            end do
        end do
        deallocate(dbase(ind_morb)%rdata)
        if (openshell) then
            ij = 0
            do i = 1, n_mos(2)
                do j = 1, n_ao
                    ij = ij + 1
                    coef_mos(i,j,2) = dbase(ind_morb+1)%rdata(ij)
                end do
            end do
            deallocate(dbase(ind_morb+1)%rdata)
        end if
        dat_morb_loaded = .true.
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
        i,  i_ab, ioff, n_ab_, nbas_lt, n_basis_, &
        lblock_ETran, &  ! Length of a data block for a given state
        NLR              ! indicate if left and right trans. matrix data stored
    real(real64) :: rval
    logical :: ok
    type(fchkparser) :: dfchk
    type(fchkdata), dimension(:), allocatable :: dbase

    dfile%error = InitError()

    dfchk = fchkparser(dfile%name)
    dbase = dfchk%read(fchk_keys)

    ok = dfchk%close()
    if (.not. ok) then
        call RaiseFileError(dfile%error, dfile%name, 'closing')
        return
    end if

    ! -- Basic information
    ! Number of basis functions
    n_basis_ = dbase(1)%idata(1)
    if (dbase(2)%dtype /= '0') then
        n_ab_ = 2
    else
        n_ab_ = 1
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
    ! structure: n_basis_, n_basis_, alpha/beta, n_states
    allocate(g2e_dens(n_basis_,n_basis_,2,n_states), &
             exc_dens(n_basis_,n_basis_,n_ab_,n_states))
    ! WARNING: We need to be careful here!
    ! Gaussian up to G16 included saved the transition densities with embedded
    ! coefficients, sqrt(2) for TD and it seems 1/sqrt(2) for CI.
    ! We need to do the inverse operation to get the correct coefficients.
    if (dfile%check_version(major='G16')) then
        g2e_dens = reshape(dbase(10)%rdata, [n_basis_,n_basis_,2,n_states]) &
        / sqrt(2.0_real64)
    else
        g2e_dens = reshape(dbase(10)%rdata, [n_basis_,n_basis_,2,n_states])
    end if

    ! Excited-state densities are stored in lower-triangular forms.
    ! We need to unpack them.
    ! They are always stored for both a and b, but in case of closed-shell.
    nbas_lt = n_basis_*(n_basis_+1)/2
    do i = 1, n_states
        ioff = (i-1)*2*nbas_lt
        do i_ab = 1, n_ab_
            call symm_tri_array(n_basis_, &
                dbase(11)%rdata(ioff+(i_ab-1)*nbas_lt+1:), linear=.True., &
                lower=.True., anti_symm=.False., arr_new=exc_dens(:,:,i_ab,i))
        end do
    end do

end subroutine build_excdata_fchk

! ======================================================================

end
