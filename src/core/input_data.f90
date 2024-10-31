submodule (input) input_data
    !! Submodule containing the definition of procedures related to the
    !! data extraction.
    use basisset, only: build_bset_DB
    use string, only: locase
    use workdata
    use exception, only: Error, InitError, RaiseArgError

    implicit none

interface

! ----------------------------------------------------------------------

module function build_mol_data_fchk(dfile) result(mol)
    class(DataFile), intent(inout) :: dfile
    !! Name of the formatted checkpoint file.
    type(MoleculeDB) :: mol
    !! Molecular specifications database.

end function build_mol_data_fchk

! ----------------------------------------------------------------------

module function build_bset_data_fchk(dfile) result(bset)
    class(DataFile), intent(inout) :: dfile
    !! Name of the formatted checkpoint file.
    type(BasisSetDB) :: bset
    !! Basis set information database.

end function build_bset_data_fchk

! ----------------------------------------------------------------------

module function build_orb_data_fchk(dfile) result(orb)
    class(DataFile), intent(inout) :: dfile
    !! Name of the formatted checkpoint file.
    type(OrbitalsDB) :: orb
    !! Orbitals information database.

end function build_orb_data_fchk

! ----------------------------------------------------------------------

module function build_exc_data_fchk(dfile) result(exc)
    class(DataFile), intent(inout) :: dfile
    !! Name of the formatted checkpoint file.
    type(ExcitationDB) :: exc
    !! Electronic excitation information.

end function build_exc_data_fchk

! ----------------------------------------------------------------------

module function build_vib_data_fchk(dfile, get_Lmat, get_Lmweig) result(vib)
    class(DataFile), intent(inout) :: dfile
    !! Name of the formatted checkpoint file
    logical, intent(in), optional :: get_Lmat
    !! Build/load dimensionless matrix of Hessian eigenvectors.
    logical, intent(in), optional :: get_Lmweig
    !! Build/load mass-weighted matrix of Hessian eigenvectors.
    type(VibrationsDB) :: vib
    !! Vibrational information.

end function build_vib_data_fchk

! ----------------------------------------------------------------------

end interface


contains

! ======================================================================
! MODULE INTERFACE
! ======================================================================

module procedure build_mol_data
    !! Main function to load molecule-related data.
    !!
    !! This function is intended to be bound to DataFile instances
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
        mol = build_mol_data_fchk(file)
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

end procedure build_mol_data

! ----------------------------------------------------------------------

module procedure build_bset_data
    !! Main function to load basis set-related data.
    !!
    !! This function is intended to be bound to DataFile instances
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
        bset = build_bset_data_fchk(file)
        if (file%error%raised() .and. .not.present(dfile)) then
            if (present(err)) then
                err = file%error
                return
            else
                print '("Error while parsing basis set data")'
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

end procedure build_bset_data

! ----------------------------------------------------------------------

module procedure build_orb_data
    !! Main function to load orbitals-related data.
    !!
    !! This function is intended to be bound to DataFile instances
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
        orb = build_orb_data_fchk(file)
        if (file%error%raised() .and. .not.present(dfile)) then
            if (present(err)) then
                err = file%error
                return
            else
                print '("Error while parsing orbital data")'
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

end procedure build_orb_data

! ----------------------------------------------------------------------

module procedure build_exc_data
    !! Main function to load and store excited-states data.
    !!
    !! This function is intended to be bound to DataFile instances
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
        exc = build_exc_data_fchk(file)
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

end procedure build_exc_data

! ----------------------------------------------------------------------

module procedure build_vib_data
    !! Main function to load and store vibrational data.
    !!
    !! This function is intended to be bound to DataFile instances
    !! but can be used directly, providing either a DataFile instance
    !! or the filename.  In the latter case, an instance of DataFile
    !! is created internally to be transmitted to the internal routines.
    !!
    !! Note that the procedure focuses on extracting core information,
    !! not necessarily data related to vibrational transitions.
    !! The eigenvectors matrix can have multiple forms and two can be
    !! typically used:
    !!
    !! - diagonalization of the mass-weighted force constants, convenient
    !!   to convert derivatives from Cartesian to mass-weighted normal
    !!   coordinates.
    !! - diagonalization of the force constants.  The resulting
    !!   eigenvectors are dimensionless, generally more convenient for
    !!   standard operations.
    !!
    !! The routine can provide either or both, doing internally all
    !! necessary conversions.
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
        vib = build_vib_data_fchk(file, get_Lmat, get_Lmweig)
        if (file%error%raised() .and. .not.present(dfile)) then
            if (present(err)) then
                err = file%error
                return
            else
                print '("Error while parsing vibration-related data")'
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

end procedure build_vib_data

! ======================================================================
! SUB-MODULE COMPONENTS
! ======================================================================



end
