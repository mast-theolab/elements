module input
    !! Input-processing module
    !!
    !! Contains data and procedures to process input data/options.
    use datatypes
    use run_env, only: CoreExecObject, ErrorHandle, run

    implicit none

    private
    public :: build_bset_data, build_exc_data, build_mol_data, &
        build_orb_data, build_vib_data, DataFile, ProgramInfo

    type, extends(CoreExecObject) :: ProgramInfo
        private
        character(len=:), allocatable :: name
        character(len=:), allocatable :: major, minor, version
    contains
        procedure, pass(prog_info) :: get_name => get_prog_name
        procedure, pass(prog_info) :: get_major => get_prog_major
        procedure, pass(prog_info) :: get_minor => get_prog_minor
        procedure, pass(prog_info) :: get_version => get_prog_version
        procedure, pass(prog_info) :: check_version => check_prog_version
    end type ProgramInfo

    type, extends(CoreExecObject) :: DataFile
        private
        character(len=:), allocatable :: name
            !! Filename associated to data file.
        character(len=:), allocatable :: type
            !! File type associated to data file.
        type(ProgramInfo) :: prog
            !! Information on program that generated file.
    contains
        procedure :: get_mol_data => build_mol_data
        procedure :: get_bset_data => build_bset_data
        procedure :: get_orb_data => build_orb_data
        procedure :: get_exc_data => build_exc_data
        procedure :: get_vib_data => build_vib_data
        procedure, private :: get_data_from_id, get_data_from_tag
        procedure :: get_filename => get_datafile_name
        procedure :: get_filetype => get_datafile_type
        procedure, pass(file_data) :: program => get_prog_name
        procedure, pass(file_data) :: version => get_prog_version
        procedure, pass(file_data) :: check_version => check_prog_version
        generic :: get_data => get_data_from_id, get_data_from_tag
    end type DataFile

    interface DataFile
        module procedure init_data_file
    end interface DataFile

    interface ProgramInfo
        module procedure init_program_version
    end interface ProgramInfo

interface

! ----------------------------------------------------------------------
! INTERFACE TO CONSTRUCTORS
! ----------------------------------------------------------------------

module function init_data_file(fname, ftype, prog_name, prog_major, &
                               prog_minor, no_prog_version_ok, exit_on_error, &
                               silent) result(dfile)
    character(len=*), intent(in) :: fname
        !! File name.
    character(len=*), intent(in), optional :: ftype
        !! File type, which overrides any detection attempt.
    character(len=*), intent(in), optional :: prog_name
        !! Program name; overrides internal search.
    character(len=*), intent(in), optional :: prog_major
        !! Major version.  If provided, overrides the automatic search.
    character(len=*), intent(in), optional :: prog_minor
        !! Minor version.  If provided, overrides the automatic search.
    logical, intent(in), optional :: no_prog_version_ok
        !! If no version is found, simply ignore, setting version to N/A.
    logical, intent(in), optional :: exit_on_error
        !! Exit if an error is encountered.  By default, `.true.`
    logical, intent(in), optional :: silent
        !! Do not print messages.  By default, `.false.`.
    type(DataFile) :: dfile
        !! DataFile instance

end function init_data_file

! ----------------------------------------------------------------------

module function init_program_version(fname, ftype, prog_name, major, minor, &
                                     no_version_ok, exit_on_error, silent &
                                     ) result(prog)
    character(len=*), intent(in) :: fname
        !! File name.
    character(len=*), intent(in) :: ftype
        !! File type.
    character(len=*), intent(in), optional :: prog_name
        !! Program name; overrides internal search.
    character(len=*), intent(in), optional :: major
        !! Major version.  If provided, overrides the automatic search.
    character(len=*), intent(in), optional :: minor
        !! Minor version.  If provided, overrides the automatic search.
    logical, intent(in), optional :: no_version_ok
        !! If no version is found, simply ignore, setting version to N/A.
    logical, intent(in), optional :: exit_on_error
        !! Exit if an error is encountered.  By default, `.true.`
    logical, intent(in), optional :: silent
        !! Do not print messages.  By default, `.false.`.
    type(ProgramInfo) :: prog
        !! ProgramInfo instance with program information.

end function init_program_version

! ----------------------------------------------------------------------
! INTERFACE TO TYPE-BOUND PROCEDURES
! ----------------------------------------------------------------------

module function build_exc_data(dfile, get_dens, fname, ftype) result(exc)
    class(DataFile), intent(inout), target, optional :: dfile
        !! DataFile instance.
    character(len=*), intent(in), optional :: fname
        !! File name containing data of interest.
    character(len=*), intent(in), optional :: ftype
        !! File type, superseeds the automatic search.
    logical, intent(in), optional :: get_dens
        !! Load electronic transition density from data file.
    type(ExcitationDB) :: exc
        !! Electronic excitation information.

end function build_exc_data

! ----------------------------------------------------------------------

module function build_bset_data(dfile, fname, ftype) result(bset)
    class(DataFile), intent(inout), target, optional :: dfile
        !! DataFile instance.
    character(len=*), intent(in), optional :: fname
        !! File name containing data of interest.
    character(len=*), intent(in), optional :: ftype
        !! File type, superseeds the automatic search.
    type(BasisSetDB) :: bset
        !! Basis set information database.

end function build_bset_data

! ----------------------------------------------------------------------

module function build_mol_data(dfile, fname, ftype, get_dens) result(mol)
    class(DataFile), intent(inout), target, optional :: dfile
        !! DataFile instance.
    character(len=*), intent(in), optional :: fname
        !! File name containing data of interest.
    character(len=*), intent(in), optional :: ftype
        !! File type, superseeds the automatic search.
    logical, intent(in), optional :: get_dens
        !! Load electronic density from data file.
    type(MoleculeDB) :: mol
        !! Molecular specifications database.

end function build_mol_data

! ----------------------------------------------------------------------

module function build_orb_data(dfile, fname, ftype) result(orb)
    class(DataFile), intent(inout), target, optional :: dfile
        !! DataFile instance.
    character(len=*), intent(in), optional :: fname
        !! File name containing data of interest.
    character(len=*), intent(in), optional :: ftype
        !! File type, superseeds the automatic search.
    type(OrbitalsDB) :: orb
        !! Orbitals information database.

end function build_orb_data

! ----------------------------------------------------------------------

module function build_vib_data(dfile, fname, ftype, get_Lmat, get_Lmweig &
                               )  result(vib)
    class(DataFile), intent(inout), target, optional :: dfile
        !! DataFile instance.
    character(len=*), intent(in), optional :: fname
        !! File name containing data of interest.
    character(len=*), intent(in), optional :: ftype
        !! File type, superseeds the automatic search.
    logical, intent(in), optional :: get_Lmat
        !! Build/load dimensionless matrix of Hessian eigenvectors.
    logical, intent(in), optional :: get_Lmweig
        !! Build/load mass-weighted matrix of Hessian eigenvectors.
    type(VibrationsDB) :: vib
        !! Vibrational information.

end function build_vib_data

! ----------------------------------------------------------------------

module function check_prog_version(major, minor, prog_info, file_data) &
        result(res)
    character(len=*), intent(in) :: major
        !! Major version.
    character(len=*), intent(in), optional :: minor
        !! Minor revision.
    class(ProgramInfo), intent(in), target, optional :: prog_info
        !! Instance of ProgramInfo.
    class(DataFile), intent(in), target, optional :: file_data
        !! Instance of DataFile.
    logical :: res
        !! Result of the query.

end function check_prog_version

! ----------------------------------------------------------------------

module function get_data_from_id(dfile, identifier, start_state, end_state, &
                                 derorder) result(prop)
    class(DataFile), intent(inout) :: dfile
        !! DataFile instance.
    integer, intent(in) :: identifier
        !! Identifier of the property of interest.
    integer, intent(in), optional :: start_state
        !! Starting or reference electronic state.
    integer, intent(in), optional :: end_state
        !! End electronic state, only for electronic transition.
    integer, intent(in), optional :: derorder
        !! Derivative order, if relevant or assumed to be 0.
    type(PropertyDB) :: prop
        !! Property information.

end function get_data_from_id

! ----------------------------------------------------------------------

module function get_data_from_tag(dfile, name, tag, start_state, end_state, &
                                  derorder) result(prop)
    class(DataFile), intent(inout) :: dfile
        !! DataFile instance.
    character(len=*), intent(in) :: name
        !! name/group name of the quantity of interest.
    character(len=*), intent(in), optional :: tag
        !! tag of the quantity within group.
    integer, intent(in), optional :: start_state
        !! Starting or reference electronic state.
    integer, intent(in), optional :: end_state
        !! End electronic state, only for electronic transition.
    integer, intent(in), optional :: derorder
        !! Derivative order, if relevant or assumed to be 0.
    type(PropertyDB) :: prop
        !! Property information.

end function get_data_from_tag

! ----------------------------------------------------------------------

module function get_datafile_name(dfile) result(name)
    class(DataFile), intent(in) :: dfile
        !! DataFile instance.
    character(len=:), allocatable :: name
        !! Name of the file associated to datafile instance.

end function get_datafile_name

! ----------------------------------------------------------------------

module function get_datafile_type(dfile) result(dtype)
    class(DataFile), intent(in) :: dfile
        !! DataFile instance.
    character(len=:), allocatable :: dtype
        !! DataFile type.

end function get_datafile_type

! ----------------------------------------------------------------------

module function get_file_type(fname, read_file, soft_check, err) result(ftype)
    character(len=*), intent(in) :: fname
        !! File name.
    logical, intent(in), optional :: read_file
        !! Read file content to guess the type.
    logical, intent(in), optional :: soft_check
        !! Perform soft check, trusting extension if unequivocal.
    type(ErrorHandle), intent(out), optional :: err
        !! Error instance.  If not present, then use runtime error handler.
        !! The error is supposed to be configured before.
    character(len=:), allocatable :: ftype
        !! File type.

end function get_file_type

! ----------------------------------------------------------------------

module function get_prog_major(prog_info, file_data) result(version)
    class(ProgramInfo), intent(in), optional :: prog_info
        !! Instance of ProgramInfo.
    class(DataFile), intent(in), optional :: file_data
        !! Instance of DataFile.
    character(len=:), allocatable :: version
        !! Major revision.

end function get_prog_major

! ----------------------------------------------------------------------

module function get_prog_minor(prog_info, file_data) result(version)
    class(ProgramInfo), intent(in), optional :: prog_info
        !! Instance of ProgramInfo.
    class(DataFile), intent(in), optional :: file_data
        !! Instance of DataFile.
    character(len=:), allocatable :: version
        !! Minor revision.

end function get_prog_minor

! ----------------------------------------------------------------------

module function get_prog_name(prog_info, file_data) result(name)
    class(ProgramInfo), intent(in), optional :: prog_info
        !! Instance of ProgramInfo.
    class(DataFile), intent(in), optional :: file_data
        !! Instance of DataFile.
    character(len=:), allocatable :: name
        !! Name of file associated to data file.

end function get_prog_name

! ----------------------------------------------------------------------

module function get_prog_version(prog_info, file_data) result(version)
    class(ProgramInfo), intent(in), optional :: prog_info
        !! Instance of ProgramInfo.
    class(DataFile), intent(in), optional :: file_data
        !! Instance of DataFile.
    character(len=:), allocatable :: version
        !! Type of file storing original data.

end function get_prog_version

! ----------------------------------------------------------------------

end interface

end module input
