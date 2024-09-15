module input
    !! Input-processing module
    !!
    !! Contains data and procedures to process input data/options.
    use iso_fortran_env, only: real64
    use exception, only: BaseException
    use workdata

    implicit none

    private
    public :: build_bset_data, build_exc_data, build_orb_data, build_vib_data

    type, public :: ProgramInfo
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

    type, public :: DataFile
        private
        character(len=:), allocatable :: name
        character(len=:), allocatable :: type
        type(ProgramInfo) :: prog
        class(BaseException), allocatable :: error
    contains
        procedure :: build_mol_data, build_bset_data, build_orb_data, &
            build_exc_data, build_vib_data
        procedure :: get_name => get_datafile_name
        procedure :: get_type => get_datafile_type
        procedure :: get_error_type => get_error_instance
        procedure, pass(file_data) :: get_program => get_prog_name
        procedure, pass(file_data) :: get_version => get_prog_version
        procedure :: has_error => check_error_status
        procedure :: get_error => get_error_msg
        procedure, pass(file_data) :: check_version => check_prog_version
    end type DataFile

    interface DataFile
        module procedure init_file
    end interface DataFile

    interface ProgramInfo
        module procedure get_program_version
    end interface ProgramInfo

interface

! ----------------------------------------------------------------------

module function init_file(fname, ftype) result(file)
    character(len=*), intent(in) :: fname
    !! File name.
    character(len=*), intent(in), optional :: ftype
    !! File type, which overrides any detection attempt.
    type(DataFile) :: file
    !! DataFile instance

end function init_file

! ----------------------------------------------------------------------

module function get_file_type(fname, read_file, err) result(ftype)
    character(len=*), intent(in) :: fname
    !! File name.
    logical, intent(in), optional :: read_file
    !! Read file content to guess the type.
    class(BaseException), allocatable :: err
    !! Error instance.
    character(len=:), allocatable :: ftype
    !! File type.

end function get_file_type

! ----------------------------------------------------------------------

module function get_program_version(fname, ftype, err) result(prog)
    character(len=*), intent(in) :: fname
    !! File name.
    character(len=*), intent(in) :: ftype
    !! File type.
    class(BaseException), allocatable, intent(out), optional :: err
    !! Error instance.
    type(ProgramInfo) :: prog
    !! ProgramInfo instance with program information.

end function get_program_version

! ----------------------------------------------------------------------

module function get_prog_name(prog_info, file_data) result(name)
    class(ProgramInfo), intent(in), optional :: prog_info
    !! Instance of ProgramInfo.
    class(DataFile), intent(in), optional :: file_data
    !! Instance of DataFile.
    character(len=:), allocatable :: name
    !! Program name.

end function get_prog_name

! ----------------------------------------------------------------------

module function get_prog_version(prog_info, file_data) result(version)
    class(ProgramInfo), intent(in), optional :: prog_info
    !! Instance of ProgramInfo.
    class(DataFile), intent(in), optional :: file_data
    !! Instance of DataFile.
    character(len=:), allocatable :: version
    !! Program version.

end function get_prog_version

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
    !! Result of the query

end function check_prog_version

! ----------------------------------------------------------------------

module function get_datafile_name(dfile) result(err)
    class(DataFile), intent(in) :: dfile
    !! DataFile instance.
    character(len=:), allocatable :: err
    !! Error instance.

end function get_datafile_name

! ----------------------------------------------------------------------

module function get_datafile_type(dfile) result(err)
    class(DataFile), intent(in) :: dfile
    !! DataFile instance.
    character(len=:), allocatable :: err
    !! Error instance.

end function get_datafile_type

! ----------------------------------------------------------------------

module function get_error_instance(dfile) result(err)
    class(DataFile), intent(in) :: dfile
    !! DataFile instance.
    class(BaseException), allocatable :: err
    !! Error instance.

end function get_error_instance

! ----------------------------------------------------------------------

module function check_error_status(dfile) result(raised)
    class(DataFile), intent(in) :: dfile
    !! DataFile instance.
    logical :: raised
    !! Status of the error

end function check_error_status

! ----------------------------------------------------------------------

module function get_error_msg(dfile) result(err)
    class(DataFile), intent(in) :: dfile
    !! DataFile instance.
    character(len=:), allocatable :: err
    !! Error instance.

end function get_error_msg

! ----------------------------------------------------------------------

module function build_mol_data(dfile, fname, ftype, err) result(mol)
    class(DataFile), intent(in), target, optional :: dfile
    !! DataFile instance.
    character(len=*), intent(in), optional :: fname
    !! File name containing data of interest.
    character(len=*), intent(in), optional :: ftype
    !! File type, superseeds the automatic search.
    class(BaseException), allocatable, intent(out), optional :: err
    !! Error instance
    type(MoleculeDB) :: mol
    !! Molecular specifications database.

end function build_mol_data

! ----------------------------------------------------------------------

module function build_bset_data(dfile, fname, ftype, err) result(bset)
    class(DataFile), intent(in), target, optional :: dfile
    !! DataFile instance.
    character(len=*), intent(in), optional :: fname
    !! File name containing data of interest.
    character(len=*), intent(in), optional :: ftype
    !! File type, superseeds the automatic search.
    class(BaseException), allocatable, intent(out), optional :: err
    !! Error instance
    type(BasisSetDB) :: bset
    !! Basis set information database.

end function build_bset_data

! ----------------------------------------------------------------------

module function build_orb_data(dfile, fname, ftype, err) result(orb)
    class(DataFile), intent(in), target, optional :: dfile
    !! DataFile instance.
    character(len=*), intent(in), optional :: fname
    !! File name containing data of interest.
    character(len=*), intent(in), optional :: ftype
    !! File type, superseeds the automatic search.
    class(BaseException), allocatable, intent(out), optional :: err
    !! Error instance
    type(OrbitalsDB) :: orb
    !! Orbitals information database.


end function build_orb_data

! ----------------------------------------------------------------------

module function build_exc_data(dfile, fname, ftype, err) result(exc)
    class(DataFile), intent(in), target, optional :: dfile
    !! DataFile instance.
    character(len=*), intent(in), optional :: fname
    !! File name containing data of interest.
    character(len=*), intent(in), optional :: ftype
    !! File type, superseeds the automatic search.
    class(BaseException), allocatable, intent(out), optional :: err
    !! Error instance
    type(ExcitationDB) :: exc
    !! Electronic excitation information.

end function build_exc_data

! ----------------------------------------------------------------------

module function build_vib_data(dfile, fname, ftype, get_Lmat, get_Lmweig, &
                               err)  result(vib)
    class(DataFile), intent(in), target, optional :: dfile
    !! DataFile instance.
    character(len=*), intent(in), optional :: fname
    !! File name containing data of interest.
    character(len=*), intent(in), optional :: ftype
    !! File type, superseeds the automatic search.
    logical, intent(in), optional :: get_Lmat
    !! Build/load dimensionless matrix of Hessian eigenvectors.
    logical, intent(in), optional :: get_Lmweig
    !! Build/load mass-weighted matrix of Hessian eigenvectors.
    class(BaseException), allocatable, intent(out), optional :: err
    !! Error instance
    type(VibrationsDB) :: vib
    !! Vibrational information.

end function build_vib_data

! ----------------------------------------------------------------------

end interface

end module input
