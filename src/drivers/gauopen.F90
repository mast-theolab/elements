module gauopen_drv
    !! Provide modern Fortran interfaces to GauOpen library.
    !!
    !! @note
    !! The source code of the GauOpen library contains two versions, depending
    !! on the storage used to store integer labels:
    !! 32 bits (I4) and 64 bits (I8)
    !! The versions are "chosen" through preprocessor directives and the code
    !! consistency with compiler directives.
    !! As a result, the standard compilation process uses the same names
    !! (symbols) for the two versions, which are then linked in separate
    !! executables or in interfaces with higher-level languages.
    !! To provide a generic interface, we need to use a generator that create
    !! distinct symbol names for the 32 and 64 bits versions.
    !! The interface here expects that the generator has been used and refers
    !! to the created files, with the original names appended with `_I4` or
    !! `_I8`.
    !! The interface recreates here the original names as proper interfaces,
    !! considering all possible versions.
    !! @endnote
    !!
    !! @warning
    !! To make the interface work with standard modern Fortran, a few
    !! hypotheses are made, that should be reasonable considering practical
    !! implementations of the language (compiler) and general programming
    !! practice, but are not required to be true by the standard.
    !!
    !! * logical use the same kinds as integers.  The standard refers to the
    !!   same storage size (NSU) as integer by default, but not that kinds must
    !!   match.  In practice, it seems most compilers use the same option to
    !!   set the representation of logical and integer, and the kind typically
    !!   matches the actual size (i.e., number of bytes).
    !! * all integers and logical are consistent: the functions/interfaces are
    !!   defined assuming that the arguments have consistent kinds.  This is
    !!   also in line with the implementation of GauOpen
    !! * `Real*8` is equivalent to the kind `real64`.  This should be the case
    !!   but the former is not standard and thus could be in theory interpreted
    !!   in a different way.
    !! @endwarning

    use iso_fortran_env, only: int32, int64, real64
    use run_env, only: run

    implicit none

    private :: Open_Read_32, Open_Read_64, Open_Write_32, Open_Write_64

    ! ===============================================
    !  GOpen Interfaces for Label-Dependent Routines
    ! ===============================================

    ! ----------------------------------------------------------------------

    interface Open_Read_GOpen
        subroutine Open_Read_I4(Name, IU, LabFil, IVers, NLab, GVers, Title, &
                                NAtoms, NBasis, NBsUse, ICharg, Multip, NE, &
                                Len12L, Len4L, IOpCl, ICGU)
            use iso_fortran_env, only: int32
            character(len=*), intent(in) :: Name
            character(len=64), intent(out) :: LabFil, GVers, Title
            integer(int32), intent(out) :: IU, IVers, NLab, NAtoms, NBasis, &
                NBsUse, ICharg, Multip, NE, Len12L, Len4L, IOpCl, ICGU
        end subroutine Open_Read_I4
        subroutine Open_Read_I8(Name, IU, LabFil, IVers, NLab, GVers, Title, &
                                NAtoms, NBasis, NBsUse, ICharg, Multip, NE, &
                                Len12L, Len4L, IOpCl, ICGU)
            use iso_fortran_env, only: int64
            character(len=*), intent(in) :: Name
            character(len=64), intent(out) :: LabFil, GVers, Title
            integer(int64), intent(out) :: IU, IVers, NLab, NAtoms, NBasis, &
                NBsUse, ICharg, Multip, NE, Len12L, Len4L, IOpCl, ICGU

        end subroutine Open_Read_I8
    end interface Open_Read_GOpen

    ! ----------------------------------------------------------------------

    interface Open_Write_GOpen
        Subroutine Open_Write_I4(Name, IU, LabFil, GVers, Title, NAtoms, &
                                 NBasis, NBsUse, ICharg, Multip, NE, IOpCl, &
                                 ICGU)
            use iso_fortran_env, only: int32
            character(len=*), intent(in) :: Name, LabFil, GVers, Title
            integer(int32), intent(in) :: NAtoms, NBasis, NBsUse, ICharg, &
                Multip, NE, IOpCl, ICGU
            integer(int32), intent(out) :: IU
        end subroutine Open_Write_I4
        subroutine Open_Write_I8(Name, IU, LabFil, GVers, Title, NAtoms, &
                                 NBasis, NBsUse, ICharg, Multip, NE, IOpCl, &
                                 ICGU)
            use iso_fortran_env, only: int64
            character(len=*), intent(in) :: Name, LabFil, GVers, Title
            integer(int64), intent(in) :: NAtoms, NBasis, NBsUse, ICharg, &
                Multip, NE, IOpCl, ICGU
            integer(int64), intent(out) :: IU
        end subroutine Open_Write_I8
    end interface Open_Write_GOpen

    ! ----------------------------------------------------------------------

    interface Rd_2E1_GOpen
        subroutine Rd_2E1_I4(IU, LR, NTot, LenBuf, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, LenBuf, LR, NTot
            real(real64), intent(out) :: RArr(LR)
        end subroutine Rd_2E1_I4
        subroutine Rd_2E1_I8(IU, LR, NTot, LenBuf, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, LenBuf, LR, NTot
            real(real64), intent(out) :: RArr(LR)
        end subroutine Rd_2E1_I8
    end interface Rd_2E1_GOpen

    ! ----------------------------------------------------------------------

    interface Rd_2EN_GOpen
#if GAUOPEN >= 3
        subroutine Rd_2EN_I4(IU, NR, LR, NTot, LenBuf, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, LenBuf, LR, NR, NTot
            real(real64), intent(out) :: RArr(NR,LR)
        end subroutine Rd_2EN_I4
        subroutine Rd_2EN_I8(IU, NR, LR, NTot, LenBuf, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, LenBuf, LR, NR, NTot
            real(real64), intent(out) :: RArr(NR,LR)
        end subroutine Rd_2EN_I8
#else
        subroutine Rd_2EN_I4(IU, NR, LR, LRNR, NTot, LenBuf, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, LenBuf, LR, LRNR, NR, NTot
            real(real64), intent(out) :: RArr(LRNR)
        end subroutine Rd_2EN_I4
        subroutine Rd_2EN_I8(IU, NR, LR, LRNR, NTot, LenBuf, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, LenBuf, LR, LRNR, NR, NTot
            real(real64), intent(out) :: RArr(LRNR)
        end subroutine Rd_2EN_I8
#endif
    end interface Rd_2EN_GOpen

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface Rd_SpA_GOpen
        subroutine Rd_SpA_I4(IU, NI, NR, NTot, LenBuf, IArr, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, LenBuf, NI, NR, NTot
            integer(int32), intent(out) :: IArr(NI*NTot)
            real(real64), intent(out) :: RArr(NR*NTot)
        end subroutine Rd_SpA_I4
        subroutine Rd_SpA_I8(IU, NI, NR, NTot, LenBuf, IArr, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, LenBuf, NI, NR, NTot
            integer(int64), intent(out) :: IArr(NI*NTot)
            real(real64), intent(out) :: RArr(NR*NTot)
        end subroutine Rd_SpA_I8
    end interface Rd_SpA_GOpen
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface Rd_SpAC_GOpen
        !! Read and return a list of values and indices for a real or complex
        !! array which is kept in sparse form.
        !! `NR` is the number of values, `NI` is the number of indices and
        !! `NTot` is the number of value sets.
        subroutine Rd_SpAC_I4(IU, NI, NR, NTot, LenBuf, IArr, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, LenBuf, NI, NR, NTot
            integer(int32), intent(out) :: IArr(NI*NTot)
            complex(real64), intent(out) :: RArr(NR*NTot)
        end subroutine Rd_SpAC_I4
        subroutine Rd_SpAC_I8(IU, NI, NR, NTot, LenBuf, IArr, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, LenBuf, NI, NR, NTot
            integer(int64), intent(out) :: IArr(NI*NTot)
            complex(real64), intent(out) :: RArr(NR*NTot)
        end subroutine Rd_SpAC_I8
    end interface Rd_SpAC_GOpen
#endif

    ! ----------------------------------------------------------------------

    interface Rd_Head_GOpen
        subroutine Rd_Head_I4(IU, NLab, NAtoms, NBasis, IAn, IAtTyp, AtmChg, &
                              C, IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, &
                              IDum9, NShlAO, NPrmAO, NShlDB, NPrmDB, NBTot)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, NLab, NAtoms, NBasis
            integer(int32), intent(out) :: IAn(NAtoms), IAtTyp(NAtoms), &
                IBfAtm(NBasis), IBfTyp(NBasis), NFC, NFV, ITran, IDum9, &
                NShlAO, NPrmAO, NShlDB, NPrmDB, NBTot
            real(real64), intent(out) :: AtmChg(NAtoms), C(3*NAtoms), &
                AtmWgt(NAtoms)
        end subroutine Rd_Head_I4
        subroutine Rd_Head_I8(IU, NLab, NAtoms, NBasis, IAn, IAtTyp, AtmChg, &
                              C, IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, &
                              IDum9, NShlAO, NPrmAO, NShlDB, NPrmDB, NBTot)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, NLab, NAtoms, NBasis
            integer(int64), intent(out) :: IAn(NAtoms), IAtTyp(NAtoms), &
                IBfAtm(NBasis), IBfTyp(NBasis), NFC, NFV, ITran, IDum9, &
                NShlAO, NPrmAO, NShlDB, NPrmDB, NBTot
            real(real64), intent(out) :: AtmChg(NAtoms), C(3*NAtoms), &
                AtmWgt(NAtoms)
        end subroutine Rd_Head_I8
    end interface Rd_Head_GOpen

    ! ----------------------------------------------------------------------

    interface Rd_HeadA_GOpen
        subroutine Rd_HeadA_I4(IU, NLab, NAtoms, NBasis, IAn, IAtTyp, AtmChg, &
                               C, IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, &
                               IDum9, Rec11)
            use iso_fortran_env, only: int32, real64
            integer, parameter :: LRec11 = 16
            integer(int32), intent(in) :: IU, NLab, NAtoms, NBasis
            integer(int32), intent(out) :: IAn(NAtoms), IAtTyp(NAtoms), &
                IBfAtm(NBasis), IBfTyp(NBasis), NFC, NFV, ITran, IDum9, &
                Rec11(LRec11)
            real(real64), intent(out) :: AtmChg(NAtoms), C(3*NAtoms), &
                AtmWgt(NAtoms)
        end subroutine Rd_HeadA_I4
        subroutine Rd_HeadA_I8(IU, NLab, NAtoms, NBasis, IAn, IAtTyp, AtmChg, &
                               C, IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, &
                               IDum9, Rec11)
            use iso_fortran_env, only: int64, real64
            integer, parameter :: LRec11 = 16
            integer(int64), intent(in) :: IU, NLab, NAtoms, NBasis
            integer(int64), intent(out) :: IAn(NAtoms), IAtTyp(NAtoms), &
                IBfAtm(NBasis), IBfTyp(NBasis), NFC, NFV, ITran, IDum9, &
                Rec11(LRec11)
            real(real64), intent(out) :: AtmChg(NAtoms), C(3*NAtoms), &
                AtmWgt(NAtoms)
        end subroutine Rd_HeadA_I8
    end interface Rd_HeadA_GOpen

    ! ----------------------------------------------------------------------

    interface Rd_IBuf_GOpen
        !! Read an integer array from the file given parameters from the header
        !! record.
        subroutine Rd_IBuf_I4(IU, LR, LenBuf, Arr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, LenBuf, LR
            integer(int32), intent(out) :: Arr(LR)
        end subroutine Rd_IBuf_I4
        subroutine Rd_IBuf_I8(IU, LR, LenBuf, Arr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, LenBuf, LR
            integer(int64), intent(out) :: Arr(LR)
        end subroutine Rd_IBuf_I8
    end interface Rd_IBuf_GOpen

    ! ----------------------------------------------------------------------

    interface Rd_Labl_GOpen
        !! Read the label record for an operator matrix.
        subroutine Rd_Labl_I4(IU, IVers, CBuf, NI, NR, NTot, LenBuf, N1, N2, &
                              N3, N4, N5, TypeA, NRI, EOF)
            use iso_fortran_env, only: int32
            integer(int32), intent(in) :: IU, IVers
            integer(int32), intent(out) :: NI, NR, NTot, LenBuf, N1, N2, N3, &
                N4, N5, TypeA, NRI
            character(len=64), intent(out) :: CBuf
            Logical(int32), intent(out) :: EOF
        end subroutine Rd_Labl_I4
        subroutine Rd_Labl_I8(IU, IVers, CBuf, NI, NR, NTot, LenBuf, N1, N2, &
                              N3, N4, N5, TypeA, NRI, EOF)
            use iso_fortran_env, only: int64
            integer(int64), intent(in) :: IU, IVers
            integer(int64), intent(out) :: NI, NR, NTot, LenBuf, N1, N2, N3, &
                N4, N5, TypeA, NRI
            character(len=64), intent(out) :: CBuf
            Logical(int64), intent(out) :: EOF
        end subroutine Rd_Labl_I8
    end interface Rd_Labl_GOpen

    ! ----------------------------------------------------------------------

    interface Rd_RInd_GOpen
        !! Read a real array stored with indices for non-zero elements.
        subroutine Rd_RInd_I4(IU, NR, LR, NRLR, NTot, LenBuf, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, NR, LR, NRLR, NTot, LenBuf
            real(real64), intent(out) :: RArr(NRLR)
        end subroutine Rd_RInd_I4
        subroutine Rd_RInd_I8(IU, NR, LR, NRLR, NTot, LenBuf, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, NR, LR, NRLR, NTot, LenBuf
            real(real64), intent(out) :: RArr(NRLR)
        end subroutine Rd_RInd_I8
    end interface Rd_RInd_GOpen

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface Wr_LChBuf_GOpen
        subroutine Wr_LChBuf_I4(IU, Label, LenBuf, X)
            use iso_fortran_env, only: int32
            character(len=*), intent(in) :: Label, X
            integer(int32), intent(in) :: IU, LenBuf
        end subroutine Wr_LChBuf_I4
        subroutine Wr_LChBuf_I8(IU, Label, LenBuf, X)
            use iso_fortran_env, only: int64
            character(len=*), intent(in) :: Label, X
            integer(int64), intent(in) :: IU, LenBuf
        end subroutine Wr_LChBuf_I8
    end interface Wr_LChBuf_GOpen
#endif

    ! ----------------------------------------------------------------------

    ! =======================================================
    !  Interfaces to ELEMENTS Routines as Wrapper of Fortran
    ! =======================================================

    ! ----------------------------------------------------------------------

    interface Open_Read
        !! Interface to Open_Read provided by `qcmatrixio`.
        !!
        !! Open the named binary array file for reading, and return
        !! the listed scalars from the initial 2 header records.
        !! `IU` receives the Fortran unit number of the open file,
        !! or `-1` if the open failed.
        !!
        !! @note
        !! 4 variants: 32/64-bits labels, 32/64-bits input integer
        !! 2 versions are provided by the generator, cross-versions need to be
        !! defined.
        !! Since two versions are not distinguishable by arguments, we add an
        !! optional argument.
        !! @endnote
        module procedure Open_Read_32, Open_Read_64

    end interface

    ! ----------------------------------------------------------------------

    interface Open_Write
        !! Interface to Open_Write provided by `qcmatrixio`.
        !!
        !!  Open the named binary array file for writing and write the
        !! named scalars to the initial 2 header records.
        !! `IU` receives the Fortran unit number of the open file,
        !! or -1 if the open failed.
        !!
        !! @note
        !! 4 variants: 32/64-bits labels, 32/64-bits input integer
        !! 2 versions are provided by the generator, cross-versions need to be
        !! defined
        !! @endnote
        module procedure Open_Write_32, Open_Write_64
    end interface

    ! ----------------------------------------------------------------------

    interface Rd_2E1
        !! Read and return an array of AO 2e integrals.
        !!
        !! Read and return an array of AO 2e integrals (a 4D array with
        !! quartets of indices and one value per index set).
        !! `NTot` is the number of non-zero values (from the header record
        !! for the object) and `LR` is the total number of elements
        !! (LenArr_I4(-N,-N,-N,N,1), where N is the number of basis functions).
        module procedure Rd_2E1_32, Rd_2E1_64
    end interface Rd_2E1

    ! ----------------------------------------------------------------------

    interface Rd_2EN
        !! Read and return an array of AO 2e integrals.
        !!
        !! Read and return an array of AO 2e integrals (a 4D array with
        !! quartets of indices and `NR` values per index set).
        !! `NTot` is the number of non-zero values (from the header record
        !! for the object) and `LR` is the total number of elements
        !! (LenArr_I4(-N,-N,-N,N,1), where N is the number of basis functions).
        module procedure Rd_2EN_32, Rd_2EN_64, Rd_2EN_1D_32, Rd_2EN_1D_64
    end interface Rd_2EN

    ! ----------------------------------------------------------------------

    interface Rd_SpA
        !! Read and return a list of values and indices for a real or complex
        !! array which is kept in sparse form.
        !! `NR` is the number of values, `NI` is the number of indices and
        !! `NTot` is the number of value sets.
        module procedure Rd_SpA_32, Rd_SpA_64
    end interface Rd_SpA

    ! ----------------------------------------------------------------------

    interface Rd_SpAC
        !! Read and return a list of values and indices for a real or complex
        !! array which is kept in sparse form.
        !! `NR` is the number of values, `NI` is the number of indices and
        !! `NTot` is the number of value sets.
        module procedure Rd_SpAC_32, Rd_SpAC_64
    end interface Rd_SpAC

    ! ----------------------------------------------------------------------

    interface Rd_Head
        !! Read the remaining header records (3 to NLab) and return in distinct
        !! arrays.
        module procedure Rd_Head_32, Rd_Head_64
    end interface Rd_Head

    ! ----------------------------------------------------------------------

    interface Rd_HeadA
        !! Read the remaining header records (3 to NLab) and return in distinct
        !! arrays.
        !! This version also returns the full array of 16 stored in record 11.
        module procedure Rd_HeadA_32, Rd_HeadA_64
    end interface Rd_HeadA

    ! ----------------------------------------------------------------------

    interface Rd_IBuf
        !! Read an integer array from the file given parameters from the header
        !! record.
        module procedure Rd_IBuf_32, Rd_IBuf_64
    end interface Rd_IBuf

    ! ----------------------------------------------------------------------

    interface Rd_Labl
        !! Read the label record for an operator matrix.
        module procedure Rd_Labl_32, Rd_Labl_64
    end interface Rd_Labl

    ! ----------------------------------------------------------------------

    interface Rd_RInd
        !! Read a real array stored with indices for non-zero elements.
        module procedure Rd_RInd_32, Rd_RInd_64
    end interface Rd_RInd

    ! ----------------------------------------------------------------------

    interface Wr_LChBuf
        !! Interface to Wr_LChBuf provided by `qcmatrix`.
        !!
        !! Write character string `X` unblocked (blanks included) to unit `IU`
        !! given the parameters from the header record.
        !!
        !! @note
        !! 4 variants: 32/64-bits labels, 32/64-bits input integer
        !! 2 versions are provided by the generator, cross-versions need to be
        !! defined
        !! @endnote
        module procedure Wr_LChBuf_32, Wr_LChBuf_64
    end interface

    ! ----------------------------------------------------------------------

    ! =========================================
    !  Interfaces to Other QCMATRIXIO Routines
    ! =========================================

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface Close_FAF
        !! Close a Matrix array file.
        !!
        !! Close the file open on Fortran unit `IU`
        subroutine Close_FAF_I4(DoWrt, IU)
            use iso_fortran_env, only: int32
            integer(int32), intent(in) :: IU
            logical(int32), intent(in) :: DoWrt
        end subroutine Close_FAF_I4
        subroutine Close_FAF_I8(DoWrt, IU)
            use iso_fortran_env, only: int64
            integer(int64), intent(in) :: IU
            logical(int64), intent(in) :: DoWrt
        end subroutine Close_FAF_I8
    end interface Close_FAF
#else
    interface Close_MatF
        !! Close a binary array file.
        !!
        !! Close the file open on Fortran unit `IU`, writing the END record
        !! if `DoWrt` is true.
        subroutine Close_MatF_I4(DoWrt, IU)
            use iso_fortran_env, only: int32
            integer(int32), intent(in) :: IU
            logical(int32), intent(in) :: DoWrt
        end subroutine Close_MatF_I4
        subroutine Close_MatF_I8(DoWrt, IU)
            use iso_fortran_env, only: int64
            integer(int64), intent(in) :: IU
            logical(int64), intent(in) :: DoWrt
        end subroutine Close_MatF_I8
    end interface Close_MatF
#endif

    ! ----------------------------------------------------------------------

    interface LenArr
        !! Return the total number of index values of an array.
        !!
        !! Returns the total number of index values of an array with the
        !! specified dimensions, accounting for possible lower-triangular
        !! indices.
        !! This does not include a possible multiple number of values per index
        !! (NI or NR in the record header/nelem in the object).
        !!
        !! @note
        !! The function acts as a wrapper to `QCM_LenArr`.
        !! @endnote
        function LenArr_I4(N1, N2, N3, N4, N5)
            use iso_fortran_env, only: int32
            integer(int32) :: LenArr_I4
            integer(int32), intent(in) :: N1, N2, N3, N4, N5
        end function LenArr_I4
        function LenArr_I8(N1, N2, N3, N4, N5)
            use iso_fortran_env, only: int64
            integer(int64) :: LenArr_I8
            integer(int64), intent(in) :: N1, N2, N3, N4, N5
        end function LenArr_I8
    end interface LenArr

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_LenArray
        !! Return the total number of index values of an array.
        !!
        !! Returns the total number of index values of an array with the
        !! specified dimensions, accounting for possible lower-triangular
        !! indices.
        !! This does not include a possible multiple number of values per index
        !! (NI or NR in the record header/nelem in the object).
        !!
        !! @note "version"
        !! This is the core version.
        !! @endnote
        function QCM_LenArray_I4(N1, N2, N3, N4, N5)
            use iso_fortran_env, only: int32
            integer(int32) :: QCM_LenArray_I4
            integer(int32), intent(in) :: N1, N2, N3, N4, N5
        end function QCM_LenArray_I4
        function QCM_LenArray_I8(N1, N2, N3, N4, N5)
            use iso_fortran_env, only: int64
            integer(int64) :: QCM_LenArray_I8
            integer(int64), intent(in) :: N1, N2, N3, N4, N5
        end function QCM_LenArray_I8
    end interface QCM_LenArray
#endif

    ! ----------------------------------------------------------------------

    interface LInd2C
        !! Return the 0-based linear index given 2 0-based indices.
        !!
        !! Returns the linear or square indexing, with indexes `I`, `J`
        !! 0-based given the C order.  The output is 0-based.
        !! If `CkArg` is true, the routine checks indices in range and the
        !! function returns `-1` if out of range.
        !! Returns the 0-based index into a linear array given 2 indices and
        !! dimensions.
        !! If `CkArg` is true, the function checks if indices are in range and
        !! returns -1 if they are out of range, otherwise the indices are
        !! assumed to be valid.
        !! `JSign` takes a +/-1 value whether the upper or lower triangle was
        !! selected (i.e., whether to apply a `JSign` flip for anti-symmetric/
        !! Hermitian matrices and/or take a complex conjugate.
        !!
        !! @note "version"
        !! Version compatible with C indexing.
        !! @endnote
        function LInd2C_I4(CkArg, N1, N2, I, J, JSign)
            use iso_fortran_env, only: int32
            integer(int32) :: LInd2C_I4
            logical(int32), intent(in) :: CkArg
            integer(int32), intent(in) :: N1, N2, I, J
            integer(int32), intent(out) :: JSign
        end function LInd2C_I4
        function LInd2C_I8(CkArg, N1, N2, I, J, JSign)
            use iso_fortran_env, only: int64
            integer(int64) :: LInd2C_I8
            logical(int64), intent(in) :: CkArg
            integer(int64), intent(in) :: N1, N2, I, J
            integer(int64), intent(out) :: JSign
        end function LInd2C_I8
    end interface LInd2C

    ! ----------------------------------------------------------------------

    interface Lind2
        !! Return the 0-based linear index given 2 1-based indices.
        !!
        !! Returns the linear or square indexing, with indexes `I`, `J`
        !! 1-based given the F order.  The output is 0-based.
        !! If `CkArg` is true, the routine checks indices in range and the
        !! function returns `-1` if out of range.
        !! Returns the 0-based index into a linear array given 2 indices and
        !! dimensions.
        !! If `CkArg` is true, the function checks if indices are in range and
        !! returns -1 if they are out of range, otherwise the indices are
        !! assumed to be valid.
        !! `JSign` takes a +/-1 value whether the upper or lower triangle was
        !! selected (i.e., whether to apply a `JSign` flip for anti-symmetric/
        !! Hermitian matrices and/or take a complex conjugate.
        !!
        !! @note "version"
        !! Version compatible with Fortran indexing.
        !! @endnote
        function Lind2_I4(CkArg, N1, N2, I, J, JSign)
            use iso_fortran_env, only: int32
            integer(int32) :: Lind2_I4
            logical(int32), intent(in) :: CkArg
            integer(int32), intent(in) :: N1, N2, I, J
            integer(int32), intent(out) :: JSign
        end function Lind2_I4
        function Lind2_I8(CkArg, N1, N2, I, J, JSign)
            use iso_fortran_env, only: int64
            integer(int64) :: Lind2_I8
            logical(int64), intent(in) :: CkArg
            integer(int64), intent(in) :: N1, N2, I, J
            integer(int64), intent(out) :: JSign
        end function Lind2_I8
    end interface Lind2

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_Lind2
        !! Return the 0-based linear index given 2 1-based indices.
        !!
        !! Returns the linear or square indexing, with indexes `I`, `J`
        !! 1-based given the F order.  The output is 0-based.
        !! If `CkArg` is true, the routine checks indices in range and the
        !! function returns `-1` if out of range.
        !! Returns the 0-based index into a linear array given 2 indices and
        !! dimensions.
        !! If `CkArg` is true, the function checks if indices are in range and
        !! returns -1 if they are out of range, otherwise the indices are
        !! assumed to be valid.
        !! `JSign` takes a +/-1 value whether the upper or lower triangle was
        !! selected (i.e., whether to apply a `JSign` flip for anti-symmetric/
        !! Hermitian matrices and/or take a complex conjugate.
        !!
        !! @note "version"
        !! This is the core version.
        !! @endnote
        function QCM_Lind2_I4(N1, N2, I, J, JSign)
            use iso_fortran_env, only: int32
            integer(int32) :: QCM_Lind2_I4
            integer(int32), intent(in) :: N1, N2, I, J
            integer(int32), intent(out) :: JSign
        end function QCM_Lind2_I4
        function QCM_Lind2_I8(N1, N2, I, J, JSign)
            use iso_fortran_env, only: int64
            integer(int64) :: QCM_Lind2_I8
            integer(int64), intent(in) :: N1, N2, I, J
            integer(int64), intent(out) :: JSign
        end function QCM_Lind2_I8
    end interface QCM_Lind2

#endif
    ! ----------------------------------------------------------------------

    interface Lind3C
        !! Return the 0-based linear index given 3 0-based indices.
        !!
        !! Returns the linear or square indexing, with indexes `I`, `J`, `K`
        !! 0-based given the C order.  The output is 0-based.
        !! If `CkArg` is true, the routine checks indices in range and the
        !! function returns `-1` if out of range.
        !! Returns the 0-based index into a linear array given 2 indices and
        !! dimensions.
        !! If `CkArg` is true, the function checks if indices are in range and
        !! returns -1 if they are out of range, otherwise the indices are
        !! assumed to be valid.
        !! `JSign` takes a +/-1 value whether the upper or lower triangle was
        !! selected (i.e., whether to apply a `JSign` flip for anti-symmetric/
        !! Hermitian matrices and/or take a complex conjugate.
        !!
        !! @note "version"
        !! Version compatible with C indexing.
        !! @endnote
        function Lind3C_I4(CkArg, N1, N2, N3, I, J, K, JSign)
            use iso_fortran_env, only: int32
            integer(int32) :: Lind3C_I4
            logical(int32), intent(in) :: CkArg
            integer(int32), intent(in) :: N1, N2, N3, I, J, K
            integer(int32), intent(out) :: JSign
        end function Lind3C_I4
        function Lind3C_I8(CkArg, N1, N2, N3, I, J, K, JSign)
            use iso_fortran_env, only: int64
            integer(int64) :: Lind3C_I8
            logical(int64), intent(in) :: CkArg
            integer(int64), intent(in) :: N1, N2, N3, I, J, K
            integer(int64), intent(out) :: JSign
        end function Lind3C_I8
    end interface Lind3C

    ! ----------------------------------------------------------------------

    interface Lind3
        !! Return the 0-based linear index given 3 1-based indices.
        !!
        !! Returns the linear or square indexing, with indexes `I`, `J`, `K`
        !! 1-based given the F order.  The output is 0-based.
        !! If `CkArg` is true, the routine checks indices in range and the
        !! function returns `-1` if out of range.
        !! Returns the 0-based index into a linear array given 2 indices and
        !! dimensions.
        !! If `CkArg` is true, the function checks if indices are in range and
        !! returns -1 if they are out of range, otherwise the indices are
        !! assumed to be valid.
        !! `JSign` takes a +/-1 value whether the upper or lower triangle was
        !! selected (i.e., whether to apply a `JSign` flip for anti-symmetric/
        !! Hermitian matrices and/or take a complex conjugate.
        !!
        !! @note "version"
        !! Version compatible with Fortran indexing.
        !! @endnote
        function Lind3_I4(CkArg, N1, N2, N3, I, J, K, JSign)
            use iso_fortran_env, only: int32
            integer(int32) :: Lind3_I4
            logical(int32), intent(in) :: CkArg
            integer(int32), intent(in) :: N1, N2, N3, I, J, K
            integer(int32), intent(out) :: JSign
        end function Lind3_I4
        function Lind3_I8(CkArg, N1, N2, N3, I, J, K, JSign)
            use iso_fortran_env, only: int64
            integer(int64) :: Lind3_I8
            logical(int64), intent(in) :: CkArg
            integer(int64), intent(in) :: N1, N2, N3, I, J, K
            integer(int64), intent(out) :: JSign
        end function Lind3_I8
    end interface Lind3

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_Lind3
        !! Return the 0-based linear index given 3 1-based indices.
        !!
        !! Returns the linear or square indexing, with indexes `I`, `J`, `K`
        !! 1-based given the F order.  The output is 0-based.
        !! If `CkArg` is true, the routine checks indices in range and the
        !! function returns `-1` if out of range.
        !! Returns the 0-based index into a linear array given 2 indices and
        !! dimensions.
        !! If `CkArg` is true, the function checks if indices are in range and
        !! returns -1 if they are out of range, otherwise the indices are
        !! assumed to be valid.
        !! `JSign` takes a +/-1 value whether the upper or lower triangle was
        !! selected (i.e., whether to apply a `JSign` flip for anti-symmetric/
        !! Hermitian matrices and/or take a complex conjugate.
        !!
        !! @note "version"
        !! This is the core version.
        !! @endnote
        function QCM_Lind3_I4(N1, N2, N3, I, J, K, JSign)
            use iso_fortran_env, only: int32
            integer(int32) :: QCM_Lind3_I4
            integer(int32), intent(in) :: N1, N2, N3, I, J, K
            integer(int32), intent(out) :: JSign
        end function QCM_Lind3_I4
        function QCM_Lind3_I8(N1, N2, N3, I, J, K, JSign)
            use iso_fortran_env, only: int64
            integer(int64) :: QCM_Lind3_I8
            integer(int64), intent(in) :: N1, N2, N3, I, J, K
            integer(int64), intent(out) :: JSign
        end function QCM_Lind3_I8
    end interface QCM_Lind3
#endif

    ! ----------------------------------------------------------------------

    interface LInd4C
        !! Return the 0-based linear index given 4 0-based indices.
        !!
        !! Returns the linear or square indexing, with indexes `I`, `J`, `K`,
        !! `L` 0-based given the C order.  The output is 0-based.
        !! If `CkArg` is true, the routine checks indices in range and the
        !! function returns `-1` if out of range.
        !! Returns the 0-based index into a linear array given 2 indices and
        !! dimensions.
        !! If `CkArg` is true, the function checks if indices are in range and
        !! returns -1 if they are out of range, otherwise the indices are
        !! assumed to be valid.
        !! `JSign` takes a +/-1 value whether the upper or lower triangle was
        !! selected (i.e., whether to apply a `JSign` flip for anti-symmetric/
        !! Hermitian matrices and/or take a complex conjugate.
        !!
        !! @note "version"
        !! Version compatible with C indexing.
        !! @endnote
        !! Linear or square indexing, I,J,K,L are 0-based, c order
        !! output is 0-based.  JSign is +/-1.
        function LInd4C_I4(CkArg, N1, N2, N3, N4, I, J, K, L, JSign)
            use iso_fortran_env, only: int32
            integer(int32) :: LInd4C_I4
            logical(int32), intent(in) :: CkArg
            integer(int32), intent(in) :: N1, N2, N3, N4, I, J, K, L
            integer(int32), intent(out) :: JSign
        end function LInd4C_I4
        function LInd4C_I8(CkArg, N1, N2, N3, N4, I, J, K, L, JSign)
            use iso_fortran_env, only: int64
            integer(int64) :: LInd4C_I8
            logical(int64), intent(in) :: CkArg
            integer(int64), intent(in) :: N1, N2, N3, N4, I, J, K, L
            integer(int64), intent(out) :: JSign
        end function LInd4C_I8
    end interface LInd4C

    ! ----------------------------------------------------------------------

    interface LInd4
        !! Return the 0-based linear index given 4 1-based indices.
        !!
        !! Returns the linear or square indexing, with indexes `I`, `J`, `K`,
        !! `L` 1-based given the F order.  The output is 0-based.
        !! If `CkArg` is true, the routine checks indices in range and the
        !! function returns `-1` if out of range.
        !! Returns the 0-based index into a linear array given 2 indices and
        !! dimensions.
        !! If `CkArg` is true, the function checks if indices are in range and
        !! returns -1 if they are out of range, otherwise the indices are
        !! assumed to be valid.
        !! `JSign` takes a +/-1 value whether the upper or lower triangle was
        !! selected (i.e., whether to apply a `JSign` flip for anti-symmetric/
        !! Hermitian matrices and/or take a complex conjugate.
        !!
        !! @note "version"
        !! Version compatible with Fortran indexing.
        !! @endnote
        !! Linear or square indexing, I,J,K,L are 1-based,
        !! output is 0-based.  JSign is +/-1.
        function LInd4_I4(CkArg, N1, N2, N3, N4, I, J, K, L, JSign)
            use iso_fortran_env, only: int32
            integer(int32) :: LInd4_I4
            logical(int32), intent(in) :: CkArg
            integer(int32), intent(in) :: N1, N2, N3, N4, I, J, K, L
            integer(int32), intent(out) :: JSign
        end function LInd4_I4
        function LInd4_I8(CkArg, N1, N2, N3, N4, I, J, K, L, JSign)
            use iso_fortran_env, only: int64
            integer(int64) :: LInd4_I8
            logical(int64), intent(in) :: CkArg
            integer(int64), intent(in) :: N1, N2, N3, N4, I, J, K, L
            integer(int64), intent(out) :: JSign
        end function LInd4_I8
    end interface LInd4

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_LInd4
        !! Return the 0-based linear index given 4 1-based indices.
        !!
        !! Returns the linear or square indexing, with indexes `I`, `J`, `K`,
        !! `L` 1-based given the F order.  The output is 0-based.
        !! If `CkArg` is true, the routine checks indices in range and the
        !! function returns `-1` if out of range.
        !! Returns the 0-based index into a linear array given 2 indices and
        !! dimensions.
        !! If `CkArg` is true, the function checks if indices are in range and
        !! returns -1 if they are out of range, otherwise the indices are
        !! assumed to be valid.
        !! `JSign` takes a +/-1 value whether the upper or lower triangle was
        !! selected (i.e., whether to apply a `JSign` flip for anti-symmetric/
        !! Hermitian matrices and/or take a complex conjugate.
        !!
        !! @note "version"
        !! This is the core version.
        !! @endnote
        function QCM_LInd4_I4(N1, N2, N3, N4, I, J, K, L, JSign)
            use iso_fortran_env, only: int32
            integer(int32) :: QCM_LInd4_I4
            integer(int32), intent(in) :: N1, N2, N3, N4, I, J, K, L
            integer(int32), intent(out) :: JSign
        end function QCM_LInd4_I4
        function QCM_LInd4_I8(N1, N2, N3, N4, I, J, K, L, JSign)
            use iso_fortran_env, only: int64
            integer(int64) :: QCM_LInd4_I8
            integer(int64), intent(in) :: N1, N2, N3, N4, I, J, K, L
            integer(int64), intent(out) :: JSign
        end function QCM_LInd4_I8
    end interface QCM_LInd4
#endif

    ! ----------------------------------------------------------------------

    interface LInd5C
        !! Return the 0-based linear index given 5 0-based indices.
        !!
        !! Returns the linear or square indexing, with indexes `I`, `J`, `K`,
        !! `L`, `M` 0-based given the C order.  The output is 0-based.
        !! If `CkArg` is true, the routine checks indices in range and the
        !! function returns `-1` if out of range.
        !! Returns the 0-based index into a linear array given 2 indices and
        !! dimensions.
        !! If `CkArg` is true, the function checks if indices are in range and
        !! returns -1 if they are out of range, otherwise the indices are
        !! assumed to be valid.
        !! `JSign` takes a +/-1 value whether the upper or lower triangle was
        !! selected (i.e., whether to apply a `JSign` flip for anti-symmetric/
        !! Hermitian matrices and/or take a complex conjugate.
        !!
        !! @note "version"
        !! Version compatible with C indexing.
        !! @endnote
        !! Linear or square indexing, I,J,K,L,M are 0-based, c order
        !! output is 0-based.  JSign is +/-1.
        function LInd5C_I4(CkArg, N1, N2, N3, N4, N5, I, J, K, L, M, JSign)
            use iso_fortran_env, only: int32
            integer(int32) :: LInd5C_I4
            logical(int32), intent(in) :: CkArg
            integer(int32), intent(in) :: N1, N2, N3, N4, N5, I, J, K, L, M
            integer(int32), intent(out) :: JSign
        end function LInd5C_I4
        function LInd5C_I8(CkArg, N1, N2, N3, N4, N5, I, J, K, L, M, JSign)
            use iso_fortran_env, only: int64
            integer(int64) :: LInd5C_I8
            logical(int64), intent(in) :: CkArg
            integer(int64), intent(in) :: N1, N2, N3, N4, N5, I, J, K, L, M
            integer(int64), intent(out) :: JSign
        end function LInd5C_I8
    end interface LInd5C

    ! ----------------------------------------------------------------------

    interface LInd5
        !! Return the 0-based linear index given 5 1-based indices.
        !!
        !! Returns the linear or square indexing, with indexes `I`, `J`, `K`,
        !! `L`, `M` 1-based given the F order.  The output is 0-based.
        !! If `CkArg` is true, the routine checks indices in range and the
        !! function returns `-1` if out of range.
        !! Returns the 0-based index into a linear array given 2 indices and
        !! dimensions.
        !! If `CkArg` is true, the function checks if indices are in range and
        !! returns -1 if they are out of range, otherwise the indices are
        !! assumed to be valid.
        !! `JSign` takes a +/-1 value whether the upper or lower triangle was
        !! selected (i.e., whether to apply a `JSign` flip for anti-symmetric/
        !! Hermitian matrices and/or take a complex conjugate.
        !!
        !! @note "version"
        !! Version compatible with Fortran indexing.
        !! @endnote
        function LInd5_I4(CkArg, N1, N2, N3, N4, N5, I, J, K, L, M, JSign)
            use iso_fortran_env, only: int32
            integer(int32) :: LInd5_I4
            logical(int32), intent(in) :: CkArg
            integer(int32), intent(in) :: N1, N2, N3, N4, N5, I, J, K, L, M
            integer(int32), intent(out) :: JSign
        end function LInd5_I4
        function LInd5_I8(CkArg, N1, N2, N3, N4, N5, I, J, K, L, M, JSign)
            use iso_fortran_env, only: int64
            integer(int64) :: LInd5_I8
            logical(int64), intent(in) :: CkArg
            integer(int64), intent(in) :: N1, N2, N3, N4, N5, I, J, K, L, M
            integer(int64), intent(out) :: JSign
        end function LInd5_I8
    end interface LInd5

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_LInd5
        !! Return the 0-based linear index given 5 1-based indices.
        !!
        !! Returns the linear or square indexing, with indexes `I`, `J`, `K`,
        !! `L`, `M` 1-based given the F order.  The output is 0-based.
        !! If `CkArg` is true, the routine checks indices in range and the
        !! function returns `-1` if out of range.
        !! Returns the 0-based index into a linear array given 2 indices and
        !! dimensions.
        !! If `CkArg` is true, the function checks if indices are in range and
        !! returns -1 if they are out of range, otherwise the indices are
        !! assumed to be valid.
        !! `JSign` takes a +/-1 value whether the upper or lower triangle was
        !! selected (i.e., whether to apply a `JSign` flip for anti-symmetric/
        !! Hermitian matrices and/or take a complex conjugate.
        !!
        !! @note "version"
        !! This is the core version.
        !! @endnote
        function QCM_LInd5_I4(N1, N2, N3, N4, N5, I, J, K, L, M, JSign)
            use iso_fortran_env, only: int32
            integer(int32) :: QCM_LInd5_I4
            integer(int32), intent(in) :: N1, N2, N3, N4, N5, I, J, K, L, M
            integer(int32), intent(out) :: JSign
        end function QCM_LInd5_I4
        function QCM_LInd5_I8(N1, N2, N3, N4, N5, I, J, K, L, M, JSign)
            use iso_fortran_env, only: int64
            integer(int64) :: QCM_LInd5_I8
            integer(int64), intent(in) :: N1, N2, N3, N4, N5, I, J, K, L, M
            integer(int64), intent(out) :: JSign
        end function QCM_LInd5_I8
    end interface QCM_LInd5
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_MaxINZR
        !! Return the (1-based) index into the last non-zero column of X,
        !! or 0 if all columns are zero.
        function QCM_MaxINZR_I4(NR, LR, X)
            use iso_fortran_env, only: int32, real64
            integer(int32) :: QCM_MaxINZR_I4
            integer(int32), intent(in) :: NR, LR
            real(real64), intent(in) :: X(NR,LR)
        end function QCM_MaxINZR_I4
        function QCM_MaxINZR_I8(NR, LR, X)
            use iso_fortran_env, only: int64, real64
            integer(int64) :: QCM_MaxINZR_I8
            integer(int64), intent(in) :: NR, LR
            real(real64), intent(in) :: X(NR,LR)
        end function QCM_MaxINZR_I8
    end interface QCM_MaxINZR
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_NumNZR
        !! Return the number of non-zero elements of X(NR,NTot) (Fortran order).
        !! Returns at least 1 so that something will be written even for zero
        !! arrays.
        function QCM_NumNZR_I4(NR, LR, X)
            use iso_fortran_env, only: int32, real64
            integer(int32) :: QCM_NumNZR_I4
            integer(int32), intent(in) :: NR, LR
            real(real64), intent(in) :: X(NR,LR)
        end function QCM_NumNZR_I4
        function QCM_NumNZR_I8(NR, LR, X)
            use iso_fortran_env, only: int64, real64
            integer(int64) :: QCM_NumNZR_I8
            integer(int64), intent(in) :: NR, LR
            real(real64), intent(in) :: X(NR,LR)
        end function QCM_NumNZR_I8
    end interface QCM_NumNZR
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN <= 2
    interface NumNZA
        !! Return the number of non-zero elements of X(NR,NTot) (Fortran order).
        !! The value can be 0.
        function NumNZA_I4(NR, NTot, X)
            use iso_fortran_env, only: int32, real64
            integer(int32) :: NumNZA_I4
            integer(int32), intent(in) :: NR, NTot
            real(real64), intent(in) :: X(NTot,NR)
        end function NumNZA_I4
        function NumNZA_I8(NR, NTot, X)
            use iso_fortran_env, only: int64, real64
            integer(int64) :: NumNZA_I8
            integer(int64), intent(in) :: NR, NTot
            real(real64), intent(in) :: X(NTot,NR)
        end function NumNZA_I8
    end interface NumNZA
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN <= 2
    interface NumNZR
        !! Return the number of non-zero elements of X(NR,NTot) (Fortran order).
        !! The value can be 0.
        function NumNZR_I4(NR, LR, X)
            use iso_fortran_env, only: int32, real64
            integer(int32) :: NumNZR_I4
            integer(int32), intent(in) :: NR, LR
            real(real64), intent(in) :: X(NR,LR)
        end function NumNZR_I4
        function NumNZR_I8(NR, LR, X)
            use iso_fortran_env, only: int64, real64
            integer(int64) :: NumNZR_I8
            integer(int64), intent(in) :: NR, LR
            real(real64), intent(in) :: X(NR,LR)
        end function NumNZR_I8
    end interface NumNZR
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_NumNZT
        !! Return the number of non-zero sets of `NR` elements in `X`.
        !!
        !! Returns the number of non-zero sets of `NR` elements in `X`,
        !! checking for any non-zero value in a set of `NR` values).
        !! The routine assumes that `NR` is the slowest running values.
        function QCM_NumNZT_I4(NR, LR, X)
            use iso_fortran_env, only: int32, real64
            integer(int32) :: QCM_NumNZT_I4
            integer(int32), intent(in) :: NR, LR
            real(real64), intent(in) :: X(LR,NR)
        end function QCM_NumNZT_I4
        function QCM_NumNZT_I8(NR, LR, X)
            use iso_fortran_env, only: int64, real64
            integer(int64) :: QCM_NumNZT_I8
            integer(int64), intent(in) :: NR, LR
            real(real64), intent(in) :: X(LR,NR)
        end function QCM_NumNZT_I8
    end interface QCM_NumNZT
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_NumNZRP
        !! Wrapper for `QCM_NumNZR` to allow `X` to be provided as a flat array
        !! (for Python interface mainly).
        function QCM_NumNZRP_I4(NR, LR, X)
            use iso_fortran_env, only: int32, real64
            integer(int32) :: QCM_NumNZRP_I4
            integer(int32), intent(in) :: NR, LR
            real(real64), intent(in) :: X(NR*LR)
        end function QCM_NumNZRP_I4
        function QCM_NumNZRP_I8(NR, LR, X)
            use iso_fortran_env, only: int64, real64
            integer(int64) :: QCM_NumNZRP_I8
            integer(int64), intent(in) :: NR, LR
            real(real64), intent(in) :: X(NR*LR)
        end function QCM_NumNZRP_I8
    end interface QCM_NumNZRP
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_NumNZTP
        !! Wrapper for `QCM_NumNZT` to allow `X` to be provided as flat array
        !! (for Python interface mainly).
        function QCM_NumNZTP_I4(NR, LR, X)
            use iso_fortran_env, only: int32, real64
            integer(int32) :: QCM_NumNZRP_I4
            integer(int32), intent(in) :: NR, LR
            real(real64), intent(in) :: X(NR*LR)
        end function QCM_NumNZTP_I4
        function QCM_NumNZTP_I8(NR, LR, X)
            use iso_fortran_env, only: int64, real64
            integer(int64) :: QCM_NumNZRP_I8
            integer(int64), intent(in) :: NR, LR
            real(real64), intent(in) :: X(NR*LR)
        end function QCM_NumNZTP_I8
    end interface QCM_NumNZTP
#endif

    ! ----------------------------------------------------------------------

    interface Rd_CBuf
        !! Read a complex array from the file given parameters from the header
        !! record.
        subroutine Rd_CBuf_I4(IU, LR, LenBuf, Arr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, LenBuf, LR
            complex(real64), intent(out) :: Arr(LR)
        end subroutine Rd_CBuf_I4
        subroutine Rd_CBuf_I8(IU, LR, LenBuf, Arr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, LenBuf, LR
            complex(real64), intent(out) :: Arr(LR)
        end subroutine Rd_CBuf_I8
    end interface Rd_CBuf

    ! ----------------------------------------------------------------------

    interface Rd_ChBuf
        !! Read a character string from the file given parameters from the
        !! header record.
#if GAUOPEN >= 3
        subroutine Rd_ChBuf_I4(IU, LR, LenBuf, CArr)
            use iso_fortran_env, only: int32
            integer(int32), intent(in) :: IU, LenBuf, LR
            character(len=*), intent(inout) :: CArr
        end subroutine Rd_ChBuf_I4
        subroutine Rd_ChBuf_I8(IU, LR, LenBuf, CArr)
            use iso_fortran_env, only: int64
            integer(int64), intent(in) :: IU, LenBuf, LR
            character(len=*), intent(inout) :: CArr
        end subroutine Rd_ChBuf_I8
#else
        module procedure Rd_ChBuf_Dum_I4, Rd_ChBuf_Dum_I8
#endif
    end interface Rd_ChBuf

    ! ----------------------------------------------------------------------

    interface Rd_ChBufI
        !! Read a character string from the file given parameters from the
        !! header record.
        !! The routine also returns the string as an array of character
        !! indexes in the default collating sequence.
#if GAUOPEN >= 3
        subroutine Rd_ChBufI_I4(IU, LR, LenBuf, CArr, IArr)
            use iso_fortran_env, only: int32
            integer(int32), intent(in) :: IU, LenBuf, LR
            character(len=*), intent(inout) :: CArr
            integer(int32), intent(out) :: IArr(LR)
        end subroutine Rd_ChBufI_I4
        subroutine Rd_ChBufI_I8(IU, LR, LenBuf, CArr, IArr)
            use iso_fortran_env, only: int64
            integer(int64), intent(in) :: IU, LenBuf, LR
            character(len=*), intent(inout) :: CArr
            integer(int64), intent(out) :: IArr(LR)
        end subroutine Rd_ChBufI_I8
#else
        module procedure Rd_ChBufI_Dum_I4, Rd_ChBufI_Dum_I8
#endif
end interface Rd_ChBufI

    ! ----------------------------------------------------------------------

    interface Rd_RBuf
        !! Read a real array from the file given parameters from the header
        !! record.
        subroutine Rd_RBuf_I4(IU, LR, LenBuf, Arr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, LenBuf, LR
            real(real64), intent(out) :: Arr(LR)
        end subroutine Rd_RBuf_I4
        subroutine Rd_RBuf_I8(IU, LR, LenBuf, Arr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, LenBuf, LR
            real(real64), intent(out) :: Arr(LR)
        end subroutine Rd_RBuf_I8
    end interface Rd_RBuf

    ! ----------------------------------------------------------------------

    interface Rd_Skip
        !! Skip the data records for an object on the file having `NTot`
        !! elements stored with `LenBuf` per record.
        subroutine Rd_Skip_I4(IU, NTot, LenBuf)
            use iso_fortran_env, only: int32
            integer(int32), intent(in) :: IU, NTot, LenBuf
        end subroutine Rd_Skip_I4
        subroutine Rd_Skip_I8(IU, NTot, LenBuf)
            use iso_fortran_env, only: int64
            integer(int64), intent(in) :: IU, NTot, LenBuf
        end subroutine Rd_Skip_I8
    end interface Rd_Skip

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface Wr_2E
        !! Write AO 2-electron integrals to Fortran unit `IU`.
        subroutine Wr_2E_I4(IU, NTot, NR, N, LR, LenBuf, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, NTot, NR, N, LR, LenBuf
            real(real64), intent(in) :: RArr(NR,LR)
        end subroutine Wr_2E_I4
        subroutine Wr_2E_I8(IU, NTot, NR, N, LR, LenBuf, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, NTot, NR, N, LR, LenBuf
            real(real64), intent(in) :: RArr(NR,LR)
        end subroutine Wr_2E_I8
    end interface Wr_2E
#else
    interface Wr_2E
        !! Write AO 2-electron integrals to Fortran unit `IU`.
        subroutine Wr_2E_I4(IU, NTot, NR, N, LR, LenBuf, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, NTot, NR, N, LR, LenBuf
            real(real64), intent(in) :: RArr(LR,NR)
        end subroutine Wr_2E_I4
        subroutine Wr_2E_I8(IU, NTot, NR, N, LR, LenBuf, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, NTot, NR, N, LR, LenBuf
            real(real64), intent(in) :: RArr(LR,NR)
        end subroutine Wr_2E_I8
    end interface Wr_2E
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_Wr_2E
        !! Write AO 2-electron integrals to Fortran unit `IU`.
        !! The array is expected as one-dimensional.
        subroutine QCM_Wr_2E_I4(IU, NTot, NR, N, LR, LenBuf, LRNR, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, NTot, NR, N, LR, LenBuf
            real(real64), intent(in) :: RArr(LR,NR)
        end subroutine QCM_Wr_2E_I4
        subroutine QCM_Wr_2E_I8(IU, NTot, NR, N, LR, LenBuf, LRNR, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, NTot, NR, N, LR, LenBuf
            real(real64), intent(in) :: RArr(LR,NR)
        end subroutine QCM_Wr_2E_I8
    end interface QCM_Wr_2E
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface Wr_SpA
        !! Write sparse matrix `RArr` to Fortran unit `IU`.
        subroutine Wr_SpA_I4(IU, NI, NR, NTot, LenBuf, IArr, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, NI, NR, NTot, LenBuf
            integer(int32), intent(in) :: IArr(NI*NTot)
            real(real64), intent(in) :: RArr(NR*NTot)
        end subroutine Wr_SpA_I4
        subroutine Wr_SpA_I8(IU, NI, NR, NTot, LenBuf, IArr, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, NI, NR, NTot, LenBuf
            integer(int64), intent(in) :: IArr(NI*NTot)
            real(real64), intent(in) :: RArr(NR*NTot)
        end subroutine Wr_SpA_I8
    end interface Wr_SpA
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface Wr_ChBuf
        !! Write character string to Fortran unit `IU`.
        subroutine Wr_ChBuf_I4(IU, LR, LenBuf, CArr)
            use iso_fortran_env, only: int32
            integer(int32), intent(in) :: IU, LR, LenBuf
            character(len=*), intent(in) :: CArr
        end subroutine Wr_ChBuf_I4
        subroutine Wr_ChBuf_I8(IU, LR, LenBuf, CArr)
            use iso_fortran_env, only: int64
            integer(int64), intent(in) :: IU, LR, LenBuf
            character(len=*), intent(in) :: CArr
        end subroutine Wr_ChBuf_I8
    end interface Wr_ChBuf
#endif

    ! ----------------------------------------------------------------------

    interface Wr_CBuf
        !! Write complex array to Fortran unit `IU`.
        subroutine Wr_CBuf_I4(IU, LR, LenBuf, Arr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, LR, LenBuf
            complex(real64), intent(in) :: Arr
        end subroutine Wr_CBuf_I4
        subroutine Wr_CBuf_I8(IU, LR, LenBuf, Arr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, LR, LenBuf
            complex(real64), intent(in) :: Arr
        end subroutine Wr_CBuf_I8
    end interface Wr_CBuf

    ! ----------------------------------------------------------------------

    interface Wr_Head
        !! Write the header records (3 to NLab) to Fortran unit `IU`.
        subroutine Wr_Head_I4(IU, NAtoms, NAt3, NBasis, IAn, IAtTyp, AtmChg, &
                              C, IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, &
                              IDum9, NShlAO, NPrmAO, NShlDB, NPrmDB, NBTot)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, NAtoms, NAt3, NBasis, NFC, NFV, &
                ITran, IDum9, NShlAO, NPrmAO, NShlDB, NPrmDB, NBTot
            integer(int32), intent(in) :: IAn(NAtoms), IAtTyp(NAtoms), &
                IBfAtm(NBasis), IBfTyp(NBasis)
            real(real64), intent(in) :: AtmChg(NAtoms), C(NAt3), AtmWgt(NAtoms)
        end subroutine Wr_Head_I4
        subroutine Wr_Head_I8(IU, NAtoms, NAt3, NBasis, IAn, IAtTyp, AtmChg, &
                              C, IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, &
                              IDum9, NShlAO, NPrmAO, NShlDB, NPrmDB, NBTot)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, NAtoms, NAt3, NBasis, NFC, NFV, &
                ITran, IDum9, NShlAO, NPrmAO, NShlDB, NPrmDB, NBTot
            integer(int64), intent(in) :: IAn(NAtoms), IAtTyp(NAtoms), &
                IBfAtm(NBasis), IBfTyp(NBasis)
            real(real64), intent(in) :: AtmChg(NAtoms), C(NAt3), AtmWgt(NAtoms)
        end subroutine Wr_Head_I8
    end interface Wr_Head

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface Wr_HeadA
        !! Write the header records (3 to NLab) to Fortran unit `IU`.
        !! The routine expects the full 16 integer array for record 11.
        subroutine Wr_HeadA_I4(IU, NAtoms, NAt3, NBasis, IAn, IAtTyp, AtmChg, &
                               C, IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, &
                               IDum9, Rec11)
            use iso_fortran_env, only: int32, real64
            integer(int32), parameter :: LRec11=16
            integer(int32), intent(in) :: IU, NAtoms, NAt3, NBasis, NFC, NFV, &
                ITran, IDum9
            integer(int32), intent(in) :: IAn(NAtoms), IAtTyp(NAtoms), &
                IBfAtm(NBasis), IBfTyp(NBasis), Rec11(LRec11)
            real(real64), intent(in) :: AtmChg(NAtoms), C(NAt3), AtmWgt(NAtoms)
        end subroutine Wr_HeadA_I4
        subroutine Wr_HeadA_I8(IU, NAtoms, NAt3, NBasis, IAn, IAtTyp, AtmChg, &
                              C, IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, &
                              IDum9, Rec11)
            use iso_fortran_env, only: int64, real64
            integer(int64), parameter :: LRec11=16
            integer(int64), intent(in) :: IU, NAtoms, NAt3, NBasis, NFC, NFV, &
                ITran, IDum9
            integer(int64), intent(in) :: IAn(NAtoms), IAtTyp(NAtoms), &
                IBfAtm(NBasis), IBfTyp(NBasis), Rec11(LRec11)
            real(real64), intent(in) :: AtmChg(NAtoms), C(NAt3), AtmWgt(NAtoms)
        end subroutine Wr_HeadA_I8
    end interface Wr_HeadA
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface Wr_DummyHead
        !! Write header with dummy information to Fortran unit `IU`.
        !! This routine is intended for files which are just used for arrays
        !! storage and do not need a molecule specification.
        subroutine Wr_DummyHead_I4(IU)
            use iso_fortran_env, only: int32
            integer(int32), intent(in) :: IU
        end subroutine Wr_DummyHead_I4
        subroutine Wr_DummyHead_I8(IU)
            use iso_fortran_env, only: int64
            integer(int64), intent(in) :: IU
        end subroutine Wr_DummyHead_I8
    end interface Wr_DummyHead
#endif

    ! ----------------------------------------------------------------------

    interface Wr_IBuf
        !! Write integer array to Fortran unit `IU`.
        subroutine Wr_IBuf_I4(IU, LR, LenBuf, Arr)
            use iso_fortran_env, only: int32
            integer(int32), intent(in) :: IU, LR, LenBuf
            integer(int32), intent(in) :: Arr(LR)
        end subroutine Wr_IBuf_I4
        subroutine Wr_IBuf_I8(IU, LR, LenBuf, Arr)
            use iso_fortran_env, only: int64
            integer(int64), intent(in) :: IU, LR, LenBuf
            integer(int64), intent(in) :: Arr(LR)
        end subroutine Wr_IBuf_I8
    end interface Wr_IBuf

    ! ----------------------------------------------------------------------

    interface Wr_Labl
        !! Write the header record for one matrix to Fortran unit `IU`.
        subroutine Wr_Labl_I4(IU, CBuf, NI, NR, NTot, LenBuf, N1, N2, N3, N4, &
                              N5, TypeA)
            use iso_fortran_env, only: int32
            integer(int32), intent(in) :: IU, NI, NR, NTot, LenBuf, N1, N2, &
                N3, N4, N5, TypeA
            character(len=*), intent(in) :: CBuf
        end subroutine Wr_Labl_I4
        subroutine Wr_Labl_I8(IU, CBuf, NI, NR, NTot, LenBuf, N1, N2, N3, N4, &
                              N5, TypeA)
            use iso_fortran_env, only: int64
            integer(int64), intent(in) :: IU, NI, NR, NTot, LenBuf, N1, N2, &
                N3, N4, N5, TypeA
            character(len=*), intent(in) :: CBuf
        end subroutine Wr_Labl_I8
    end interface Wr_Labl

    ! ----------------------------------------------------------------------

    interface Wr_RBuf
        !! Write real array to Fortran unit `IU`.
        subroutine Wr_RBuf_I4(IU, LR, LenBuf, Arr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, LR, LenBuf
            real(real64), intent(in) :: Arr(LR)
        end subroutine Wr_RBuf_I4
        subroutine Wr_RBuf_I8(IU, LR, LenBuf, Arr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, LR, LenBuf
            real(real64), intent(in) :: Arr(LR)
        end subroutine Wr_RBuf_I8
    end interface Wr_RBuf

    ! ----------------------------------------------------------------------

    interface Wr_RInd
        !! Write non-zero elements as sparse matrix in Fortran unit `IU`.
        !! The indexes of non-zero elements in full array are also stored.
        subroutine Wr_RInd_I4(IU, NR, LR, NTot, LenBuf, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, NR, LR, NTot, LenBuf
            real(real64), intent(in) :: RArr(NR,LR)
        end subroutine Wr_RInd_I4
        subroutine Wr_RInd_I8(IU, NR, LR, NTot, LenBuf, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, NR, LR, NTot, LenBuf
            real(real64), intent(in) :: RArr(NR,LR)
        end subroutine Wr_RInd_I8
    end interface Wr_RInd


    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_Wr_RInd
        !! Write non-zero elements as sparse matrix in Fortran unit `IU`.
        !! The indexes of non-zero elements in full array are also stored.
        !! In this version, the array is one-dimensional.
        subroutine QCM_Wr_RInd_I4(IU, NR, LR, NTot, LenBuf, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IU, NR, LR, NTot, LenBuf
            real(real64) :: RArr(NR*LR)
        end subroutine QCM_Wr_RInd_I4
        subroutine QCM_Wr_RInd_I8(IU, NR, LR, NTot, LenBuf, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IU, NR, LR, NTot, LenBuf
            real(real64) :: RArr(NR*LR)
        end subroutine QCM_Wr_RInd_I8
    end interface QCM_Wr_RInd
#endif

    ! ----------------------------------------------------------------------

    interface ExpAO1
        !! Expand AO matrix stored in compressed form to full 4D shape.
        subroutine ExpAO1_I4(N, LR, RI, RO)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: N, LR
            real(real64), intent(in) :: RI(LR)
            real(real64), intent(out) :: RO(N,N,N,N)
        end subroutine ExpAO1_I4
        subroutine ExpAO1_I8(N, LR, RI, RO)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: N, LR
            real(real64), intent(in) :: RI(LR)
            real(real64), intent(out) :: RO(N,N,N,N)
        end subroutine ExpAO1_I8
    end interface ExpAO1

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface ExpAO1C
        !! Expand AO matrix stored in compressed form to full 4D shape.
        subroutine ExpAO1C_I4(N, RI, RO)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: N
            real(real64), intent(in) :: RI(2,*)
            real(real64), intent(out) :: RO(2,N,N,N,N)
        end subroutine ExpAO1C_I4
        subroutine ExpAO1C_I8(N, RI, RO)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: N
            real(real64), intent(in) :: RI(2,*)
            real(real64), intent(out) :: RO(2,N,N,N,N)
        end subroutine ExpAO1C_I8
    end interface ExpAO1C
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface ExpAO2
        !! Expand AO matrix stored in compressed form to full 4D shape.
        subroutine ExpAO2_I4(N, LR, RI, RO)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: N
            real(real64), intent(in) :: RI(2,LR)
            real(real64), intent(out) :: RO(2,N,N,N,N)
        end subroutine ExpAO2_I4
        subroutine ExpAO2_I8(N, LR, RI, RO)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: N
            real(real64), intent(in) :: RI(2,LR)
            real(real64), intent(out) :: RO(2,N,N,N,N)
        end subroutine ExpAO2_I8
    end interface ExpAO2
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface ExpAO3
        !! Expand AO matrix stored in compressed form to full 4D shape.
        subroutine ExpAO3_I4(N, LR, RI, RO)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: N
            real(real64), intent(in) :: RI(3,LR)
            real(real64), intent(out) :: RO(3,N,N,N,N)
        end subroutine ExpAO3_I4
        subroutine ExpAO3_I8(N, LR, RI, RO)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: N
            real(real64), intent(in) :: RI(3,LR)
            real(real64), intent(out) :: RO(3,N,N,N,N)
        end subroutine ExpAO3_I8
    end interface ExpAO3
#endif

    ! ----------------------------------------------------------------------

    interface ExpAON
        !! Expand AO matrix with `NE` elements per index, stored in compressed
        !! form to full 5D shape.
        subroutine ExpAON_I4(NE, N, LR, RI, RO)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: NE, N, LR
            real(real64), intent(in) :: RI(NE,LR)
            real(real64), intent(out) :: RO(NE,N,N,N,N)
        end subroutine ExpAON_I4
        subroutine ExpAON_I8(NE, N, LR, RI, RO)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: NE, N, LR
            real(real64), intent(in) :: RI(NE,LR)
            real(real64), intent(out) :: RO(NE,N,N,N,N)
        end subroutine ExpAON_I8
    end interface ExpAON

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_ExpAO1C
        !! Expand AO matrix with `NE` elements per index, stored in compressed
        !! form to full 5D shape.
        subroutine QCM_ExpAO1C_I4(N, LR, RI, RO)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: N, LR
            complex(real64), intent(in) :: RI(LR)
            complex(real64), intent(out) :: RO(N,N,N,N)
        end subroutine QCM_ExpAO1C_I4
        subroutine QCM_ExpAO1C_I8(N, LR, RI, RO)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: N, LR
            complex(real64), intent(in) :: RI(LR)
            complex(real64), intent(out) :: RO(N,N,N,N)
        end subroutine QCM_ExpAO1C_I8
    end interface QCM_ExpAO1C
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_ExpAO2
        !! Expand AO matrix with `NE` elements per index, stored in compressed
        !! form to full 5D shape.
        !! This routine returns the final array as one-dimensional.
        subroutine QCM_ExpAO2_I4(N, LR, NELR, NOut, RI, RO)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: N, LR, NELR, NOut
            real(real64), intent(in) :: RI(NELR)
            real(real64), intent(out) :: RO(NOut)
        end subroutine QCM_ExpAO2_I4
        subroutine QCM_ExpAO2_I8(N, LR, NELR, NOut, RI, RO)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: N, LR, NELR, NOut
            real(real64), intent(in) :: RI(NELR)
            real(real64), intent(out) :: RO(NOut)
        end subroutine QCM_ExpAO2_I8
    end interface QCM_ExpAO2
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_ExpAO3
        !! Expand AO matrix with `NE` elements per index, stored in compressed
        !! form to full 5D shape.
        !! This routine returns the final array as one-dimensional.
        subroutine QCM_ExpAO3_I4(N, LR, NELR, NOut, RI, RO)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: N, LR, NELR, NOut
            real(real64), intent(in) :: RI(NELR)
            real(real64), intent(out) :: RO(NOut)
        end subroutine QCM_ExpAO3_I4
        subroutine QCM_ExpAO3_I8(N, LR, NELR, NOut, RI, RO)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: N, LR, NELR, NOut
            real(real64), intent(in) :: RI(NELR)
            real(real64), intent(out) :: RO(NOut)
        end subroutine QCM_ExpAO3_I8
    end interface QCM_ExpAO3
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_ExpAON
        !! Expand AO matrix with `NE` elements per index, stored in compressed
        !! form to full 5D shape.
        !! This routine returns the final array as one-dimensional.
        subroutine QCM_ExpAON_I4(NE, N, LR, NELR, NOut, RI, RO)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: NE, N, LR, NELR, NOut
            real(real64), intent(in) :: RI(NELR)
            real(real64), intent(out) :: RO(NOut)
        end subroutine QCM_ExpAON_I4
        subroutine QCM_ExpAON_I8(NE, N, LR, NELR, NOut, RI, RO)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: NE, N, LR, NELR, NOut
            real(real64), intent(in) :: RI(NELR)
            real(real64), intent(out) :: RO(NOut)
        end subroutine QCM_ExpAON_I8
    end interface QCM_ExpAON
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface UnPack2E
        subroutine UnPack2E_I4(NR, N4, NTot, AI, AV, AO)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: NR, N4, NTot, AI(4,NTot)
            real(real64), intent(in) :: AV(NR,NTot)
            real(real64), intent(out) :: AO(NR,N4)
        end subroutine UnPack2E_I4
        subroutine UnPack2E_I8(NR, N4, NTot, AI, AV, AO)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: NR, N4, NTot, AI(4,NTot)
            real(real64), intent(in) :: AV(NR,NTot)
            real(real64), intent(out) :: AO(NR,N4)
        end subroutine UnPack2E_I8
    end interface UnPack2E
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface QCM_UnPack2E
        subroutine QCM_UnPack2E_I4(NR, N4, NTot, AI, AV, AO)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: NR, N4, NTot, AI(4,NTot)
            real(real64), intent(in) :: AV(NR,NTot)
            real(real64), intent(out) :: AO(NR,N4)
        end subroutine QCM_UnPack2E_I4
        subroutine QCM_UnPack2E_I8(NR, N4, NTot, AI, AV, AO)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: NR, N4, NTot, AI(4,NTot)
            real(real64), intent(in) :: AV(NR,NTot)
            real(real64), intent(out) :: AO(NR,N4)
        end subroutine QCM_UnPack2E_I8
    end interface QCM_UnPack2E
#endif

    ! ----------------------------------------------------------------------

    interface AClear
        !! Clear N first elements in array A.
        subroutine AClear_I4(N, A)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: N
            real(real64), intent(out) :: A(N)
        end subroutine AClear_I4
        subroutine AClear_I8(N, A)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: N
            real(real64), intent(out) :: A(N)
        end subroutine AClear_I8
    end interface AClear

    ! ----------------------------------------------------------------------

    interface IClear
        !! Clear N first elements in array A.
        subroutine IClear_I4(N, IA)
            use iso_fortran_env, only: int32
            integer(int32), intent(in) :: N
            integer(int32), intent(out) :: IA(N)
        end subroutine IClear_I4
        subroutine IClear_I8(N, IA)
            use iso_fortran_env, only: int64
            integer(int64), intent(in) :: N
            integer(int64), intent(out) :: IA(N)
        end subroutine IClear_I8
    end interface IClear

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface Clean_Label
        !! Canonicalize the label in CBuf by removing leading and trailing
        !! spaces, collapsing multiple spaces, and captializing.
        subroutine Clean_Label_I4(CBuf)
            character(len=*), intent(inout) :: CBuf
        end subroutine Clean_Label_I4
    end interface Clean_Label
#endif

    ! ----------------------------------------------------------------------

    ! =======================================
    !  Interfaces to Other QCMATRIX Routines
    ! =======================================

    ! ----------------------------------------------------------------------

    interface Wr_LIBuf
        !! Write integer array `IX` unblocked (zeros included) to unit `IU`
        !! given the parameters from the header record.
        subroutine Wr_LIBuf_I4(IU, Label, NI, LenBuf, N1, N2, N3, N4, N5, &
                               TypeA, IX)
            use iso_fortran_env, only: int32
            character(len=*), intent(in) :: Label
            integer(int32), intent(in) :: IU, NI, LenBuf, N1, N2, N3, N4, N5, &
                TypeA
            integer(int32), intent(in) :: IX(*)
        end subroutine Wr_LIBuf_I4
        subroutine Wr_LIBuf_I8(IU, Label, NI, LenBuf, N1, N2, N3, N4, N5, &
                               TypeA, IX)
            use iso_fortran_env, only: int64
            character(len=*), intent(in) :: Label
            integer(int64), intent(in) :: IU, NI, LenBuf, N1, N2, N3, N4, N5, &
                TypeA
            integer(int64), intent(in) :: IX(*)
        end subroutine Wr_LIBuf_I8
    end interface Wr_LIBuf

    ! ----------------------------------------------------------------------

    interface Wr_LRBuf
        !! Write real array `X` unblocked (zeros included) to unit `IU`
        !! given the parameters from the header record.
        subroutine Wr_LRBuf_I4(IU, Label, NR, LenBuf, N1, N2, N3, N4, N5, &
                               TypeA, X)
            use iso_fortran_env, only: int32, real64
            character(len=*), intent(in) :: Label
            integer(int32), intent(in) :: IU, NR, LenBuf, N1, N2, N3, N4, N5, &
                TypeA
            real(real64), intent(in) :: X(*)
        end subroutine Wr_LRBuf_I4
        subroutine Wr_LRBuf_I8(IU, Label, NR, LenBuf, N1, N2, N3, N4, N5, &
                               TypeA, X)
            use iso_fortran_env, only: int64, real64
            character(len=*), intent(in) :: Label
            integer(int64), intent(in) :: IU, NR, LenBuf, N1, N2, N3, N4, N5, &
                TypeA
            real(real64), intent(in) :: X(*)
        end subroutine Wr_LRBuf_I8
    end interface Wr_LRBuf

    ! ----------------------------------------------------------------------

    interface Wr_LCBuf
        !! Write complex array `X` unblocked (zeros included) to unit `IU`
        !! given the parameters from the header record.
        subroutine Wr_LCBuf_I4(IU, Label, NR, LenBuf, N1, N2, N3, N4, N5, &
                               TypeA, X)
            use iso_fortran_env, only: int32, real64
            character(len=*), intent(in) :: Label
            integer(int32), intent(in) :: IU, NR, LenBuf, N1, N2, N3, N4, N5, &
                TypeA
            complex(real64), intent(in) :: X(*)
        end subroutine Wr_LCBuf_I4
        subroutine Wr_LCBuf_I8(IU, Label, NR, LenBuf, N1, N2, N3, N4, N5, &
                               TypeA, X)
            use iso_fortran_env, only: int64, real64
            character(len=*), intent(in) :: Label
            integer(int64), intent(in) :: IU, NR, LenBuf, N1, N2, N3, N4, N5, &
                TypeA
            complex(real64), intent(in) :: X(*)
        end subroutine Wr_LCBuf_I8
    end interface Wr_LCBuf

    ! ----------------------------------------------------------------------

    interface Wr_LRInd
        !! Write a real array compressed, with indices for non-zero
        !! sets of (NR) elements.
        subroutine Wr_LRInd_I4(IU, Label, NR, LenBuf, N1, N2, N3, N4, N5, &
                               TypeA, X)
            use iso_fortran_env, only: int32, real64
            character(len=*), intent(in) :: Label
            integer(int32), intent(in) :: IU, NR, LenBuf, N1, N2, N3, N4, N5, &
                TypeA
            real(real64), intent(in) :: X(NR,*)
        end subroutine Wr_LRInd_I4
        subroutine Wr_LRInd_I8(IU, Label, NR, LenBuf, N1, N2, N3, N4, N5, &
                               TypeA, X)
            use iso_fortran_env, only: int64, real64
            character(len=*), intent(in) :: Label
            integer(int64), intent(in) :: IU, NR, LenBuf, N1, N2, N3, N4, N5, &
                TypeA
            real(real64), intent(in) :: X(NR,*)
        end subroutine Wr_LRInd_I8
    end interface Wr_LRInd

    ! ----------------------------------------------------------------------

    interface Wr_LAO2E
        subroutine Wr_LAO2E_I4(IU, Label, NR, LenBuf, N, RInt)
            use iso_fortran_env, only: int32, real64
            character(len=*), intent(in) :: Label
            integer(int32), intent(in) :: IU, NR, LenBuf, N
            real(real64), intent(in) :: RInt(*)
        end subroutine Wr_LAO2E_I4
        subroutine Wr_LAO2E_I8(IU, Label, NR, LenBuf, N, RInt)
            use iso_fortran_env, only: int64, real64
            character(len=*), intent(in) :: Label
            integer(int64), intent(in) :: IU, NR, LenBuf, N
            real(real64), intent(in) :: RInt(*)
        end subroutine Wr_LAO2E_I8
    end interface Wr_LAO2E

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface ExpD2E
        !! Expand D2E in compressed form to full dimensions.
        subroutine ExpD2E_I4(NAt, N, NTot, RD2E, ID2E, D2E)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: NAt, N, NTot
            integer(int32), intent(in) :: ID2E(5,NTot)
            real(real64), intent(in) :: RD2E(3,NTot)
            real(real64), intent(out) :: D2E(N,N,N,N,3,NAt)
        end subroutine ExpD2E_I4
        subroutine ExpD2E_I8(NAt, N, NTot, RD2E, ID2E, D2E)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: NAt, N, NTot
            integer(int64), intent(in) :: ID2E(5,NTot)
            real(real64), intent(in) :: RD2E(3,NTot)
            real(real64), intent(out) :: D2E(N,N,N,N,3,NAt)
        end subroutine ExpD2E_I8
    end interface ExpD2E
#endif

    ! ----------------------------------------------------------------------

#if GAUOPEN >= 3
    interface Pr_SpA
        !! Print a pair of matrices with arbitrary dimensions.
        subroutine Pr_SpA_I4(IOut ,NI, NRI, NR, NTot, IArr, RArr)
            use iso_fortran_env, only: int32, real64
            integer(int32), intent(in) :: IOut ,NI, NRI, NR, NTot
            integer(int32), intent(in) :: IArr(NI,NTot)
            real(real64), intent(in) :: RArr(NRI,NR,NTot)
        end subroutine Pr_SpA_I4
        subroutine Pr_SpA_I8(IOut ,NI, NRI, NR, NTot, IArr, RArr)
            use iso_fortran_env, only: int64, real64
            integer(int64), intent(in) :: IOut ,NI, NRI, NR, NTot
            integer(int64), intent(in) :: IArr(NI,NTot)
            real(real64), intent(in) :: RArr(NRI,NR,NTot)
        end subroutine Pr_SpA_I8
    end interface Pr_SpA
#endif

    ! ----------------------------------------------------------------------

contains

! ======================================================================

! =========================================================
!  Interfaces Routines to Label-Dependent GauOpen Routines
! =========================================================

! ======================================================================

subroutine Open_Read_32(Name, IU, LabFil, IVers, NLab, GVers, Title, &
                        NAtoms, NBasis, NBsUse, ICharg, Multip, NE, &
                        Len12L, Len4L, IOpCl, ICGU, LabVer)
    !! 32-bit version of `Open_Read` from `qcmatrixio`.
    !!
    !! This version can process files written with 32 and 64-bits labels.
    character(len=*), intent(in) :: Name
    character(len=64), intent(out) :: LabFil, GVers, Title
    integer(int32), intent(out) :: IU, IVers, NLab, NAtoms, NBasis, &
        NBsUse, ICharg, Multip, NE, Len12L, Len4L, IOpCl, ICGU
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IU64, IVers64, NLab64, NAtoms64, NBasis64, &
        NBsUse64, ICharg64, Multip64, NE64, Len12L64, Len4L64, IOpCl64, ICGU64
    integer(int32) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

    select case (isize)
    case (32)
        call Open_Read_GOpen(Name, IU, LabFil, IVers, NLab, GVers, Title, &
                             NAtoms, NBasis, NBsUse, ICharg, Multip, NE, &
                             Len12L, Len4L, IOpCl, ICGU)
    case (64)
        call Open_Read_GOpen(Name, IU64, LabFil, IVers64, NLab64, GVers, &
                             Title, NAtoms64, NBasis64, NBsUse64, ICharg64, &
                             Multip64, NE64, Len12L64, Len4L64, IOpCl64, &
                             ICGU64)
        IU = int(IU64, kind=int32)
        IVers = int(IVers64, kind=int32)
        NLab = int(NLab64, kind=int32)
        NAtoms = int(NAtoms64, kind=int32)
        NBasis = int(NBasis64, kind=int32)
        NBsUse = int(NBsUse64, kind=int32)
        ICharg = int(ICharg64, kind=int32)
        Multip = int(Multip64, kind=int32)
        NE = int(NE64, kind=int32)
        Len12L = int(Len12L64, kind=int32)
        Len4L = int(Len4L64, kind=int32)
        IOpCl = int(IOpCl64, kind=int32)
        ICGU = int(ICGU64, kind=int32)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Open_Read')
    end select

end subroutine Open_Read_32

! ======================================================================

subroutine Open_Read_64(Name, IU, LabFil, IVers, NLab, GVers, Title, &
                        NAtoms, NBasis, NBsUse, ICharg, Multip, NE, &
                        Len12L, Len4L, IOpCl, ICGU, LabVer)
    !! 64-bit version of `Open_Read` from `qcmatrixio`.
    !!
    !! This version can process files written with 32 and 64-bits labels.
    character(len=*), intent(in) :: Name
    character(len=64), intent(out) :: LabFil, GVers, Title
    integer(int64), intent(out) :: IU, IVers, NLab, NAtoms, NBasis, &
        NBsUse, ICharg, Multip, NE, Len12L, Len4L, IOpCl, ICGU
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IU32, IVers32, NLab32, NAtoms32, NBasis32, &
        NBsUse32, ICharg32, Multip32, NE32, Len12L32, Len4L32, IOpCl32, ICGU32
    integer(int64) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

    select case (isize)
    case (64)
        call Open_Read_GOpen(Name, IU, LabFil, IVers, NLab, GVers, Title, &
                             NAtoms, NBasis, NBsUse, ICharg, Multip, NE, &
                             Len12L, Len4L, IOpCl, ICGU)
    case (32)
        call Open_Read_GOpen(Name, IU32, LabFil, IVers32, NLab32, GVers, &
                             Title, NAtoms32, NBasis32, NBsUse32, ICharg32, &
                             Multip32, NE32, Len12L32, Len4L32, IOpCl32, ICGU32)
        IU = int(IU32, kind=int64)
        IVers = int(IVers32, kind=int64)
        NLab = int(NLab32, kind=int64)
        NAtoms = int(NAtoms32, kind=int64)
        NBasis = int(NBasis32, kind=int64)
        NBsUse = int(NBsUse32, kind=int64)
        ICharg = int(ICharg32, kind=int64)
        Multip = int(Multip32, kind=int64)
        NE = int(NE32, kind=int64)
        Len12L = int(Len12L32, kind=int64)
        Len4L = int(Len4L32, kind=int64)
        IOpCl = int(IOpCl32, kind=int64)
        ICGU = int(ICGU32, kind=int64)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Open_Read')
    end select

end subroutine Open_Read_64

! ======================================================================

subroutine Open_Write_32(Name, IU, LabFil, GVers, Title, NAtoms, NBasis, &
                         NBsUse, ICharg, Multip, NE, IOpCl, ICGU, LabVer)
    !! 32-bit version of `Open_Write` from `qcmatrixio`.
    !!
    !! This version can write files using 32 and 64-bits labels.
    character(len=*), intent(in) :: Name, LabFil, GVers, Title
    integer(int32), intent(in) :: NAtoms, NBasis, NBsUse, ICharg, &
        Multip, NE, IOpCl, ICGU
    integer(int32), intent(out) :: IU
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IU64, NAtoms64, NBasis64, NBsUse64, ICharg64, Multip64, &
        NE64, IOpCl64, ICGU64
    integer(int32) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

    select case (isize)
    case (32)
        call Open_Write_GOpen(Name, IU, LabFil, GVers, Title, NAtoms, NBasis, &
                              NBsUse, ICharg, Multip, NE, IOpCl, ICGU)
    case (64)
        NAtoms64 = int(NAtoms, kind=int64)
        NBasis64 = int(NBasis, kind=int64)
        NBsUse64 = int(NBsUse, kind=int64)
        ICharg64 = int(ICharg, kind=int64)
        Multip64 = int(Multip, kind=int64)
        NE64 = int(NE, kind=int64)
        IOpCl64 = int(IOpCl, kind=int64)
        ICGU64 = int(ICGU, kind=int64)
        call Open_Write_GOpen(Name, IU64, LabFil, GVers, Title, NAtoms64, &
                              NBasis64, NBsUse64, ICharg64, Multip64, NE64, &
                              IOpCl64, ICGU64)
        IU = int(IU64, kind=int32)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Open_Write')
    end select

end subroutine Open_Write_32

! ======================================================================

subroutine Open_Write_64(Name, IU, LabFil, GVers, Title, NAtoms, NBasis, &
                         NBsUse, ICharg, Multip, NE, IOpCl, ICGU, LabVer)
    !! 64-bit version of `Open_Write` from `qcmatrixio`.
    !!
    !! This version can write files using 32 and 64-bits labels.
    character(len=*), intent(in) :: Name, LabFil, GVers, Title
    integer(int64), intent(in) :: NAtoms, NBasis, NBsUse, ICharg, &
        Multip, NE, IOpCl, ICGU
    integer(int64), intent(out) :: IU
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IU32, NAtoms32, NBasis32, NBsUse32, ICharg32, Multip32, &
        NE32, IOpCl32, ICGU32
    integer(int64) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

    select case (isize)
    case (64)
        call Open_Write_GOpen(Name, IU, LabFil, GVers, Title, NAtoms, NBasis, &
                              NBsUse, ICharg, Multip, NE, IOpCl, ICGU)
    case (32)
        NAtoms32 = int(NAtoms, kind=int32)
        NBasis32 = int(NBasis, kind=int32)
        NBsUse32 = int(NBsUse, kind=int32)
        ICharg32 = int(ICharg, kind=int32)
        Multip32 = int(Multip, kind=int32)
        NE32 = int(NE, kind=int32)
        IOpCl32 = int(IOpCl, kind=int32)
        ICGU32 = int(ICGU, kind=int32)
        call Open_Write_GOpen(Name, IU32, LabFil, GVers, Title, NAtoms32, &
                              NBasis32, NBsUse32, ICharg32, Multip32, NE32, &
                              IOpCl32, ICGU32)
        IU = int(IU32, kind=int64)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Open_Write')
    end select

end subroutine Open_Write_64

! ======================================================================

subroutine Rd_2E1_32(IU, LR, NTot, LenBuf, RArr, LabVer)
    !! 32-bit version of `Rd_2E1` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int32), intent(in) :: IU, LenBuf, LR, NTot
    real(real64), intent(out) :: RArr(LR)
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IU64, LenBuf64, LR64, NTot64
    integer(int32) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

    select case (isize)
    case (32)
        call Rd_2E1_GOpen(IU, LR, NTot, LenBuf, RArr)
    case (64)
        IU64 = int(IU, kind=int64)
        LR64 = int(LR, kind=int64)
        NTot64 = int(NTot, kind=int64)
        LenBuf64 = int(LenBuf, kind=int64)
        call Rd_2E1_GOpen(IU64, LR64, NTot64, LenBuf64, RArr)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_2E1')
    end select

end subroutine Rd_2E1_32

! ======================================================================

subroutine Rd_2E1_64(IU, LR, NTot, LenBuf, RArr, LabVer)
    !! 64-bit version of `Rd_2E1` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int64), intent(in) :: IU, LenBuf, LR, NTot
    real(real64), intent(out) :: RArr(LR)
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IU32, LenBuf32, LR32, NTot32
    integer(int64) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

    select case (isize)
    case (64)
        call Rd_2E1_GOpen(IU, LR, NTot, LenBuf, RArr)
    case (32)
        IU32 = int(IU, kind=int32)
        LR32 = int(LR, kind=int32)
        NTot32 = int(NTot, kind=int32)
        LenBuf32 = int(LenBuf, kind=int32)
        call Rd_2E1_GOpen(IU32, LR32, NTot32, LenBuf32, RArr)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_2E1')
    end select

end subroutine Rd_2E1_64

! ======================================================================

subroutine Rd_2EN_32(IU, NR, LR, NTot, LenBuf, RArr, LabVer)
    !! 32-bit version of `Rd_2EN` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int32), intent(in) :: IU, LenBuf, LR, NR, NTot
    real(real64), intent(out), target :: RArr(NR,LR)
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IU64, LenBuf64, LR64, NR64, NTot64
    integer(int32) :: isize
#if GAUOPEN <= 2
    integer(int32) :: LRNR
    integer(int64) :: LRNR64
    real(real64), dimension(:), pointer :: Arr
#endif

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

#if GAUOPEN >= 3
    select case (isize)
    case (32)
        call Rd_2EN_GOpen(IU, NR, LR, NTot, LenBuf, RArr)
    case (64)
        IU64 = int(IU, kind=int64)
        NR64 = int(NR, kind=int64)
        LR64 = int(LR, kind=int64)
        NTot64 = int(NTot, kind=int64)
        LenBuf64 = int(LenBuf, kind=int64)
        call Rd_2EN_GOpen(IU64, NR64, LR64, NTot64, LenBuf64, RArr)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_2EN')
    end select
#else
    LRNR = LR*NR
    Arr(1:LRNR) => RArr
    select case (isize)
    case (32)
        call Rd_2EN_GOpen(IU, NR, LR, LRNR, NTot, LenBuf, Arr)
    case (64)
        IU64 = int(IU, kind=int64)
        NR64 = int(NR, kind=int64)
        LR64 = int(LR, kind=int64)
        LRNR64 = int(LRNR, kind=int64)
        NTot64 = int(NTot, kind=int64)
        LenBuf64 = int(LenBuf, kind=int64)
        call Rd_2EN_GOpen(IU64, NR64, LR64, LRNR64, NTot64, LenBuf64, Arr)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_2EN')
    end select
#endif

end subroutine Rd_2EN_32

! ======================================================================

subroutine Rd_2EN_64(IU, NR, LR, NTot, LenBuf, RArr, LabVer)
    !! 64-bit version of `Rd_2EN` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int64), intent(in) :: IU, LenBuf, LR, NR, NTot
    real(real64), intent(out), target :: RArr(NR,LR)
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IU32, LenBuf32, LR32, NR32, NTot32
    integer(int64) :: isize

#if GAUOPEN <= 2
    integer(int64) :: LRNR
    integer(int32) :: LRNR32
    real(real64), dimension(:), pointer :: Arr
#endif

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

#if GAUOPEN >= 3
    select case (isize)
    case (64)
        call Rd_2EN_GOpen(IU, NR, LR, NTot, LenBuf, RArr)
    case (32)
        IU32 = int(IU, kind=int32)
        NR32 = int(NR, kind=int32)
        LR32 = int(LR, kind=int32)
        NTot32 = int(NTot, kind=int32)
        LenBuf32 = int(LenBuf, kind=int32)
        call Rd_2EN_GOpen(IU32, NR32, LR32, NTot32, LenBuf32, RArr)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_2EN')
    end select
#else
    LRNR = LR*NR
    Arr(1:LRNR) => RArr
    select case (isize)
    case (64)
        call Rd_2EN_GOpen(IU, NR, LR, LRNR, NTot, LenBuf, Arr)
    case (32)
        IU32 = int(IU, kind=int32)
        NR32 = int(NR, kind=int32)
        LR32 = int(LR, kind=int32)
        LRNR32 = int(LRNR, kind=int32)
        NTot32 = int(NTot, kind=int32)
        LenBuf32 = int(LenBuf, kind=int32)
        call Rd_2EN_GOpen(IU32, NR32, LR32, LRNR32, NTot32, LenBuf32, Arr)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_2EN')
    end select
#endif

end subroutine Rd_2EN_64

! ======================================================================

subroutine Rd_2EN_1D_32(IU, NR, LR, NTot, LenBuf, RArr, LabVer)
    !! 32-bit version of `Rd_2EN` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int32), intent(in) :: IU, LenBuf, LR, NR, NTot
    real(real64), intent(out), target :: RArr(NR*LR)
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IU64, LenBuf64, LR64, NR64, NTot64
    integer(int32) :: isize
    integer(int32) :: LRNR
    integer(int64) :: LRNR64

#if GAUOPEN >= 3
    real(real64), dimension(:,:), pointer :: Arr
#endif

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

#if GAUOPEN >= 3
    Arr(1:NR,1:LR) => RArr
    select case (isize)
    case (32)
        call Rd_2EN_GOpen(IU, NR, LR, NTot, LenBuf, Arr)
    case (64)
        IU64 = int(IU, kind=int64)
        NR64 = int(NR, kind=int64)
        LR64 = int(LR, kind=int64)
        NTot64 = int(NTot, kind=int64)
        LenBuf64 = int(LenBuf, kind=int64)
        call Rd_2EN_GOpen(IU64, NR64, LR64, NTot64, LenBuf64, Arr)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_2EN')
    end select
#else
    LRNR = LR*NR
    select case (isize)
    case (32)
        call Rd_2EN_GOpen(IU, NR, LR, LRNR, NTot, LenBuf, RArr)
    case (64)
        IU64 = int(IU, kind=int64)
        NR64 = int(NR, kind=int64)
        LR64 = int(LR, kind=int64)
        LRNR64 = int(LRNR, kind=int64)
        NTot64 = int(NTot, kind=int64)
        LenBuf64 = int(LenBuf, kind=int64)
        call Rd_2EN_GOpen(IU64, NR64, LR64, LRNR64, NTot64, LenBuf64, RArr)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_2EN')
    end select
#endif

end subroutine Rd_2EN_1D_32

! ======================================================================

subroutine Rd_2EN_1D_64(IU, NR, LR, NTot, LenBuf, RArr, LabVer)
    !! 64-bit version of `Rd_2EN` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int64), intent(in) :: IU, LenBuf, LR, NR, NTot
    real(real64), intent(out), target :: RArr(NR*LR)
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IU32, LenBuf32, LR32, NR32, NTot32
    integer(int64) :: isize
    integer(int64) :: LRNR
    integer(int32) :: LRNR32
#if GAUOPEN >= 3
    real(real64), dimension(:,:), pointer :: Arr
#endif

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

#if GAUOPEN >= 3
    Arr(1:NR,1:LR) => RArr
    select case (isize)
    case (64)
        call Rd_2EN_GOpen(IU, NR, LR, NTot, LenBuf, Arr)
    case (32)
        IU32 = int(IU, kind=int32)
        NR32 = int(NR, kind=int32)
        LR32 = int(LR, kind=int32)
        NTot32 = int(NTot, kind=int32)
        LenBuf32 = int(LenBuf, kind=int32)
        call Rd_2EN_GOpen(IU32, NR32, LR32, NTot32, LenBuf32, Arr)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_2EN')
    end select
#else
    LRNR = LR*NR
    select case (isize)
    case (64)
        call Rd_2EN_GOpen(IU, NR, LR, LRNR, NTot, LenBuf, RArr)
    case (32)
        IU32 = int(IU, kind=int32)
        NR32 = int(NR, kind=int32)
        LR32 = int(LR, kind=int32)
        LRNR32 = int(LRNR, kind=int32)
        NTot32 = int(NTot, kind=int32)
        LenBuf32 = int(LenBuf, kind=int32)
        call Rd_2EN_GOpen(IU32, NR32, LR32, LRNR32, NTot32, LenBuf32, RArr)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_2EN')
    end select
#endif

end subroutine Rd_2EN_1D_64

! ======================================================================

subroutine Rd_SpA_32(IU, NI, NR, NTot, LenBuf, IArr, RArr, LabVer)
    !! 32-bit version of `Rd_SpA` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int32), intent(in) :: IU, LenBuf, NI, NR, NTot
    integer(int32), intent(out) :: IArr(NI*NTot)
    real(real64), intent(out) :: RArr(NR*NTot)
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IU64, LenBuf64, NI64, NR64, NTot64
    integer(int64), dimension(:), allocatable :: IArr64
    integer(int32) :: isize

#if GAUOPEN >= 3
    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

    select case (isize)
    case (32)
        call Rd_SpA_GOpen(IU, NI, NR, NTot, LenBuf, IArr, RArr)
    case (64)
        IU64 = int(IU, kind=int64)
        NI64 = int(NI, kind=int64)
        NTot64 = int(NTot, kind=int64)
        LenBuf64 = int(LenBuf, kind=int64)
        allocate(IArr64(NI*NTot))
        call Rd_SpA_GOpen(IU64, NI64, NR64, NTot64, LenBuf64, IArr64, RArr)
        IArr = int(IArr64, kind=int32)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_SpA')
    end select
#else
    call run%error%raise_deverror('version', &
        'Rd_SpA is not supported by this version of GauOpen')
#endif

end subroutine Rd_SpA_32

! ======================================================================

subroutine Rd_SpA_64(IU, NI, NR, NTot, LenBuf, IArr, RArr, LabVer)
    !! 64-bit version of `Rd_SpA` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int64), intent(in) :: IU, LenBuf, NI, NR, NTot
    integer(int64), intent(out) :: IArr(NI*NTot)
    real(real64), intent(out) :: RArr(NR*NTot)
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IU32, LenBuf32, NI32, NR32, NTot32
    integer(int32), dimension(:), allocatable :: IArr32
    integer(int64) :: isize

#if GAUOPEN >= 3
    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

    select case (isize)
    case (64)
        call Rd_SpA_GOpen(IU, NI, NR, NTot, LenBuf, IArr, RArr)
    case (32)
        IU32 = int(IU, kind=int32)
        NI32 = int(NI, kind=int32)
        NTot32 = int(NTot, kind=int32)
        LenBuf32 = int(LenBuf, kind=int32)
        allocate(IArr32(NI*NTot))
        call Rd_SpA_GOpen(IU32, NI32, NR32, NTot32, LenBuf32, IArr32, RArr)
        IArr = int(IArr32, kind=int64)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_SpA')
    end select
#else
    call run%error%raise_deverror('version', &
        'Rd_SpA is not supported by this version of GauOpen')
#endif

end subroutine Rd_SpA_64

! ======================================================================

subroutine Rd_SpAC_32(IU, NI, NR, NTot, LenBuf, IArr, RArr, LabVer)
    !! 32-bit version of `Rd_SpAC` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int32), intent(in) :: IU, LenBuf, NI, NR, NTot
    integer(int32), intent(out) :: IArr(NI*NTot)
    complex(real64), intent(out) :: RArr(NR*NTot)
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IU64, LenBuf64, NI64, NR64, NTot64
    integer(int64), dimension(:), allocatable :: IArr64
    integer(int32) :: isize

#if GAUOPEN >= 3
    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

    select case (isize)
    case (32)
        call Rd_SpAC_GOpen(IU, NI, NR, NTot, LenBuf, IArr, RArr)
    case (64)
        IU64 = int(IU, kind=int64)
        NI64 = int(NI, kind=int64)
        NTot64 = int(NTot, kind=int64)
        LenBuf64 = int(LenBuf, kind=int64)
        allocate(IArr64(NI*NTot))
        call Rd_SpAC_GOpen(IU64, NI64, NR64, NTot64, LenBuf64, IArr64, RArr)
        IArr = int(IArr64, kind=int32)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_SpAC')
    end select
#else
    call run%error%raise_deverror('version', &
        'Rd_SpAC is not supported by this version of GauOpen')
#endif

end subroutine Rd_SpAC_32

! ======================================================================

subroutine Rd_SpAC_64(IU, NI, NR, NTot, LenBuf, IArr, RArr, LabVer)
    !! 64-bit version of `Rd_SpAC` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int64), intent(in) :: IU, LenBuf, NI, NR, NTot
    integer(int64), intent(out) :: IArr(NI*NTot)
    complex(real64), intent(out) :: RArr(NR*NTot)
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IU32, LenBuf32, NI32, NR32, NTot32
    integer(int32), dimension(:), allocatable :: IArr32
    integer(int64) :: isize

#if GAUOPEN >= 3
    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

    select case (isize)
    case (64)
        call Rd_SpAC_GOpen(IU, NI, NR, NTot, LenBuf, IArr, RArr)
    case (32)
        IU32 = int(IU, kind=int32)
        NI32 = int(NI, kind=int32)
        NTot32 = int(NTot, kind=int32)
        LenBuf32 = int(LenBuf, kind=int32)
        allocate(IArr32(NI*NTot))
        call Rd_SpAC_GOpen(IU32, NI32, NR32, NTot32, LenBuf32, IArr32, RArr)
        IArr = int(IArr32, kind=int64)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_SpAC')
    end select
#else
    call run%error%raise_deverror('version', &
        'Rd_SpAC is not supported by this version of GauOpen')
#endif

end subroutine Rd_SpAC_64

! ======================================================================

subroutine Rd_Head_32(IU, NLab, NAtoms, NBasis, IAn, IAtTyp, AtmChg, C, &
                      IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, IDum9, NShlAO, &
                      NPrmAO, NShlDB, NPrmDB, NBTot, LabVer)
    !! 32-bit version of `Rd_Head` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int32), intent(in) :: IU, NLab, NAtoms, NBasis
    integer(int32), intent(out) :: IAn(NAtoms), IAtTyp(NAtoms), &
        IBfAtm(NBasis), IBfTyp(NBasis), NFC, NFV, ITran, IDum9, &
        NShlAO, NPrmAO, NShlDB, NPrmDB, NBTot
    real(real64), intent(out) :: AtmChg(NAtoms), C(3*NAtoms), &
        AtmWgt(NAtoms)
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IDum9_64, ITran64, IU64, NAtoms64, NBasis64, NBTot64, &
        NFC64, NFV64, NLab64, NPrmAO64, NPrmDB64, NShlAO64, NShlDB64
    integer(int64), dimension(:), allocatable :: IAn64, IAtTyp64, IBfAtm64, &
        IBfTyp64
    integer(int32) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

    select case (isize)
    case (32)
        call Rd_Head_GOpen(IU, NLab, NAtoms, NBasis, IAn, IAtTyp, AtmChg, &
                           C, IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, &
                           IDum9, NShlAO, NPrmAO, NShlDB, NPrmDB, NBTot)
    case (64)
        IU64 = int(IU, kind=int64)
        NLab64 = int(NLab, kind=int64)
        NAtoms64 = int(NAtoms, kind=int64)
        NBasis64 = int(NBasis, kind=int64)
        allocate(IAn64(NAtoms), IAtTyp64(NAtoms), IBfAtm64(NBasis), &
                 IBfTyp64(NBasis))
        call Rd_Head_GOpen(IU64, NLab64, NAtoms64, NBasis64, IAn64, IAtTyp64, &
                           AtmChg, C, IBfAtm64, IBfTyp64, AtmWgt, NFC64, &
                           NFV64, ITran64, IDum9_64, NShlAO64, NPrmAO64, &
                           NShlDB64, NPrmDB64, NBTot64)
        IAn = int(IAn64, kind=int32)
        IAtTyp = int(IAtTyp64, kind=int32)
        IBfAtm = int(IBfAtm64, kind=int32)
        IBfTyp = int(IBfTyp64, kind=int32)
        NFC = int(NFC64, kind=int32)
        NFV = int(NFV64, kind=int32)
        ITran = int(ITran64, kind=int32)
        IDum9 = int(IDum9_64, kind=int32)
        NShlAO = int(NShlAO64, kind=int32)
        NPrmAO = int(NPrmAO64, kind=int32)
        NShlDB = int(NShlDB64, kind=int32)
        NPrmDB = int(NPrmDB64, kind=int32)
        NBTot = int(NBTot64, kind=int32)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_Head')
    end select

end subroutine Rd_Head_32

! ======================================================================

subroutine Rd_Head_64(IU, NLab, NAtoms, NBasis, IAn, IAtTyp, AtmChg, C, &
                      IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, IDum9, NShlAO, &
                      NPrmAO, NShlDB, NPrmDB, NBTot, LabVer)
    !! 64-bit version of `Rd_Head` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int64), intent(in) :: IU, NLab, NAtoms, NBasis
    integer(int64), intent(out) :: IAn(NAtoms), IAtTyp(NAtoms), &
        IBfAtm(NBasis), IBfTyp(NBasis), NFC, NFV, ITran, IDum9, &
        NShlAO, NPrmAO, NShlDB, NPrmDB, NBTot
    real(real64), intent(out) :: AtmChg(NAtoms), C(3*NAtoms), &
        AtmWgt(NAtoms)
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IDum9_32, ITran32, IU32, NAtoms32, NBasis32, NBTot32, &
        NFC32, NFV32, NLab32, NPrmAO32, NPrmDB32, NShlAO32, NShlDB32
    integer(int32), dimension(:), allocatable :: IAn32, IAtTyp32, IBfAtm32, &
        IBfTyp32
    integer(int64) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

    select case (isize)
    case (64)
        call Rd_Head_GOpen(IU, NLab, NAtoms, NBasis, IAn, IAtTyp, AtmChg, &
                           C, IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, &
                           IDum9, NShlAO, NPrmAO, NShlDB, NPrmDB, NBTot)
    case (32)
        IU32 = int(IU, kind=int32)
        NLab32 = int(NLab, kind=int32)
        NAtoms32 = int(NAtoms, kind=int32)
        NBasis32 = int(NBasis, kind=int32)
        allocate(IAn32(NAtoms), IAtTyp32(NAtoms), IBfAtm32(NBasis), &
                 IBfTyp32(NBasis))
        call Rd_Head_GOpen(IU32, NLab32, NAtoms32, NBasis32, IAn32, IAtTyp32, &
                           AtmChg, C, IBfAtm32, IBfTyp32, AtmWgt, NFC32, &
                           NFV32, ITran32, IDum9_32, NShlAO32, NPrmAO32, &
                           NShlDB32, NPrmDB32, NBTot32)
        IAn = int(IAn32, kind=int64)
        IAtTyp = int(IAtTyp32, kind=int64)
        IBfAtm = int(IBfAtm32, kind=int64)
        IBfTyp = int(IBfTyp32, kind=int64)
        NFC = int(NFC32, kind=int64)
        NFV = int(NFV32, kind=int64)
        ITran = int(ITran32, kind=int64)
        IDum9 = int(IDum9_32, kind=int64)
        NShlAO = int(NShlAO32, kind=int64)
        NPrmAO = int(NPrmAO32, kind=int64)
        NShlDB = int(NShlDB32, kind=int64)
        NPrmDB = int(NPrmDB32, kind=int64)
        NBTot = int(NBTot32, kind=int64)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_Head')
    end select

end subroutine Rd_Head_64

! ======================================================================

subroutine Rd_HeadA_32(IU, NLab, NAtoms, NBasis, IAn, IAtTyp, AtmChg, C, &
                       IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, IDum9, Rec11, &
                       LabVer)
    !! 32-bit version of `Rd_HeadA` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer, parameter :: LRec11 = 16
    integer(int32), intent(in) :: IU, NLab, NAtoms, NBasis
    integer(int32), intent(out) :: IAn(NAtoms), IAtTyp(NAtoms), &
        IBfAtm(NBasis), IBfTyp(NBasis), NFC, NFV, ITran, IDum9, &
        Rec11(LRec11)
    real(real64), intent(out) :: AtmChg(NAtoms), C(3*NAtoms), &
        AtmWgt(NAtoms)
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IDum9_64, ITran64, IU64, NAtoms64, NBasis64, NFC64, &
        NFV64, NLab64
    integer(int64), dimension(:), allocatable :: IAn64, IAtTyp64, IBfAtm64, &
        IBfTyp64, Rec11_64
    integer(int32) :: isize

#if GAUOPEN >= 3
    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

    select case (isize)
    case (32)
        call Rd_HeadA_GOpen(IU, NLab, NAtoms, NBasis, IAn, IAtTyp, AtmChg, &
                            C, IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, &
                            IDum9, Rec11)
    case (64)
        IU64 = int(IU, kind=int64)
        NLab64 = int(NLab, kind=int64)
        NAtoms64 = int(NAtoms, kind=int64)
        NBasis64 = int(NBasis, kind=int64)
        allocate(IAn64(NAtoms), IAtTyp64(NAtoms), IBfAtm64(NBasis), &
                 IBfTyp64(NBasis), Rec11_64(LRec11))
        call Rd_HeadA_GOpen(IU64, NLab64, NAtoms64, NBasis64, IAn64, &
                            IAtTyp64, AtmChg, C, IBfAtm64, IBfTyp64, AtmWgt, &
                            NFC64, NFV64, ITran64, IDum9_64, Rec11_64)
        IAn = int(IAn64, kind=int32)
        IAtTyp = int(IAtTyp64, kind=int32)
        IBfAtm = int(IBfAtm64, kind=int32)
        IBfTyp = int(IBfTyp64, kind=int32)
        NFC = int(NFC64, kind=int32)
        NFV = int(NFV64, kind=int32)
        ITran = int(ITran64, kind=int32)
        IDum9 = int(IDum9_64, kind=int32)
        Rec11 = int(Rec11_64, kind=int32)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_HeadA')
    end select
#else
    call run%error%raise_deverror('version', &
        'Rd_HeadA is not supported by this version of GauOpen')
#endif

end subroutine Rd_HeadA_32

! ======================================================================

subroutine Rd_HeadA_64(IU, NLab, NAtoms, NBasis, IAn, IAtTyp, AtmChg, C, &
                       IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, IDum9, Rec11, &
                       LabVer)
    !! 64-bit version of `Rd_HeadA` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer, parameter :: LRec11 = 16
    integer(int64), intent(in) :: IU, NLab, NAtoms, NBasis
    integer(int64), intent(out) :: IAn(NAtoms), IAtTyp(NAtoms), &
        IBfAtm(NBasis), IBfTyp(NBasis), NFC, NFV, ITran, IDum9, &
        Rec11(LRec11)
    real(real64), intent(out) :: AtmChg(NAtoms), C(3*NAtoms), &
        AtmWgt(NAtoms)
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IDum9_32, ITran32, IU32, NAtoms32, NBasis32, NFC32, &
        NFV32, NLab32
    integer(int32), dimension(:), allocatable :: IAn32, IAtTyp32, IBfAtm32, &
        IBfTyp32, Rec11_32
    integer(int64) :: isize

#if GAUOPEN >= 3
    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

    select case (isize)
    case (64)
        call Rd_HeadA_GOpen(IU, NLab, NAtoms, NBasis, IAn, IAtTyp, AtmChg, &
                            C, IBfAtm, IBfTyp, AtmWgt, NFC, NFV, ITran, &
                            IDum9, Rec11)
    case (32)
        IU32 = int(IU, kind=int32)
        NLab32 = int(NLab, kind=int32)
        NAtoms32 = int(NAtoms, kind=int32)
        NBasis32 = int(NBasis, kind=int32)
        allocate(IAn32(NAtoms), IAtTyp32(NAtoms), IBfAtm32(NBasis), &
                 IBfTyp32(NBasis), Rec11_32(LRec11))
        call Rd_HeadA_GOpen(IU32, NLab32, NAtoms32, NBasis32, IAn32, &
                            IAtTyp32, AtmChg, C, IBfAtm32, IBfTyp32, AtmWgt, &
                            NFC32, NFV32, ITran32, IDum9_32, Rec11_32)
        IAn = int(IAn32, kind=int64)
        IAtTyp = int(IAtTyp32, kind=int64)
        IBfAtm = int(IBfAtm32, kind=int64)
        IBfTyp = int(IBfTyp32, kind=int64)
        NFC = int(NFC32, kind=int64)
        NFV = int(NFV32, kind=int64)
        ITran = int(ITran32, kind=int64)
        IDum9 = int(IDum9_32, kind=int64)
        Rec11 = int(Rec11_32, kind=int64)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_HeadA')
    end select
#else
    call run%error%raise_deverror('version', &
        'Rd_HeadA is not supported by this version of GauOpen')
#endif

end subroutine Rd_HeadA_64

! ======================================================================

subroutine Rd_IBuf_32(IU, LR, LenBuf, Arr, LabVer)
    !! 32-bit version of ` Rd_IBuf` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int32), intent(in) :: IU, LenBuf, LR
    integer(int32), intent(out) :: Arr(LR)
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IU64, LenBuf64, LR64
    integer(int64), dimension(:), allocatable :: Arr64
    integer(int32) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

    select case (isize)
    case (32)
        call Rd_IBuf_GOpen(IU, LR, LenBuf, Arr)
    case (64)
        IU64 = int(IU, kind=int64)
        LR64 = int(LR, kind=int64)
        LenBuf64 = int(LenBuf, kind=int64)
        allocate(Arr64(LR))
        call Rd_IBuf_GOpen(IU64, LR64, LenBuf64, Arr64)
        Arr = int(Arr64, kind=int32)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_IBuf')
    end select

end subroutine Rd_IBuf_32

! ======================================================================

subroutine Rd_IBuf_64(IU, LR, LenBuf, Arr, LabVer)
    !! 64-bit version of ` Rd_IBuf` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int64), intent(in) :: IU, LenBuf, LR
    integer(int64), intent(out) :: Arr(LR)
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IU32, LenBuf32, LR32
    integer(int32), dimension(:), allocatable :: Arr32
    integer(int64) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

    select case (isize)
    case (64)
        call Rd_IBuf_GOpen(IU, LR, LenBuf, Arr)
    case (32)
        IU32 = int(IU, kind=int32)
        LR32 = int(LR, kind=int32)
        LenBuf32 = int(LenBuf, kind=int32)
        allocate(Arr32(LR))
        call Rd_IBuf_GOpen(IU32, LR32, LenBuf32, Arr32)
        Arr = int(Arr32, kind=int64)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_IBuf')
    end select

end subroutine Rd_IBuf_64

! ======================================================================

subroutine Rd_Labl_32(IU, IVers, CBuf, NI, NR, NTot, LenBuf, N1, N2, N3, N4, &
                      N5, TypeA, NRI, EOF, LabVer)
    !! 32-bit version of ` Rd_Lab` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int32), intent(in) :: IU, IVers
    integer(int32), intent(out) :: NI, NR, NTot, LenBuf, N1, N2, N3, N4, N5, &
        TypeA, NRI
    character(len=64), intent(out) :: CBuf
    Logical(int32), intent(out) :: EOF
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IU64, IVers64, LenBuf64, N1_64, N2_64, N3_64, N4_64, &
        N5_64, NI64, NR64, NRI64, NTot64, TypeA64
    logical(int64) :: EOF64
    integer(int32) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

    select case (isize)
    case (32)
        call Rd_Labl_GOpen(IU, IVers, CBuf, NI, NR, NTot, LenBuf, N1, N2, N3, &
                           N4, N5, TypeA, NRI, EOF)
    case (64)
        IU64 = int(IU, kind=int64)
        IVers64 = int(IVers, kind=int64)
        call Rd_Labl_GOpen(IU64, IVers64, CBuf, NI64, NR64, NTot64, LenBuf64, &
                           N1_64, N2_64, N3_64, N4_64, N5_64, TypeA64, NRI64, &
                           EOF64)
        NI = int(NI64, kind=int32)
        NR = int(NR64, kind=int32)
        NTot = int(NTot64, kind=int32)
        LenBuf = int(LenBuf64, kind=int32)
        N1 = int(N1_64, kind=int32)
        N2 = int(N2_64, kind=int32)
        N3 = int(N3_64, kind=int32)
        N4 = int(N4_64, kind=int32)
        N5 = int(N5_64, kind=int32)
        TypeA = int(TypeA64, kind=int32)
        NRI = int(NRI64, kind=int32)
        EOF = EOF64
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_Labl')
    end select

end subroutine Rd_Labl_32

! ======================================================================

subroutine Rd_Labl_64(IU, IVers, CBuf, NI, NR, NTot, LenBuf, N1, N2, N3, N4, &
                      N5, TypeA, NRI, EOF, LabVer)
    !! 64-bit version of ` Rd_Lab` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int64), intent(in) :: IU, IVers
    integer(int64), intent(out) :: NI, NR, NTot, LenBuf, N1, N2, N3, N4, N5, &
        TypeA, NRI
    character(len=64), intent(out) :: CBuf
    Logical(int64), intent(out) :: EOF
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IU32, IVers32, LenBuf32, N1_32, N2_32, N3_32, N4_32, &
        N5_32, NI32, NR32, NRI32, NTot32, TypeA32
    logical(int32) :: EOF32
    integer(int64) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

    select case (isize)
    case (64)
        call Rd_Labl_GOpen(IU, IVers, CBuf, NI, NR, NTot, LenBuf, N1, N2, N3, &
                           N4, N5, TypeA, NRI, EOF)
    case (32)
        IU32 = int(IU, kind=int32)
        IVers32 = int(IVers, kind=int32)
        call Rd_Labl_GOpen(IU32, IVers32, CBuf, NI32, NR32, NTot32, LenBuf32, &
                           N1_32, N2_32, N3_32, N4_32, N5_32, TypeA32, NRI32, &
                           EOF32)
        NI = int(NI32, kind=int64)
        NR = int(NR32, kind=int64)
        NTot = int(NTot32, kind=int64)
        LenBuf = int(LenBuf32, kind=int64)
        N1 = int(N1_32, kind=int64)
        N2 = int(N2_32, kind=int64)
        N3 = int(N3_32, kind=int64)
        N4 = int(N4_32, kind=int64)
        N5 = int(N5_32, kind=int64)
        TypeA = int(TypeA32, kind=int64)
        NRI = int(NRI32, kind=int64)
        EOF = EOF32
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_Labl')
    end select

end subroutine Rd_Labl_64

! ======================================================================

subroutine Rd_RInd_32(IU, NR, LR, NRLR, NTot, LenBuf, RArr, LabVer)
    !! 32-bit version of ` Rd_Lab` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int32), intent(in) :: IU, NR, LR, NRLR, NTot, LenBuf
    real(real64), intent(out) :: RArr(NRLR)
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IU64, LenBuf64, LR64, NR64, NRLR64, NTot64
    integer(int32) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

    select case (isize)
    case (32)
        call Rd_RInd_GOpen(IU, NR, LR, NRLR, NTot, LenBuf, RArr)
    case (64)
        IU64 = int(IU, kind=int64)
        LenBuf64 = int(LenBuf, kind=int64)
        LR64 = int(LR, kind=int64)
        NR64 = int(NR, kind=int64)
        NRLR64 = int(NRLR, kind=int64)
        NTot64 = int(NTot, kind=int64)
        call Rd_RInd_GOpen(IU64, NR64, LR64, NRLR64, NTot64, LenBuf64, RArr)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_RInd')
    end select

end subroutine Rd_RInd_32

! ======================================================================

subroutine Rd_RInd_64(IU, NR, LR, NRLR, NTot, LenBuf, RArr, LabVer)
    !! 64-bit version of ` Rd_Lab` from `qcmatrixio`.
    !!
    !! This version can read data with 32 and 64-bits labels.
    integer(int64), intent(in) :: IU, NR, LR, NRLR, NTot, LenBuf
    real(real64), intent(out) :: RArr(NRLR)
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IU32, LenBuf32, LR32, NR32, NRLR32, NTot32
    integer(int64) :: isize

    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

    select case (isize)
    case (64)
        call Rd_RInd_GOpen(IU, NR, LR, NRLR, NTot, LenBuf, RArr)
    case (32)
        IU32 = int(IU, kind=int32)
        LenBuf32 = int(LenBuf, kind=int32)
        LR32 = int(LR, kind=int32)
        NR32 = int(NR, kind=int32)
        NRLR32 = int(NRLR, kind=int32)
        NTot32 = int(NTot, kind=int32)
        call Rd_RInd_GOpen(IU32, NR32, LR32, NRLR32, NTot32, LenBuf32, RArr)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Rd_RInd')
    end select

end subroutine Rd_RInd_64

! ======================================================================

subroutine Wr_LChBuf_32(IU, Label, LenBuf, X, LabVer)
    !! 32-bit version of `Wr_LChBuf` from `qcmatrix`.
    !!
    !! This version can write 32 and 64-bits labels.
    character(len=*), intent(in) :: Label, X
    integer(int32), intent(in) :: IU, LenBuf
    integer(int32), intent(in), optional :: LabVer

    integer(int64) :: IU64, LenBuf64
    integer(int32) :: isize

#if GAUOPEN >= 3
    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 32
    end if

    select case (isize)
    case (32)
        call Wr_LChBuf_GOpen(IU, Label, LenBuf, X)
    case (64)
        IU64 = int(IU, kind=int64)
        LenBuf64 = int(LenBuf, kind=int64)
        call Wr_LChBuf_GOpen(IU64, Label, LenBuf64, X)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Wr_LChBuf')
    end select
#else
    call run%error%raise_deverror('version', &
        'Wr_LChBuf is not supported by this version of GauOpen')
#endif

end subroutine Wr_LChBuf_32

! ======================================================================

subroutine Wr_LChBuf_64(IU, Label, LenBuf, X, LabVer)
    !! 64-bit version of `Wr_LChBuf` from `qcmatrix`.
    !!
    !! This version can write 32 and 64-bits labels.
    character(len=*), intent(in) :: Label, X
    integer(int64), intent(in) :: IU, LenBuf
    integer(int64), intent(in), optional :: LabVer

    integer(int32) :: IU32, LenBuf32
    integer(int64) :: isize

#if GAUOPEN >= 3
    if (present(LabVer)) then
        isize = LabVer
    else
        isize = 64
    end if

    select case (isize)
    case (64)
        call Wr_LChBuf_GOpen(IU, Label, LenBuf, X)
    case (32)
        IU32 = int(IU, kind=int32)
        LenBuf32 = int(LenBuf, kind=int32)
        call Wr_LChBuf_GOpen(IU32, Label, LenBuf32, X)
    case default
        call run%error%raise_deverror('argval', &
            'Unsupported label version for Wr_LChBuf')
    end select
#else
    call run%error%raise_deverror('version', &
        'Wr_LChBuf is not supported by this version of GauOpen')
#endif

end subroutine Wr_LChBuf_64

! =======================================
!  Wrapper Routines to Logical Functions
! =======================================

! Some functions cannot be distinguished since they do not involve any
! integers as arguments.
! We do a simple wrapper by calling one of the version (I4, that should use
! the default compiler parameters, which should be the default compilation of
! ELEMENTS in general).
! By doing this call, we ensure there is in theory no problem of conversion.

! ======================================================================

function AOInts(CBuf) result(res)
    logical :: res
    character(len=*), intent(in) :: CBuf
    interface
        function AOInts_I4(CBuf0)
            use iso_fortran_env, only: int32
            logical(int32) :: AOInts_I4
            character(len=*), intent(in) :: CBuf0
        end function AOInts_I4
    end interface
    res = logical(AOInts_I4(CBuf), kind=kind(res))
end function AOInts

! ======================================================================

function DAOInts(CBuf) result(res)
    logical :: res
    character(len=*), intent(in) :: CBuf
#if GAUOPEN >= 3
    interface
        function DAOInts_I4(CBuf0)
            use iso_fortran_env, only: int32
            logical(int32) :: DAOInts_I4
            character(len=*), intent(in) :: CBuf0
        end function DAOInts_I4
    end interface
    res = logical(DAOInts_I4(CBuf), kind=kind(res))
#else
    call run%error%raise_deverror('version', &
        'DAOInts is not supported by this version of GauOpen')
#endif

end function DAOInts

! ======================================================================

#if GAUOPEN >= 3
function QCM_AOInts(CBuf) result(res)
    logical :: res
    character(len=*), intent(in) :: CBuf
    interface
        function QCM_AOInts_I4(CBuf0)
            use iso_fortran_env, only: int32
            logical(int32) :: QCM_AOInts_I4
            character(len=*), intent(in) :: CBuf0
        end function QCM_AOInts_I4
    end interface
    res = logical(QCM_AOInts_I4(CBuf), kind=kind(res))
end function QCM_AOInts
#endif

! ======================================================================

#if GAUOPEN >= 3
function QCM_DAOInts(CBuf) result(res)
    logical :: res
    character(len=*), intent(in) :: CBuf
    interface
        function QCM_DAOInts_I4(CBuf0)
            use iso_fortran_env, only: int32
            logical(int32) :: QCM_DAOInts_I4
            character(len=*), intent(in) :: CBuf0
        end function QCM_DAOInts_I4
    end interface
    res = logical(QCM_DAOInts_I4(CBuf), kind=kind(res))
end function QCM_DAOInts
#endif

! ======================================================================

! ===============================================
!  Dummy Routines for Interface with GFAF Module
! ===============================================

! Older version of GauOpen did not provide some basic features used by the
! gfaf_io module.  The following routines provides a wrapper for it.

! ======================================================================

#if GAUOPEN <= 2
subroutine Rd_ChBuf_Dum_I4(IU, LR, LenBuf, CArr)
    integer(int32), intent(in) :: IU, LenBuf, LR
    character(len=*), intent(inout) :: CArr

    call run%error%raise_deverror('version', &
        'Rd_ChBuf is not supported by this version of GauOpen')
end subroutine Rd_ChBuf_Dum_I4
! ======================================================================
subroutine Rd_ChBuf_Dum_I8(IU, LR, LenBuf, CArr)
    integer(int64), intent(in) :: IU, LenBuf, LR
    character(len=*), intent(inout) :: CArr

    call run%error%raise_deverror('version', &
        'Rd_ChBuf is not supported by this version of GauOpen')
end subroutine Rd_ChBuf_Dum_I8
#endif


! ======================================================================

#if GAUOPEN <= 2
subroutine Rd_ChBufI_Dum_I4(IU, LR, LenBuf, CArr, IArr)
    integer(int32), intent(in) :: IU, LenBuf, LR
    character(len=*), intent(inout) :: CArr
    integer(int32), intent(out) :: IArr(LR)

    call run%error%raise_deverror('version', &
        'Rd_ChBufI is not supported by this version of GauOpen')
end subroutine Rd_ChBufI_Dum_I4
! ======================================================================
subroutine Rd_ChBufI_Dum_I8(IU, LR, LenBuf, CArr, IArr)
    integer(int64), intent(in) :: IU, LenBuf, LR
    character(len=*), intent(inout) :: CArr
    integer(int64), intent(out) :: IArr(LR)

    call run%error%raise_deverror('version', &
        'Rd_ChBufI is not supported by this version of GauOpen')
end subroutine Rd_ChBufI_Dum_I8
#endif

! ======================================================================

end module gauopen_drv