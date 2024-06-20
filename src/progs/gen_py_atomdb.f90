program gen_py_atomdb
    !! Program to generate a source file for Python with the atomic DB.
    !!
    !! The program generates a source file compatible with Python,
    !! building the atomic DB as a list of dictionary.
    !! By default, the program creates a function atomic_data, that can
    !! return a dictionary of dictionaries containing the relevant data,
    !! to be used by ESTAMPES (base.data.atom)

    use iso_fortran_env, only: real64
    use atominfo, only: atdata
    use physics, only: bohr => bohr_radius
    use parse_cmdline, only: CmdArgDB
    use string, only: locase

    implicit none

    integer :: iout, iunit
    character(len=1024) :: arg
    type(CmdArgDB) :: opts
    
    opts = CmdArgDB(progname='gen_py_atomdb')

    if (opts%has_error()) then
        write(*, '(a)') 'Failed to initialize the command-line parser'
        write(*, '(a)') opts%get_error()
        stop 1
    end if

    call opts%add_arg_char( &
        'string', label='outfile', &
        help='Output file to write the Python source.')
    call opts%add_arg_char( &
        'string', shortname='-u', longname='--unit', &
        help='Unit for distance: au, bohr or ang, angstrom')
    
    call opts%parse_args()
    if (opts%has_error()) then
        write(*, '(a)') 'Error found while parsing'
        write(*, '(a)') opts%get_error()
        stop 1
    end if

    call opts%get_value('outfile', arg)
    open(newunit=iout, file=arg)

    if (opts%is_user_set('unit')) then
        call opts%get_value('unit', arg)
        select case(locase(arg))
        case('au', 'bohr')
            iunit = 0
        case('ang', 'angstrom')
            iunit = 1
        case default
            write(*, '(a)') 'ERROR: Unsupported type of unit'
            stop 1
        end select
    else
        iunit = 1
    end if

    call write_atomic_data_head(iout)
    call write_atomic_data_entries(iout, iunit==1)

    close(iout)

contains

subroutine write_atomic_data_head(out)
    !! Write the header of the function `atomic_data`.
    !!
    !! Writes the header part (declaration, header doc,
    !! initialization) in `out`.
    integer, intent(in) :: out
    !! Identifier of the output file.

    write(out, '(a)') 'def atomic_data(*atoms: TypeAtLab) -> TypeAtData:'
    write(out, '(a)') '    """Generates atomic data.'
    write(out, '(a)') ''
    write(out, '(a)') '    Generates a dictionary containing atomic data for each atom given in'
    write(out, '(a)') '    argument.  Each item of the returned dictionary contains the'
    write(out, '(a)') '    following elements:'
    write(out, '(a)') ''
    write(out, '(a)') '    `name`'
    write(out, '(a)') '        Full name of the atom.'
    write(out, '(a)') '    `num`'
    write(out, '(a)') '        Atomic number.'
    write(out, '(a)') '    `mass`'
    write(out, '(a)') '        Atomic mass (in amu).'
    write(out, '(a)') '    `rcov`'
    write(out, '(a)') '        Covalent radius (in Ang), as a tuple (single, double, triple).'
    write(out, '(a)') '    `rvdw`'
    write(out, '(a)') '        Van der Waals radius (in Ang), as dictionary, with keys:'
    write(out, '(a)') '        "bondi64"'
    write(out, '(a)') '          A. Bondi, J. Phys. Chem. A 1964 (68) 441'
    write(out, '(a)') '          https://doi.org/10.1021/j100785a001'
    write(out, '(a)') '        "truhlar09:"'
    write(out, '(a)') '          M. Mantina, A.C. Chamberlin, R. Valero, C.J. Cramer,'
    write(out, '(a)') '          D.G. Truhlar, J. Phys. Chem. A 2009 (113) 5809.'
    write(out, '(a)') '          https://doi.org/10.1021/jp8111556'
    write(out, '(a)') '          Extension of Bondi''s set with some additional atoms from'
    write(out, '(a)') '          main group, but some values from Bondi were ignored.'
    write(out, '(a)') '        "alvarez13"'
    write(out, '(a)') '          S. Alvarez, Dalt. Trans. 2013 (42) 8617'
    write(out, '(a)') '          https://dx.doi.org/10.1039/c3dt50599e'
    write(out, '(a)') '          Statistical analysis from Cambridge Structure Database'
    write(out, '(a)') '        "rahm16"'
    write(out, '(a)') '          M. Rahm, R. Hoffmann, N.W. Ashcroft,'
    write(out, '(a)') '          Chem. Eur. J. 2016 (22) 14625.'
    write(out, '(a)') '          https://doi.org/10.1002/chem.201602949'
    write(out, '(a)') '          Radii built by considering a threshold in density of'
    write(out, '(a)') '          0.001 e.bohr^-3, computed at the PBE0/ANO-RCC level'
    write(out, '(a)') '        "truhlar_ext"'
    write(out, '(a)') '          Extended basis set considering values of bondi64 if not'
    write(out, '(a)') '          provided in original work of Trular and coworkers.'
    write(out, '(a)') '    `rvis`'
    write(out, '(a)') '        Visualization-related radius (in Ang).'
    write(out, '(a)') '    `rgb`'
    write(out, '(a)') '        Color of the atom (as tuple of integer 0-255).'
    write(out, '(a)') ''
    write(out, '(a)') '    Parameters'
    write(out, '(a)') '    ----------'
    write(out, '(a)') '    *atoms'
    write(out, '(a)') '        List of atomic symbols (or numbers).'
    write(out, '(a)') ''
    write(out, '(a)') '    Returns'
    write(out, '(a)') '    -------'
    write(out, '(a)') '    dict'
    write(out, '(a)') '        Dictionary of atomic data, containing the elements listed above'
    write(out, '(a)') '        grouped by atomic symbols.'
    write(out, '(a)') ''
    write(out, '(a)') '    Raises'
    write(out, '(a)') '    ------'
    write(out, '(a)') '    KeyError'
    write(out, '(a)') '        Unrecognized atomic symbol.'
    write(out, '(a)') ''
    write(out, '(a)') '    Notes'
    write(out, '(a)') '    -----'
    write(out, '(a)') '    * The function has a very basic support of symbols and atomic'
    write(out, '(a)') '      numbers.  A more robust procedure should rely on a first'
    write(out, '(a)') '      conversion by estampes.tools.convert_labsymb'
    write(out, '(a)') '    """'
    write(out, '(a)') '    at_data = {}'
    write(out, '(a)') '    for atom in set(atoms):'
    write(out, '(a)') '        try:'
    write(out, '(a)') '            at_symb = atom.title()'
    write(out, '(a)') '            at_idx = at_symb'
    write(out, '(a)') '        except AttributeError:'
    write(out, '(a)') '            at_symb = ELEMENTS[atom]'
    write(out, '(a)') '            at_idx = atom'

end subroutine write_atomic_data_head


subroutine write_atomic_data_entries(out, to_ang)
    !! Write the list of entries for function `atomic_data`.
    !!
    !! Writes the list of entries for the function `atomic_data
    integer, intent(in) :: out
    !! Identifier of the output file.
    logical, intent(in) :: to_ang
    !! Convert to angstrom any distance (default).

    1000 format(8x,"if at_symb == '",a,"':")
    1001 format(8x,"elif at_symb == '",a,"':")
    1002 format(12x,"at_data[at_idx] = {")
    1003 format(12x,"}")

    integer :: iat

    write(out, 1000) trim(atdata(1)%symbol)
    write(out, 1002)
    call write_dict_entry(1, out, 12, to_ang)
    write(out, 1003)

    do iat = 2, size(atdata)
        write(out, 1001) trim(atdata(iat)%symbol)
        write(out, 1002)
        call write_dict_entry(iat, out, 12, to_ang)
        write(out, 1003)
    end do

    write(out, '(8x,a)') "else:"
    write(out, '(8x,a)') "    raise KeyError('Atomic symbol not recognized')"
    write(out, '(4x,a)') "return at_data"

end subroutine write_atomic_data_entries


subroutine write_dict_entry(iat, out, indent, to_ang)
    !! Write a dictionary entry with atomic data for a single atom.
    !!
    !! Write a dictionary entry containing the data for atom `iat`.
    !! it is possible to specify a basic indentation to shift the entry.
    !!
    !! Note: the routine does not include the enclosing {} to give more
    !!       flexibility to the calling procedure.
    integer, intent(in) :: iat
    !! Atomic index/number.
    integer, intent(in) :: out
    !! Identifier of the output file.
    integer, intent(in), optional :: indent
    !! Base indentation for the printing.
    logical, intent(in), optional :: to_ang
    !! Convert to angstrom any distance (default).

    integer :: i
    real(real64) :: fact_r
    character(len=40) :: fmt_key, key, val, vals(3)

    1000 format(a,": '",a,"',")   ! string
    1001 format(a,": ",a,",")     ! keyword
    1002 format(a,": (",a,", ",a,", ",a,"),")  ! tuple of keywords
    1003 format(a,": ",a)         ! no comma

    1010 format(i0)     ! integer
    1011 format(f0.5)   ! float with 5 digits (mass)
    1012 format(f6.4)   ! float with 4 digits (rcov, rvdw, rvis)

    1020 format("(",i0,"x,a)")  ! format for key

    if (present(indent)) then
        write(fmt_key, 1020) 4 + abs(indent)
    else
        ! default directly ready for function atomic data
        write(fmt_key, 1020) 16
    end if

    if (present(to_ang)) then
        if (to_ang) then
            fact_r = bohr
        else
            fact_r = 1.0_real64
        end if
    else
        fact_r = bohr
    end if

    if (iat > 0 .and. iat <= size(atdata)) then
        write(key, fmt_key) "'symb'"
        write(out, 1000) trim(key), trim(atdata(iat)%symbol)

        write(key, fmt_key) "'name'"
        write(out, 1000) trim(key), trim(atdata(iat)%name)

        write(key, fmt_key) "'num'"
        write(val, 1010) atdata(iat)%number
        write(out, 1001) trim(key), trim(val)

        write(key, fmt_key) "'mass'"
        if (atdata(iat)%mass == 0.0_real64) then
            val = 'None'
        else
            write(val, 1011) atdata(iat)%mass
        end if
        write(out, 1001) trim(key), trim(val)

        write(key, fmt_key) "'rcov'"
        do i = 1, 3
            ! We use epsilon here since radii can be computed/transformed
            if (abs(atdata(iat)%rcov(i)) < epsilon(fact_r)) then
                vals(i) = 'None'
            else
                write(vals(i), 1012) atdata(iat)%rcov(i)*fact_r
            end if
        end do
        write(out, 1002) trim(key), (trim(vals(i)), i=1,3)

        write(key, fmt_key) "'rvdw'"
        ! We construct a dictionary here for the different values.
        write(out, 1003) trim(key), '{'
        write(key, fmt_key) "    'bondi64'"
        if (abs(atdata(iat)%rvdw(1)) < epsilon(fact_r)) then
            val = 'None'
        else
            write(val, 1012) atdata(iat)%rvdw(1)*fact_r
        end if
        write(out, 1001) trim(key), trim(val)
        write(key, fmt_key) "    'truhlar09'"
        if (abs(atdata(iat)%rvdw(2)) < epsilon(fact_r)) then
            val = 'None'
        else
            write(val, 1012) atdata(iat)%rvdw(2)*fact_r
        end if
        write(out, 1001) trim(key), trim(val)
        write(key, fmt_key) "    'alvarez13'"
        if (abs(atdata(iat)%rvdw(3)) < epsilon(fact_r)) then
            val = 'None'
        else
            write(val, 1012) atdata(iat)%rvdw(3)*fact_r
        end if
        write(out, 1001) trim(key), trim(val)
        write(key, fmt_key) "    'rahm16'"
        if (abs(atdata(iat)%rvdw(4)) < epsilon(fact_r)) then
            val = 'None'
        else
            write(val, 1012) atdata(iat)%rvdw(4)*fact_r
        end if
        write(key, fmt_key) "    'truhlar_ext'"
        if (abs(atdata(iat)%get_rvdw('truhlar_ext')) < epsilon(fact_r)) then
            val = 'None'
        else
            write(val, 1012) atdata(iat)%get_rvdw('truhlar_ext')*fact_r
        end if
        write(out, 1001) trim(key), trim(val)
        
        write(out, fmt_key) '}'

        write(key, fmt_key) "'rvis'"
        if (abs(atdata(iat)%rvis) < epsilon(fact_r)) then
            val = 'None'
        else
            write(val, 1012) atdata(iat)%rvis*fact_r
        end if
        write(out, 1001) trim(key), trim(val)

        write(key, fmt_key) "'rgb'"
        if (all(atdata(iat)%color == 0)) then
            write(out, 1001) trim(key), 'None'
        else
            do i = 1, 3
                write(vals(i), 1010) atdata(iat)%color(i)
            end do
            write(out, 1002) trim(key), (trim(vals(i)), i=1,3)
        end if
    else
        print *, 'Unsupported value of iat'
        stop 1
    end if

end subroutine write_dict_entry

end program gen_py_atomdb
