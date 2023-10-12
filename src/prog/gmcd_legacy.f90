module gmcd_legacy
    !! Module with legacy subroutines from SOS/GUVCDE
    !!
    !! New or "modernized" procedures, used to preserve some of the
    !!   original behavior and facilitate the migration.
    use iso_fortran_env, only: real64
    use exception, only: ArgumentError, BaseException, InitError, RaiseError
    use output, only: iu_out
    use electronic, only: convert_AO2MO
    use moldata
    use physic, only: PhysFact
    use gmcd_output, only: shell_lmxyz

contains

! ======================================================================

subroutine build_MOs(n_ab, n_ao, n_mos, c_ia, txt_fmt)
    !! Build properties in MOs
    !!
    !! Reads properties stored on scratch files and expressed in atomic
    !!   orbitals and converts them to molecular orbitals.
    !! New scratch files are created to store them.
    !! In addition, txt files can be created.
    implicit none
    integer, intent(in) :: n_ab
    !! Number of unique MO sets (1 for closed-shell, 2 for open-shell)
    integer, intent(in) :: n_ao
    !! Number of atomic orbitals
    integer, dimension(:), intent(in) :: n_mos
    !! Number of alpha/beta molecular orbitals
    real(real64), dimension(:,:,:), intent(in) :: c_ia
    !! Coefficients of MOs in AOs basis
    logical, optional, intent(in) :: txt_fmt
    !! if present and True, save a formatted version of the file as well

    integer(int32), parameter :: nfiles = 8, ncols_txt = 5
    integer(int32), dimension(nfiles), parameter :: &
        LPx = [3, 6, 3, 3, 3, 3, 3, 3]
        ! LPx = [3, 6, 3, 3, 3, 3, 3, 3, 3, 3]
    integer :: i, iab, ifile, iorb, irec0, iu_fi, iu_fo, iu_txt, ix, llab, &
        lrec, N1, N2
    real(real64), dimension(n_ao,n_ao) :: q_ao
    real(real64), dimension(:, :), allocatable :: q_mo
    real(real64), dimension(:, :), allocatable :: tmp_arr
    character(len=2) :: qlab
    character(len=1), dimension(3), target :: lab_xyz = ['X', 'Y', 'Z']
    character(len=2), dimension(6), target :: &
        lab_2Dxyz = ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']
    character(len=:), dimension(:), pointer :: lab => null()
    character(len=5), dimension(nfiles), parameter :: fnames_ao = [ &
        'P.SCR', &  !  1 P  dipole
        'Q.SCR', &  !  2 Q  quadrupole
        'L.SCR', &  !  3 L  magnetic r x grad / 2
        'V.SCR', &  !  4 V  gradient
        'W.SCR', &  !  5 W  magnetic (r-rj) x grad / 2
        'T.SCR', &  !  6 T  magnetic (r-rij) x grad / 2
        ! 'U.SCR', &  !  7 U  magnetic (r-rj) x grad / 2 * "pop" ab
        ! 'C.SCR', &  !  8 C  magnetic (r-rij) x grad / 2* "pop" ab
        'D.SCR', &  !  9 D  GIAO like r x rj / 2
        'E.SCR'  &  ! 10 E  GIAO like r x ri / 2
        ]
    character(len=7) :: new_file
    character(len=80) :: fmt_h, fmt_d

    if (txt_fmt) then
        write(fmt_h, '("(6x,",i0,"(i14,:))")') ncols_txt
        write(fmt_d, '("(i6,",i0,"(e14.6,:))")') ncols_txt
    end if

    allocate(tmp_arr(n_ao,maxval(n_mos)))
    do iab = 1, n_ab
        n_mo = n_mos(iab)
        allocate(q_mo(n_mo,n_mo))
        do ifile = 1, nfiles
            if (iab == 1) then
                llab = 1
                qlab = fnames_ao(ifile)(1:1)
            else
                llab = 2
                qlab = fnames_ao(ifile)(1:1) // 'B'
            end if
            new_file = qlab(:llab) // 'MO' // fnames_ao(ifile)(2:)
            write(iu_out, '(" Transforming ",a," to ",a,".")') &
                fnames_ao(ifile), new_file
            inquire(iolength=lrec) q_ao(:,1)
            open(newunit=iu_fi, file=fnames_ao(ifile), access='direct', &
                action='read', recl=lrec)
            inquire(iolength=lrec) q_mo(:,1)
            open(newunit=iu_fo, file=new_file, access='direct', &
                action='write', recl=lrec)
            if (LPx(ifile) == 3) then
                lab => lab_xyz
            else
                lab => lab_2Dxyz
            end if

            do ix = 1, LPx(ifile)
                irec0 = (ix-1)*n_ao
                do iorb = 1, n_ao
                    read(iu_fi, rec=irec0+iorb) (q_ao(i,iorb),i=1,n_ao)
                end do
                call convert_AO2MO(n_ao, n_mo, c_ia(:,:,iab), q_ao, q_mo, &
                                   tmp_arr)
                irec0 = (ix-1)*n_mos(iab)
                do iorb = 1, n_mo
                    write(iu_fo, rec=irec0+iorb) (q_mo(i,iorb),i=1,n_mo)
                end do
                if (txt_fmt) then
                    open(newunit=iu_txt, &
                         file=qlab(:llab)//lab(ix)//'.MO.SCR.txt')
                    write(iu_txt, *) qlab(:llab)//lab(ix)//' MO'
                    write(iu_txt, *) ' '
                    N1 = 1
                    do
                        N2 = min(N1+ncols_txt-1, n_mo)
                        write(iu_txt, fmt_h) (i, i=N1,N2)
                        do iorb = 1, n_mo
                            write(iu_txt, fmt_d) iorb, (q_mo(iorb,i),i=N1,N2)
                        end do
                        if (N2 == n_mo) exit
                        N1 = N1 + 5
                    end do
                end if
            end do
            close(iu_fi)
            close(iu_fo)
        end do
        deallocate(q_mo)
    end do
end subroutine build_MOs

! ======================================================================

subroutine write_control(iout, err)
    !! SOS: writes control output
    !!
    !! This subroutine reproduces the behavior of SOS [writecon]
    !! with ndiff = 0
    implicit none
    integer, intent(in) :: iout
    !! unit for output
    class(BaseException), allocatable, intent(out) :: err
    !! Error instance

    integer :: ix, N, nc, ncols
    integer :: i, ia, iao, iprim, j, k, n_mo_1e, n_mo_2e
    integer, dimension(:), allocatable :: atlist
    character(len=:), dimension(:), allocatable :: fcomp
    character(len=10), dimension(:), allocatable :: fcoords
    type(PhysFact) :: phys
    class(BaseException), allocatable :: suberr

    err = InitError()

    ! compute the number of doubly occupied orbitals
    ! number of singly occupied orbitals = (multip-1)
    n_mo_1e = multip - 1
    n_mo_2e = (n_el - n_mo_1e)/2
    if (n_mo_2e*2 + n_mo_1e /= n_el) then
        call RaiseError(err, &
                        'Error in defining the number of occupied orbitals')
        return
    end if
    write(iout, *)  'ROAAI control output'
    write(iout, *)  'MOLECULAR PARAMETERS'
    write(iout, *)  '------------------------------------------'
    write(iout, *) 'Number of basis functions           :    ', n_basis
    write(iout, *) 'Number of electrons                 :    ', n_el
    write(iout, *) 'Charge                              :    ', charge
    write(iout, *) 'Multiplicity                        :    ', multip
    write(iout, *) 'Number of doubly occupied orbitals  :    ', n_mo_2e
    write(iout, *) 'Number of singly occupied a-orbitals:    ', n_els(1)-n_mo_2e
    write(iout, *) 'Number of singly occupied b-orbitals:    ', n_els(n_ab)-n_mo_2e
    write(iout, *) 'Number of atoms                     :    ', n_at
    write(iout, *) 'Sum of atomic charges               :    ', sum(int(at_chg))
    write(iout, *)
    write(iout, *)

    1000 format(i3,' atom, ',i4,' primitive functions:',/, &
                '  shell prim          alpha      expansion    P-expansion', &
                '   AO range')
    1002 format(1x,i5,a2,i5,3f15.6,:,i4,' - ',i4)
    1003 format(1x,i5,a2,i5,2f15.6,:,15x,i4,' - ',i4)
    ! The list are for a control on the shells
    allocate(atlist(n_ao), fcoords(n_ao))

    iao = 1
    do ia = 1, n_at
        write(iout, 1000) ia, nprim_per_at(ia)
        do iprim = 1, nprim_per_at(ia)
            associate (bs => bset_info(ia,iprim))
            if (bs%shell_last) then
                if(allocated(fcomp)) deallocate(fcomp)
                fcomp = shell_lmxyz(bs%shelltype, bs%pure, suberr)
                if (suberr%raised()) then
                    select type(suberr)
                        class is (ArgumentError)
                            call RaiseError(err, 'Unsupported shell type')
                            return
                        class default
                            call RaiseError(err, &
                                            'Unknown error from shell_lmxyz')
                            return
                    end select
                end if
                N = size(fcomp)
                fcoords(iao:iao+N-1) = fcomp
                atlist(iao:iao+N-1) = ia
                if (bs%shelltype /= 'SP') then
                    write(iout, 1003) bs%shellid, bs%shelltype, iprim, &
                        bs%alpha, bs%coeff(0), iao, iao+N-1
                else
                    write(iout, 1002) bs%shellid, bs%shelltype, iprim, &
                        bs%alpha, bs%coeff(0), bs%coeff(-1), iao, iao+N-1
                end if
                iao = iao + N
            else
                if (bs%shelltype /= 'SP') then
                    write(iout, 1003) bs%shellid, bs%shelltype, iprim, &
                        bs%alpha, bs%coeff(0)
                else
                    write(iout, 1002) bs%shellid, bs%shelltype, iprim, &
                        bs%alpha, bs%coeff(0), bs%coeff(-1)
                end if
            end if
            end associate
        end do
    end do

    write(iout, *) ' AO    type       atom'
    do iao = 1, n_mos(1)
        write(iout, '(i6,1x,a10,i6)') iao, fcoords(iao), atlist(iao)
    end do
    deallocate(fcoords, atlist)

    ncols = 5
    1010 format(12x,5(:,f13.6))
    1011 format(12x,5(:,i13))
    1012 format(i12,5(:,f13.6))
    write(iout, *) 'Eigenvalues'
    write(iout, *) 'Alpha'
    do i = 1, n_mos(1), ncols
        nc = min(i+ncols-1,n_mos(1))
        write(iout, 1010) (en_mos(j,1), j=i,nc)
        if (.false.) then ! Add coefficients
            write(iout, 1011) (j, j=i,nc)
            do j = 1, n_mos(1)
                write(iout, 1012) j, (coef_mos(k,j,1), k=i,nc)
            end do
        end if
    end do
    write(iout, *) 'Beta'
    do i = 1, n_mos(n_ab), ncols
        nc = min(i+ncols-1,n_mos(n_ab))
        write(iout, 1010) (en_mos(j,n_ab), j=i,nc)
        if (.false.) then ! Add coefficients
            write(iout, 1011) (j, j=i,nc)
            do j = 1, n_mos(n_ab)
                write(iout, 1012) j, (coef_mos(k,j,n_ab), k=i,nc)
            end do
        end if
    end do
    
    1020 format(1x,a2,i6,3f16.6,f10.5)
    write(iout, *) 'Geometry (Ang):'
    do ia = 1, n_at
        write(iout, 1020) at_lab(ia), ia, &
            (phys%bohr2Ang(at_crd(ix,ia)),ix=1,3), at_chg(ia)
    end do
    
    return
end subroutine write_control

! ======================================================================

end module gmcd_legacy
