module exc_sos
    !! Module fot the excited (EXC) sum-over-state scheme.
    !!
    !! Module providing the necessary tools to compute transition
    !!   moments using the EXC sum-over-states scheme.
    !!
    !! @note
    !!   References:
    !!
    !!   1. P. Bour, Chem. Phys. Lett. 1997, 265, 65-70
    !!   2. P. Bour, Chem. Phys. Lett. 1998, 288, 363-370
    !! @endnote
    !!
    use iso_fortran_env, only: real64

    use exception, only: ArgumentError, BaseException, InitError, &
        RaiseArgError, RaiseError

    implicit none

contains

! ======================================================================

function sos_eiOg(ldim, n_ab, n_mos, n_els, t_mo, O_ij, forbid, model) &
    result(res)
    !! Compute < e_i | O | g > using the SOS formalism
    !!
    !! Computes the ground-to-excited integral < e_i | O | g > using the
    !!   summation over molecular orbitals, where O is an arbitrary
    !!   operator.
    !! The function expects all quantities wrt molecular orbitals.
    !! If forbid is true, the transition is forbidden in standard
    !!   theory, for instance singlet-triplet transitions.
    integer, intent(in) :: ldim
    !! Leading dimension / number of components of O
    integer, intent(in) :: n_ab
    !! Number of unique MO sets (1 for closed-shell, 2 for open-shell)
    integer, dimension(:), intent(in) :: n_mos
    !! Number of molecular orbitals.
    integer, dimension(:), intent(in) :: n_els
    !! Number of electrons in the molecular orbitals.
    real(real64), dimension(:,:,:), intent(in) :: t_mo
    !! Electronic transition amplitudes
    real(real64), dimension(:,:,:,:), intent(in) :: O_ij
    !! MO-integrals for the quantity/property of interest
    logical, intent(in), optional :: forbid
    !! Transition is forbidden (transition moment = 0).
    integer, intent(in), optional :: model
    !! SOS model
    !! 1: use pure TD-DFT amplitudes
    !! 2: use pure Slater determinants
    !! 3: hybrid, all TD-DFT amplitudes + Slater permutations
    real(real64), dimension(:), allocatable :: res
    !! Integral < e_j | O | g >

    integer :: ia, iab, ib

    allocate(res(ldim))

    res = 0.0_real64
    if (present(forbid)) then
        if (forbid) return
    end if
    if (model == 2) then
        do iab = 1, n_ab
            do ia = 1, n_els(iab)
                do ib = n_els(iab)+1, n_mos(iab)
                    ! res = res + t_mo(ia,ib,iab)*O_ij(:,ib,ia,iab)
                    res = res + t_mo(ia,ib,iab)*O_ij(:,ia,ib,iab)
                end do
            end do
        end do
    else
        do iab = 1, n_ab
            do ia = 1, n_mos(iab)
                do ib = 1, n_mos(iab)
                    ! res = res + t_mo(ia,ib,iab)*O_ij(:,ib,ia,iab)
                    res = res + t_mo(ia,ib,iab)*O_ij(:,ia,ib,iab)
                end do
            end do
        end do
    end if
    if (n_ab == 1) res = res*sqrt(2.0_real64)

end function sos_eiOg

! ======================================================================

function sos_ejOei(ldim, n_ab, n_mos, n_els, t_mo, prefac, forbid, model) &
    result(res)
    !! Compute < e_j | O | e_i > using the SOS formalism
    !!
    !! Computes the excited-to-excited integral < e_j | O | e_i > using
    !!   the summation over molecular orbitals.
    !! The function expects all quantities wrt molecular orbitals.
    !! If forbid is true, the transition is forbidden in standard
    !!   theory, for instance singlet-triplet transitions.
    integer, intent(in) :: ldim
    !! Leading dimension / number of components of O
    integer, intent(in) :: n_ab
    !! Number of unique MO sets (1 for closed-shell, 2 for open-shell)
    integer, dimension(:), intent(in) :: n_mos
    !! Number of molecular orbitals.
    integer, dimension(:), intent(in) :: n_els
    !! Number of electrons in the molecular orbitals.
    real(real64), dimension(:,:,:), intent(in) :: t_mo
    !! Electronic transition amplitudes
    real(real64), dimension(:,:,:,:), intent(in) :: prefac
    !! Prefactor, computed by [sos_prefac_ejOei]
    logical, intent(in), optional :: forbid
    !! Transition is forbidden (transition moment = 0).
    integer, intent(in), optional :: model
    !! SOS model
    !! 1: use pure TD-DFT amplitudes
    !! 2: use pure Slater determinants
    !! 3: hybrid, all TD-DFT amplitudes + Slater permutations
    real(real64), dimension(:), allocatable :: res
    !! Integral < e_j | O | e_i >

    integer :: ia, iab, ib

    allocate(res(ldim))

    res = 0.0_real64
    if (present(forbid)) then
        if (forbid) return
    end if
    if (model == 2) then
        do iab = 1, n_ab
            do ia = 1, n_els(iab)
                do ib = n_els(iab)+1, n_mos(iab)
                    res = res + t_mo(ia,ib,iab)*prefac(:,ia,ib,iab)
                end do
            end do
        end do
    else
        do iab = 1, n_ab
            do ia = 1, n_mos(iab)
                do ib = 1, n_mos(iab)
                    res = res + t_mo(ia,ib,iab)*prefac(:,ia,ib,iab)
                end do
            end do
        end do
    end if

end function sos_ejOei

! ======================================================================

function sos_prefac_ejOei(ldim, n_ab, n_mos, n_els, t_mo, O_gg, O_ij, &
                          model, add_ijaa) result(prefac)
    !! Prefactor to < e_j | O | e_i > for the SOS formalism
    !!
    !! Computes the prefactor term for < e_j | O | e_i > in the SOS
    !!   formalism to accelerate calculations.
    !! The function expects all quantities wrt molecular orbitals.
    integer, intent(in) :: ldim
    !! Leading dimension / number of components of O.
    integer, intent(in) :: n_ab
    !! Number of unique MO sets (1 for closed-shell, 2 for open-shell).
    integer, dimension(:), intent(in) :: n_mos
    !! Number of molecular orbitals.
    integer, dimension(:), intent(in) :: n_els
    !! Number of electrons in the molecular orbitals.
    real(real64), dimension(:,:,:), intent(in) :: t_mo
    !! Electronic transition amplitudes.
    real(real64), dimension(:,:), intent(in) :: O_gg
    !! Ground-state moment of the quantity/property.
    real(real64), dimension(:,:,:,:), intent(in) :: O_ij
    !! MO-integrals for the quantity/property of interest.
    integer, intent(in), optional :: model
    !! SOS model
    !! 1: use pure TD-DFT amplitudes
    !! 2: use pure Slater determinants
    !! 3: hybrid, all TD-DFT amplitudes + Slater permutations
    logical, intent(in) :: add_ijaa
    !! Include <ia|ja> terms.
    real(real64), dimension(:,:,:,:), allocatable :: prefac
    !! Returned prefactor.

    integer :: ia, iab, ib, ic, n_mo
    real(real64), dimension(:), allocatable :: x, u_gg

    n_mo = maxval(n_mos(:n_ab))
    allocate(prefac(ldim,n_mo,n_mo,n_ab), x(ldim))

    prefac = 0.0_real64

    if (n_ab == 1) then
        u_gg = O_gg(:,1)
    else
        u_gg = O_gg(:,1) + O_gg(:,2)
    end if

    if (model == 1) then
        !! Use pure TD-DFT transition amplitudes
        do iab = 1, n_ab
            !$omp parallel private(x)
            !$omp do collapse(2)
            do ia = 1, n_mos(iab)
                do ib = 1, n_mos(iab)
                    prefac(:,ia,ib,iab) = t_mo(ia,ib,iab) &
                        * (u_gg - O_ij(:,ia,ia,iab) + O_ij(:,ib,ib,iab))
                end do
            end do
            !$omp end do
            !$omp end parallel
        end do
    else if (model == 2) then
        !! Use pure Slater determinants
        if (add_ijaa) then
            do iab = 1, n_ab
                !$omp parallel private(x)
                !$omp do collapse(2)
                do ia = 1, n_els(iab)
                    do ib = n_els(iab)+1, n_mos(iab)
                        x = t_mo(ia,ib,iab) &
                            * (u_gg - O_ij(:,ia,ia,iab) + O_ij(:,ib,ib,iab))
                        do ic = 1, n_els(iab)
                            if (ic /= ia) &
                                x = x + t_mo(ic,ib,iab)*O_ij(:,ic,ia,iab)
                        end do
                        do ic = n_els(iab)+1, n_mos(iab)
                            if (ic /= ib) &
                                x = x + t_mo(ia,ic,iab)*O_ij(:,ic,ib,iab)
                        end do
                        prefac(:,ia,ib,iab) = x
                    end do
                end do
                !$omp end do
                !$omp end parallel
            end do
        else
            do iab = 1, n_ab
                !$omp parallel private(x)
                !$omp do collapse(2)
                do ia = 1, n_els(iab)
                    do ib = n_els(iab)+1, n_mos(iab)
                        x = t_mo(ia,ib,iab) &
                            * (u_gg - O_ij(:,ia,ia,iab) + O_ij(:,ib,ib,iab))
                        do ic = n_els(iab)+1, n_mos(iab)
                            if (ic /= ib) &
                                x = x + t_mo(ia,ic,iab)*O_ij(:,ic,ib,iab)
                        end do
                        prefac(:,ia,ib,iab) = x
                    end do
                end do
                !$omp end do
                !$omp end parallel
            end do
        end if
    else if (model == 3) then
        if (add_ijaa) then
            do iab = 1, n_ab
                !$omp parallel private(x)
                !$omp do collapse(2)
                do ia = 1, n_mos(iab)
                    do ib = 1, n_mos(iab)
                        x = t_mo(ia,ib,iab) &
                            * (u_gg - O_ij(:,ia,ia,iab) + O_ij(:,ib,ib,iab))
                        do ic = 1, n_els(iab)
                            if (ic /= ia) &
                                x = x + t_mo(ic,ib,iab)*O_ij(:,ic,ia,iab)
                        end do
                        do ic = n_els(iab)+1, n_mos(iab)
                            if (ic /= ib) &
                                x = x + t_mo(ia,ic,iab)*O_ij(:,ic,ib,iab)
                        end do
                        prefac(:,ia,ib,iab) = x
                    end do
                end do
                !$omp end do
                !$omp end parallel
            end do
        else
            do iab = 1, n_ab
                !$omp parallel private(x)
                !$omp do collapse(2)
                do ia = 1, n_mos(iab)
                    do ib = 1, n_mos(iab)
                        x = t_mo(ia,ib,iab) &
                            * (u_gg - O_ij(:,ia,ia,iab) + O_ij(:,ib,ib,iab))
                        do ic = n_els(iab)+1, n_mos(iab)
                            if (ic /= ib) &
                                x = x + t_mo(ia,ic,iab)*O_ij(:,ic,ib,iab)
                        end do
                        prefac(:,ia,ib,iab) = x
                    end do
                end do
                !$omp end do
                !$omp end parallel
            end do
        end if
    end if
    ! if (n_ab == 1) then
    !     ! closed shell
    !     !$omp parallel private(x)
    !     !$omp do collapse(2)
    !     do ia = 1, n_mo
    !         do ib = 1, n_mo
    !             x = t_mo(ia,ib,1) &
    !                 * (O_gg(:,1) - O_ij(:,ia,ia,1) + O_ij(:,ib,ib,1))
    !             do ic = n_els(1)+1, n_mo
    !                 if (ic /= ib) &
    !                     x = x + t_mo(ia,ic,1)*O_ij(:,ic,ib,1)
    !             end do
    !             prefac(:,ia,ib,1) = x
    !         end do
    !     end do
    !     !$omp end do
    !     !$omp end parallel
    ! else
    !     do iab = 1, n_ab
    !         !$omp parallel private(x)
    !         !$omp do collapse(2)
    !         do ia = 1, n_mos(iab)
    !             do ib = 1, n_mos(iab)
    !                 ! x = t_mo(ia,ib,iab) &
    !                 !     * (O_gg(:,1) + O_gg(:,2) - O_ij(:,ia,ia,iab) &
    !                 !        + O_ij(:,ib,ib,iab))
    !                 x = t_mo(ia,ib,iab) &
    !                     * (O_gg(:,iab) - O_ij(:,ia,ia,iab) &
    !                        + O_ij(:,ib,ib,iab))
    !                 do ic = 1, n_els(iab)
    !                     if (ic /= ia) &
    !                         x = x + t_mo(ic,ib,iab)*O_ij(:,ic,ia,iab)
    !                 end do
    !                 do ic = n_els(iab)+1, n_mos(iab)
    !                     if (ic /= ib) &
    !                         x = x + t_mo(ia,ic,iab)*O_ij(:,ic,ib,iab)
    !                 end do
    !                 prefac(:,ia,ib,iab) = x
    !             end do
    !         end do
    !         !$omp end do
    !         !$omp end parallel
    !     end do

    ! end if

end function sos_prefac_ejOei

! ======================================================================

subroutine sos_MCD_tensor_LORG_corr(n_states, n_els, fstate, kstate, r_gg, &
                                    r_ij, p_ij, ov_eieg, ov_eiej, G_gf)
    !! Apply a LORG-like correction to SOS MCD tensor for a given state
    !!
    !! The subroutine applies a LORG-like[^ref:LORG] correction to
    !!   remove the origin dependence of the MCD tensor G computed with
    !!   the SOS approach as described in Ref. [^ref:MCD_SOS_GIAO].
    !! The correction is applied for a specific state in the summation.
    !! More specifically, compute the G components in t2, t3, t5, t6 in
    !!   Eq. 4 of Ref. [^ref:MCD_SOS_GIAO].
    !!
    !! \begin{align}
    !!   T_2
    !!   &=
    !!   \mu_{fk} \frac{1}{2} \left [
    !!      \frac{1}{2 N_e} \left(
    !!          \sum_{l \neq g}
    !!          \frac{\mu_{kl} \times \nabla_{lg}}{E_l - E_g}
    !!      \right)
    !!   \right ]
    !!   \\ %%%%%%%%%%
    !!   T_3
    !!   &=
    !!   \mu_{fk} \frac{1}{2} \left [
    !!      \frac{1}{2 N_e} \left(
    !!          \sum_{l \neq k}
    !!          \frac{\mu_{lg} \times \nabla_{kl}}{E_k - E_l}
    !!      \right)
    !!   \right ]
    !!   \\ %%%%%%%%%%
    !!   T_5
    !!   &=
    !!   \mu_{kg} \frac{1}{2} \left [
    !!      \frac{1}{2 N_e} \left(
    !!          \sum_{l \neq k}
    !!          \frac{\mu_{fl} \times \nabla_{lk}}{E_k - E_l}
    !!      \right)
    !!   \right ]
    !!   \\ %%%%%%%%%%
    !!   T_6
    !!   &=
    !!   \mu_{kg} \frac{1}{2} \left [
    !!      \frac{1}{2 N_e} \left(
    !!          \sum_{l \neq f}
    !!          \frac{\mu_{lk} \times \nabla_{fl}}{E_l - E_f}
    !!      \right)
    !!   \right ]
    !! \end{align}
    !!
    !! @note
    !!
    !!   * The procedure must be called for each intermediate state
    !!     `kstate`.
    !!   * The initial state is assumed to be the ground state.
    !! @endnote
    !!
    !! @todo
    !!
    !!    why l=g excluded from T5 and T6 in original implementation?
    !!
    !! @endtodo
    !!
    !! [^ref:LORG]: A.E. Hansen, T.D. Bouman, _J. Chem. Phys._ 1985,
    !!         **82**, 5035 (https://doi.org/10.1063/1.448625);
    !!         A.E. Hansen, T.D. Bouman, _J. Chem. Phys._ 1986,
    !!         **84**, 2433 (https://doi.org/10.1063/1.450857).
    !! [^ref:MCD_SOS_GIAO]: P. Štětepánek, P. Bouř, _J. Comput. Chem._
    !!         2015, **36**, 723 (https://doi.org/10.1002/jcc.23845)
    use math, only: cross

    integer, intent(in) :: n_states
    !! Number of excited electronic states
    integer, dimension(2), intent(in) :: n_els
    !! Total number of electrons
    integer, intent(in) :: fstate
    !! Identifier of the final state
    integer, intent(in) :: kstate
    !! Identifier of the intermediate state
    real(real64), dimension(:), intent(in) :: r_gg
    !! Ground-state moment < g | r | g >
    real(real64), dimension(:,:,0:), intent(in) :: r_ij
    !! Integral between electronic states (0: ground) < e_i | r | e_j >
    real(real64), dimension(:,:,0:), intent(in) :: p_ij
    !! Integral between electronic states (0: ground) < e_i | p | e_j >
    real(real64), dimension(:), intent(in) :: ov_eieg
    !! 1/(e_i - e_g), inverse of electronic transition energy, in a.u.
    real(real64), dimension(0:,0:), intent(in) :: ov_eiej
    !! 1/(e_i - e_j), inverse of electronic energy difference, in a.u.
    real(real64), intent(inout) :: G_gf(3,3)
    !! MCD G tensor

    integer :: ix, jx, lstate
    real(real64) :: ov_2Ne
    real(real64), dimension(3) :: t2, t3, t5, t6

    t2 = 0.0_real64
    t3 = 0.0_real64
    t5 = 0.0_real64
    t6 = 0.0_real64

    ! Contrib. to T3, case l=g: <g|r|g> x <k|p|g> / (e_k - e_g)
    if (kstate /= 0) then
        t3 = t3 + cross(r_gg, p_ij(:,kstate,0)) * ov_eieg(kstate)
    end if

    do lstate = 1, n_states
        ! Contrib. to T2: <k|r|l> x <l|p|g> / (e_l - e_g)
        if (kstate /= 0) then
            t2 = t2 + cross(r_ij(:,lstate,kstate), p_ij(:,lstate,0)) &
                * ov_eieg(lstate)
        end if
        if (lstate /= kstate) then
            ! Contrib. to T3: <l|r|g> x <k|p|l> / (e_k - e_l)
            if (kstate /= 0) then
                t3 = t3 + cross(r_ij(:,lstate,0), p_ij(:,kstate,lstate)) &
                    * ov_eiej(kstate,lstate)
            end if
            ! Contrib. to T5: <f|r|l> x <l|p|k> / (e_k - e_l)
            if (kstate /= fstate) then
                t5 = t5 + cross(r_ij(:,fstate,lstate), p_ij(:,lstate,kstate)) &
                    * ov_eiej(kstate,lstate)
            end if
        end if
        ! Contrib. to T6: <l|r|k> x <f|p|l> / (e_l - e_f)
        if (lstate /= fstate) then
            t6 = t6 + cross(r_ij(:,lstate,kstate), p_ij(:,fstate,lstate)) &
                * ov_eiej(lstate,fstate)
        end if
    end do

    ov_2Ne = 1.0_real64/(2.0_real64*(n_els(1)+n_els(2)))
    if (kstate == 0) then
        do ix = 1, 3
            do jx = 1, 3
                G_gf(jx,ix) = G_gf(jx,ix) - ov_2Ne*( &
                    + (t5(ix)-t6(ix))*r_gg(jx))
            end do
        end do
    else if (kstate == fstate) then
        do ix = 1, 3
            do jx = 1, 3
                G_gf(jx,ix) = G_gf(jx,ix) - ov_2Ne*( &
                    (t2(ix)-t3(ix))*r_ij(jx,fstate,fstate))
            end do
        end do
    else
        do ix = 1, 3
            do jx = 1, 3
                G_gf(jx,ix) = G_gf(jx,ix) - ov_2Ne*( &
                    (t2(ix)-t3(ix))*r_ij(jx,fstate,kstate) &
                    + (t5(ix)-t6(ix))*r_ij(jx,kstate,0))
            end do
        end do
    end if

end subroutine sos_MCD_tensor_LORG_corr

! ======================================================================

end module exc_sos
