module vibronic_td

    use numeric, only: realwp, f0, f1, f2, fhalf, fquart, f8th, f16th

    implicit none

contains

! ======================================================================

function get_a_final(n_vib, time, temperature, omega_f) result(a_f)
    !! Construct the final-state **a** matrix for vibronic TD.
    !!
    !! Constructs the diagonal elements of the final-state **a** matrix,
    !! as defined in
    !! [DOI:10.1021/ct400450k](https://dx.doi.org/10.1021/ct400450k).
    !!
    !! The formula for \(a_f\) is:
    !! \[ a_f = \frac{\omega_f}{\sinh(\omega_f \tau)} \]
    !! where \(\tau = it\)

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    real(realwp), intent(in) :: time
    !! Time variable (in atomic units).
    real(realwp), intent(in) :: omega_f(:)
    !! Vector of final-state vibrational frequencies (in atomic units).
    real(realwp), intent(in) :: temperature
    !! Temperature in Kelvin

    ! Local
    complex(realwp), dimension(n_vib) :: a_f
    !! Diagonal elements of the a matrix for the final state.

    if (temperature /= f0) stop 'Temperature is not zero'

    a_f = omega_f / sinh(omega_f * cmplx(f0, time, realwp))

end function get_a_final

! ======================================================================

function get_a_initial(n_vib, time, temperature, omega_i) result(a_i)
    !! Construct the initial-state **a** matrix for vibronic TD.
    !!
    !! Constructs the diagonal elements of the initial-state **a** matrix,
    !! as defined in
    !! [DOI:10.1021/ct400450k](https://dx.doi.org/10.1021/ct400450k).
    !!
    !! The formula for \(a_i\) is:
    !! \[ a_i = 2 \omega_i \]
    !!
    !! @todo
    !! Temperature effects are currently not supported.
    !! @endtodo

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    real(realwp), intent(in) :: time
    !! Time variable (in atomic units).
    real(realwp), intent(in) :: omega_i(:)
    !! Vector of initial-state vibrational frequencies (in atomic units).
    real(realwp), intent(in) :: temperature
    !! Temperature in Kelvin.
    complex(realwp), dimension(n_vib) :: a_i
    !! Diagonal elements of the **a** matrix for the initial state.

    if (temperature /= f0) stop 'Temperature is not zero'

    a_i = f2 * omega_i

end function get_a_initial

! ======================================================================

function get_chi00(n_vib, time, omega_i, omega_f, K_vec, a_i, a_f, C_inv, &
                    D_inv, v, eta, phase) result(chi00)
    !! Construct the vibronic TD \(\chi_{00}\).
    !!
    !! Constructs \(\chi_{00}\) as:
    !! \[
    !!      \chi_{00}
    !!      =
    !!      \sqrt{\det(C^{-1}D^{-1} (a_i a_f))} \exp\left(
    !!          \sum\limits_{k} \frac{1}{2}v_k\eta_k
    !!          - \sum\limits_{k} K_k \omega_{i_k} K_k
    !!          + \frac{i}{2} \sum\limits_{k} \omega_{f_k} t
    !!          \right)
    !! \]
    !!
    !! @note
    !! - requires the phase of the previous time step to ensure continuity
    !! - the phase is updated
    !! - the zero point energy of the final state in the last term is included
    !! @endnote
    use math, only: det
    use numeric, only: pi

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    real(realwp), intent(in) :: time
    !! Time variable (in atomic units).
    real(realwp), intent(in) :: omega_i(:)
    !! Vector of initial-state vibrational frequencies (in atomic units).
    real(realwp), intent(in) :: omega_f(:)
    !! Vector of final-state vibrational frequencies (in atomic units).
    real(realwp), intent(in) :: K_vec(:)
    !! Duschinsky shift vector ***K*** (in atomic units).
    complex(realwp), intent(in) :: v(:)
    !! Vibronic TD ***v*** vector.
    complex(realwp), intent(in) :: eta(:)
    !! Vibronic TD ***&eta;*** vector.
    complex(realwp), intent(in) :: a_i(:)
    !! Diagonal of the initial-state **a** matrix.
    complex(realwp), intent(in) :: a_f(:)
    !! Diagonal of the final-state **a** matrix.
    complex(realwp), intent(in) :: C_inv(:, :)
    !! Inverse of the vibronic TD **C** matrix.
    complex(realwp), intent(in) :: D_inv(:, :)
    !! Inverse of the vibronic TD **D** matrix.
    real(realwp), intent(inout) :: phase
    !! Phase of &chi;00 in the previous time step.
    complex(realwp) :: chi00
    !! Calculated value of &chi;00.

    ! Local
    integer :: i
    real(realwp) :: new_phase
    complex(realwp) :: tmp
    complex(realwp), dimension(n_vib, n_vib) :: tmp_mat

    tmp = - sum(K_vec(:) * omega_i(:) * K_vec(:))
    tmp = tmp + fhalf * sum(v(:) * eta(:))
    tmp_mat = matmul(C_inv, D_inv)

    do i = 1, n_vib
        tmp_mat(:, i) = a_i(:) * a_f(:) * tmp_mat(:, i)
    end do

    chi00 = det(tmp_mat)
    chi00 = chi00 * exp(f2 * tmp + cmplx(f0, sum(omega_f) * time, realwp))
    chi00 = sqrt(chi00)
    chi00 = chi00 * cmplx(f0, f1, realwp)

    ! to make sure the phase is continuous
    new_phase = atan2(chi00%im, chi00%re)
    do while(abs(new_phase - phase) > fhalf * pi)
        chi00 = -chi00
        if (new_phase - phase > f0) then
            new_phase = new_phase - pi
        else
            new_phase = new_phase + pi
        end if
    end do
    phase = new_phase

end function get_chi00

! ======================================================================

function get_chi10(n_vib, chi00, eta) result(chi10)
    !! Construct the vibronic TD \(\chi_{10}\).
    !!
    !! Constructs \(\chi_{10}\) as:
    !! \[ \chi_{10} = -0.5 \eta \chi_{00} \]

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    complex(realwp) :: chi00
    !! &chi;00 value.
    complex(realwp) :: eta(:)
    !! ***&eta;*** vector.
    complex(realwp) :: chi10(n_vib)
    !! &chi;10 vector.

    chi10 = -fhalf * chi00 * eta

end function get_chi10

! ======================================================================

function get_chi01(n_vib, chi00, eta) result(chi01)
    !! Construct the vibronic TD \(\chi_{01}\).
    !!
    !! Constructs \(\chi_{01}\) as:
    !! \[ \chi_{10} = -0.5 \eta \chi_{00} \]

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    complex(realwp) :: chi00
    !! &chi;00 value.
    complex(realwp) :: eta(:)
    !! Vibronic TD ***&eta;*** vector.
    complex(realwp) :: chi01(n_vib)
    !! &chi;01 vector.

    chi01 = -fhalf * chi00 * eta

end function get_chi01

! ======================================================================

function get_chi11(n_vib, chi00, eta, C_inv, D_inv) result(chi11)
    !! Construct the vibronic TD \(\chi_{11}\).
    !!
    !! Constructs \(\chi_{11}\) as:
    !! \[
    !!      \chi_{11}
    !!      =
    !!      (-0.5 C^{-1}_{ij} + 0.5 D^{-1}_{ij} + 0.25 \eta_i \eta_j )
    !!      \chi_{00}
    !! \]

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    complex(realwp) :: chi00
    !! &chi;00 value.
    complex(realwp) :: eta(:)
    !! Vibronic TD ***&eta;*** vector.
    complex(realwp) :: C_inv(:, :)
    !! Inverse of the vibronic TD **C** matrix.
    complex(realwp) :: D_inv(:, :)
    !! Inverse of the vibronic TD **D** matrix.
    complex(realwp) :: chi11(n_vib, n_vib)
    !! &chi;11 matrix.

    ! Local
    integer :: i, j

    chi11 = (-fhalf * C_inv + fhalf * D_inv) * chi00
    do i = 1, n_vib
        do j = 1, n_vib
            chi11(i, j) = chi11(i, j) + fquart * eta(i) * eta(j) * chi00
        end do
    end do

end function get_chi11

! ======================================================================

function get_chi20(n_vib, chi00, eta, C_inv, D_inv) result(chi20)
    !! Construct the vibronic TD \(\chi_{20}\).
    !!
    !! Constructs \(\chi_{20}\) as:
    !! \[
    !!      \chi_{20}
    !!      =
    !!      (0.5 C^{-1}_{ij} + 0.5 D^{-1}_{ij} + 0.25 \eta_i \eta_j )
    !!      \chi_{00}
    !! \]

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    complex(realwp) :: chi00
    !! &chi;00 value.
    complex(realwp) :: eta(:)
    !! Vibronic TD ***&eta;*** vector.
    complex(realwp) :: C_inv(:, :)
    !! Inverse of the vibronic TD **C** matrix.
    complex(realwp) :: D_inv(:, :)
    !! Inverse of the vibronic TD **D** matrix.
    complex(realwp) :: chi20(n_vib, n_vib)
    !! &chi;20 matrix.

    ! Local
    integer :: i, j

    chi20 = (fhalf * C_inv + fhalf * D_inv) * chi00

    do i = 1, n_vib
        do j = 1, n_vib
            chi20(i, j) = chi20(i, j) + fquart * eta(i) * eta(j) * chi00
        end do
    end do

end function get_chi20

! ======================================================================

function get_chi21(n_vib, chi00, eta, C_inv, D_inv) result(chi21)
    !! Construct the vibronic TD \(\chi_{21}\).
    !!
    !! Constructs \(\chi_{21}\) as:
    !! \[
    !!      \chi_{21}
    !!      =
    !!      (-\frac{1}{8}\eta_i\eta_j\eta_k
    !!       - \frac{1}{4}P(\eta D^{-1})
    !!       + \frac{1}{4}P(\eta C^{-1}))
    !!      \chi_{00}
    !! \]
    !! where \(P(X)\) permute the indices.

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    complex(realwp) :: chi00
    !! &chi;00 value.
    complex(realwp) :: eta(:)
    !! Vibronic TD ***&eta;*** vector.
    complex(realwp) :: C_inv(:, :)
    !! Inverse of the vibronic TD **C** matrix.
    complex(realwp) :: D_inv(:, :)
    !! Inverse of the vibronic TD **D** matrix.
    complex(realwp) :: chi21(n_vib, n_vib, n_vib)
    !! &chi;21 tensor.

    ! Local
    integer :: i, j, k

    do i = 1, n_vib
        do j = 1, n_vib
            do k = 1, n_vib
                chi21(i, j, k) = chi00 * ( &
                    - f8th * eta(i) * eta(j) * eta(k) &
                    + fquart * ( &
                        - eta(i) * D_inv(j, k) - eta(j) * D_inv(i, k) &
                        - eta(k) * D_inv(i, j) + eta(i) * C_inv(j, k) &
                        + eta(j) * C_inv(i, k) - eta(k) * C_inv(i, j) &
                    ))
            end do
        end do
    end do

end function get_chi21

! ======================================================================

function get_chi30(n_vib, chi00, eta, C_inv, D_inv) result(chi30)
    !! Construct the vibronic TD \(\chi_{30}\).
    !!
    !! Constructs \(\chi_{30}\) as:
    !! \[
    !!      \chi_{30}
    !!      =
    !!      (-\frac{1}{8}\eta_i\eta_j\eta_k
    !!       - \frac{1}{4}P(\eta D^{-1})
    !!       - \frac{1}{4}P(\eta C^{-1}))
    !!      \chi_{00}
    !! \]
    !! where \(P(X)\) permute the indices.

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    complex(realwp) :: chi00
    !! &chi;00 value.
    complex(realwp) :: eta(:)
    !! Vibronic TD ***&eta;*** vector.
    complex(realwp) :: C_inv(:, :)
    !! Inverse of the vibronic TD **C** matrix.
    complex(realwp) :: D_inv(:, :)
    !! Inverse of the vibronic TD **D** matrix.
    complex(realwp) :: chi30(n_vib, n_vib, n_vib)
    !! &chi;30 tensor.

    ! Local
    integer :: i, j, k

    do i = 1, n_vib
        do j = 1, n_vib
            do k = 1, n_vib
                chi30(i, j, k) = chi00 * ( &
                    - f8th * eta(i) * eta(j) * eta(k) + fquart * ( &
                        - eta(i) * D_inv(j, k) - eta(j) * D_inv(i, k) &
                        - eta(k) * D_inv(i, j) - eta(i) * C_inv(j, k) &
                        - eta(j) * C_inv(i, k) - eta(k) * C_inv(i, j) &
                    ))
            end do
        end do
    end do

end function get_chi30

! ======================================================================

function get_chi31(n_vib, chi00, eta, C_inv, D_inv) result(chi31)
    !! Construct the vibronic TD \(\chi_{31}\).
    !!
    !! Constructs \(\chi_{31}\) as:
    !! \[
    !!      \chi_{31}
    !!      =
    !!      (\frac{1}{16}\eta_i\eta_j\eta_k\eta_l
    !!       + \frac{1}{4}D^{-1}_{ij}D^{-1}_{kl}
    !!       + \frac{1}{4}D^{-1}_{ik}D^{-1}_{lj}
    !!       + \frac{1}{4}D^{-1}_{il}D^{-1}_{jk}
    !!       - \frac{1}{4}C^{-1}_{ij}C^{-1}_{kl}
    !!       - \frac{1}{4}C^{-1}_{ik}C^{-1}_{lj}
    !!       - \frac{1}{4}C^{-1}_{il}C^{-1}_{jk}
    !!       + \frac{1}{4}P(D^{-1}C^{-1})
    !!       - \frac{1}{8}P(\eta \eta C^{-1})
    !!       + \frac{1}{8}P(\eta \eta D^{-1}))
    !!      \chi_{00}
    !! \]
    !! where \(P(X)\) permute the indices.

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    complex(realwp) :: chi00
    !! &chi;00 value.
    complex(realwp) :: eta(:)
    !! Vibronic TD ***&eta;*** vector.
    complex(realwp) :: C_inv(:, :)
    !! Inverse of the vibronic TD **C** matrix.
    complex(realwp) :: D_inv(:, :)
    !! Inverse of the vibronic TD **D** matrix.
    complex(realwp) :: chi31(n_vib, n_vib, n_vib, n_vib)
    !! &chi;31 tensor.

    ! Local
    integer :: i, j, k, l

    do i = 1, n_vib
        do j = 1, n_vib
            do k = 1, n_vib
                do l = 1, n_vib
                    ! eta^4
                    chi31(i, j, k, l) = chi00 * ( &
                        f16th * eta(i) * eta(j) * eta(k) * eta(l) &
                        ! Z_b^4
                        + fquart * D_inv(i, j) * D_inv(k, l) &
                        + fquart * D_inv(i, k) * D_inv(l, j) &
                        + fquart * D_inv(i, l) * D_inv(j, k) &
                        ! U^4
                        - fquart * C_inv(i, j) * C_inv(k, l) &
                        - fquart * C_inv(i, k) * C_inv(l, j) &
                        - fquart * C_inv(i, l) * C_inv(j, k) &
                        ! Zb^2 U^2
                        - fquart * D_inv(i, j) * C_inv(k, l) &
                        - fquart * D_inv(i, k) * C_inv(j, l) &
                        - fquart * D_inv(k, j) * C_inv(i, l) &
                        + fquart * D_inv(i, l) * C_inv(k, j) &
                        + fquart * D_inv(l, j) * C_inv(k, i) &
                        + fquart * D_inv(k, l) * C_inv(i, j) &
                        ! eta^2 U^2
                        - f8th * eta(i) * eta(j) * C_inv(k, l) &
                        - f8th * eta(i) * eta(k) * C_inv(j, l) &
                        - f8th * eta(k) * eta(j) * C_inv(i, l) &
                        + f8th * eta(i) * eta(l) * C_inv(k, j) &
                        + f8th * eta(l) * eta(j) * C_inv(k, i) &
                        + f8th * eta(k) * eta(l) * C_inv(i, j) &
                        ! eta^2 Z_b^2
                        + f8th * eta(i) * eta(j) * D_inv(k, l) &
                        + f8th * eta(i) * eta(k) * D_inv(j, l) &
                        + f8th * eta(k) * eta(j) * D_inv(i, l) &
                        + f8th * eta(i) * eta(l) * D_inv(k, j) &
                        + f8th * eta(l) * eta(j) * D_inv(k, i) &
                        + f8th * eta(k) * eta(l) * D_inv(i, j) &
                        )
                end do
            end do
        end do
    end do

end function get_chi31

! ======================================================================

function get_C_inv(n_vib, time, temperature, omega_i, omega_f, J_mat) &
    result(C_inv)
    !! Construct the inverse of vibronic TD **C** matrix.
    !!
    !! Constructs the inverse of the C matrix,
    !! \[ C = c_f + J^T c_i J \]
    !! where \(c = \frac{\omega}{\tanh(\frac{1}{2}\omega \tau)}\)
    !! and \(\tau = it\)
    !!
    !! @todo
    !! Temperature effects are currently not supported.
    !! @endtodo
    use math, only: inv_mat

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    real(realwp), intent(in) :: time
    !! Time variable (in atomic units).
    real(realwp), intent(in) :: omega_i(:)
    !! Vector of initial-state vibrational frequencies (in atomic units).
    real(realwp), intent(in) :: omega_f(:)
    !! Vector of final-state vibrational frequencies (in atomic units).
    real(realwp), intent(in) :: temperature
    !! Temperature in Kelvin.
    real(realwp), intent(in) :: J_mat(:, :)
    !! Duschinsky matrix **J**.
    complex(realwp), dimension(n_vib, n_vib) :: C_inv
    !! Inverse of the vibronic TD **C** matrix.

    ! Local
    integer :: i, j
    complex(realwp) :: tmp
    complex(realwp), dimension(n_vib, n_vib) :: C

    if (temperature /= f0) stop 'Temperature is not zero'

    C = (f0, f0)
    do i = 1, n_vib
        tmp = tanh(fhalf * omega_f(i) * cmplx(f0, time, kind=realwp))
        C(i, i) = omega_f(i) / tmp
    end do
    do i = 1, n_vib
        do j = 1, n_vib
            tmp = sum(J_mat(:, i) * omega_i(:) * J_mat(:, j))
            C(i, j) = C(i, j) + tmp
        end do
    end do

    call inv_mat(C, C_inv)

end function get_C_inv

! ======================================================================

function get_D_inv(n_vib, time, temperature, omega_i, omega_f, J_mat) result(D_inv)
    !! Construct the inverse of vibronic TD **D** matrix.
    !!
    !! Constructs the inverse of the D matrix,
    !! \[ D = d_f + J^T d_i J \]
    !! where \(d = \omega \tanh(\frac{1}{2}\omega \tau)\)
    !! and \(\tau = it\)
    !!
    !! @todo
    !! Temperature effects are currently not supported.
    !! @endtodo
    use math, only: inv_mat

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    real(realwp), intent(in) :: time
    !! Time variable (in atomic units).
    real(realwp), intent(in) :: temperature
    !! Temperature in Kelvin.
    real(realwp), intent(in) :: omega_i(:)
    !! Vector of initial-state vibrational frequencies (in atomic units).
    real(realwp), intent(in) :: omega_f(:)
    !! Vector of final-state vibrational frequencies (in atomic units).
    real(realwp), intent(in) :: J_mat(:, :)
    !! Duschinsky matrix **J**.
    complex(realwp), dimension(n_vib, n_vib) :: D_inv
    !! Inverse of the vibronic TD **D** matrix.

    ! Local
    integer :: i, j
    complex(realwp) :: tmp
    complex(realwp), dimension(n_vib, n_vib) :: D

    if (temperature /= f0) stop 'Temperature is not zero'

    D = (f0, f0)
    do i = 1, n_vib
        tmp = tanh(fhalf * omega_f(i) * cmplx(f0, time, kind=realwp))
        D(i, i) = omega_f(i) * tmp
    end do

    do i = 1, n_vib
        do j = 1, n_vib
            tmp = sum(J_mat(:, i) * omega_i(:) * J_mat(:, j))
            D(i, j) = D(i, j) + tmp
        end do
    end do

    call inv_mat(D, D_inv)

end function get_D_inv

! ======================================================================

function get_eta(n_vib, v, D_inv) result(eta)
    !! Construct the vibronic TD **&eta;** vector.
    !!
    !! Constructs the \(\eta\) vector, defined as:
    !! \[ \eta = 2 D^{-1} v \]

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes.
    complex(realwp), intent(in) :: D_inv(:, :)
    !! Inverse of the vibronic TD **D** matrix.
    complex(realwp), intent(in) :: v(:)
    !! Vibronic TD ***v*** vector.
    complex(realwp), dimension(n_vib) :: eta
    !! Vibronic TD ***&eta;*** vector.

    eta = f2 * matmul(D_inv, v)

end function get_eta

! ======================================================================

function get_v(n_vib, omega_i, K_vec, J_mat) result(v)
    !! Construct the vibronic TD **v** vector.
    !!
    !! Constructs the \(v\) vector, defined as:
    !! \[ v = J^T \omega_i K \]

    ! Arguments
    integer, intent(in) :: n_vib
    !! Number of vibrational modes
    real(realwp), intent(in) :: omega_i(:)
    !! Vector of initial-state vibrational frequencies (in atomic units).
    real(realwp), intent(in) :: K_vec(:)
    !! Duschinsky shift vector K (in atomic units).
    real(realwp), intent(in) :: J_mat(:, :)
    !! Duschinsky matrix J.
    complex(realwp), dimension(n_vib) :: v
    !! Vibronic TD ***v*** vector.

    ! Local
    integer :: i

    do i = 1, n_vib
        v(i) = sum(J_mat(:, i) * omega_i(:) * K_vec(:))
    end do

end function get_v

! ======================================================================

end module vibronic_td
