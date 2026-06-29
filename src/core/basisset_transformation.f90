module basisset_transformation
    !! Basis sets-related modules
    !!
    !! Procedure related to basis sets and their definition
    use iso_fortran_env, only: real64
    use datatypes, only: PrimitiveFunction
    use basisset, only: transfo_cart2pure
    use blas_drv, only: xgemm
    use numeric
    
    implicit none

contains

! ======================================================================

    subroutine convert_matrix_from_pure_to_cart(bset_info, nprim_per_atom, matrix_pure, matrix_cart)
        !! Convert a given matrix from pure to cartesian basis.
        implicit none

        type(PrimitiveFunction), dimension(:,:), intent(in) :: bset_info
        integer, dimension(:), intent(in) :: nprim_per_atom
        real(real64), dimension(:,:), intent(in) :: matrix_pure
        real(real64), dimension(:,:), intent(out) :: matrix_cart
        
        real(real64), dimension(:,:), allocatable :: mu_t, nu_t
        real(real64), dimension(:,:), allocatable :: tmp
        
        integer :: center_mu, shell_mu, center_nu, shell_nu
        integer :: l_mu, l_nu
        integer :: mu_cart_length, nu_cart_length
        integer :: mu_pure_length, nu_pure_length
        integer :: mu_idx_pure, nu_idx_pure
        integer :: mu_idx_cart, nu_idx_cart

        integer, dimension(:), allocatable :: len_shell_per_atom_mu, len_shell_per_atom_nu
        integer :: shell_per_atom_mu, shell_per_atom_nu
        integer :: offset_mu, offset_nu

        allocate(tmp(size(matrix_cart, 1), size(matrix_pure, 1)))

        mu_idx_pure = 1
        mu_idx_cart = 1
        
        do center_mu = 1, size(bset_info, 1)
            
            offset_mu = 1
            shell_per_atom_mu = calculate_shell_per_atom(bset_info(center_mu, :),nprim_per_atom,center_mu)

            allocate(len_shell_per_atom_mu(shell_per_atom_mu))    
            len_shell_per_atom_mu = calculate_len_shell_per_atom(bset_info(center_mu, :),nprim_per_atom,&
                                                                 center_mu, shell_per_atom_mu)
            
            do shell_mu = 1, shell_per_atom_mu

                l_mu = bset_info(center_mu, offset_mu)%l
                
                if (l_mu == -1) then
                    mu_cart_length = 4
                    mu_pure_length = 4
                else
                    mu_cart_length = (l_mu + 1) * (l_mu + 2) / 2
                    mu_pure_length = 2 * l_mu + 1
                endif

                mu_t = transfo_cart2pure(l_mu)

                offset_mu = offset_mu + len_shell_per_atom_mu(shell_mu)
                nu_idx_pure = 1
                nu_idx_cart = 1
               
                do center_nu = 1, size(bset_info, 1)

                    offset_nu = 1
                    shell_per_atom_nu = calculate_shell_per_atom(bset_info(center_nu, :),nprim_per_atom,center_nu)

                    allocate(len_shell_per_atom_nu(shell_per_atom_nu))    
                    len_shell_per_atom_nu = calculate_len_shell_per_atom(bset_info(center_nu, :),nprim_per_atom,&
                                                                         center_nu, shell_per_atom_nu)

                    do shell_nu = 1, shell_per_atom_nu

                        l_nu = bset_info(center_nu, offset_nu)%l

                        if (l_nu == -1) then
                            nu_cart_length = 4
                            nu_pure_length = 4
                        else
                            nu_cart_length = (l_nu + 1) * (l_nu + 2) / 2
                            nu_pure_length = 2 * l_nu + 1
                        endif

                        nu_t = transfo_cart2pure(l_nu)

                        call xgemm("T", "N",                                                    &
                                   nu_cart_length,                                              &
                                   mu_pure_length,                                              & 
                                   nu_pure_length,                                              &
                                   f1,                                                          &
                                   nu_t,                                                        &
                                   nu_pure_length,                                              &
                                   matrix_pure(nu_idx_pure : nu_idx_pure + nu_pure_length - 1,  &
                                               mu_idx_pure : mu_idx_pure + mu_pure_length - 1), &
                                   nu_pure_length,                                              &
                                   f0,                                                          &
                                   tmp,                                                         &
                                   nu_cart_length)

                        call xgemm("N", "N",                                                    &
                                   nu_cart_length,                                              &
                                   mu_cart_length,                                              &
                                   mu_pure_length,                                              &
                                   f1,                                                          &
                                   tmp,                                                         &
                                   nu_cart_length,                                              &
                                   mu_t,                                                        &
                                   mu_pure_length,                                              &
                                   f0,                                                          &
                                   matrix_cart(nu_idx_cart : nu_idx_cart + nu_cart_length - 1,  &
                                               mu_idx_cart : mu_idx_cart + mu_cart_length - 1), &
                                   nu_cart_length)

                        offset_nu = offset_nu + len_shell_per_atom_nu(shell_nu)

                        nu_idx_pure = nu_idx_pure + nu_pure_length
                        nu_idx_cart = nu_idx_cart + nu_cart_length 

                    end do

                    deallocate(len_shell_per_atom_nu)

                end do
            
                mu_idx_pure = mu_idx_pure + mu_pure_length
                mu_idx_cart = mu_idx_cart + mu_cart_length

            end do

            deallocate(len_shell_per_atom_mu)

        end do

        deallocate(tmp)

    end subroutine convert_matrix_from_pure_to_cart
        
    function calculate_shell_per_atom(bset_info, nprim_per_atom, ia) result (shell_per_atom)

        implicit none
        
        type(PrimitiveFunction), dimension(:), intent(in) :: bset_info
        integer, dimension(:), intent(in) :: nprim_per_atom
        integer, intent(in) :: ia
        integer :: shell_per_atom

        integer :: i

        shell_per_atom = 0

        do i = 1, nprim_per_atom(ia)

            if (bset_info(i)%shell_first .eqv. .true.) then
                shell_per_atom = shell_per_atom + 1
            endif

        enddo

    end function calculate_shell_per_atom

    function calculate_len_shell_per_atom(bset_info, nprim_per_atom, ia, n_shells_per_atom) result (len_shell_per_atom)

        implicit none
        
        type(PrimitiveFunction), dimension(:), intent(in) :: bset_info
        integer, dimension(:), intent(in) :: nprim_per_atom
        integer, intent(in) :: ia
        integer, intent(in) :: n_shells_per_atom
        
        integer, dimension(n_shells_per_atom) :: len_shell_per_atom

        integer :: i, counter, offset

        counter = 1
        offset = 1

        do i = 1, nprim_per_atom(ia)

            if (bset_info(i)%shell_last .eqv. .false.) then
               counter = counter + 1
            else
               len_shell_per_atom(offset) = counter
               counter = 1
               offset = offset + 1
            endif

        enddo

    end function calculate_len_shell_per_atom

    function calculate_n_ao_cart(bset, nprim_per_atom) result(n_ao_cartesian)

        implicit none
        
        type(PrimitiveFunction), dimension(:,:), intent(in) :: bset
        integer, dimension(:), intent(in) :: nprim_per_atom
        integer :: n_ao_cartesian, n_d, n_f

        integer :: i, center

        n_ao_cartesian = 0
        n_d = 0
        n_f = 0

        do center = 1, size(bset, 1)
            do i = 1, nprim_per_atom(center)
                
                if (bset(center, i)%shell_first .eqv. .true.) then

                    if (bset(center, i)%shelltype .eq. 'SP' ) then

                        n_ao_cartesian = n_ao_cartesian + 4

                    elseif (bset(center, i)%shelltype .eq. 'S' ) then

                        n_ao_cartesian = n_ao_cartesian + 1

                    elseif (bset(center, i)%shelltype .eq. 'P' ) then

                        n_ao_cartesian = n_ao_cartesian + 3

                    elseif (bset(center, i)%shelltype .eq. 'D' ) then

                        n_ao_cartesian = n_ao_cartesian + 6

                    elseif (bset(center, i)%shelltype .eq. 'F' ) then

                        n_ao_cartesian = n_ao_cartesian + 10

                    endif
                endif

            enddo
        enddo

    end function calculate_n_ao_cart

    subroutine construct_bset_info_cart(bset_info, nprim_per_atom, bset_info_cart)

        implicit none
        
        type(PrimitiveFunction), dimension(:,:), intent(in) :: bset_info
        integer, dimension(:), intent(in) :: nprim_per_atom
        type(PrimitiveFunction), dimension(:,:), intent(out) :: bset_info_cart

        integer :: center, iprim

        do center = 1, size(bset_info, 1)
            do iprim = 1, nprim_per_atom(center)

                    bset_info_cart(center, iprim)%l = bset_info(center, iprim)%l
                    bset_info_cart(center, iprim)%shelltype = bset_info(center, iprim)%shelltype
                    bset_info_cart(center, iprim)%shellid = bset_info(center, iprim)%shellid
                    bset_info_cart(center, iprim)%pure  = .false.

                    bset_info_cart(center, iprim)%shell_first = bset_info(center, iprim)%shell_first
                    bset_info_cart(center, iprim)%shell_last = bset_info(center, iprim)%shell_last
                    bset_info_cart(center, iprim)%alpha = bset_info(center, iprim)%alpha
                    bset_info_cart(center, iprim)%coeff = bset_info(center, iprim)%coeff

                if (bset_info(center, iprim)%shelltype .eq. 'D') then

                    bset_info_cart(center, iprim)%ndim = 6

                else if (bset_info(center, iprim)%shelltype .eq. 'F') then

                    bset_info_cart(center, iprim)%ndim = 10

                else
                    bset_info_cart(center, iprim)%ndim = bset_info(center, iprim)%ndim
                endif

            enddo
        enddo

    end subroutine construct_bset_info_cart

! ======================================================================

end module basisset_transformation
