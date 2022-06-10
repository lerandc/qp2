BEGIN_PROVIDER [double precision, greens_omega, (greens_omega_N)]
    implicit none 
    integer                 :: i

    ! linearly spaced data for now
    do i = 1, greens_omega_N
        greens_omega(i) = (greens_omega_min*(greens_omega_N - i) &
                        +  greens_omega_max*(i-1) ) / (greens_omega_N-1)
    end do

END_PROVIDER


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   Real Implementations   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

BEGIN_PROVIDER [complex*16, greens_A, (greens_omega_N)]
    implicit none

    double precision        :: alpha(lanczos_N), beta(lanczos_N), epsilon
    double precision        :: psi_coef_excited(N_det, N_states), E0
    integer(bit_kind)       :: det_copy(N_int,2,N_det)
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N), bbeta(lanczos_N), cfraction_c
    integer                 :: i, i_ok

    call build_A_wavefunction(elec_alpha_num+1,1,psi_coef_excited,det_copy)

    call lanczos_tridiag_reortho_r(h_matrix_all_dets, psi_coef_excited,&
                                   alpha, beta,&
                                   lanczos_N, N_det)

    bbeta(1) = (1.0, 0.0)
    do i = 2, lanczos_N
        bbeta(i) = beta(i)**2.0
    end do

    
    epsilon = 0.001 ! small, a limit is taken here
    E0 = psi_energy(1)
    z = E0 + (greens_omega + (0.0, 1.0)*epsilon)
    
    ! this could be embarrisingly split across processes
    ! TODO: refresh on OMP directives
    do i = 1, greens_omega_N
        zalpha = z(i) - alpha
        greens_A(i) = cfraction_c((0.0, 0.0), bbeta, zalpha, lanczos_N)
    end do


END_PROVIDER

BEGIN_PROVIDER [complex*16, greens_R, (greens_omega_N)]
    implicit none

    double precision        :: alpha(lanczos_N), beta(lanczos_N), epsilon
    double precision        :: psi_coef_excited(N_det, N_states), E0
    integer(bit_kind)       :: det_copy(N_int,2,N_det)
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N), bbeta(lanczos_N), cfraction_c
    integer                 :: i

    call build_R_wavefunction(1,1,psi_coef_excited,det_copy)

    call lanczos_tridiag_reortho_r(h_matrix_all_dets, psi_coef_excited,&
                                   alpha, beta,&
                                   lanczos_N, N_det)

    bbeta(1) = (1.0, 0.0)
    do i = 2, lanczos_N
        bbeta(i) = beta(i)**2.0
    end do

    epsilon = 0.001 ! small, a limit is taken here
    E0 = psi_energy(1)
    z = E0 - (greens_omega + (0.0, 1.0)*epsilon)

    ! this could be embarrisingly split across processes
    ! TODO: refresh on OMP directives
    do i = 1, greens_omega_N
        zalpha = z(i) - alpha
        greens_R(i) = -1.0*cfraction_c((0.0, 0.0), bbeta, zalpha, lanczos_N)
    end do

END_PROVIDER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Complex Implementations   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

BEGIN_PROVIDER [complex*16, greens_A_complex, (greens_omega_N)]
    implicit none

    double precision        :: alpha(lanczos_N), beta(lanczos_N), epsilon, E0
    complex*16              :: psi_coef_excited(N_det, N_states) 
    integer(bit_kind)       :: det_copy(N_int,2,N_det)
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N), bbeta(lanczos_N), cfraction_c
    integer                 :: i

    call build_A_wavefunction_complex(elec_alpha_num+1,1,psi_coef_excited,det_copy)

    call lanczos_tridiag_reortho_c(h_matrix_all_dets_complex, psi_coef_excited,&
                                   alpha, beta,&
                                   lanczos_N, N_det)

    bbeta(1) = 1.0
    do i = 2, lanczos_N
        bbeta(i) = beta(i)**2.0
    end do

    epsilon = 0.001 ! small, a limit is taken here
    E0 = psi_energy(1)
    z = E0 + (greens_omega + (0.0, 1.0)*epsilon)

    ! this could be embarrisingly split across processes
    ! TODO: refresh on OMP directives
    do i = 1, greens_omega_N
        zalpha = z(i) - alpha
        greens_A_complex(i) = cfraction_c((0.0, 0.0), bbeta, zalpha, lanczos_N)
    end do

END_PROVIDER

BEGIN_PROVIDER [complex*16, greens_R_complex, (greens_omega_N)]
    implicit none

    double precision        :: alpha(lanczos_N), beta(lanczos_N), epsilon, E0
    complex*16              :: psi_coef_excited(N_det, N_states) 
    integer(bit_kind)       :: det_copy(N_int,2,N_det)
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N), bbeta(lanczos_N), cfraction_c
    integer                 :: i

    call build_R_wavefunction_complex(1,1,psi_coef_excited,det_copy)

    call lanczos_tridiag_reortho_c(h_matrix_all_dets_complex, psi_coef_excited,&
                                   alpha, beta,&
                                   lanczos_N, N_det)

    bbeta(1) = 1.0
    do i = 2, lanczos_N
        bbeta(i) = beta(i)**2.0
    end do

    epsilon = 0.001 ! small, a limit is taken here
    E0 = psi_energy(1)
    z = E0 - (greens_omega + (0.0, 1.0)*epsilon)

    
    ! this could be embarrisingly split across processes
    ! TODO: refresh on OMP directives
    do i = 1, greens_omega_N
        zalpha = z(i) - alpha
        greens_R_complex(i) = -1.0*cfraction_c((0.0, 0.0), bbeta, zalpha, lanczos_N)
    end do

END_PROVIDER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  UTILITIES !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! reduced versions of excitation routines below
subroutine build_A_wavefunction(i_particle, ispin, coef_out, det_out)
    implicit none
    BEGIN_DOC
    ! Applies the creation operator : a^{dager}_(i_particle) of
    ! spin = ispin to the current wave function (psi_det, psi_coef)
    END_DOC
    integer, intent(in)            :: i_particle,ispin
    integer(bit_kind), intent(out) :: det_out(N_int,2,N_det)
    double precision, intent(out)  :: coef_out(N_det,N_states)
  
    integer :: k
    integer :: i_ok
    double precision :: phase
    do k=1,N_det
      coef_out(k,:) = psi_coef(k,:)
      det_out(:,:,k) = psi_det(:,:,k)
      call add_electron(det_out(1,1,k),i_particle,ispin,i_ok)
      if (i_ok == 1) then
        call get_phase_ca(psi_det(1,1,k),i_particle,ispin,phase)
        coef_out(k,:) = phase * coef_out(k,:)
      else
        coef_out(k,:) = 0.d0
        det_out(:,:,k) = psi_det(:,:,k)
      endif
    enddo
end

subroutine build_R_wavefunction(i_hole, ispin, coef_out, det_out)
    implicit none
    BEGIN_DOC
    ! Applies the annihilation operator: a_(i_hole) of
    ! spin = ispin to the current wave function (psi_det, psi_coef)
    END_DOC
    integer, intent(in)            :: i_hole,ispin
    integer(bit_kind), intent(out) :: det_out(N_int,2,N_det)
    double precision, intent(out)  :: coef_out(N_det,N_states)
  
    integer :: k
    integer :: i_ok
    double precision :: phase
    do k=1,N_det
      coef_out(k,:) = psi_coef(k,:)
      det_out(:,:,k) = psi_det(:,:,k)
      call remove_electron(det_out(1,1,k),i_hole,ispin,i_ok)
      if (i_ok == 1) then
        call get_phase_ca(psi_det(1,1,k),i_hole,ispin,phase)
        coef_out(k,:) = phase * coef_out(k,:)
      else
        coef_out(k,:) = 0.d0
        det_out(:,:,k) = psi_det(:,:,k)
      endif
    enddo
end

subroutine build_A_wavefunction_complex(i_particle, ispin, coef_out, det_out)
    implicit none
    BEGIN_DOC
    ! Applies the creation operator : a^{dager}_(i_particle) of
    ! spin = ispin to the current wave function (psi_det, psi_coef_complex)
    END_DOC
    integer, intent(in)            :: i_particle,ispin
    integer(bit_kind), intent(out) :: det_out(N_int,2,N_det)
    complex*16, intent(out)       :: coef_out(N_det,N_states)
  
    integer :: k
    integer :: i_ok
    double precision :: phase
    do k=1,N_det
      coef_out(k,:) = psi_coef_complex(k,:)
      det_out(:,:,k) = psi_det(:,:,k)
      call add_electron(det_out(1,1,k),i_particle,ispin,i_ok)
      if (i_ok == 1) then
        call get_phase_ca(psi_det(1,1,k),i_particle,ispin,phase)
        coef_out(k,:) = phase * coef_out(k,:)
      else
        coef_out(k,:) = (0.d0, 0.d0)
        det_out(:,:,k) = psi_det(:,:,k)
      endif
    enddo
end

subroutine build_R_wavefunction_complex(i_hole, ispin, coef_out, det_out)
    implicit none
    BEGIN_DOC
    ! Applies the annihilation operator: a_(i_hole) of
    ! spin = ispin to the current wave function (psi_det, psi_coef_complex)
    END_DOC
    integer, intent(in)            :: i_hole,ispin
    integer(bit_kind), intent(out) :: det_out(N_int,2,N_det)
    complex*16, intent(out)       :: coef_out(N_det,N_states)
  
    integer :: k
    integer :: i_ok
    double precision :: phase
    do k=1,N_det
      coef_out(k,:) = psi_coef_complex(k,:)
      det_out(:,:,k) = psi_det(:,:,k)
      call remove_electron(det_out(1,1,k),i_hole,ispin,i_ok)
      if (i_ok == 1) then
        call get_phase_ca(psi_det(1,1,k),i_hole,ispin,phase)
        coef_out(k,:) = phase * coef_out(k,:)
      else
        coef_out(k,:) = (0.d0, 0.d0)
        det_out(:,:,k) = psi_det(:,:,k)
      endif
    enddo
end

subroutine add_electron(key_in,i_particle,ispin,i_ok)
    implicit none
    BEGIN_DOC
    ! Apply the creation : a^{dager}_(i_particle) of spin = ispin
    ! on key_in
    ! ispin = 1  == alpha
    ! ispin = 2  == beta
    ! i_ok = 1  == the excitation is possible
    ! i_ok = -1 == the excitation is not possible
    END_DOC
    integer, intent(in)            :: i_particle,ispin
    integer(bit_kind), intent(inout) :: key_in(N_int,2)
    integer, intent(out)           :: i_ok
    integer                        :: k,j,i
    integer(bit_kind)              :: mask
    use bitmasks

    ASSERT (i_particle <= mo_num)
    i_ok = 1
    
    ! particle
    k = shiftr(i_particle-1,bit_kind_shift)+1
    j = i_particle-shiftl(k-1,bit_kind_shift)-1
    mask = ibset(0_bit_kind,j)
    if (iand(key_in(k,ispin),mask) == 0_bit_kind) then
        key_in(k,ispin) = ibset(key_in(k,ispin),j)
    else 
        i_ok= -1
        return
    end if

end

subroutine remove_electron(key_in,i_hole,ispin,i_ok)
    implicit none
    BEGIN_DOC
    ! Apply the annihilation operator : a_(i_hole) of spin = ispin
    ! on key_in
    ! ispin = 1  == alpha
    ! ispin = 2  == beta
    ! i_ok = 1  == the excitation is possible
    ! i_ok = -1 == the excitation is not possible
    END_DOC
    integer, intent(in)            :: i_hole,ispin
    integer(bit_kind), intent(inout) :: key_in(N_int,2)
    integer, intent(out)           :: i_ok
    integer                        :: k,j,i
    integer(bit_kind)              :: mask
    use bitmasks
    ASSERT (i_hole > 0 )
    i_ok = 1
    
    ! hole
    k = shiftr(i_hole-1,bit_kind_shift)+1
    j = i_hole-shiftl(k-1,bit_kind_shift)-1
    mask = ibset(0_bit_kind,j)
    ! check whether position j is occupied
    if (iand(key_in(k,ispin),mask) /= 0_bit_kind) then
        key_in(k,ispin) = ibclr(key_in(k,ispin),j)
    else 
        i_ok= -1
        return
    end if
  
  end


subroutine get_phase_ca(det, iorb, ispin, phase)
    implicit none
    BEGIN_DOC
    ! Calculate relative phase of determinant after creation/annihilation
    ! of ispin at iorb
    END_DOC

    integer, intent(in)            :: iorb, ispin
    integer                        :: N_perm, i,j,k
    integer(bit_kind), intent(in)  :: det(N_int,2)
    integer(bit_kind)              :: mask
    double precision, intent(out)  :: phase
    double precision, parameter    :: phase_dble(0:1) = (/ 1.d0, -1.d0 /)
    ! If we add an alpha electron, parity is not affected by beta part of determinant
    ! (only need number of alpha occupied orbs below iorb)
    
    ! If we add a beta electron, the parity is affected by alpha part
    ! (need total number of occupied alpha orbs (all of which come before beta)
    ! and total number of beta occupied orbs below iorb)

    mask = 2**(iorb - 1) - 1
    k = shiftr(iorb-1,bit_kind_shift)+1

    ! top of the order
    N_perm = popcnt(iand(mask, det(k,ispin)))

    ! all necessary orbitals occupied below
    do i = 1, ispin
        do j = 1, k-1
            N_perm = N_perm + popcnt(det(j,i))
        end do
    end do

    phase = phase_dble(iand(N_perm, 1))
end