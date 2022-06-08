! get H in det space
! get wave function with added/removed electron
! calculate lanczos tridiag
! square beta vector
! with a supplied omega range and ground state energy/epsilon, calculate z
! 


BEGIN_PROVIDER [double precision, greens_A_real, (greens_omega_N)]
    implicit none

    double precision        :: alpha(lanczos_N), bbeta(lanczos_N), epsilon
    double precision        :: psi_coef_excited(N_det, N_states)
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N)
    integer                 :: i

    ! psi_coef, h_matrix_all_dets, n_det, lanczos_N, psi_energy

    ! TODO: find the actual a^{dagger} operator
    call build_A_wavefunction(1,1,psi_det,psi_coef_excited)

    call lanczos_tridiag_reortho_r(h_matrix_all_dets, psi_coef_excited&
                                   alpha, bbeta,&
                                   lanczos_N, N_det)

    bbeta(1) = 1.0
    do i = 2, lanczos_N
        bbeta(i) = bbeta(i)**2.0
    end do

    epsilon = 0.001 ! small, a limit is taken here
    z = psi_energy(1) + omega + (0.0, 1.0)*epsilon

    ! this could be embarrisingly split across processes
    ! TODO: refresh on OMP directives
    do i = 1, greens_omega_N
        zalpha = z(i) - alpha
        greens_A_real(i) = cfraction_c((0.0, 0.0), bbeta, zalpha, lanczos_N)
    end do

END_PROVIDER

subroutine build_A_wavefunction(i_particle, ispin, coef_out, det_out)
    implicit none
    BEGIN_DOC
    ! Applies the single excitation operator : a^{dager}_(i_particle) a_(i_hole) of
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
        call get_phase(psi_det(1,1,k),det_out(1,1,k),phase,N_int)
        coef_out(k,:) = phase * coef_out(k,:)
      else
        coef_out(k,:) = 0.d0
        det_out(:,:,k) = psi_det(:,:,k)
      endif
    enddo
end

subroutine add_electron(key_in,,i_particle,ispin,i_ok)
    implicit none
    BEGIN_DOC
    ! Apply the single excitation operator : a^{dager}_(i_particle) a_(i_hole) of spin = ispin
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

! subroutine remove_electron()

! end