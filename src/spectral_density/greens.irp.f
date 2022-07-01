BEGIN_PROVIDER [integer, greens_omega_N]
    implicit none
    greens_omega_N = ceiling((greens_omega_max - greens_omega_min)/ greens_omega_resolution) + 1
END_PROVIDER

BEGIN_PROVIDER [double precision, greens_omega, (greens_omega_N)]
    implicit none 
    integer                 :: i

    ! linearly spaced data for now
    do i = 0, greens_omega_N-1
        greens_omega(i+1) = greens_omega_min + i * greens_omega_resolution
    end do
END_PROVIDER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!   Real Implementations   !!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

BEGIN_PROVIDER [complex*16, greens_A, (greens_omega_N, n_iorb_A, ns_dets)]
&BEGIN_PROVIDER[double precision, lanczos_alpha_A, (lanczos_N, n_iorb_A,ns_dets)]
&BEGIN_PROVIDER[double precision, lanczos_beta_A, (lanczos_N, n_iorb_A,ns_dets)]
    implicit none

    double precision        :: alpha(lanczos_N), beta(lanczos_N), epsilon
    double precision        :: E0, norm, dnrm2
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N), bbeta(lanczos_N), cfraction_c
    integer                 :: i, iorb, nnz, s_max_sze, i_n_det, N_det_l
    integer, allocatable    :: H_c(:), t_H_c(:), H_p(:)
    double precision, allocatable  ::  H_v(:), t_H_v(:), v(:)
    double precision, allocatable  :: psi_coef_excited(:,:) 
    integer(bit_kind), allocatable :: det_excited(:,:,:)

    greens_A = 0.0
    s_max_sze = max_row_sze_factor*nnz_max_per_row

    do i_n_det = 1, ns_dets

        N_det_l = n_det_sequence(i_n_det)
        allocate(psi_coef_excited(N_det_l, N_states), det_excited(N_int,2,N_det_l))

        do iorb = 1, n_iorb_A
            ! reset bit masks
            call set_ref_bitmask(iorb_A(iorb), 1, .false.)

            ! prepare input vector
            ! add electron to LUMO
            print *, '####### Building excited wavefunction for iorb ', iorb_A(iorb)
            call build_A_wavefunction(iorb_A(iorb),1,psi_coef_excited,det_excited,N_det_l)
            norm = dnrm2(N_det_l, psi_coef_excited(:,1), 1)
            psi_coef_excited = psi_coef_excited / norm

            ! prepare sparse Hamiltonian arrays
            print *, '####### Forming sparse arrays'
            allocate(H_c(s_max_sze), H_v(s_max_sze), H_p(N_det_l+1))
            call form_sparse_dH(H_p, H_c, H_v, s_max_sze, det_excited,&
                                iorb_A(iorb), 1, .false., N_det_l)

            nnz = H_p(N_det_l+1)-1

            print *, '### Moving allocations'
            allocate(t_H_c(nnz), t_H_v(nnz))
            t_H_c = H_c(:nnz)
            t_H_v = H_v(:nnz)

            call move_alloc(t_H_c, H_c)
            call move_alloc(t_H_v, H_v)

            print *, '####### Performing Lanczos tridiagonalization'
            call lanczos_tridiag_sparse_reortho_r(H_v, H_c, H_p, psi_coef_excited(:,1),&
                                        alpha, beta,&
                                        lanczos_N, nnz, N_det_l)

            
            print *, '####### Calculating Greens function'
            bbeta(1) = (1.d0, 0.d0)
            do i = 2, lanczos_N
                bbeta(i) = -1.d0*beta(i)**2.0
            end do

            lanczos_alpha_A(:, iorb, i_n_det) = alpha
            lanczos_beta_A(:, iorb, i_n_det) = real(bbeta)

            epsilon = greens_epsilon ! small, a limit is taken here
            E0 = psi_energy(1)
            z = E0 + (greens_omega + (0.d0, 1.d0)*epsilon)
            
            !$OMP PARALLEL DO PRIVATE(i, zalpha) SHARED(alpha, bbeta, lanczos_N, greens_A)&
            !$OMP SCHEDULE(GUIDED)
            do i = 1, greens_omega_N
                zalpha = z(i) - alpha
                greens_A(i, iorb, i_n_det) = cfraction_c((0.d0, 0.d0), bbeta, zalpha,lanczos_N)
            end do
            !$OMP END PARALLEL DO

            deallocate(H_c, H_v, H_p)
        end do

        deallocate(psi_coef_excited, det_excited)
    end do

END_PROVIDER

BEGIN_PROVIDER [complex*16, greens_R, (greens_omega_N, n_iorb_R, ns_dets)]
&BEGIN_PROVIDER[double precision, lanczos_alpha_R, (lanczos_N, n_iorb_R,ns_dets)]
&BEGIN_PROVIDER[double precision, lanczos_beta_R, (lanczos_N, n_iorb_R,ns_dets)]
    implicit none

    double precision        :: alpha(lanczos_N), beta(lanczos_N), epsilon
    double precision        :: E0, norm, dnrm2
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N), bbeta(lanczos_N), cfraction_c
    integer                 :: i, iorb, nnz, s_max_sze, i_n_det, N_det_l
    integer, allocatable    :: H_c(:), H_p(:), t_H_c(:)
    double precision , allocatable  ::  H_v(:), t_H_v(:)
    double precision, allocatable  :: psi_coef_excited(:,:) 
    integer(bit_kind), allocatable :: det_excited(:,:,:)

    greens_R = 0.0
    s_max_sze = max_row_sze_factor*nnz_max_per_row

    do i_n_det = 1, ns_dets

        N_det_l = n_det_sequence(i_n_det)
        allocate(psi_coef_excited(N_det_l, N_states), det_excited(N_int,2,N_det_l))

        do iorb = 1, n_iorb_R
            ! reset bit masks
            call set_ref_bitmask(iorb_R(iorb), 1, .true.)

            ! prepare input vector
            ! remove electron from HOMO
            print *, '####### Building excited wavefunction for iorb ', iorb_R(iorb)
            call build_R_wavefunction(iorb_R(iorb),1,psi_coef_excited,det_excited, N_det_l)
            norm = dnrm2(N_det_l, psi_coef_excited(:,1), 1)
            psi_coef_excited = psi_coef_excited / norm

            ! prepare sparse Hamiltonian arrays
            print *, '####### Forming sparse arrays'
            allocate(H_c(s_max_sze), H_v(s_max_sze), H_p(N_det_l+1))
            call form_sparse_dH(H_p, H_c, H_v, s_max_sze, det_excited,&
                                iorb_R(iorb), 1, .true., N_det_l)

            nnz = H_p(N_det_l+1)-1

            print *, '### Moving allocations'
            allocate(t_H_c(nnz), t_H_v(nnz))
            t_H_c = H_c(:nnz)
            t_H_v = H_v(:nnz)

            call move_alloc(t_H_c, H_c)
            call move_alloc(t_H_v, H_v)

            print *, '####### Performing Lanczos tridiagonalization'
            call lanczos_tridiag_sparse_reortho_r(H_v, H_c, H_p, psi_coef_excited(:,1),&
                                            alpha, beta,&
                                            lanczos_N, nnz, N_det_l)

            print *, '####### Calculating Greens function'
            bbeta(1) = (1.d0, 0.d0)
            do i = 2, lanczos_N
                bbeta(i) = -1.0*beta(i)**2.0
            end do

            lanczos_alpha_R(:, iorb, i_n_det) = alpha
            lanczos_beta_R(:, iorb, i_n_det) = real(bbeta)

            epsilon = greens_epsilon ! small, a limit is taken here
            E0 = psi_energy(1)
            z = E0 - (-1.d0*greens_omega + (0.d0, 1.d0)*epsilon) ! omega is abs. energy value, hole part is only positive in negative energy range

            !$OMP PARALLEL DO PRIVATE(i, zalpha) SHARED(alpha, bbeta, lanczos_N, greens_R)&
            !$OMP SCHEDULE(GUIDED)
            do i = 1, greens_omega_N
                zalpha = z(i) - alpha
                greens_R(i, iorb, i_n_det) = -1.0*cfraction_c((0.d0, 0.d0), bbeta, zalpha, lanczos_N)
            end do
            !$OMP END PARALLEL DO

            deallocate(H_c, H_v, H_p)
        end do

        deallocate(psi_coef_excited, det_excited)
    end do

END_PROVIDER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Complex Implementations   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

BEGIN_PROVIDER [complex*16, greens_A_complex, (greens_omega_N, n_iorb_A)]
    implicit none

    double precision        :: alpha(lanczos_N), beta(lanczos_N), epsilon, E0, norm, dznrm2
    complex*16              :: psi_coef_excited(N_det, N_states) 
    integer(bit_kind)       :: det_excited(N_int,2,N_det)
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N), bbeta(lanczos_N), cfraction_c
    integer                 :: i, iorb, nnz, s_max_sze
    integer, allocatable    :: H_c(:), H_p(:), t_H_c(:)
    complex*16 , allocatable  ::  H_v(:), t_H_v(:)


    s_max_sze =1000*nnz_max_per_row
    ! prepare sparse Hamiltonian arrays
    print *, '####### Forming sparse arrays'
    allocate(H_c(s_max_sze), H_v(s_max_sze), H_p(N_det+1))
    call form_sparse_zH(H_p, H_c, H_v, s_max_sze)

    nnz = H_p(N_det+1)-1
    
    print *, '### Moving allocations'
    allocate(t_H_c(nnz), t_H_v(nnz))
    t_H_c = H_c(:nnz)
    t_H_v = H_v(:nnz)
    call move_alloc(t_H_c, H_c)
    call move_alloc(t_H_v, H_v)

    ! prepare input vector
    ! add electron to LUMO
    do iorb = 1, n_iorb_A
        print *, '####### Building excited wavefunction for iorb ', iorb_A(iorb)
        call build_A_wavefunction_complex(iorb_A(iorb),1,psi_coef_excited,det_excited)

        norm = dznrm2(N_det, psi_coef_excited(:,1), 1)
        psi_coef_excited = psi_coef_excited / norm

        print *, '####### Performing Lanczos tridiagonalization'
        call lanczos_tridiag_sparse_reortho_c(H_v, H_c, H_p, psi_coef_excited(:,1),&
                                    alpha, beta,&
                                    lanczos_N, nnz, N_det)

        print *, '####### Calculating Greens function'
        bbeta(1) = 1.0
        do i = 2, lanczos_N
            bbeta(i) = -1.0*beta(i)**2.0
        end do

        epsilon = greens_epsilon ! small, a limit is taken here
        E0 = psi_energy(1)
        z = E0 + (greens_omega + (0.0, 1.0)*epsilon)

        !$OMP PARALLEL DO PRIVATE(i, zalpha) SHARED(alpha, bbeta, lanczos_N, greens_A_complex)&
        !$OMP SCHEDULE(GUIDED)
        do i = 1, greens_omega_N
            zalpha = z(i) - alpha
            greens_A_complex(i, iorb) = cfraction_c((0.0, 0.0), bbeta, zalpha, lanczos_N)
        end do
        !$OMP END PARALLEL DO
    end do

END_PROVIDER

BEGIN_PROVIDER [complex*16, greens_R_complex, (greens_omega_N, n_iorb_R)]
    implicit none

    double precision        :: alpha(lanczos_N), beta(lanczos_N), epsilon, E0, norm, dznrm2
    complex*16              :: psi_coef_excited(N_det, N_states) 
    integer(bit_kind)       :: det_excited(N_int,2,N_det)
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N), bbeta(lanczos_N), cfraction_c
    integer                 :: i, iorb, nnz, s_max_sze
    integer, allocatable    :: H_c(:), H_p(:), t_H_c(:)
    complex*16 , allocatable  ::  H_v(:), t_H_v(:)

    s_max_sze =1000*nnz_max_per_row
    ! prepare sparse Hamiltonian arrays
    print *, '####### Forming sparse arrays'
    allocate(H_c(s_max_sze), H_v(s_max_sze), H_p(N_det+1))
    call form_sparse_zH(H_p, H_c, H_v, s_max_sze)

    nnz = H_p(N_det+1)-1

    print *, '### Moving allocations'
    allocate(t_H_c(nnz), t_H_v(nnz))
    t_H_c = H_c(:nnz)
    t_H_v = H_v(:nnz)

    call move_alloc(t_H_c, H_c)
    call move_alloc(t_H_v, H_v)

    ! prepare input vector
    ! remove electron from HOMO
    do iorb = 1, n_iorb_R
        print *, '####### Building excited wavefunction for iorb ', iorb_R(iorb)
        call build_R_wavefunction_complex(iorb_R(iorb),1,psi_coef_excited,det_excited)

        norm = dznrm2(N_det, psi_coef_excited(:,1), 1)
        psi_coef_excited = psi_coef_excited / norm

        print *, '####### Performing Lanczos tridiagonalization'
        call lanczos_tridiag_sparse_reortho_c(H_v, H_c, H_p, psi_coef_excited(:,1),&
                                    alpha, beta,&
                                    lanczos_N, nnz, N_det)

        print *, '####### Calculating Greens function'
        bbeta(1) = 1.0
        do i = 2, lanczos_N
            bbeta(i) = -1.0*beta(i)**2.0
        end do

        epsilon = greens_epsilon ! small, a limit is taken here
        E0 = psi_energy(1)
        z = E0 - (greens_omega + (0.0, 1.0)*epsilon)

        
        !$OMP PARALLEL DO PRIVATE(i, zalpha) SHARED(alpha, bbeta, lanczos_N, greens_R_complex)&
        !$OMP SCHEDULE(GUIDED)
        do i = 1, greens_omega_N
            zalpha = z(i) - alpha
            greens_R_complex(i, iorb) = -1.0*cfraction_c((0.0, 0.0),bbeta, zalpha,lanczos_N)
        end do
        !$OMP END PARALLEL DO
    end do

END_PROVIDER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  UTILITIES !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! reduced versions of excitation routines below
subroutine build_A_wavefunction(i_particle, ispin, coef_out, det_out, N_det_l)
    implicit none
    BEGIN_DOC
    ! Applies the creation operator : a^{dager}_(i_particle) of
    ! spin = ispin to the current wave function (psi_det, psi_coef)
    END_DOC
    integer, intent(in)            :: i_particle,ispin,N_det_l
    integer(bit_kind), intent(out) :: det_out(N_int,2,N_det_l)
    double precision, intent(out)  :: coef_out(N_det_l,N_states)
  
    integer :: k
    integer :: i_ok
    double precision :: phase
    do k=1,N_det_l
      coef_out(k,:) = psi_coef(k,:)
      det_out(:,:,k) = psi_det(:,:,k)
      call add_electron(det_out(:,:,k),i_particle,ispin,i_ok)
      if (i_ok == 1) then
        call get_phase_ca(psi_det(:,:,k),i_particle,ispin,phase)
        coef_out(k,:) = phase * coef_out(k,:)
      else
        coef_out(k,:) = 0.d0
        det_out(:,:,k) = psi_det(:,:,k)
      endif
    enddo
end

subroutine build_R_wavefunction(i_hole, ispin, coef_out, det_out, N_det_l)
    implicit none
    BEGIN_DOC
    ! Applies the annihilation operator: a_(i_hole) of
    ! spin = ispin to the current wave function (psi_det, psi_coef)
    END_DOC
    integer, intent(in)            :: i_hole,ispin,N_det_l
    integer(bit_kind), intent(out) :: det_out(N_int,2,N_det)
    double precision, intent(out)  :: coef_out(N_det,N_states)
  
    integer :: k
    integer :: i_ok
    double precision :: phase
    do k=1,N_det_l
      coef_out(k,:) = psi_coef(k,:)
      det_out(:,:,k) = psi_det(:,:,k)
      call remove_electron(det_out(:,:,k),i_hole,ispin,i_ok)
      if (i_ok == 1) then
        call get_phase_ca(psi_det(:,:,k),i_hole,ispin,phase)
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
      call add_electron(det_out(:,:,k),i_particle,ispin,i_ok)
      if (i_ok == 1) then
        call get_phase_ca(psi_det(:,:,k),i_particle,ispin,phase)
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
      call remove_electron(det_out(:,:,k),i_hole,ispin,i_ok)
      if (i_ok == 1) then
        call get_phase_ca(psi_det(:,:,k),i_hole,ispin,phase)
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

subroutine set_ref_bitmask(iorb, ispin, ac_type)
    ! TODO: add complex implementation

    implicit none
    BEGIN_DOC
    ! Reset the bitmask and bitmask energy for added/removed electron
    ! ac_type == F if adding electron, T if removing
    END_DOC

    integer :: i, j, a, b, occ(N_int*bit_kind_size,2)
    integer, intent(in)  :: iorb, ispin
    logical, intent(in)  :: ac_type
    integer, allocatable :: t_occ(:)

    !! Resetting bitmasks
    ! adjust alpha/beta sizes based on annihilation vs creation and spin type
    a = elec_alpha_num + merge(merge(-1, 1, ac_type),0,ispin==1)
    b = elec_beta_num  + merge(merge(-1, 1, ac_type),0,ispin==2)

    elec_num_tab(1) = a
    elec_num_tab(2) = b
    elec_num = a + b

    allocate(t_occ(max(a,b)))

    do i = 1, max(a,b)
        t_occ(i) = i
    end do

    call list_to_bitstring(ref_bitmask(1,1), t_occ(:a), a, N_int)
    call list_to_bitstring(ref_bitmask(1,2), t_occ(:b), b, N_int)

    deallocate(t_occ)

    do i = 1, N_int
        ref_closed_shell_bitmask(i,1) = ref_bitmask(i,1)
        ref_closed_shell_bitmask(i,2) = ref_bitmask(i,2)
    end do

    ! clear out extra orbitals
    if (a < b) then
        do i = a+1, b
            call clear_bit_to_integer(i, ref_closed_shell_bitmask(1,2),N_int)
        end do
    else if (b < a) then
        do i = b+1, a
            call clear_bit_to_integer(i, ref_closed_shell_bitmask(1,1),N_int)
        end do
    end if

    !! Resetting bitmask energy
    call bitstring_to_list(ref_bitmask(1,1), occ(1,1), i, N_int)
    call bitstring_to_list(ref_bitmask(1,2), occ(1,2), i, N_int)
    ref_bitmask_energy = 0.d0

    do i = 1, a
        ref_bitmask_energy += mo_one_e_integrals_diag(occ(i,1))
    end do

    do i = 1, b
        ref_bitmask_energy += mo_one_e_integrals_diag(occ(i,2))
    end do 

    do j = 1, a
        do i = j+1, a
            ref_bitmask_energy += mo_two_e_integrals_jj_anti(occ(i,1),occ(j,1))
        end do
    end do

    do j = 1, b
        do i = j+1, b
            ref_bitmask_energy += mo_two_e_integrals_jj_anti(occ(i,2),occ(j,2))
        end do

        do i = 1, a
            ref_bitmask_energy += mo_two_e_integrals_jj(occ(i,1),occ(j,2))
        end do
    end do

    !! Force provides for single excitation things
    PROVIDE fock_op_cshell_ref_bitmask
end