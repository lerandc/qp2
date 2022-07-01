program spectral_density
    implicit none
    BEGIN_DOC
    ! Program that calculates the spectral density.
    END_DOC

    logical                         :: force_reads
    integer                         :: N_det_read
    integer(bit_kind), allocatable  :: psi_det_read(:,:,:)
    double precision , allocatable  :: psi_coef_read(:,:)
    complex*16       , allocatable  :: psi_coef_complex_read(:,:)

    read_wf = .true. 

    call ezfio_get_determinants_n_det(N_det_read)
    N_det = N_det_read
    if (N_det == N_det_read) then 
        allocate(psi_det_read(N_int, 2, N_det))
        call ezfio_get_determinants_psi_det(psi_det_read)
        psi_det = psi_det_read
        deallocate(psi_det_read)
    end if

    if (is_complex) then 
        ! allocate(psi_coef_complex_read(N_det, N_states))
        ! call ezfio_get_determinants_psi_coef_complex(psi_coef_complex_read)
        ! psi_coef_complex = psi_coef_complex_read
        ! deallocate(psi_coef_complex_read)

        ! force_reads = size(psi_coef_complex, 1) == N_det_read .and.&
        ! size(psi_det, 3)  == N_det_read .and. read_wf
        ! if(force_reads) then
        !     if(spectral_density_calc_A) then
        !         call ezfio_set_spectral_density_spectral_density_A_complex(spectral_density_A_complex)
        !         if(write_greens_f) then
        !             call ezfio_set_spectral_density_greens_A_complex(greens_A_complex)
        !         end if

        !         if(write_lanczos_ab) then
        !             call ezfio_set_spectral_density_lanczos_alpha_A_complex(lanczos_alpha_A_complex)
        !             call ezfio_set_spectral_density_lanczos_beta_A_complex(lanczos_beta_A_complex)
        !         end if
        !     end if 
        !     if(spectral_density_calc_R) then
        !         call ezfio_set_spectral_density_spectral_density_R_complex(spectral_density_R_complex)
        !         if(write_greens_f) then
        !             call ezfio_set_spectral_density_greens_R_complex(greens_R_complex)
        !         end if

        !         if(write_lanczos_ab) then
        !             call ezfio_set_spectral_density_lanczos_alpha_R_complex(lanczos_alpha_R_complex)
        !             call ezfio_set_spectral_density_lanczos_beta_R_complex(lanczos_beta_R_complex)
        !         end if
        !     end if
        ! end if
        print *, "waiting a bit to finish this"
    else
        allocate(psi_coef_read(N_det, N_states))
        call ezfio_get_determinants_psi_coef(psi_coef_read)
        psi_coef = psi_coef_read
        deallocate(psi_coef_read)

        force_reads = size(psi_coef, 1) == N_det_read .and.&
        size(psi_det, 3)  == N_det_read .and. read_wf
        
        if(force_reads) then
            if(spectral_density_calc_A) then
                call ezfio_set_spectral_density_spectral_density_A(spectral_density_A)
                if(write_greens_f) then
                    call ezfio_set_spectral_density_greens_A(greens_A)
                end if

                if(write_lanczos_ab) then
                    call ezfio_set_spectral_density_lanczos_alpha_A(lanczos_alpha_A)
                    call ezfio_set_spectral_density_lanczos_beta_A(lanczos_beta_A)
                end if
                
            end if 
            if(spectral_density_calc_R) then
                call ezfio_set_spectral_density_spectral_density_R(spectral_density_R)
                if(write_greens_f) then
                    call ezfio_set_spectral_density_greens_R(greens_R)
                end if

                if(write_lanczos_ab) then
                    call ezfio_set_spectral_density_lanczos_alpha_R(lanczos_alpha_R)
                    call ezfio_set_spectral_density_lanczos_beta_R(lanczos_beta_R)
                end if
            end if
        end if
    end if

    ! if you are writing out alpha/beta to be able to study convergence,
    ! you will also need ground state energy
    if(write_lanczos_ab) call ezfio_set_davidson_psi_energy(psi_energy)

    
end

BEGIN_PROVIDER [double precision, spectral_density_A, (greens_omega_N, n_iorb_A, ns_dets)]
    implicit none

    double precision :: pi

    pi = acos(-1.0)
    spectral_density_A = (-1.0/pi) * aimag(greens_A)

END_PROVIDER

BEGIN_PROVIDER [double precision, spectral_density_R, (greens_omega_N, n_iorb_R, ns_dets)]
    implicit none

    double precision :: pi

    pi = acos(-1.0)
    spectral_density_R = (-1.0/pi) * aimag(greens_R)

END_PROVIDER

! BEGIN_PROVIDER [double precision, spectral_density_A_complex, (greens_omega_N, n_iorb_A, ns_dets)]
!     implicit none

!     double precision :: pi

!     pi = acos(-1.0)
!     spectral_density_A_complex = (-1.0/pi) * aimag(greens_A_complex)

! END_PROVIDER

! BEGIN_PROVIDER [double precision, spectral_density_R_complex, (greens_omega_N, n_iorb_R, ns_dets)]
!     implicit none

!     double precision :: pi

!     pi = acos(-1.0)
!     spectral_density_R_complex = (-1.0/pi) * aimag(greens_R_complex)

! END_PROVIDER