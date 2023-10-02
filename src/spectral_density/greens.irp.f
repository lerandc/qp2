BEGIN_PROVIDER [integer, greens_omega_N]
    implicit none
    greens_omega_N = ceiling((greens_omega_max - greens_omega_min)/ greens_omega_resolution) + 1
    call ezfio_set_spectral_density_greens_omega_N(greens_omega_N)
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
    double precision        :: E0, norm, dnrm2, pi
    double precision        :: t0, t1
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N), bbeta(lanczos_N), cfraction_c
    integer                 :: i, iorb, nnz, nnz_l, i_n_det, N_det_l
    integer(kind=8)              :: s_max_sze
    integer, allocatable    :: H_c(:), t_H_c(:), H_p(:), H_c_all(:), H_p_all(:), I_k(:)
    double precision, allocatable  ::  H_v(:), t_H_v(:), psi_coef_excited(:,:)
    integer(bit_kind), allocatable :: det_excited(:,:,:)
    character(len=72)       :: filename

    ! calculate the maximum size of the sparse arrays with some overflow protection
    ! could still be much improved
    if (nnz_max_per_row > 0) then
        s_max_sze = max_row_sze_factor*nnz_max_per_row
    else 
        s_max_sze = max_row_sze_factor*1000000 
    end if

    if (s_max_sze < 0) then
        print *, "Desired max row size is hitting integer overflow. Setting max size to 2^32"
        s_max_sze = 2**32
    end if
    
    print *, "Calculating spectral densities for added electrons in orbitals: "
    
    do i = 1, n_iorb_A
        if (modulo(i, 5) == 0) then
            write(*, "(I6)"), iorb_A(i)
        else
            write(*, "(I6)", advance="no"), iorb_A(i)
        end if
    end do

    if (modulo(n_iorb_A, 5) >0) write(*,*)," "
    
    print *, " with N dets of "
    
    do i = 1, ns_dets
        if (modulo(i, 5) == 0) then
            write(*, "(I10)"), n_det_sequence(i)
        else
            write(*, "(I10)", advance="no"), n_det_sequence(i)
        end if
    end do
    
    if (modulo(ns_dets,5) > 0) write(*,*)," "
    

    greens_A = (0.d0, 0.d0)
    ! loop over number of determinants 
    do i_n_det = 1, ns_dets
        
        N_det_l = n_det_sequence(i_n_det)
        allocate(psi_coef_excited(N_det_l, N_states), det_excited(N_int,2,N_det_l))

        ! Calculate full sparsity structure for N_det_l
        call wall_time(t0)
        allocate(H_c_all(s_max_sze), H_p_all(N_det_l+1), I_k(N_det_l))
        call get_sparsity_structure(H_p_all, H_c_all, s_max_sze, N_det_l)
        
        nnz = H_p_all(N_det_l+1)-1

        ! move vectors to smaller allocations
        allocate(t_H_c(nnz))
        t_H_c = H_c_all(:nnz)
        call move_alloc(t_H_c, H_c_all)

        call wall_time(t1)
        write(*, "(A33, F8.2, A10)"), "Sparsity structure calculated in ", t1-t0, " seconds"

        ! loop over orbitals
        do iorb = 1, n_iorb_A

            write(*, "(A34, I6, A5, I11, A13)"),&
                '#### Calculating density for iorb ', iorb_A(iorb),&
                ' and ', N_det_l, ' determinants'

            ! reset bit masks
            call set_ref_bitmask(1, .false.)

            ! prepare orthonormal wavefunction in space of N+1 determinants
            ! add electron to orbital
            call build_A_wavefunction(iorb_A(iorb),1,psi_coef_excited,det_excited,N_det_l,I_k)
            norm = dnrm2(N_det_l, psi_coef_excited(:,1), 1)
            psi_coef_excited = psi_coef_excited / norm

            ! prepare sparse Hamiltonian arrays
            call wall_time(t0)
            allocate(H_c(nnz), H_p(N_det_l+1))

            ! reduce sparse Hamiltonian to set of valid determinants for this excitation
            ! call sparse_csr_MM(H_c_all,H_p_all, I_k, H_c, H_p, N_det_l, nnz)
            call sparse_csr_MM_row_part(H_c_all,H_p_all, I_k, H_c, H_p, N_det_l, nnz)
            call wall_time(t1)
            write(*, "(A33, F8.2, A10)"), "Intrinic space reduced in ", t1-t0, " seconds"

            nnz_l = H_p(N_det_l+1)-1

            ! move vectors to smaller allocations
            allocate(t_H_c(nnz_l), H_v(nnz_l))
            t_H_c = H_c(:nnz_l)
            call move_alloc(t_H_c, H_c)

            ! calculate entries of Hamiltonian
            call calc_sparse_dH(H_p, H_c, H_v, N_det_l, nnz_l, det_excited)
            ! call form_sparse_dH(H_p, H_c, H_v, s_max_sze, det_excited,&
            !                     iorb_A(iorb), 1, .false., N_det_l)

            call wall_time(t1)
            write(*, "(A33, F8.2, A10)"), "Sparse Hamiltonian calculated in ", t1-t0, " seconds"

            call wall_time(t0)
            ! calculate tridiagonalization of Hamiltoninan
            call lanczos_tridiag_sparse_reortho_r(H_v, H_c, H_p, psi_coef_excited(:,1),&
                                        alpha, beta,&
                                        lanczos_N, nnz_l, N_det_l)
            call wall_time(t1)
            write(*, "(A21, F8.2, A10)"), "Lanczos iteration in ", t1-t0, " seconds"
            
            ! prepare beta array for continued fractions
            bbeta(1) = (1.d0, 0.d0)
            do i = 2, lanczos_N
                bbeta(i) = -1.d0*beta(i)**2.0
            end do

            lanczos_alpha_A(:, iorb, i_n_det) = alpha
            lanczos_beta_A(:, iorb, i_n_det) = real(bbeta)

            epsilon = greens_epsilon ! broadening factor
            E0 = psi_energy(1)
            z = E0 + (greens_omega + (0.d0, 1.d0)*epsilon)
            
            ! calculate greens functions
            !$OMP PARALLEL DO PRIVATE(i, zalpha) SHARED(alpha, bbeta, lanczos_N, greens_A)&
            !$OMP SCHEDULE(GUIDED)
            do i = 1, greens_omega_N
                zalpha = z(i) - alpha
                greens_A(i, iorb, i_n_det) = cfraction_c((0.d0, 0.d0), bbeta, zalpha,lanczos_N)
            end do
            !$OMP END PARALLEL DO

            if (dump_intermediate_output) then
                print *, "Dumping intermediate outputs to: "
                ! write out abcissa
                write(filename, *), "omega_A.out"
                print *, filename
                call dump_array_real(filename, greens_omega, greens_omega_N)

                ! if requested, write out greens functions
                if (write_greens_f) then
                    write(filename, "(A, A, I4.4, A, I10.10, A)"), "greens_A",&
                                                        "_iorb_", iorb_A(iorb), &
                                                        "_ndets_", N_det_l,&
                                                        ".out"
                    print *, filename
                    call dump_array_complex(filename, greens_A(:, iorb, i_n_det), greens_omega_N)
                end if

                ! write out spectral density
                write(filename, "(A, A, I4.4, A, I10.10, A)"), "spectral_density_A",&
                                                    "_iorb_", iorb_A(iorb), &
                                                    "_ndets_", N_det_l,&
                                                    ".out"
                print *, filename

                pi = acos(-1.d0)
                call dump_array_real(filename, (-1.d0/pi) * aimag(greens_A(:, iorb, i_n_det)), greens_omega_N)

                ! if requested, write out lanczos vectors

                if (write_lanczos_ab) then
                    write(filename, "(A, A, I4.4, A, I10.10, A, I6.6, A)"), "lanczos_alpha_A",&
                                                        "_iorb_", iorb_A(iorb),&
                                                        "_ndets_", N_det_l,&
                                                        "_nvecs_", lanczos_N,&
                                                        ".out"
                    print *, filename
                    call dump_array_real(filename, lanczos_alpha_A(:, iorb, i_n_det), lanczos_N)

                    write(filename, "(A, A, I4.4, A, I10.10, A, I6.6, A)"), "lanczos_beta_A",&
                                                        "_iorb_", iorb_A(iorb),&
                                                        "_ndets_", N_det_l,&
                                                        "_nvecs_", lanczos_N,&
                                                        ".out"
                    print *, filename
                    call dump_array_real(filename, lanczos_beta_A(:, iorb, i_n_det), lanczos_N)
                end if
            end if

            deallocate(H_c, H_v, H_p)

            call write_time(0)
            call print_memory_usage()
        end do

        deallocate(psi_coef_excited, det_excited, H_c_all, H_p_all, I_k)
    end do

END_PROVIDER

BEGIN_PROVIDER [complex*16, greens_R, (greens_omega_N, n_iorb_R, ns_dets)]
&BEGIN_PROVIDER[double precision, lanczos_alpha_R, (lanczos_N, n_iorb_R,ns_dets)]
&BEGIN_PROVIDER[double precision, lanczos_beta_R, (lanczos_N, n_iorb_R,ns_dets)]
    implicit none

    double precision        :: alpha(lanczos_N), beta(lanczos_N), epsilon
    double precision        :: E0, norm, dnrm2, pi
    double precision        :: t0, t1
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N), bbeta(lanczos_N), cfraction_c
    integer                 :: i, iorb, i_n_det, N_det_l
    integer(kind=8)         :: nnz, nnz_l, s_max_sze
    integer, allocatable    :: H_c(:), t_H_c(:), H_c_all(:), I_k(:)
    integer(kind=8), allocatable :: H_p(:), H_p_all(:)
    double precision , allocatable ::  H_v(:), t_H_v(:), psi_coef_excited(:,:) 
    integer(bit_kind), allocatable :: det_excited(:,:,:)
    character(len=72)       :: filename


    ! variables specifically used for orbital coupling procedures
    logical                 :: orbital_coupling, hash_success
    integer                 :: hash_table_size, n_coupled_dets, hash_prime, ref_count, ref_degree_a, ref_degree_b, test_degree_a, test_degree_b, j, k
    integer(kind=8)         :: nnz_max_u
    integer, allocatable    :: I_cut_k(:,:), hash_vals(:), H_e_all(:), t_H_e(:)
    integer, allocatable    :: uH_c(:), t_uH_c(:), I_det_ind(:,:)
    integer(kind=8), allocatable :: uH_p(:)
    double precision, allocatable ::  uH_v(:)
    integer(bit_kind), allocatable :: I_det_k(:,:,:,:), hash_alpha(:,:), hash_beta(:,:), det_basis(:,:,:), t_det_basis(:,:,:)
    double precision, allocatable  :: psi_coef_intrinsic_excited(:,:,:), psi_coef_coupled_excited(:,:)
    integer :: hash_value, test_col, mem_err, intrinsic_rank, max_union_rank
    double precision :: test_val

    ! calculate the maximum size of the sparse arrays with some overflow protection
    ! could still be much improved
    s_max_sze = 0
    print *, s_max_sze
    if (nnz_max_per_row > 0) then
        s_max_sze = max_row_sze_factor*nnz_max_per_row
    else 
        s_max_sze = max_row_sze_factor*1000000 
    end if
    
    print *, s_max_sze
    if (s_max_sze <= 0) then
        print *, "Desired max row size is hitting integer overflow. Setting max size to 2^32"
        s_max_sze = 2_8**32
    end if
    print *, s_max_sze
    
    print *, "Calculating spectral densities for removed electrons in orbitals: "
    
    do i = 1, n_iorb_R
        if (modulo(i, 5) == 0) then
            write(*, "(I6)"), iorb_R(i)
        else
            write(*, "(I6)", advance="no"), iorb_R(i)
        end if
    end do
    
    if (modulo(n_iorb_R, 5) >0) write(*,*)," "
    
    print *, " with N dets of "
    
    do i = 1, ns_dets
        if (modulo(i, 5) == 0) then
            write(*, "(I10)"), n_det_sequence(i)
        else
            write(*, "(I10)", advance="no"), n_det_sequence(i)
        end if
    end do
    
    if (modulo(ns_dets, 5) > 0) write(*,*)," "
    
    ! greens_R = (0.d0, 0.d0)

    if (spectral_density_orbital_coupling == .true.) then

        call set_ref_bitmask(1, .true.)
        ! work on just the terminal number of determinants for now
        N_det_l = n_det_sequence(ns_dets) 
        
        t0 = 0.0
        t1 = 0.0
        call wall_time(t0)

        !! get all sets of I_cut, I_det, new wavefunctions
        allocate(I_cut_k(N_det_l, n_iorb_R), I_det_k(N_int, 2, N_det_l, n_iorb_R),&
                    psi_coef_intrinsic_excited(N_det_l, N_states, n_iorb_R), stat=mem_err)

        ! if(mem_err /= 0) then
        !     print *, "Allocation failed at before build_R_wavefunction"
        ! end if
        max_union_rank = 0
        do iorb = 1, n_iorb_R

            call build_R_wavefunction(iorb_R(iorb), 1, &
                                    psi_coef_intrinsic_excited(:,:,iorb), &
                                    I_det_k(:,:,:,iorb), &
                                    N_det_l, &
                                    I_cut_k(:,iorb),&
                                    intrinsic_rank,&
                                    )
            max_union_rank += intrinsic_rank
        end do
                                
                                
        !! construct the hash table
        ! eventually, need a better heuristic for the hash table size to balance size in memory
        ! at larger intrinsic space sizes, this works better since there are more generated determinants 
        ! falling into intersecting spaces
        hash_table_size =N_det_l * n_iorb_R * 32
        hash_prime = 8191

        allocate(hash_vals(hash_table_size), hash_alpha(N_int, hash_table_size), hash_beta(N_int, hash_table_size), stat=mem_err)
        allocate(I_det_ind(N_det_l, n_iorb_r))
        allocate(det_basis(N_int, 2, hash_table_size), stat=mem_err)

        hash_alpha = 0
        hash_beta = 0
        hash_vals = 0
        det_basis = 0

        call build_hash_table(hash_alpha, hash_beta, hash_vals, hash_prime, hash_table_size,&
                              I_cut_k, I_det_k, I_det_ind, n_iorb_R, N_det_l, det_basis, n_coupled_dets)
    
        allocate(t_det_basis(N_int, 2, n_coupled_dets), stat=mem_err)
        t_det_basis(:,:,:n_coupled_dets) = det_basis(:,:,:n_coupled_dets)
        call move_alloc(t_det_basis, det_basis, stat=mem_err)

        call wall_time(t1)
        write(*, "(A42, F8.2, A10)"), "Union space basis and hash table built in ", t1-t0, " seconds"
        write(*, "(A30, I10, A1, I10, A1, I10)"), "Size of union space: ", n_coupled_dets, "/", max_union_rank, "/", N_det_l*n_iorb_R

        print *, "Building ground state Hamiltonian"
        call wall_time(t0)
        print *, "Max allowable size:", s_max_sze
        allocate(H_p_all(N_det_l+1), H_c_all(s_max_sze), H_e_all(s_max_sze))

        !!! get ground state sparsity structure with triples and values corresponding to excitation degree
        call get_sparsity_structure_with_triples(H_p_all, H_c_all, H_e_all, s_max_sze, N_det_l)
        ! shrink allocations
        nnz = H_p_all(N_det_l+1)-1
        
        allocate(t_H_c(nnz), t_H_e(nnz), stat=mem_err)
        t_H_c(:nnz) = H_c_all(:nnz)
        t_H_e(:nnz) = H_e_all(:nnz)
        call move_alloc(t_H_c, H_c_all, stat=mem_err)
        call move_alloc(t_H_e, H_e_all, stat=mem_err)
        
        call wall_time(t1)
        write(*, "(A58, F8.2, A10)"), "Ground state sparsity structure with triples built in ", t1-t0, " seconds"
        
        !!! construct sparse hamiltonian on coupled excitation space
        ! heuristic maximum number of nonzero elements on the union space should be nnz with triples
        print *, "Constructing coupled matrix"
        call wall_time(t0)
        nnz_max_u = 2_8*nnz !n_iorb_R*n_iorb_R*nnz
        print *, nnz_max_u
        allocate(uH_p(n_coupled_dets+1), uH_c(nnz_max_u), stat=mem_err)
        if(mem_err /= 0) then
            print *, "Allocation failed for union state matrix", mem_err
        end if
        call uH_structure_from_gH(H_p_all, H_c_all, H_e_all, uH_p, uH_c,&
                            N_det_l, n_coupled_dets, nnz, nnz_max_u, n_iorb_R,&
                            hash_alpha, hash_beta, hash_vals, hash_prime, hash_table_size,&
                            I_cut_k, I_det_k, I_det_ind)

        nnz_max_u = uH_p(n_coupled_dets+1)-1
        allocate(t_uH_c(nnz_max_u), uH_v(nnz_max_u), stat=mem_err)
        t_uH_c(:nnz_max_u) = uH_c(:nnz_max_u)
        call move_alloc(t_uH_c, uH_c, stat=mem_err)
        
        call wall_time(t1)
        write(*, "(A50, F8.2, A15, I12, A8)"), "Coupled state H generated from ground state in ", t1-t0, " seconds with ", nnz_max_u, " entries"
        
        !!! compute an averaged wave function
        call wall_time(t0)
        print *, "Calculating matrix entries"
        call calc_sparse_dH(uH_p, uH_c, uH_v, n_coupled_dets, nnz_max_u, det_basis)
        
        print *, "Building coupled wavefunction"
        allocate(psi_coef_coupled_excited(n_coupled_dets, N_states), stat=mem_err)
        psi_coef_coupled_excited = 0.d0
        call build_coupled_wavefunction(psi_coef_intrinsic_excited, psi_coef_coupled_excited, &
                                        I_cut_k, I_det_k, I_det_ind, N_det_l, n_coupled_dets, n_iorb_R, hash_alpha,&
                                        hash_beta, hash_vals, hash_prime, hash_table_size)
        
        call wall_time(t1)
        write(*, "(A60, F8.2, A10)"), "Coupled state entries and excited wavefunction calculated in ", t1-t0, " seconds"

        print *, "Beginning Lanczos iteration"
        call wall_time(t0)
        !!! proceed with lanczos iteration as normal
        call lanczos_tridiag_sparse_reortho_row_part_r(uH_v, uH_c, uH_p, psi_coef_coupled_excited(:,1),&
                                                alpha, beta,&
                                                lanczos_N, nnz_max_u, n_coupled_dets)

        call wall_time(t1)
        write(*, "(A21, F8.2, A10)"), "Lanczos iteration in ", t1-t0, " seconds"

        lanczos_alpha_R = 0.d0
        lanczos_beta_R = 0.d0
        greens_R = (0.d0, 0.d0)

        ! prepare beta array for continued fractions
        bbeta(1) = (1.d0, 0.d0)
        do i = 2, lanczos_N
            bbeta(i) = -1.d0*beta(i)**2.0
        end do

        ! save only to first index for now
        lanczos_alpha_R(:, 1, 1) = alpha
        lanczos_beta_R(:, 1, 1) = real(bbeta)

        epsilon = greens_epsilon ! broadening factor
        E0 = psi_energy(1)
        z = E0 - (-1.d0*greens_omega + (0.d0, 1.d0)*epsilon) ! omega is abs. energy value

        ! calculate greens functions
        !$OMP PARALLEL DO PRIVATE(i, zalpha) SHARED(alpha, bbeta, lanczos_N, greens_R)&
        !$OMP SCHEDULE(GUIDED)
        do i = 1, greens_omega_N
            zalpha = z(i) - alpha
            greens_R(i, 1, 1) = -1.d0*cfraction_c( (0.d0, 0.d0), bbeta, zalpha, lanczos_N)
        end do
        !$OMP END PARALLEL DO
                
        ! clean up
        deallocate(H_p_all, H_c_all, H_e_all, uH_p, uH_c, uH_v)
        deallocate(hash_vals, hash_alpha, hash_beta)
        deallocate(det_basis, I_cut_k, I_det_k, I_det_ind, psi_coef_intrinsic_excited, psi_coef_coupled_excited)
    else
        ! loop over number of determinants 
        do i_n_det = 1, ns_dets

            N_det_l = n_det_sequence(i_n_det)
            allocate(psi_coef_excited(N_det_l, N_states), det_excited(N_int,2,N_det_l))

            ! Calculate full sparsity structure for N_det_l
            call wall_time(t0)
            allocate(H_c_all(s_max_sze), H_p_all(N_det_l+1), I_k(N_det_l))
            call get_sparsity_structure(H_p_all, H_c_all, s_max_sze, N_det_l)
            
            nnz = H_p_all(N_det_l+1)-1

            ! move vectors to smaller allocations
            allocate(t_H_c(nnz))
            t_H_c = H_c_all(:nnz)
            call move_alloc(t_H_c, H_c_all)

            call wall_time(t1)        
            write(*, "(A33, F8.2, A10)"), "Sparsity structure calculated in ", t1-t0, " seconds"

            ! loop over orbitals
            do iorb = 1, n_iorb_R

                write(*, "(A34, I6, A5, I11, A13)"),&
                    '#### Calculating density for iorb ', iorb_R(iorb),&
                    ' and ', N_det_l, ' determinants'

                ! reset bit masks
                call set_ref_bitmask(1, .true.)

                ! prepare orthonormal wavefunction in space of N-1 determinants
                ! remove electron from orbital
                call build_R_wavefunction(iorb_R(iorb),1,psi_coef_excited,det_excited, N_det_l, I_k, intrinsic_rank)
                norm = dnrm2(N_det_l, psi_coef_excited(:,1), 1)
                psi_coef_excited = psi_coef_excited / norm

                ! prepare sparse Hamiltonian arrays
                call wall_time(t0)
                allocate(H_c(nnz), H_p(N_det_l+1))
                H_c = 0
                H_p = 0

                ! reduce sparse Hamiltonian to set of valid determinants for this excitation
                call sparse_csr_MM_row_part(H_c_all,H_p_all, I_k, H_c, H_p, N_det_l, nnz)
                call wall_time(t1)
                write(*, "(A33, F8.2, A10)"), "Intrinic space reduced in ", t1-t0, " seconds"         
                
                nnz_l = H_p(N_det_l+1)-1

                ! move vectors to smaller allocations
                allocate(t_H_c(nnz_l), H_v(nnz_l))
                t_H_c = H_c(:nnz_l)
                call move_alloc(t_H_c, H_c)

                ! calculate entries of Hamiltonian
                call calc_sparse_dH(H_p, H_c, H_v, N_det_l, nnz_l, det_excited)
                ! call form_sparse_dH(H_p, H_c, H_v, s_max_sze, det_excited,&
                !                     iorb_R(iorb), 1, .true., N_det_l)

                call wall_time(t1)
                write(*, "(A33, F8.2, A10)"), "Sparse Hamiltonian calculated in ", t1-t0, " seconds"

                call wall_time(t0)
                ! calculate tridiagonalization of Hamiltoninan
                call lanczos_tridiag_sparse_reortho_r(H_v, H_c, H_p, psi_coef_excited(:,1),&
                                                alpha, beta,&
                                                lanczos_N, nnz, N_det_l)

                call wall_time(t1)
                write(*, "(A21, F8.2, A10)"), "Lanczos iteration in ", t1-t0, " seconds"
                ! prepare beta array for continued fractions
                bbeta(1) = (1.d0, 0.d0)
                do i = 2, lanczos_N
                    bbeta(i) = -1.d0*beta(i)**2.0
                end do

                lanczos_alpha_R(:, iorb, i_n_det) = alpha
                lanczos_beta_R(:, iorb, i_n_det) = real(bbeta)

                epsilon = greens_epsilon ! broadening factor
                E0 = psi_energy(1)
                z = E0 - (-1.d0*greens_omega + (0.d0, 1.d0)*epsilon) ! omega is abs. energy value

                ! calculate greens functions
                !$OMP PARALLEL DO PRIVATE(i, zalpha) SHARED(alpha, bbeta, lanczos_N, greens_R)&
                !$OMP SCHEDULE(GUIDED)
                do i = 1, greens_omega_N
                    zalpha = z(i) - alpha
                    greens_R(i, iorb, i_n_det) = -1.d0*cfraction_c( (0.d0, 0.d0), bbeta, zalpha, lanczos_N)
                end do
                !$OMP END PARALLEL DO

                if (dump_intermediate_output) then
                    print *, "Dumping intermediate outputs to: "
                    ! write out abcissa
                    write(filename, *), "omega_R.out"
                    print *, filename
                    call dump_array_real(filename, -1.d0*greens_omega, greens_omega_N)

                    ! if requested, write out greens functions
                    if (write_greens_f) then
                        write(filename, "(A, A, I4.4, A, I10.10, A)"), "greens_R",&
                                                            "_iorb_", iorb_R(iorb), &
                                                            "_ndets_", N_det_l,&
                                                            ".out"
                        print *, filename
                        call dump_array_complex(filename, greens_R(:, iorb, i_n_det), greens_omega_N)
                    end if

                    ! write out spectral density
                    write(filename, "(A, A, I4.4, A, I10.10, A)"), "spectral_density_R",&
                                                        "_iorb_", iorb_R(iorb), &
                                                        "_ndets_", N_det_l,&
                                                        ".out"
                    print *, filename

                    pi = acos(-1.d0)
                    call dump_array_real(filename, (-1.d0/pi) * aimag(greens_R(:, iorb, i_n_det)), greens_omega_N)

                    ! if requested, write out lanczos vectors

                    if (write_lanczos_ab) then
                        write(filename, "(A, A, I4.4, A, I10.10, A, I6.6, A)"), "lanczos_alpha_R",&
                                                            "_iorb_", iorb_R(iorb),&
                                                            "_ndets_", N_det_l,&
                                                            "_nvecs_", lanczos_N,&
                                                            ".out"
                        print *, filename
                        call dump_array_real( filename, lanczos_alpha_R(:,iorb,i_n_det),lanczos_N)

                        write(filename, "(A, A, I4.4, A, I10.10, A, I6.6, A)"), "lanczos_beta_R",&
                                                            "_iorb_", iorb_R(iorb),&
                                                            "_ndets_", N_det_l,&
                                                            "_nvecs_", lanczos_N,&
                                                            ".out"
                        print *, filename
                        call dump_array_real(filename, lanczos_beta_R(:,iorb,i_n_det),lanczos_N)
                    end if
                end if

                deallocate(H_c, H_v, H_p)

                call write_time(0)
                call print_memory_usage()
            end do

            deallocate(psi_coef_excited, det_excited, H_c_all, H_p_all, I_k)
        end do
    end if

END_PROVIDER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   Complex Implementations   !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

BEGIN_PROVIDER [complex*16, greens_A_complex, (greens_omega_N, n_iorb_A,ns_dets)]
&BEGIN_PROVIDER[double precision,lanczos_alpha_A_complex, (lanczos_N, n_iorb_A,ns_dets)]
&BEGIN_PROVIDER[double precision, lanczos_beta_A_complex, (lanczos_N, n_iorb_A,ns_dets)]   
    implicit none

    double precision        :: alpha(lanczos_N), beta(lanczos_N), epsilon, E0, norm, dznrm2, pi, t0, t1
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N), bbeta(lanczos_N), cfraction_c
    integer                 :: i, iorb, nnz, nnz_l, i_n_det, N_det_l
    integer(kind=8)         :: s_max_sze
    integer, allocatable    :: H_c(:), H_p(:), t_H_c(:), H_c_all(:), H_p_all(:), I_k(:)
    integer(bit_kind), allocatable :: det_excited(:,:,:)
    complex*16 , allocatable  ::  H_v(:), t_H_v(:), psi_coef_excited(:,:)
    character(len=72)       :: filename

    ! calculate the maximum size of the sparse arrays with some overflow protection
    ! could still be much improved
    if (nnz_max_per_row > 0) then
        s_max_sze = max_row_sze_factor*nnz_max_per_row
    else 
        s_max_sze = max_row_sze_factor*1000000 
    end if

    if (s_max_sze < 0) then
        print *, "Desired max row size is hitting integer overflow. Setting max size to 2^32"
        s_max_sze = 2**32
    end if
    
    print *, "Calculating spectral densities for added electrons in orbitals: "
    
    do i = 1, n_iorb_A
        if (modulo(i, 5) == 0) then
            write(*, "(I6)"), iorb_A(i)
        else
            write(*, "(I6)", advance="no"), iorb_A(i)
        end if
    end do
    
    if (modulo(n_iorb_A, 5) >0) write(*,*)," "
    
    print *, " with N dets of "
    
    do i = 1, ns_dets
        if (modulo(i, 5) == 0) then
            write(*, "(I10)"), n_det_sequence(i)
        else
            write(*, "(I10)", advance="no"), n_det_sequence(i)
        end if
    end do
    
    if (modulo(ns_dets, 5) > 0) write(*,*)," "
    
    
    greens_A_complex = (0.d0, 0.d0)
    ! loop over number of determinants 
    do i_n_det = 1, ns_dets
        N_det_l = n_det_sequence(i_n_det)
        allocate(psi_coef_excited(N_det_l, N_states), det_excited(N_int,2,N_det_l))

        ! Calculate full sparsity structure for N_det_l
        call wall_time(t0)
        allocate(H_c_all(s_max_sze), H_p_all(N_det_l+1), I_k(N_det_l))
        call get_sparsity_structure(H_p_all, H_c_all, s_max_sze, N_det_l)

        if (.false.) then
            print *, "Early stopping"
            call dump_array_int("H_columns.out", H_c_all, H_p_all(N_det_l+1)-1)
            call dump_array_int("H_rows.out", H_p_all, N_det_l+1)
            stop
        end if
        
        nnz = H_p_all(N_det_l+1)-1

        ! move vectors to smaller allocations
        allocate(t_H_c(nnz))
        t_H_c = H_c_all(:nnz)
        call move_alloc(t_H_c, H_c_all)

        call wall_time(t1)
        write(*, "(A33, F8.2, A10)"), "Sparsity structure calculated in ", t1-t0, " seconds"

        ! loop over orbitals
        do iorb = 1, n_iorb_A

            write(*, "(A34, I6, A5, I11, A13)"),&
                '#### Calculating density for iorb ', iorb_A(iorb),&
                ' and ', N_det_l, ' determinants'

            call set_ref_bitmask_complex(1, .false.)
        
            ! prepare orthonormal wavefunction in space of N+1 determinants
            ! add electron to orbital
            call build_A_wavefunction_complex(iorb_A(iorb),1,psi_coef_excited,det_excited,N_det_l,I_k)
            norm = dznrm2(N_det_l, psi_coef_excited(:,1), 1)
            psi_coef_excited = psi_coef_excited / norm

             ! prepare sparse Hamiltonian arrays
            call wall_time(t0)
            allocate(H_c(nnz), H_p(N_det_l+1))

            ! reduce sparse Hamiltonian to set of valid determinants for this excitation
            call sparse_csr_MM_row_part(H_c_all,H_p_all, I_k, H_c, H_p, N_det_l, nnz)

            nnz_l = H_p(N_det_l+1)-1

            ! move vectors to smaller allocations
            allocate(t_H_c(nnz_l), H_v(nnz_l))
            t_H_c = H_c(:nnz_l)
            call move_alloc(t_H_c, H_c)

            ! calculate entries of Hamiltonian
            call calc_sparse_zH(H_p, H_c, H_v, N_det_l, nnz_l, det_excited)
            ! call form_sparse_zH(H_p, H_c, H_v, s_max_sze,det_excited,&
            !                     iorb_A(iorb), 1, .false., N_det_l)
            call wall_time(t1)
            write(*, "(A33, F8.2, A10)"), "Sparse Hamiltonian calculated in ", t1-t0, " seconds"

        
            call wall_time(t0)
            ! calculate tridiagonalization of Hamiltoninan
            call lanczos_tridiag_sparse_reortho_c(H_v, H_c, H_p, psi_coef_excited(:,1),&
                                        alpha, beta,&
                                        lanczos_N, nnz, N_det_l)
            call wall_time(t1)
            write(*, "(A21, F8.2, A10)"), "Lanczos iteration in ", t1-t0, " seconds"

            ! prepare beta array for continued fractions
            bbeta(1) = (1.d0, 0.d0)
            do i = 2, lanczos_N
                bbeta(i) = -1.d0*beta(i)**2.0
            end do

            lanczos_alpha_A_complex(:, iorb, i_n_det) = alpha
            lanczos_beta_A_complex(:, iorb, i_n_det) = real(bbeta)

            epsilon = greens_epsilon ! boradening factor
            E0 = psi_energy(1)
            z = E0 + (greens_omega + (0.d0, 1.d0)*epsilon)

            ! calculate greens functions
            !$OMP PARALLEL DO PRIVATE(i, zalpha) SHARED(alpha, bbeta, lanczos_N, greens_A_complex)&
            !$OMP SCHEDULE(GUIDED)
            do i = 1, greens_omega_N
                zalpha = z(i) - alpha
                greens_A_complex(i, iorb, i_n_det) = cfraction_c((0.d0, 0.d0), bbeta, zalpha, lanczos_N)
            end do
            !$OMP END PARALLEL DO

            if (dump_intermediate_output) then
                print *, "Dumping intermediate outputs to: "
                ! write out abcissa
                write(filename, *), "omega_A_complex.out"
                print *, filename
                call dump_array_real(filename, greens_omega, greens_omega_N)

                ! if requested, write out greens functions
                if (write_greens_f) then
                    write(filename, "(A, A, I4.4, A, I10.10, A)"), "greens_A_complex",&
                                                        "_iorb_", iorb_A(iorb), &
                                                        "_ndets_", N_det_l,&
                                                        ".out"
                    print *, filename
                    call dump_array_complex(filename,&
                                            greens_A_complex(:,iorb,i_n_det), greens_omega_N)
                end if

                ! write out spectral density
                write(filename, "(A, A, I4.4, A, I10.10, A)"), "spectral_density_A_complex",&
                                                    "_iorb_", iorb_A(iorb), &
                                                    "_ndets_", N_det_l,&
                                                    ".out"
                print *, filename

                pi = acos(-1.d0)
                call dump_array_real(filename, (-1.d0/pi) * aimag(greens_A_complex(:, iorb, i_n_det)), greens_omega_N)

                ! if requested, write out lanczos vectors

                if (write_lanczos_ab) then
                    write(filename, "(A, A, I4.4, A, I10.10, A, I6.6, A)"), "lanczos_alpha_A_complex",&
                                                        "_iorb_", iorb_A(iorb),&
                                                        "_ndets_", N_det_l,&
                                                        "_nvecs_", lanczos_N,&
                                                        ".out"
                    print *, filename
                    call dump_array_real(filename, lanczos_alpha_A_complex(:, iorb, i_n_det), lanczos_N)

                    write(filename, "(A, A, I4.4, A, I10.10, A, I6.6, A)"), "lanczos_beta_A_complex",&
                                                        "_iorb_", iorb_A(iorb),&
                                                        "_ndets_", N_det_l,&
                                                        "_nvecs_", lanczos_N,&
                                                        ".out"
                    print *, filename
                    call dump_array_real(filename, lanczos_beta_A_complex(:, iorb, i_n_det), lanczos_N)
                end if
            end if

            deallocate(H_c, H_v, H_p)

            call write_time(0)
            call print_memory_usage()
        end do

        deallocate(psi_coef_excited, det_excited, H_c_all, H_p_all, I_k)
    end do
END_PROVIDER

BEGIN_PROVIDER [complex*16, greens_R_complex, (greens_omega_N, n_iorb_R,ns_dets)]
&BEGIN_PROVIDER[double precision,lanczos_alpha_R_complex, (lanczos_N, n_iorb_R,ns_dets)]
&BEGIN_PROVIDER[double precision, lanczos_beta_R_complex, (lanczos_N, n_iorb_R,ns_dets)]   
    implicit none

    double precision        :: alpha(lanczos_N), beta(lanczos_N), epsilon, E0, norm, dznrm2, pi, t0, t1
    complex*16              :: z(greens_omega_N), zalpha(lanczos_N), bbeta(lanczos_N), cfraction_c
    integer                 :: i, iorb, nnz, nnz_l, i_n_det, N_det_l
    integer(kind=8)         :: s_max_sze
    integer, allocatable    :: H_c(:), H_p(:), t_H_c(:),  H_c_all(:), H_p_all(:), I_k(:)
    integer(bit_kind), allocatable :: det_excited(:,:,:)
    complex*16 , allocatable  ::  H_v(:), t_H_v(:), psi_coef_excited(:,:)
    character(len=72)       :: filename


    ! calculate the maximum size of the sparse arrays with some overflow protection
    ! could still be much improved
    if (nnz_max_per_row > 0) then
        s_max_sze = max_row_sze_factor*nnz_max_per_row
    else 
        s_max_sze = max_row_sze_factor*1000000 
    end if

    if (s_max_sze < 0) then
        print *, "Desired max row size is hitting integer overflow. Setting max size to 2^32"
        s_max_sze = 2**32
    end if

    print *, "Calculating spectral densities for removed electrons in orbitals: "
    
    do i = 1, n_iorb_R
        if (modulo(i, 5) == 0) then
            write(*, "(I6)"), iorb_R(i)
        else
            write(*, "(I6)", advance="no"), iorb_R(i)
        end if
    end do
    
    if (modulo(n_iorb_R, 5) >0) write(*,*)," "
    
    print *, " with N dets of "
    
    do i = 1, ns_dets
        if (modulo(i, 5) == 0) then
            write(*, "(I10)"), n_det_sequence(i)
        else
            write(*, "(I10)", advance="no"), n_det_sequence(i)
        end if
    end do
    
    if (modulo(ns_dets, 5) > 0) write(*,*)," "
    
    greens_R_complex = (0.d0, 0.d0)
    do i_n_det = 1, ns_dets
        N_det_l = n_det_sequence(i_n_det)
        allocate(psi_coef_excited(N_det_l, N_states), det_excited(N_int,2,N_det_l))

        ! Calculate full sparsity structure for N_det_l
        call wall_time(t0)
        allocate(H_c_all(s_max_sze), H_p_all(N_det_l+1), I_k(N_det_l))
        call get_sparsity_structure(H_p_all, H_c_all, s_max_sze, N_det_l)
        
        nnz = H_p_all(N_det_l+1)-1

        ! move vectors to smaller allocations
        allocate(t_H_c(nnz))
        t_H_c = H_c_all(:nnz)
        call move_alloc(t_H_c, H_c_all)

        call wall_time(t1)
        write(*, "(A33, F8.2, A10)"), "Sparsity structure calculated in ", t1-t0, " seconds"

        ! loop over orbitals
        do iorb = 1, n_iorb_R

            write(*, "(A34, I6, A5, I11, A13)"),&
                '#### Calculating density for iorb ', iorb_R(iorb),&
                ' and ', N_det_l, ' determinants'

            call set_ref_bitmask_complex(1, .true.)
        
            ! prepare orthonormal wavefunction in space of N-1 determinants
            ! remove electron from orbital
            call build_R_wavefunction_complex(iorb_R(iorb),1,psi_coef_excited,det_excited,N_det_l,I_k)
            norm = dznrm2(N_det_l, psi_coef_excited(:,1), 1)
            psi_coef_excited = psi_coef_excited / norm

            ! prepare sparse Hamiltonian arrays
            call wall_time(t0)
            allocate(H_c(nnz), H_p(N_det_l+1))

            ! reduce sparse Hamiltonian to set of valid determinants for this excitation
            call sparse_csr_MM_row_part(H_c_all,H_p_all, I_k, H_c, H_p, N_det_l, nnz)

            nnz_l = H_p(N_det_l+1)-1

            ! move vectors to smaller allocations
            allocate(t_H_c(nnz_l), H_v(nnz_l))
            t_H_c = H_c(:nnz_l)
            call move_alloc(t_H_c, H_c)

            ! calculate entries of Hamiltonian
            call calc_sparse_zH(H_p, H_c, H_v, N_det_l, nnz_l, det_excited)
            ! call form_sparse_zH(H_p, H_c, H_v, s_max_sze,det_excited,&
            !                     iorb_R(iorb), 1, .true., N_det_l)
            call wall_time(t1)
            write(*, "(A33, F8.2, A10)"), "Sparse Hamiltonian calculated in ", t1-t0, " seconds"

            call wall_time(t0)
            ! calculate tridiagonalization of Hamiltoninan
            call lanczos_tridiag_sparse_reortho_c(H_v, H_c, H_p, psi_coef_excited(:,1),&
                                        alpha, beta,&
                                        lanczos_N, nnz, N_det_l)

            call wall_time(t1)
            write(*, "(A21, F8.2, A10)"), "Lanczos iteration in ", t1-t0, " seconds"
            ! prepare beta array for continued fractions
            bbeta(1) = (1.d0, 0.d0)
            do i = 2, lanczos_N
                bbeta(i) = -1.d0*beta(i)**2.0
            end do

            lanczos_alpha_R_complex(:, iorb, i_n_det) = alpha
            lanczos_beta_R_complex(:, iorb, i_n_det) = real(bbeta)

            epsilon = greens_epsilon ! small, a limit is taken here
            E0 = psi_energy(1)
            z = E0 - (-1.d0*greens_omega + (0.d0, 1.d0)*epsilon) ! omega is abs. energy value

            ! calculate greens functions
            !$OMP PARALLEL DO PRIVATE(i, zalpha) SHARED(alpha, bbeta, lanczos_N, greens_R_complex)&
            !$OMP SCHEDULE(GUIDED)
            do i = 1, greens_omega_N
                zalpha = z(i) - alpha
                greens_R_complex(i, iorb, i_n_det) = -1.d0*cfraction_c((0.d0, 0.d0), bbeta, zalpha, lanczos_N)
            end do
            !$OMP END PARALLEL DO

            if (dump_intermediate_output) then
                print *, "Dumping intermediate outputs to: "
                ! write out abcissa
                write(filename, *), "omega_R_complex.out"
                print *, filename
                call dump_array_real(filename, greens_omega, greens_omega_N)

                ! if requested, write out greens functions
                if (write_greens_f) then
                    write(filename, "(A, A, I4.4, A, I10.10, A)"), "greens_R_complex",&
                                                        "_iorb_", iorb_R(iorb), &
                                                        "_ndets_", N_det_l,&
                                                        ".out"
                    print *, filename
                    call dump_array_complex(filename,&
                                            greens_R_complex(:,iorb,i_n_det), greens_omega_N)
                end if

                ! write out spectral density
                write(filename, "(A, A, I4.4, A, I10.10, A)"), "spectral_density_R_complex",&
                                                    "_iorb_", iorb_R(iorb), &
                                                    "_ndets_", N_det_l,&
                                                    ".out"
                print *, filename

                pi = acos(-1.d0)
                call dump_array_real(filename, (-1.d0/pi) * aimag(greens_R_complex(:, iorb, i_n_det)), greens_omega_N)

                ! if requested, write out lanczos vectors

                if (write_lanczos_ab) then
                    write(filename, "(A, A, I4.4, A, I10.10, A, I6.6, A)"), "lanczos_alpha_R_complex",&
                                                        "_iorb_", iorb_R(iorb),&
                                                        "_ndets_", N_det_l,&
                                                        "_nvecs_", lanczos_N,&
                                                        ".out"
                    print *, filename
                    call dump_array_real(filename, lanczos_alpha_R_complex(:, iorb, i_n_det), lanczos_N)

                    write(filename, "(A, A, I4.4, A, I10.10, A, I6.6, A)"), "lanczos_beta_R_complex",&
                                                        "_iorb_", iorb_R(iorb),&
                                                        "_ndets_", N_det_l,&
                                                        "_nvecs_", lanczos_N,&
                                                        ".out"
                    print *, filename
                    call dump_array_real(filename, lanczos_beta_R_complex(:, iorb, i_n_det), lanczos_N)
                end if
            end if

            deallocate(H_c, H_v, H_p)

            call write_time(0)
            call print_memory_usage()
        end do

        deallocate(psi_coef_excited, det_excited,  H_c_all, H_p_all, I_k)
    end do
END_PROVIDER

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!  UTILITIES !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!! reduced versions of excitation routines below
subroutine build_A_wavefunction(i_particle, ispin, coef_out, det_out, N_det_l, I_k)
    implicit none
    BEGIN_DOC
    ! Applies the creation operator : a^{dager}_(i_particle) of
    ! spin = ispin to the current wave function (psi_det, psi_coef)
    END_DOC
    integer, intent(in)            :: i_particle,ispin,N_det_l
    integer, intent(out)           :: I_k(N_det_l)
    integer(bit_kind), intent(out) :: det_out(N_int,2,N_det_l)
    double precision, intent(out)  :: coef_out(N_det_l,N_states)
  
    integer :: k, rank
    integer :: i_ok
    double precision :: phase

    I_k = 0
    rank = 0
    do k=1,N_det_l
      coef_out(k,:) = psi_coef(k,:)
      det_out(:,:,k) = psi_det(:,:,k)
      call add_electron(det_out(:,:,k),i_particle,ispin,i_ok)
      if (i_ok == 1) then
        call get_phase_ca(psi_det(:,:,k),i_particle,ispin,phase)
        coef_out(k,:) = phase * coef_out(k,:)
        I_k(k) = 1
        rank += 1
      else
        coef_out(k,:) = 0.d0
        det_out(:,:,k) = psi_det(:,:,k)
      endif
    enddo

    write(*, "(A20, I8, A1, I8)"), "Rank of N+1 space: ", rank, "/", N_det_l
end

subroutine build_R_wavefunction(i_hole, ispin, coef_out, det_out, N_det_l, I_k, out_rank)
    implicit none
    BEGIN_DOC
    ! Applies the annihilation operator: a_(i_hole) of
    ! spin = ispin to the current wave function (psi_det, psi_coef)
    END_DOC
    integer, intent(in)            :: i_hole,ispin,N_det_l
    integer, intent(out), optional :: out_rank
    integer, intent(out)           :: I_k(N_det_l)
    integer(bit_kind), intent(out) :: det_out(N_int,2,N_det_l)
    double precision, intent(out)  :: coef_out(N_det_l,N_states)
  
    integer :: k, rank
    integer :: i_ok
    double precision :: phase

    I_k = 0
    rank = 0
    do k=1,N_det_l
      coef_out(k,:) = psi_coef(k,:)
      det_out(:,:,k) = psi_det(:,:,k)
      call remove_electron(det_out(:,:,k),i_hole,ispin,i_ok)
      if (i_ok == 1) then
        call get_phase_ca(psi_det(:,:,k),i_hole,ispin,phase)
        coef_out(k,:) = phase * coef_out(k,:)
        I_k(k) = 1
        rank += 1
      else
        coef_out(k,:) = 0.d0
        det_out(:,:,k) = psi_det(:,:,k)
      endif
    enddo

    write(*, "(A20, I8, A1, I8)"), "Rank of N-1 space: ", rank, "/", N_det_l
    out_rank = rank
end

subroutine build_A_wavefunction_complex(i_particle, ispin, coef_out, det_out, N_det_l, I_k)
    implicit none
    BEGIN_DOC
    ! Applies the creation operator : a^{dager}_(i_particle) of
    ! spin = ispin to the current wave function (psi_det, psi_coef_complex)
    END_DOC
    integer, intent(in)            :: i_particle,ispin,N_det_l
    integer, intent(out)           :: I_k(N_det_l)
    integer(bit_kind), intent(out) :: det_out(N_int,2,N_det_l)
    complex*16, intent(out)       :: coef_out(N_det_l,N_states)
  
    integer :: k, rank
    integer :: i_ok
    double precision :: phase

    I_k = 0
    rank = 0
    do k=1,N_det_l
      coef_out(k,:) = psi_coef_complex(k,:)
      det_out(:,:,k) = psi_det(:,:,k)
      call add_electron(det_out(:,:,k),i_particle,ispin,i_ok)
      if (i_ok == 1) then
        call get_phase_ca(psi_det(:,:,k),i_particle,ispin,phase)
        coef_out(k,:) = phase * coef_out(k,:)
        I_k(k) = 1
        rank += 1
      else
        coef_out(k,:) = (0.d0, 0.d0)
        det_out(:,:,k) = psi_det(:,:,k)
      endif
    enddo

    write(*, "(A20, I8, A1, I8)"), "Rank of N+1 space: ", rank, "/", N_det_l
end

subroutine build_R_wavefunction_complex(i_hole,ispin,coef_out, det_out, N_det_l, I_k)
    implicit none
    BEGIN_DOC
    ! Applies the annihilation operator: a_(i_hole) of
    ! spin = ispin to the current wave function (psi_det, psi_coef_complex)
    END_DOC
    integer, intent(in)            :: i_hole,ispin,N_det_l
    integer, intent(out)           :: I_k(N_det_l)
    integer(bit_kind), intent(out) :: det_out(N_int,2,N_det_l)
    complex*16, intent(out)       :: coef_out(N_det_l,N_states)
  
    integer :: k, rank
    integer :: i_ok
    double precision :: phase

    I_k = 0
    rank = 0
    do k=1,N_det_l
      coef_out(k,:) = psi_coef_complex(k,:)
      det_out(:,:,k) = psi_det(:,:,k)
      call remove_electron(det_out(:,:,k),i_hole,ispin,i_ok)
      if (i_ok == 1) then
        call get_phase_ca(psi_det(:,:,k),i_hole,ispin,phase)
        coef_out(k,:) = phase * coef_out(k,:)
        I_k(k) = 1
        rank += 1
      else
        coef_out(k,:) = (0.d0, 0.d0)
        det_out(:,:,k) = psi_det(:,:,k)
      endif
    enddo
    write(*, "(A20, I8, A1, I8)"), "Rank of N-1 space: ", rank, "/", N_det_l
end

subroutine build_coupled_wavefunction(intrinsic_psi_coef, coupled_psi_coef, I_cut_k, I_det_k, I_det_ind, N_det_s, N_det_t, n_iorb, hash_alpha, hash_beta, hash_vals, hash_prime, ht_size)
    implicit none
    BEGIN_DOC

    END_DOC

    !! routine arguments
    integer, intent(in)           :: N_det_s, N_det_t, n_iorb, ht_size, hash_prime, hash_vals(ht_size), I_cut_k(N_det_s, n_iorb), I_det_ind(N_det_s, n_iorb)
    integer(bit_kind), intent(in) :: hash_alpha(N_int, ht_size), hash_beta(N_int, ht_size), I_det_k(N_int, 2, N_det_s, n_iorb)
    double precision, intent(in)  :: intrinsic_psi_coef(N_det_s, N_states, n_iorb)

    double precision, intent(out) :: coupled_psi_coef(N_det_t, N_states)

    !! internal variables
    integer           :: i, j, k, target_row
    integer(bit_kind) :: target_row_det(N_int, 2)
    double precision  :: norm
    
    !! function calls
    integer           :: hash_value
    double precision  :: dnrm2

    ! initialize output wavefunction
    coupled_psi_coef = 0.0

    ! iterate over excitations, source determinants and sum coefficients by target generated deteterminant
    do i = 1, n_iorb
        do j = 1, N_det_s
            if (I_cut_k(j, i) == 1) then
                ! target_row_det = I_det_k(:, :, j, i)
                ! target_row = hash_value(hash_alpha, hash_beta, hash_vals, hash_prime, ht_size, n_iorb, &
                !             target_row_det(:,1), target_row_det(:,2))

                target_row = I_det_ind(j, i)

                do k = 1, N_states
                    coupled_psi_coef(target_row, k) += intrinsic_psi_coef(j, k, i)
                end do
            end if
        end do
    end do

    ! normalize wavefunction
    do k = 1, N_states
        norm = dnrm2(N_det_t, coupled_psi_coef(:,k), 1)
        coupled_psi_coef(:, k) = coupled_psi_coef(:, k) / norm
    end do

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

subroutine set_ref_bitmask(ispin, ac_type)
    implicit none
    BEGIN_DOC
    ! Reset the bitmask and bitmask energy for added/removed electron
    ! iorb is orbital index
    ! ispin is spin of hole/particle
    ! ac_type == F if adding electron, T if removing
    END_DOC

    integer :: i, j, a, b, occ(N_int*bit_kind_size,2)
    integer, intent(in)  :: ispin
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
    ! unlike provider, no call to bitstring to list since occ is straightforward
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
    ! this should be the same for real/complex 
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
    call provide_fock_op_cshell_ref_bitmask
    ! PROVIDE fock_op_cshell_ref_bitmask
end


subroutine set_ref_bitmask_complex(ispin, ac_type)
    implicit none
    BEGIN_DOC
    ! Reset the bitmask and bitmask energy for added/removed electron
    ! iorb is orbital index
    ! ispin is spin of hole/particle
    ! ac_type == F if adding electron, T if removing
    END_DOC

    integer :: i, j, k, a, b, occ(N_int*bit_kind_size,2), kpt, korb
    integer :: n_occ_ab(2)
    integer, intent(in)  :: ispin
    logical, intent(in)  :: ac_type
    integer, allocatable :: t_occ(:)

    !! Resetting bitmasks
    ! adjust alpha/beta sizes based on annihilation vs creation and spin type
    ! should this be divided by number of kpts?
    a = elec_alpha_num + merge(merge(-1, 1, ac_type),0,ispin==1)
    b = elec_beta_num  + merge(merge(-1, 1, ac_type),0,ispin==2)

    elec_num_tab(1) = a
    elec_num_tab(2) = b
    elec_num = a + b

    allocate(t_occ(max(a,b)))

    kpt = 1
    korb = 1
    do i = 1, max(a,b)
        t_occ(i) = korb + (kpt - 1) * mo_num_per_kpt
        kpt += 1
        if (kpt > kpt_num) then
            kpt = 1
            korb += 1
        end if
    end do

    call list_to_bitstring(ref_bitmask(1,1), t_occ(:a), a, N_int)
    call list_to_bitstring(ref_bitmask(1,2), t_occ(:b), b, N_int)

    deallocate(t_occ)

    do k = 1, kpt_num
        do i = 1, N_int
            ref_bitmask_kpts(i,1,k) = iand(ref_bitmask(i,1), kpts_bitmask(i,k))
            ref_bitmask_kpts(i,2,k) = iand(ref_bitmask(i,2), kpts_bitmask(i,k))
        end do
    end do

    ! reset elec numbers
    elec_alpha_num_kpts = 0
    elec_beta_num_kpts = 0
    kpt = 1
    if (a < b) then
        do i = 1, a
            elec_alpha_num_kpts(kpt) += 1
            elec_beta_num_kpts(kpt) += 1
            kpt += 1
            if (kpt > kpt_num) kpt = 1
        end do
        do i = a+1, b
            elec_beta_num_kpts(kpt) += 1
            kpt += 1
            if (kpt > kpt_num) kpt = 1
        end do
    else if ( b < a) then
        do i = 1, b
            elec_alpha_num_kpts(kpt) += 1
            elec_beta_num_kpts(kpt) += 1
            kpt += 1
            if (kpt > kpt_num) kpt = 1
        end do
        
        do i = b+1, a
            elec_alpha_num_kpts(kpt) += 1
            kpt += 1
            if (kpt > kpt_num) kpt = 1
        end do
    end if

    ! initialize closed shells
    do i = 1, N_int
        ref_closed_shell_bitmask(i,1) = ref_bitmask(i,1)
        ref_closed_shell_bitmask(i,2) = ref_bitmask(i,2)
    end do

    ! clear out extra orbitals
    occ = 0
    if (a < b) then
        do k = 1, kpt_num
            call bitstring_to_list_ab(ref_bitmask_kpts(1,1,k),occ, n_occ_ab,N_int)
            do j = elec_alpha_num_kpts(k)+1, elec_beta_num_kpts(k) 
                i = occ(j,2)
                call clear_bit_to_integer(i, ref_closed_shell_bitmask(1,2),N_int)
            end do
        end do
    else if (b < a) then
        do k = 1, kpt_num
            call bitstring_to_list_ab(ref_bitmask_kpts(1,1,k),occ, n_occ_ab,N_int)
            do j = elec_beta_num_kpts(k)+1, elec_alpha_num_kpts(k) 
                i = occ(j,1)
                call clear_bit_to_integer(i, ref_closed_shell_bitmask(1,1),N_int)
            end do
        end do
    end if

    ! set kpt version
    do k = 1, kpt_num
        do i = 1, N_int
            ref_closed_shell_bitmask_kpts(i,1,k) = iand(ref_closed_shell_bitmask(i,1),kpts_bitmask(i,k))
            ref_closed_shell_bitmask_kpts(i,2,k) = iand(ref_closed_shell_bitmask(i,2),kpts_bitmask(i,k))
        end do
    end do

    !! Reset the energy
    ! this should be the same for real/complex 
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
    PROVIDE fock_op_cshell_ref_bitmask_kpts

end

subroutine dump_array_int(fname, arr, sze)
    implicit none

    integer, intent(in) :: sze
    integer, intent(in) :: arr(sze)
    character(len=72), intent(in) :: fname


    integer :: i

    open(1, file=fname, action="WRITE")
    do i = 1, sze
        write(1, "(I8)"), arr(i)
    end do

    close(1)

end

subroutine dump_array_real(fname, arr, sze)
    implicit none

    integer, intent(in) :: sze
    double precision, intent(in) :: arr(sze)
    character(len=72), intent(in) :: fname


    integer :: i

    open(1, file=fname, action="WRITE")
    do i = 1, sze
        write(1, "(E16.8)"), arr(i)
    end do

    close(1)

end

subroutine dump_array_complex(fname, arr, sze)
    implicit none

    integer, intent(in) :: sze
    complex*16, intent(in) :: arr(sze)
    character(len=72), intent(in) :: fname


    integer :: i

    open(1, file=fname, action="WRITE")
    do i = 1, sze
        write(1, "(E16.8, E16.8)"), real(arr(i)), aimag(arr(i))
    end do

    close(1)

end