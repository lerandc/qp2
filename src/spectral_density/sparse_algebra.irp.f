subroutine measure_sparsity(finished)
    implicit none

    logical, intent(out) :: finished
    integer              :: n, i, j, nnz, nnz_tot
    integer, allocatable :: n_list(:), n_test(:)
    double precision     :: hij

    allocate(n_list(2))
    n_list = (/1000, 25000/)!, 50000, 100000, 250000, 500000, N_det/)
    
    allocate(n_test(3))
    n_test(:2) = n_list
    n_test(3) = 1000

    call move_alloc(n_test, n_list)

    call i_H_j(psi_det(1,1,1), psi_det(1,1,1), N_int, hij) ! force provide early
    do n = 1, size(n_list, 1)
        
        !$OMP PARALLEL PRIVATE(i,j,hij,nnz) SHARED(n, psi_det, n_list, nnz_tot) 
        nnz_tot = 0
        nnz = 0
        !$OMP DO SCHEDULE(GUIDED)
        do i = 1, n_list(n)
            do j = i, n_list(n)
                call i_H_j(psi_det(1,1,i), psi_det(1,1,j),N_int, hij)
                if (abs(hij) > 0) then
                    nnz += 1
                end if
            end do
        end do
        !$OMP END DO

        !$OMP CRITICAL
        nnz_tot = nnz_tot + nnz
        !$OMP END CRITICAL
        !$OMP END PARALLEL

        print *, "N_det: ", n_list(n), "nnz: ", 2*nnz_tot, "frac: ", (nnz_tot*2.d0)/(n_list(n)**2.0)
    end do

    finished = .TRUE.
end


! to convert to provider, provide all 4 array for intel MKL CSR representation
subroutine form_sparse_H(finished)
    ! use MKL_SPBLAS

    implicit none
    BEGIN_DOC
    ! Form a compressed sparse row matrix representation of the Hamiltonian
    ! in the space of the determinants
    END_DOC

    logical              :: finished
    integer              :: n, i, j, k
    integer              :: nnz, nnz_tot
    integer              :: n_vals, n_threads, ID
    integer, allocatable :: nnz_arr(:), coo_r(:), coo_c(:), k_arr(:)
    integer, allocatable :: coo_r_all(:), coo_c_all(:), coo_r_t(:), coo_c_t(:)
    integer :: OMP_get_num_threads, OMP_get_thread_num
    double precision     :: hij, frac
    double precision, allocatable :: coo_v(:), coo_v_t(:), coo_v_all(:), h_ref(:,:)

    ! force provide early so that threads don't each try to provide
    call i_H_j(psi_det(1,1,1), psi_det(1,1,1), N_int, hij) 
    
    !$OMP PARALLEL SHARED(n, nnz_tot, nnz_arr, n_threads, h_ref) &
    !$OMP PRIVATE(i,j,k, hij, nnz, ID, coo_r, coo_c, coo_v, coo_r_t, coo_c_t, coo_v_t) 
    !$ n_threads = OMP_get_num_threads()
    !$ ID = OMP_get_thread_num() + 1
    n = 8
    frac = 0.2
    n_vals = max(nint(n*n*frac/n_threads), 128)

    allocate(coo_r(n_vals))
    allocate(coo_c(n_vals))
    allocate(coo_v(n_vals))

    !$OMP SINGLE
    allocate(nnz_arr(n_threads), h_ref(n,n))
    h_ref = 0.0
    nnz_tot = 0
    !$OMP END SINGLE
    nnz = 0

    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, n
        do j = i, n
            call i_H_j(psi_det(1,1,i), psi_det(1,1,j),N_int, hij)
            if (abs(hij) > 0) then
                nnz += 1

                if (nnz <= size(coo_r, 1)) then
                  coo_r(nnz) = i
                  coo_c(nnz) = j
                  coo_v(nnz) = hij
                  h_ref(i,j) = hij
                  h_ref(j,i) = hij
                else
                  ! dynamically increase array size on thread
                  print *, 'allocating'
                  allocate(coo_r_t(nnz + 1024))
                  allocate(coo_c_t(nnz + 1024))
                  allocate(coo_v_t(nnz + 1024))
                  
                  coo_r_t(:nnz-1) = coo_r
                  coo_c_t(:nnz-1) = coo_c
                  coo_v_t(:nnz-1) = coo_v
                  
                  coo_r_t(nnz) = i
                  coo_c_t(nnz) = j
                  coo_v_t(nnz) = hij
                  h_ref(i,j) = hij
                  h_ref(j,i) = hij

                  call move_alloc(coo_r_t, coo_r)
                  call move_alloc(coo_c_t, coo_c)
                  call move_alloc(coo_v_t, coo_v)
                end if

            end if
        end do
    end do
    !$OMP END DO

    nnz_arr(ID) = nnz

    !$OMP CRITICAl
    nnz_tot = nnz_tot + nnz
    !$OMP END CRITICAl

    !$OMP BARRIER

    !$OMP SINGLE
    print *, nnz_arr
    allocate(coo_r_all(nnz_tot))
    allocate(coo_c_all(nnz_tot))
    allocate(coo_v_all(nnz_tot))
    !$OMP END SINGLE
    
    k = 0
    do i = 1, ID-1
        k += nnz_arr(i)
    end do

    do i = k+1, k + nnz_arr(ID)
        coo_r_all(i) = coo_r(i-k)
        coo_c_all(i) = coo_c(i-k)
        coo_v_all(i) = coo_v(i-k)
    end do
    !$OMP END PARALLEL

    print *, "n: ", n, "nnz: ", 2*nnz_tot-n

    ! do i = 1, nnz_tot
    !     print *, coo_r_all(i), coo_c_all(i), coo_v_all(i)
    ! end do


    ! ! convert to sparse matrix in CSR format
    ! type(SPARSE_MATRIX_T) :: sH
    ! integer               :: status, mkl_d_create_coo

    ! status = mkl_d_create_coo(sH, SPARSE_INDEX_BASE_ONE, &
    !                           n, n, nnz_tot,&
    !                           coo_r_all, coo_c_all, coo_v_all)

    ! print *, "Status: ", status

    ! Assumption: all diagonal entries are non-zero
    ! make CSR arrays from COO arrays
    integer              :: nnz_csr, ii, scn_a, kk
    integer, allocatable :: coo_s(:), coo_n(:), csr_s(:), csr_c(:)
    double precision, allocatable :: csr_v(:)

    !$OMP PARALLEL SHARED(n, nnz_csr, coo_r_all, coo_c_all, coo_v_all, csr_s, csr_c, csr_v, coo_s, coo_n, k_arr, n_threads)&
    !$OMP PRIVATE(i,j,ii, scn_a, ID)
    !$ n_threads = OMP_get_num_threads()
    !$ ID = OMP_get_thread_num() + 1

    nnz_csr = 2*nnz_tot-n! COO was taken from upper triangle, account for double counting of diagonal

    !$OMP SINGLE
    allocate(coo_s(n), coo_n(n), csr_s(n+1), csr_c(nnz_csr), csr_v(nnz_csr))
    !$OMP END SINGLE

    !$OMP BARRIER
    !$OMP SINGLE
    print *, "here"
    !$OMP END SINGLE

    ! calculate the starting index of each row in COO, since they were not sorted
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, n
        j = 1
        do while(coo_r_all(j) .ne. i)
            j += 1
        end do
        coo_s(i) = j
        
        k = 0
        do while(coo_r_all(j) ==  i)
            k +=1
            j += 1
        end do

        coo_n(i) = k
    end do
    !$OMP END DO

    
    ! calculate CSR pointer ranges
    ! row i data goes from csr_s(i) to csr_s(i+1) - 1
    ! first, count all entries in parallel, then perform scan to set pointer ranges
    ! then reduce with inclsuive scan

    !$OMP SINGLE
    csr_s(1) = 1
    csr_s(n+1) = 1 
    !$OMP END SINGLE

    !$OMP BARRIER
    !$OMP SINGLE
    print *, "here"
    !$OMP END SINGLE
    
    !$OMP DO SCHEDULE(GUIDED) 
    do i = 2, n+1
        k = 0 ! count all entries with column i-1
        do j = 1, nnz_tot
            k = k + merge(1,0,coo_c_all(j)==i-1)
        end do
        csr_s(i) = coo_n(i-1) + k - 1 ! account for double counting
    end do
    !$OMP END DO

    !$OMP BARRIER
    !$OMP SINGLE
    print *, "here"
    print *, csr_s
    scn_a = 0
    !$OMP SIMD REDUCTION(inscan, +:scn_a)
    do i = 1, n+1
        scn_a = scn_a + csr_s(i)
        !$OMP SCAN INCLUSIVE(scn_a)
        csr_s(i) = scn_a
    end do
    !$OMP END SINGLE

    
    !$OMP BARRIER
    !$OMP SINGLE
    print *, "here"
    print *, csr_s
    !$OMP END SINGLE

    ! loop through rows and construct CSR matrix

    !$OMP SINGLE
    allocate(k_arr(n_threads))
    k_arr = 0
    !$OMP END SINGLE

    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, n
        k_arr(ID) = 0
        print *, "ID ", ID, " has ", i
        do j = 1, i-1 !TODO: possible to only loop over number of entries?
            ! get pre-diagonal entries
            ! search jth row for column i
            do ii = coo_s(j)+1, coo_s(j) + coo_n(j) -1
                if (coo_c_all(ii) == i) then
                    csr_v(csr_s(i) + k_arr(ID)) = coo_v_all(ii)
                    csr_c(csr_s(i) + k_arr(ID)) = j !coo_c_all(ii)
                    k_arr(ID) += 1
                    exit
                end if
            end do
        ! get diagonal+ entries
        end do
        
        do ii = coo_s(i), coo_s(i) + coo_n(i) - 1
            csr_v(csr_s(i) + k_arr(ID)) = coo_v_all(ii)
            csr_c(csr_s(i) + k_arr(ID)) = coo_c_all(ii)
            k_arr(ID) += 1
        end do
    end do
    !$OMP END DO

    !$OMP BARRIER


    !$OMP END PARALLEL
    print *, "-------------- Reference H"
    do i = 1, n
        do j = 1, n
            print *, i, j, h_ref(i,j)
        end do
    end do
    
    print *, "-------------- COO representation"
    do i = 1, nnz_tot
        print *, coo_r_all(i), coo_c_all(i), coo_v_all(i)
    end do

    print *, "-------------- CSR representation"
    do i = 1, nnz_csr
        print *,  i, csr_c(i), csr_v(i)
    end do
    
    print *, "----- CSR pointers"
    do i = 1, n+1
        print *, csr_s(i)
    end do

    finished = .TRUE.
end

! subroutine sparse_mat_vec
