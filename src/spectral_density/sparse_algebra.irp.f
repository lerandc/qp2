BEGIN_PROVIDER [integer, nnz_max_per_row]
&BEGIN_PROVIDER [integer, nnz_max_tot]

    BEGIN_DOC
    ! Total possible number of entries in each column/row in full determinant space
    ! Vastly overestimates actual number needed, but a useful base heuristic.
    END_DOC
    implicit none

    integer*8     :: s_a, s_b, d_a, d_b

    s_a = elec_alpha_num * (mo_num - elec_alpha_num)
    s_b = elec_beta_num * (mo_num - elec_beta_num)
    
    d_a = ( elec_alpha_num * (elec_alpha_num - 1) / 2) * &
    ( (mo_num - elec_alpha_num) * (mo_num - elec_alpha_num - 1) / 2)
    
    d_b = ( elec_beta_num * (elec_beta_num - 1) / 2) * &
    ( (mo_num - elec_beta_num) * (mo_num - elec_beta_num - 1) / 2)
    
    nnz_max_per_row = 1 + s_a + s_b + s_a*s_b + d_a + d_b
    nnz_max_tot = N_det * nnz_max_per_row

END_PROVIDER

subroutine sparse_csr_dmv(A_v, A_c, A_p, x, y, sze, nnz)
    implicit none
    BEGIN_DOC
    ! Compute the matrix vector product y = Ax with sparse A in CSR format
    ! Not general. Assumes A is sze x sze and symmetric; x, y are sze
    !
    ! A_v is matrix values
    ! A_c is matrix column indices
    ! A_p are the row pointers
    ! x is an input vector of length sze
    ! y is an output vector of length sze
    ! nnz is the total number of nonzero entries in A
    END_DOC

    integer,          intent(in)  :: sze, nnz, A_c(nnz), A_p(sze+1)
    double precision, intent(in)  :: A_v(nnz), x(sze)
    double precision, intent(out) :: y(sze)
    double precision, allocatable :: y_t(:)
    integer                       :: i, j

    !$OMP PARALLEL PRIVATE(i, j, y_t) SHARED(y, x, A_c, A_p, A_v, sze)
    allocate(y_t(sze))
    y_t = 0.d0
    !$OMP BARRIER
    ! loop over rows
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, sze
        
        
        ! make sure row actually is in reduced determinant space
        if (A_p(i+1) - A_p(i) > 0) then
            ! calculate diagonal separately to avoid double counting
            ! first entry per column is guaranteed to be diagonal since all diagonal
            ! elements of H are nonzero
            j = A_p(i)
            y_t(i) = y_t(i) + A_v(j) * x(A_c(j))
        end if
        
        ! loop over columns
        do j = A_p(i)+1, A_p(i+1)-1
            ! calculate element of owned row
            y_t(i) = y_t(i) + A_v(j) * x(A_c(j))

            ! calculate element of owned column
            y_t(A_c(j)) = y_t(A_c(j)) + A_v(j) * x(i)
        end do
    end do
    !$OMP END DO

    ! perform reduction of final answer
    !$OMP CRITICAL
    y = y + y_t
    !$OMP END CRITICAL

    deallocate(y_t)
    !$OMP END PARALLEL
end

subroutine sparse_csr_zmv(A_v, A_c, A_p, x, y, sze, nnz)
    implicit none
    BEGIN_DOC
    ! Compute the matrix vector product y = Ax with sparse A in CSR format
    ! Not general. Assumes A is sze x sze; x, y are sze
    END_DOC

    integer,          intent(in)  :: sze, nnz, A_c(nnz), A_p(sze+1)
    complex*16,       intent(in)  :: A_v(nnz), x(sze)
    complex*16,       intent(out) :: y(sze)
    complex*16,      allocatable  :: y_t(:)
    integer                       :: i, j

    
    !$OMP PARALLEL PRIVATE(i, j, y_t) SHARED(y, x, A_c, A_p, A_v, sze)
    allocate(y_t(sze))
    y_t = (0.d0, 0.d0)
    !$OMP BARRIER
    ! loop over rows
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, sze
        
        
        ! make sure row actually is in reduced determinant space
        if (A_p(i+1) - A_p(i) > 0) then
            ! calculate diagonal separately to avoid double counting
            ! first entry per column is guaranteed to be diagonal since all diagonal
            ! elements of H are nonzero
            j = A_p(i)
            y_t(i) = y_t(i) + A_v(j) * x(A_c(j))
        end if
        
        ! loop over columns
        do j = A_p(i)+1, A_p(i+1)-1
            ! calculate element of owned row
            y_t(i) = y_t(i) + A_v(j) * x(A_c(j))

            ! calculate element of owned column
            y_t(A_c(j)) = y_t(A_c(j)) + conjg(A_v(j)) * x(i)
        end do
    end do
    !$OMP END DO

    ! perform reduction of final answer
    !$OMP CRITICAL
    y = y + y_t
    !$OMP END CRITICAL

    deallocate(y_t)
    !$OMP END PARALLEL
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Utilties for building and reutilizing sparse matrices !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine sparse_csr_MM(A_c, A_p, I_k, B_c, B_p, sze, nnz_in)
    BEGIN_DOC
    ! Matrix matrix multiplication routine for sparse matrices
    ! NOT general
    ! Performs A I = B, returning B in CSR format
    ! Performs the multiplication between a symmetric, CSR matrix of 1s/0s
    ! against a low-rank identity matrix
    ! Used for transferring the sparsity calculation of the full system to the
    ! reduced pattern when a hole or particle is introduced!
    ! 
    ! A_c are the column indices of the full Hamiltonian
    ! A_p are the row pointers of the full Hamiltonian
    ! I_k is a length sze vector representing the low-rank identity matrix,
    !       whose zeros are the determinants (nodes) removed from the full space
    ! B_c are the output column indices of the reduced Hamiltonian
    ! B_p are the output row pointers of the reduced Hamiltonian
    ! sze is the full rank of A
    ! nnz_in is the number of nonzeros in the full Hamiltonian
    END_DOC
    implicit none

    integer, intent(in) :: sze, nnz_in, A_c(nnz_in), A_p(sze+1), I_k(sze)
    integer, intent(out):: B_c(nnz_in), B_p(sze+1)
    integer :: i, j, k, ii, kk, lcol, nnz_cnt, old_row, nnz_tot
    integer              :: n_vals, n_threads, ID, scn_a, nnz
    integer :: OMP_get_num_threads, OMP_get_thread_num
    integer, allocatable :: nnz_arr(:), coo_r(:), coo_c(:)
    integer, allocatable :: coo_r_all(:), coo_c_all(:), coo_r_t(:), coo_c_t(:)
    integer, allocatable:: coo_s(:), csr_n(:)
    double precision     :: frac


    !$OMP PARALLEL SHARED(nnz_tot, nnz_arr, n_threads,&
    !$OMP                 coo_r_all, coo_c_all, coo_s, I_k, A_c, A_p, B_c, B_p)& 
    !$OMP PRIVATE(i,j,old_row, k,ii,kk, scn_a, ID, nnz, nnz_cnt, coo_r, coo_c, coo_r_t, coo_c_t, lcol) 


    !$ n_threads = OMP_get_num_threads()
    !$ ID = OMP_get_thread_num() + 1

    !$OMP SINGLE
    allocate(coo_s(sze))
    B_p = 0
    !$OMP END SINGLE
    !$OMP BARRIER

    frac = 0.2
    n_vals = max(nint(sze*sze*frac/n_threads), 128)

    allocate(coo_r(n_vals))
    allocate(coo_c(n_vals))

    !$OMP SINGLE
    allocate(nnz_arr(n_threads))
    nnz_tot = 0
    !$OMP END SINGLE


    ! loop over rows
    nnz_cnt = 0
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, sze

        if (I_k(i) == 0) cycle ! node is emoved from graph

        nnz = A_p(i+1)-1 - A_p(i) ! maximum reallocation size

        if (nnz_cnt + nnz > size(coo_r, 1)) then
            allocate(coo_r_t(nnz_cnt + 10*nnz))
            allocate(coo_c_t(nnz_cnt + 10*nnz))
            
            coo_r_t(:nnz_cnt) = coo_r
            coo_c_t(:nnz_cnt) = coo_c
            
            call move_alloc(coo_r_t, coo_r)
            call move_alloc(coo_c_t, coo_c)
        end if

        ! loop over columns in row
        do j = A_p(i), A_p(i+1)-1
            lcol = A_c(j)
            if (I_k(lcol) == 1) then ! current column still in space
                B_p(i+1) += 1 ! increase pointer spacing
                nnz_cnt += 1
                coo_r(nnz_cnt) = i ! store rows/columns in COO buffers
                coo_c(nnz_cnt) = lcol
            end if

        end do

    end do
    !$OMP END DO

    nnz_arr(ID) = nnz_cnt

    !$OMP CRITICAl
    nnz_tot = nnz_tot + nnz_cnt
    !$OMP END CRITICAl
    !$OMP BARRIER
    
    !$OMP SINGLE
    print *, "Total non-zero entries in Hamiltonian: ", nnz_tot, " max size:", nnz_in
    !$OMP END SINGLE

    !$OMP SINGLE
    allocate(coo_r_all(nnz_tot))
    allocate(coo_c_all(nnz_tot))
    !$OMP END SINGLE

    k = 0
    do i = 1, ID-1
        k += nnz_arr(i)
    end do

    ! coo_*_all are in shared memory; since threads have disjoint work, this is safe
    coo_r_all(k+1:k+nnz_arr(ID)) = coo_r
    coo_c_all(k+1:k+nnz_arr(ID)) = coo_c
    !$OMP BARRIER

    ! calculate the starting index of each row in COO, since there is no sorting guarantee
    ! iterate over range, keep track of current row; if row changes, record the row start
    ii = (nnz_tot / n_threads)*(ID-1) + 1 ! start of range
    kk = (nnz_tot / n_threads) * ID + 3 ! end of range, slight overlap to catch boundary pointers
    if (ID == n_threads) kk = nnz_tot
    old_row = coo_r_all(ii)

    do i = ii, kk 
        if(coo_r_all(i) .ne. old_row) then 
            coo_s(coo_r_all(i)) = i
            old_row = coo_r_all(i)
        end if 
    end do

    ! make sure first row in COO is accounted for
    !$OMP SINGLE
    coo_s(coo_r_all(1)) = 1
    !$OMP END SINGLE

    ! calculate CSR pointer ranges
    ! row i data goes from csr_s(i) to csr_s(i+1) - 1
    ! use inclusive scan to reduce counts of rows in B_p to set pointer ranges
    !$OMP SINGLE
    B_p(1) = 1
    !$OMP END SINGLE
    
    ! need to strongly enforce synchronization here
    !$OMP BARRIER
    !$OMP SINGLE
    scn_a = 0
    !$OMP SIMD REDUCTION(inscan, +:scn_a)
    do i = 1, sze+1
        scn_a = scn_a + B_p(i)
        !$OMP SCAN INCLUSIVE(scn_a)
        B_p(i) = scn_a
    end do
    !$OMP END SINGLE
    !$OMP BARRIER

    ! loop through rows and construct CSR matrix from COO temp arrays
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, sze
        B_c(B_p(i):B_p(i+1)-1) = coo_c_all(coo_s(i):coo_s(i+1)-1)
    end do
    !$OMP END DO

    deallocate(coo_r, coo_c)

    !$OMP SINGLE
    deallocate(coo_s, nnz_arr, coo_r_all, coo_c_all)
    !$OMP END SINGLE
    !$OMP END PARALLEL
    
end

subroutine sparse_csr_MM_row_part(A_c, A_p, I_k, B_c, B_p, sze, nnz_in)
    BEGIN_DOC
    ! Matrix matrix multiplication routine for sparse matrices
    ! With row partitioning during worksharing
    ! NOT general
    ! Performs A I = B, returning B in CSR format
    ! Performs the multiplication between a symmetric, CSR matrix of 1s/0s
    ! against a low-rank identity matrix
    ! Used for transferring the sparsity calculation of the full system to the
    ! reduced pattern when a hole or particle is introduced!
    ! 
    ! A_c are the column indices of the full Hamiltonian
    ! A_p are the row pointers of the full Hamiltonian
    ! I_k is a length sze vector representing the low-rank identity matrix,
    !       whose zeros are the determinants (nodes) removed from the full space
    ! B_c are the output column indices of the reduced Hamiltonian
    ! B_p are the output row pointers of the reduced Hamiltonian
    ! sze is the full rank of A
    ! nnz_in is the number of nonzeros in the full Hamiltonian
    END_DOC
    implicit none

    integer, intent(in) :: sze, nnz_in, A_c(nnz_in), A_p(sze+1), I_k(sze)
    integer, intent(out):: B_c(nnz_in), B_p(sze+1)
    integer :: i, j, k, ii, kk, lcol, nnz_cnt, old_row, nnz_tot
    integer              :: n_vals, n_threads, ID, scn_a, nnz
    integer :: OMP_get_num_threads, OMP_get_thread_num
    integer, allocatable :: nnz_arr(:), coo_r(:), coo_c(:)
    integer, allocatable :: coo_r_all(:), coo_c_all(:), coo_r_t(:), coo_c_t(:)
    integer, allocatable:: coo_s(:), csr_n(:)
    double precision     :: frac

    integer :: target_N, row_start, row_end
    integer, allocatable :: pointer_blocks(:), row_starts(:)

    !$OMP PARALLEL SHARED(nnz_tot, nnz_arr, n_threads, target_N, nnz_in, pointer_blocks, row_starts,&
    !$OMP                 coo_r_all, coo_c_all, coo_s, I_k, A_c, A_p, B_c, B_p)& 
    !$OMP PRIVATE(i,j,old_row, k,ii,kk, scn_a, ID, n_vals, nnz, nnz_cnt, coo_r, coo_c, coo_r_t, coo_c_t, lcol, row_start, row_end) 


    !$ n_threads = OMP_get_num_threads()
    !$ ID = OMP_get_thread_num() + 1

    !$OMP SINGLE
    allocate(coo_s(sze), pointer_blocks(sze+1), row_starts(n_threads+1))
    pointer_blocks = 0
    row_starts = 0
    B_p = 0
    target_N = nnz_in / n_threads
    !$OMP END SINGLE
    !$OMP BARRIER

    ! split the work so that different groups of contiguous rows have roughly equal entries
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, sze+1
        pointer_blocks(i) = A_p(i) / target_N + 1
    end do
    !$OMP END DO

    ii = (sze / n_threads)*(ID-1) + 1 ! start of range
    kk = (sze / n_threads) * ID + 3 ! end of range, slight overlap to catch boundary pointers
    if (ID == n_threads) kk = sze
    old_row = pointer_blocks(ii)
    
    do i = ii, kk 
        if(pointer_blocks(i) .ne. old_row) then 
            row_starts(pointer_blocks(i)) = i
            old_row = pointer_blocks(i)
        end if 
    end do
    
    !$OMP SINGLE
    row_starts(1) = 1
    row_starts(n_threads+1) = sze+1
    !$OMP END SINGLE
    
    row_start = row_starts(ID)
    row_end = row_starts(ID+1) - 1
    if (ID == n_threads) row_end = sze

    ! allocate arrays and start main work loops
    n_vals = A_p(row_end+1)-A_p(row_start)

    allocate(coo_r(n_vals))
    allocate(coo_c(n_vals))
    coo_c = 0
    coo_r = 0

    !$OMP SINGLE
    allocate(nnz_arr(n_threads))
    nnz_tot = 0
    !$OMP END SINGLE


    ! loop over rows
    nnz_cnt = 0
    do i = row_start, row_end

        if (I_k(i) == 0) cycle ! node is removed from graph

        ! nnz = A_p(i+1)-1 - A_p(i) ! maximum reallocation size

        ! loop over columns in row
        do j = A_p(i), A_p(i+1)-1
            lcol = A_c(j)
            if (I_k(lcol) == 1) then ! current column still in space
                B_p(i+1) += 1 ! increase pointer spacing
                nnz_cnt += 1
                coo_r(nnz_cnt) = i ! store rows/columns in COO buffers
                coo_c(nnz_cnt) = lcol
            end if

        end do

    end do

    nnz_arr(ID) = nnz_cnt

    !$OMP CRITICAl
    nnz_tot = nnz_tot + nnz_cnt
    !$OMP END CRITICAl
    !$OMP BARRIER
    
    !$OMP SINGLE
    print *, "Total non-zero entries in Hamiltonian: ", nnz_tot, " max size:", nnz_in
    !$OMP END SINGLE
    
    !$OMP SINGLE
    allocate(coo_r_all(nnz_tot))
    allocate(coo_c_all(nnz_tot))
    coo_c_all = 0
    coo_r_all = 0
    !$OMP END SINGLE
    !$OMP BARRIER
    
    k = 0
    do i = 1, ID-1
        k += nnz_arr(i)
    end do
    
    ! coo_*_all are in shared memory; since threads have disjoint work, this is safe
    coo_r_all(k+1:k+nnz_arr(ID)) = coo_r(1:nnz_cnt)
    coo_c_all(k+1:k+nnz_arr(ID)) = coo_c(1:nnz_cnt)
    !$OMP BARRIER

    ! calculate the starting index of each row in COO, since there is no sorting guarantee
    ! iterate over range, keep track of current row; if row changes, record the row start
    ii = (nnz_tot / n_threads)*(ID-1) + 1 ! start of range
    kk = (nnz_tot / n_threads) * ID + 3 ! end of range, slight overlap to catch boundary pointers
    if (ID == n_threads) kk = nnz_tot
    old_row = coo_r_all(ii)
    
    do i = ii, kk 
        if(coo_r_all(i) .ne. old_row) then 
            coo_s(coo_r_all(i)) = i
            old_row = coo_r_all(i)
        end if 
    end do

    ! make sure first row in COO is accounted for
    !$OMP SINGLE
    coo_s(coo_r_all(1)) = 1
    !$OMP END SINGLE
    
    ! calculate CSR pointer ranges
    ! row i data goes from csr_s(i) to csr_s(i+1) - 1
    ! use inclusive scan to reduce counts of rows in B_p to set pointer ranges
    !$OMP SINGLE
    B_p(1) = 1
    !$OMP END SINGLE
    
    ! need to strongly enforce synchronization here
    !$OMP BARRIER
    !$OMP SINGLE
    scn_a = 0
    !$OMP SIMD REDUCTION(inscan, +:scn_a)
    do i = 1, sze+1
        scn_a = scn_a + B_p(i)
        !$OMP SCAN INCLUSIVE(scn_a)
        B_p(i) = scn_a
    end do
    !$OMP END SINGLE
    !$OMP BARRIER
    
    ! loop through rows and construct CSR matrix from COO temp arrays
    do i = row_start, row_end
        B_c(B_p(i):B_p(i+1)-1) = coo_c_all(coo_s(i):coo_s(i+1)-1)
    end do
    
    deallocate(coo_r, coo_c)
    
    !$OMP BARRIER
    !$OMP SINGLE
    deallocate(coo_s, nnz_arr, coo_r_all, coo_c_all, pointer_blocks, row_starts)
    !$OMP END SINGLE
    !$OMP END PARALLEL
    
end

subroutine get_sparsity_structure(csr_s, csr_c, sze, N_det_l)
    implicit none
    BEGIN_DOC
    ! Form a compressed sparse row matrix representation of the Hamiltonian
    ! connectivity in the space of the determinants. 

    ! csr_s are the row pointers
    ! csr_c are the column indices
    ! sze is the maximum possible number of nonzero entries
    ! N_det_l is the largest index of the determinants to include (when sorted by energy) 
    END_DOC

    integer, intent(in)           :: N_det_l
    integer(kind=8), intent(in)        :: sze
    integer, intent(out)          :: csr_s(N_det_l+1), csr_c(sze)

    integer              :: n, i, j, k, l_row, old_row
    integer              :: nnz, nnz_cnt, nnz_tot
    integer              :: n_vals, n_vals_row, n_threads, ID
    integer              :: nnz_csr, ii, scn_a, kk
    integer :: OMP_get_num_threads, OMP_get_thread_num
    integer, allocatable :: nnz_arr(:), coo_r(:), coo_c(:), l_cols(:)
    integer, allocatable :: coo_r_all(:), coo_c_all(:), coo_r_t(:), coo_c_t(:)
    integer, allocatable:: coo_s(:), coo_n(:)
    double precision     :: frac
    

    !$OMP PARALLEL SHARED(nnz_tot, nnz_arr, nnz_csr, n_threads, N_det, N_det_l, nnz_max_per_row, n_vals_row,&
    !$OMP                 coo_r_all, coo_c_all, csr_s, csr_c, coo_s, coo_n)& 
    !$OMP PRIVATE(i,j,old_row, k,ii,kk, scn_a, ID, nnz, nnz_cnt, coo_r, coo_c, coo_r_t, coo_c_t, l_cols, l_row) 

    !$ n_threads = OMP_get_num_threads()
    !$ ID = OMP_get_thread_num() + 1

    !$OMP SINGLE
    allocate(coo_s(N_det_l), coo_n(N_det_l))
    coo_n = 0
    !$OMP END SINGLE
    !$OMP BARRIER

    ! initial allocation sizes for vectors
    frac = 0.2
    n_vals = max(nint(N_det_l*N_det_l*frac/n_threads), 128)
    n_vals_row = nnz_max_per_row

    allocate(coo_r(n_vals))
    allocate(coo_c(n_vals))
    allocate(l_cols(n_vals_row))
    
    !$OMP SINGLE
    allocate(nnz_arr(n_threads))
    nnz_tot = 0
    !$OMP END SINGLE
    
    nnz_cnt = 0
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det ! this loop needs to go over all the determinants, since this loop is not in determinant order but rather k_a order
        nnz = 0
        l_cols = 0

        ! check if row in N+1/N-1 determinant space
        ! if so, grab all the indices of all nonzero columns for that row in the upper half of H
        call get_all_sparse_columns(i, l_cols, l_row, nnz, n_vals_row, N_det_l)

        if (nnz == 0) cycle

        ! reallocate arrays if necessary
        if (nnz_cnt + nnz > size(coo_r, 1)) then
            allocate(coo_r_t(nnz_cnt + 10*nnz))
            allocate(coo_c_t(nnz_cnt + 10*nnz))
            
            coo_r_t(:nnz_cnt) = coo_r
            coo_c_t(:nnz_cnt) = coo_c
            
            call move_alloc(coo_r_t, coo_r)
            call move_alloc(coo_c_t, coo_c)
        end if
        
        
        ! Calculate nonzero entries and temperorarily store in COO format
        coo_n(l_row) = nnz ! store number of entries in row for later
        do j = 1, nnz
            nnz_cnt += 1
            coo_r(nnz_cnt) = l_row
            coo_c(nnz_cnt) = l_cols(j)
        end do
    end do
    !$OMP END DO

    nnz_arr(ID) = nnz_cnt

    !$OMP CRITICAl
    nnz_tot = nnz_tot + nnz_cnt
    !$OMP END CRITICAl
    !$OMP BARRIER
    
    !$OMP SINGLE
    print *, "Total non-zero entries in full Hamiltonian: ", nnz_tot, " max size:", sze
    !$OMP END SINGLE
    
    !$OMP SINGLE
    allocate(coo_r_all(nnz_tot))
    allocate(coo_c_all(nnz_tot))
    !$OMP END SINGLE
    !$OMP BARRIER
    
    k = 0
    do i = 1, ID-1
        k += nnz_arr(i)
    end do

    ! coo_*_all are in shared memory; since threads have disjoint work, this is safe
    coo_r_all(k+1:k+nnz_arr(ID)) = coo_r
    coo_c_all(k+1:k+nnz_arr(ID)) = coo_c
    !$OMP BARRIER

    ! calculate the starting index of each row in COO, since there is no sorting guarantee
    ! iterate over range, keep track of current row; if row changes, record the row start
    ii = (nnz_tot / n_threads)*(ID-1) + 1 ! start of range
    kk = (nnz_tot / n_threads) * ID + 3 ! end of range, slight overlap to catch boundary pointers
    if (ID == n_threads) kk = nnz_tot
    old_row = coo_r_all(ii)

    do i = ii, kk 
        if(coo_r_all(i) .ne. old_row) then 
            coo_s(coo_r_all(i)) = i
            old_row = coo_r_all(i)
        end if 
    end do

    ! make sure first row in COO is accounted for
    !$OMP SINGLE
    coo_s(coo_r_all(1)) = 1
    !$OMP END SINGLE

    ! calculate CSR pointer ranges
    ! row i data goes from csr_s(i) to csr_s(i+1) - 1
    ! use inclusive scan to reduce counts of rows in coo_n to set pointer ranges

    !$OMP SINGLE
    csr_s(1) = 1
    csr_s(2:N_det_l+1) = coo_n
    !$OMP END SINGLE
    
    ! need to strongly enforce synchronization here
    !$OMP BARRIER
    !$OMP SINGLE
    scn_a = 0
    !$OMP SIMD REDUCTION(inscan, +:scn_a)
    do i = 1, N_det_l+1
        scn_a = scn_a + csr_s(i)
        !$OMP SCAN INCLUSIVE(scn_a)
        csr_s(i) = scn_a
    end do
    !$OMP END SINGLE
    !$OMP BARRIER

    ! loop through rows and construct CSR matrix from COO temp arrays
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det_l
        csr_c(csr_s(i):csr_s(i+1)-1) = coo_c_all(coo_s(i):coo_s(i+1)-1)
    end do
    !$OMP END DO

    deallocate(coo_r, coo_c, l_cols)

    !$OMP SINGLE
    deallocate(coo_s, coo_n, nnz_arr, coo_r_all, coo_c_all)
    !$OMP END SINGLE
    !$OMP END PARALLEL

end

subroutine get_sparsity_structure_with_triples(csr_s,csr_c,csr_e,sze,N_det_l)
    implicit none
    BEGIN_DOC
    ! Form a compressed sparse row matrix representation of the Hamiltonian
    ! connectivity in the space of the determinants. 

    ! csr_s are the row pointers
    ! csr_c are the column indices
    ! csr_e are the excitatio degrees between the row and column determinant
    ! sze is the maximum possible number of nonzero entries
    ! N_det_l is the largest index of the determinants to include (when sorted by energy) 
    END_DOC

    integer, intent(in)           :: N_det_l
    integer(kind=8), intent(in)        :: sze
    integer, intent(out)          :: csr_s(N_det_l+1), csr_c(sze), csr_e(sze)

    integer              :: n, i, j, k, l_row, old_row
    integer              :: nnz, nnz_cnt, nnz_tot
    integer              :: n_vals, n_vals_row, n_threads, ID
    integer              :: nnz_csr, ii, scn_a, kk
    integer :: OMP_get_num_threads, OMP_get_thread_num
    integer, allocatable :: nnz_arr(:), coo_r(:), coo_c(:), coo_e(:), l_cols(:), l_exc(:)
    integer, allocatable :: coo_r_all(:), coo_c_all(:), coo_e_all(:), coo_r_t(:), coo_c_t(:), coo_e_t(:)
    integer, allocatable:: coo_s(:), coo_n(:)
    double precision     :: frac
    

    !$OMP PARALLEL SHARED(nnz_tot, nnz_arr, nnz_csr, n_threads, N_det, N_det_l, nnz_max_per_row, n_vals_row,&
    !$OMP                 coo_r_all, coo_c_all, coo_e_all, csr_s, csr_c, csr_e, coo_s, coo_n)& 
    !$OMP PRIVATE(i,j,old_row, k,ii,kk, scn_a, ID, nnz, nnz_cnt, coo_r, coo_c, coo_e, coo_r_t, coo_c_t, coo_e_t, l_cols, l_exc, l_row) 

    !$ n_threads = OMP_get_num_threads()
    !$ ID = OMP_get_thread_num() + 1

    !$OMP SINGLE
    allocate(coo_s(N_det_l), coo_n(N_det_l))
    allocate(nnz_arr(n_threads))
    coo_n = 0
    nnz_tot = 0
    !$OMP END SINGLE
    !$OMP BARRIER

    ! initial allocation sizes for vectors
    frac = 0.2
    n_vals = max(nint(N_det_l*N_det_l*frac/n_threads), 128)
    n_vals_row = nnz_max_per_row

    allocate(coo_r(n_vals), coo_c(n_vals), coo_e(n_vals), l_cols(n_vals_row), l_exc(n_vals_row))
    
    nnz_cnt = 0
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det ! this loop needs to go over all the determinants, since this loop is not in determinant order but rather k_a order
        nnz = 0
        l_cols = 0

        ! check if row in N+1/N-1 determinant space
        ! if so, grab all the indices of all nonzero columns for that row in the upper half of H
        call get_all_sparse_columns_with_triples(i,l_cols,l_exc,l_row,nnz,n_vals_row,N_det_l)

        if (nnz == 0) cycle

        ! reallocate arrays if necessary
        if (nnz_cnt + nnz > size(coo_r, 1)) then
            allocate(coo_r_t(nnz_cnt + 10*nnz))
            allocate(coo_c_t(nnz_cnt + 10*nnz))
            allocate(coo_e_t(nnz_cnt + 10*nnz))
            
            coo_r_t(:nnz_cnt) = coo_r
            coo_c_t(:nnz_cnt) = coo_c
            coo_e_t(:nnz_cnt) = coo_e
            
            call move_alloc(coo_r_t, coo_r)
            call move_alloc(coo_c_t, coo_c)
            call move_alloc(coo_e_t, coo_e)
        end if
        
        
        ! Calculate nonzero entries and temperorarily store in COO format
        coo_n(l_row) = nnz ! store number of entries in row for later
        do j = 1, nnz
            nnz_cnt += 1
            coo_r(nnz_cnt) = l_row
            coo_c(nnz_cnt) = l_cols(j)
            coo_e(nnz_cnt) = l_exc(j)
        end do
    end do
    !$OMP END DO

    nnz_arr(ID) = nnz_cnt

    !$OMP CRITICAl
    nnz_tot = nnz_tot + nnz_cnt
    !$OMP END CRITICAl
    !$OMP BARRIER
    
    !$OMP SINGLE
    print *, "Total non-zero entries in full Hamiltonian: ", nnz_tot, " max size:", sze
    !$OMP END SINGLE
    
    !$OMP SINGLE
    allocate(coo_r_all(nnz_tot))
    allocate(coo_c_all(nnz_tot))
    allocate(coo_e_all(nnz_tot))
    !$OMP END SINGLE
    !$OMP BARRIER
    
    k = 0
    do i = 1, ID-1
        k += nnz_arr(i)
    end do

    ! coo_*_all are in shared memory; since threads have disjoint work, this is safe
    coo_r_all(k+1:k+nnz_arr(ID)) = coo_r
    coo_c_all(k+1:k+nnz_arr(ID)) = coo_c
    coo_c_all(k+1:k+nnz_arr(ID)) = coo_e
    !$OMP BARRIER

    ! calculate the starting index of each row in COO, since there is no sorting guarantee
    ! iterate over range, keep track of current row; if row changes, record the row start
    ii = (nnz_tot / n_threads)*(ID-1) + 1 ! start of range
    kk = (nnz_tot / n_threads) * ID + 3 ! end of range, slight overlap to catch boundary pointers
    if (ID == n_threads) kk = nnz_tot
    old_row = coo_r_all(ii)

    do i = ii, kk 
        if(coo_r_all(i) .ne. old_row) then 
            coo_s(coo_r_all(i)) = i
            old_row = coo_r_all(i)
        end if 
    end do

    ! make sure first row in COO is accounted for
    !$OMP SINGLE
    coo_s(coo_r_all(1)) = 1
    !$OMP END SINGLE

    ! calculate CSR pointer ranges
    ! row i data goes from csr_s(i) to csr_s(i+1) - 1
    ! use inclusive scan to reduce counts of rows in coo_n to set pointer ranges

    !$OMP SINGLE
    csr_s(1) = 1
    csr_s(2:N_det_l+1) = coo_n
    !$OMP END SINGLE
    
    ! need to strongly enforce synchronization here
    !$OMP BARRIER
    !$OMP SINGLE
    scn_a = 0
    !$OMP SIMD REDUCTION(inscan, +:scn_a)
    do i = 1, N_det_l+1
        scn_a = scn_a + csr_s(i)
        !$OMP SCAN INCLUSIVE(scn_a)
        csr_s(i) = scn_a
    end do
    !$OMP END SINGLE
    !$OMP BARRIER

    ! loop through rows and construct CSR matrix from COO temp arrays
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det_l
        csr_c(csr_s(i):csr_s(i+1)-1) = coo_c_all(coo_s(i):coo_s(i+1)-1)
        csr_e(csr_s(i):csr_s(i+1)-1) = coo_e_all(coo_s(i):coo_s(i+1)-1)
    end do
    !$OMP END DO

    deallocate(coo_r, coo_c, coo_e, l_cols, l_exc)

    !$OMP SINGLE
    deallocate(coo_s, coo_n, nnz_arr, coo_r_all, coo_c_all, coo_e_all)
    !$OMP END SINGLE
    !$OMP END PARALLEL

end



subroutine uH_structure_from_gH(H_p, H_c, H_e, uH_p, uH_c, N_det_g, N_det_u,nnz_g,nnz_max,N_exc,&
                      hash_alpha, hash_beta, hash_value, hash_prime, hash_table_size,&
                      I_cut_k, I_det_k)
    BEGIN_DOC
    ! Calculate the N+1/N-1 Hamiltonian from the ground state Hamiltonian with triples in CSR format
    ! N+1/N-1 hamiltonian constructed on the union of the excitation spaces generated by a set of 
    ! creation/annhiliation operators
    END_DOC

    implicit none
    
    integer, intent(in) :: N_det_g, N_det_u, nnz_g, nnz_max, N_exc
    integer, intent(in) :: H_p(N_det_g+1), H_c(nnz_g), H_e(nnz_g)
    
    integer, intent(in) :: hash_table_size, hash_prime
    integer, intent(in) :: hash_alpha(hash_table_size), hash_beta(hash_table_size), hash_value(hash_table_size)
    logical             :: hash_success
    
    integer, intent(in) :: I_cut_k(N_det_g, N_exc)
    integer(bit_kind), intent(in)   :: I_det_k(N_int, 2, N_det_u, N_exc)
    
    integer, intent(out):: uH_p(N_det_u+1), uH_c(nnz_max)

    
    integer :: i, j, row, col, target_row, target_col
    
end

subroutine get_all_sparse_columns(k_a, columns, row, nnz, nnz_max, N_det_l)
    BEGIN_DOC
    ! Find all nonzero column indices belonging to row k_a in the bilinear matrix ordering
    !
    ! k_a is current row index, in bilinear matrix order
    ! columns are the columns belong to to k_a
    ! row is the index in the energy ordering corresponding to k_a
    ! nnz is the number of entries in columns
    ! nnz_max is the maximum number of entries allowed in column
    ! N_det_l is the largest index of the determinants to include (when sorted by energy) 
    END_DOC

    implicit none

    integer, intent(in)      :: k_a, nnz_max, N_det_l
    integer, intent(out)     :: nnz, columns(nnz_max), row
    integer :: i, j, k, k_b
    integer :: krow, kcol, lrow, lcol
    integer :: lidx, kidx, tidx, n_buffer
    integer :: n_singles_a, n_singles_b, n_doubles_aa, n_doubles_bb, n_doubles_ab
    integer :: n_buffer_old, n_doubles_ab_tot
    integer :: n_off_diagonal, n_offset, tdegree

    integer(bit_kind) :: ref_det(N_int, 2), tmp_det(N_int, 2), tmp_det2(N_int, 2), sdet_a(N_int), sdet_b(N_int)
    integer(bit_kind), allocatable    :: buffer(:,:)
    integer, allocatable              :: singles_a(:), singles_b(:)
    integer, allocatable              :: doubles_aa(:), doubles_bb(:), doubles_ab(:)
    integer, allocatable              :: idx(:), all_idx(:), srt_idx(:)
    

    ! initialize arrays
    allocate(buffer(N_int, N_det_l), idx(N_det_l))

    n_singles_a = 0
    n_singles_b = 0
    n_doubles_aa = 0
    n_doubles_ab = 0
    n_doubles_bb = 0
    
    ! store index and make sure it isn't out of range of the reduced determinant space
    kidx = psi_bilinear_matrix_order(k_a)

    if (kidx > N_det_l) return ! determinant not included in this subset

    ! get pointer indices and a ref determinant
    krow = psi_bilinear_matrix_rows(k_a)
    kcol = psi_bilinear_matrix_columns(k_a)
    ref_det(:,1) = psi_det_alpha_unique(:, krow)
    ref_det(:,2) = psi_det_beta_unique (:, kcol)

    ! Finding (1,0) and (2,0) excitations
    ! loop over same beta different alpha
    n_buffer = 0
    do i = psi_bilinear_matrix_columns_loc(kcol), psi_bilinear_matrix_columns_loc(kcol+1)-1
        lidx = psi_bilinear_matrix_order(i)
        
        ! check if determinant is in upper half of reduced Hamiltonian matrix
        if ((lidx <= kidx) .or. (lidx > N_det_l)) cycle

        lcol = psi_bilinear_matrix_columns(i)
        lrow = psi_bilinear_matrix_rows(i)
        
        tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
        tmp_det(:,2) = psi_det_beta_unique (:, lcol)

        ! add determinant to buffer
        ! buffer is list of alpha spin determinants
        n_buffer += 1
        buffer(:,n_buffer) = tmp_det(:,1)
        idx(n_buffer) = lidx
    end do

    allocate(singles_a(n_buffer), doubles_aa(n_buffer))

    ! grab indices of all determinants in buffer related to ref_det by (1,0) or (2,0) excitations 
    sdet_a = ref_det(:,1)
    call get_all_spin_singles_and_doubles(buffer, idx, sdet_a, &
                        N_int, n_buffer, singles_a, doubles_aa, n_singles_a, n_doubles_aa)

    deallocate(buffer, idx)
    allocate(buffer(N_int, N_det), idx(N_det))

    ! Finding (0,1) and (0,2) excitations 
    k_b = psi_bilinear_matrix_order_transp_reverse(k_a)
    krow = psi_bilinear_matrix_transp_rows(k_b) !this is unnecessary, technically
    kcol = psi_bilinear_matrix_transp_columns(k_b)
    
    ! loop over same alpha different beta
    n_buffer = 0
    do i = psi_bilinear_matrix_transp_rows_loc(krow), psi_bilinear_matrix_transp_rows_loc(krow+1)-1
        tidx = psi_bilinear_matrix_transp_order(i)
        lidx = psi_bilinear_matrix_order(tidx)

        ! check if determinant is in upper half of reduced Hamiltonian matrix
        if ((lidx <= kidx) .or. (lidx > N_det_l)) cycle
        
        lcol = psi_bilinear_matrix_transp_columns(i)
        lrow = psi_bilinear_matrix_transp_rows(i)
        
        tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
        tmp_det(:,2) = psi_det_beta_unique (:, lcol)

        ! add determinant to buffer
        ! buffer is list of beta spin determinants
        n_buffer += 1
        buffer(:,n_buffer) = tmp_det(:,2)
        idx(n_buffer) = lidx
    end do
    
    allocate(singles_b(n_buffer), doubles_bb(n_buffer))

    ! grab indices of all determinants in buffer related to ref_det by (0,1) or (0,2) excitations 
    sdet_b = ref_det(:,2)
    call get_all_spin_singles_and_doubles(buffer, idx, sdet_b, &
                        N_int, n_buffer, singles_b, doubles_bb, n_singles_b, n_doubles_bb)
                        
    deallocate(buffer, idx)

    ! Finding (1,1) excitations
    allocate(buffer(N_int, N_det), idx(N_det))
    allocate(doubles_ab(N_det))

    n_buffer = 0
    n_buffer_old = 0
    n_doubles_ab_tot = 0

    ! starting from list of beta singles does not give all the (1,1) excitations
    ! so we need to search over either all beta or all alpha at some point
    ! TODO: perhaps better to change outer loop to alpha, since we always excite in alpha channel? would it be significantly different?
    ! start from (0,1) to excite to (1,1)
    do j = 1, n_det_beta_unique

        kcol = psi_bilinear_matrix_columns(k)

        tmp_det2(:,2) = psi_det_beta_unique (:, j)

        ! check if a single excitation different
        tdegree = 0
        do i = 1, N_int
            tdegree += popcnt(ieor(tmp_det2(i,2), ref_det(i,2)))
        end do

        if (tdegree /= 2) then
            cycle
        end if

        ! loop over same beta different alpha
        do i = psi_bilinear_matrix_columns_loc(j), psi_bilinear_matrix_columns_loc(j+1)-1
            lidx = psi_bilinear_matrix_order(i)

            ! check if determinant is in upper half of reduced Hamiltonian matrix
            if ((lidx <= kidx) .or. (lidx > N_det_l)) cycle

            lcol = psi_bilinear_matrix_columns(i)
            lrow = psi_bilinear_matrix_rows(i)

            tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
            tmp_det(:,2) = psi_det_beta_unique (:, lcol)

            ! add determinant to buffer
            ! buffer is list of alpha spin determinants
            n_buffer += 1
            buffer(:,n_buffer) = tmp_det(:,1)
            idx(n_buffer) = lidx
        end do

        sdet_a = ref_det(:,1)

        ! all determinants are (X,1) excitations from ref_det
        ! so we just need to check alpha channel now
        ! grab indices of all determinants in buffer related to ref_det by (1,1) excitations 
        call get_all_spin_singles(buffer(:,n_buffer_old+1:n_buffer), idx(n_buffer_old+1:n_buffer),&
                                sdet_a, N_int, n_buffer-n_buffer_old,&
                                doubles_ab(n_doubles_ab_tot+1:n_doubles_ab_tot+n_buffer-n_buffer_old),&
                                n_doubles_ab)


        n_buffer_old = n_buffer
        n_doubles_ab_tot += n_doubles_ab
    end do

    ! add all indices into list and sort
    ! number of off-diagonal terms needed to caclulate to fill this row in sparse Hamlitonian
    n_off_diagonal = n_singles_a + n_singles_b + n_doubles_aa + n_doubles_ab_tot + n_doubles_bb

    allocate(all_idx(n_off_diagonal+1), srt_idx(n_off_diagonal+1))

    do i = 1, n_off_diagonal+1
        srt_idx(i) = i
    end do

    n_offset = 0
    all_idx(n_offset+1:n_offset+1)              = kidx
                            n_offset           += 1
    all_idx(n_offset+1:n_offset+n_singles_a)    = singles_a(:n_singles_a)
                            n_offset           += n_singles_a
    all_idx(n_offset+1:n_offset+n_singles_b)    = singles_b(:n_singles_b)
                            n_offset           += n_singles_b
    all_idx(n_offset+1:n_offset+n_doubles_aa)   = doubles_aa(:n_doubles_aa)
                            n_offset           += n_doubles_aa
    all_idx(n_offset+1:n_offset+n_doubles_ab_tot) = doubles_ab(:n_doubles_ab_tot)
                            n_offset           += n_doubles_ab_tot
    all_idx(n_offset+1:n_offset+n_doubles_bb)   = doubles_bb(:n_doubles_bb)
                            n_offset           += n_doubles_bb


    call insertion_isort(all_idx, srt_idx, n_off_diagonal+1)

    columns(:n_off_diagonal+1) = all_idx
    nnz = n_off_diagonal + 1
    row = kidx
    
    deallocate(buffer, idx, singles_a, doubles_aa,&
                singles_b, doubles_bb, doubles_ab,&
                srt_idx, all_idx)
end

subroutine get_all_sparse_columns_with_triples( k_a, columns,exc_degree,row,nnz,nnz_max,N_det_l )
    BEGIN_DOC
    ! Find all nonzero column indices belonging to row k_a in the bilinear matrix ordering
    !
    ! k_a is current row index, in bilinear matrix order
    ! columns are the columns belong to to k_a
    ! row is the index in the energy ordering corresponding to k_a
    ! nnz is the number of entries in columns
    ! nnz_max is the maximum number of entries allowed in column
    ! N_det_l is the largest index of the determinants to include (when sorted by energy) 
    END_DOC

    implicit none

    integer, intent(in)      :: k_a, nnz_max, N_det_l
    integer, intent(out)     :: nnz, columns(nnz_max), row, exc_degree(nnz_max)
    integer :: i, j, k, k_b
    integer :: krow, kcol, lrow, lcol
    integer :: lidx, kidx, tidx, n_buffer
    integer :: n_singles_a, n_singles_b, n_doubles_aa, n_doubles_bb, n_doubles_ab
    integer :: n_buffer_old, n_doubles_ab_tot, n_triples_tot
    integer :: n_off_diagonal, n_offset, tdegree_alpha, tdegree_beta

    integer(bit_kind) :: ref_det(N_int, 2), tmp_det(N_int, 2), tmp_det2(N_int, 2), sdet_a(N_int), sdet_b(N_int)
    integer(bit_kind), allocatable    :: buffer(:,:)
    integer, allocatable              :: singles_a(:), singles_b(:)
    integer, allocatable              :: doubles_aa(:), doubles_bb(:), doubles_ab(:), triples(:), triples_case(:)
    integer, allocatable              :: idx(:), all_idx(:), srt_idx(:), all_degree(:)
    

    ! initialize arrays
    allocate(buffer(N_int, N_det_l), idx(N_det_l))

    n_singles_a = 0
    n_singles_b = 0
    n_doubles_aa = 0
    n_doubles_ab = 0
    n_doubles_bb = 0
    
    ! store index and make sure it isn't out of range of the reduced determinant space
    kidx = psi_bilinear_matrix_order(k_a)

    if (kidx > N_det_l) return ! determinant not included in this subset

    ! get pointer indices and a ref determinant
    krow = psi_bilinear_matrix_rows(k_a)
    kcol = psi_bilinear_matrix_columns(k_a)
    ref_det(:,1) = psi_det_alpha_unique(:, krow)
    ref_det(:,2) = psi_det_beta_unique (:, kcol)

    ! Finding (1,0) and (2,0) excitations
    ! loop over same beta different alpha
    n_buffer = 0
    do i = psi_bilinear_matrix_columns_loc(kcol), psi_bilinear_matrix_columns_loc(kcol+1)-1
        lidx = psi_bilinear_matrix_order(i)
        
        ! check if determinant is in upper half of reduced Hamiltonian matrix
        if ((lidx <= kidx) .or. (lidx > N_det_l)) cycle

        lcol = psi_bilinear_matrix_columns(i)
        lrow = psi_bilinear_matrix_rows(i)
        
        tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
        tmp_det(:,2) = psi_det_beta_unique (:, lcol)

        ! add determinant to buffer
        ! buffer is list of alpha spin determinants
        n_buffer += 1
        buffer(:,n_buffer) = tmp_det(:,1)
        idx(n_buffer) = lidx
    end do

    allocate(singles_a(n_buffer), doubles_aa(n_buffer))

    ! grab indices of all determinants in buffer related to ref_det by (1,0) or (2,0) excitations 
    sdet_a = ref_det(:,1)
    call get_all_spin_singles_and_doubles(buffer, idx, sdet_a, &
                        N_int, n_buffer, singles_a, doubles_aa, n_singles_a, n_doubles_aa)

    deallocate(buffer, idx)
    allocate(buffer(N_int, N_det_l), idx(N_det_l))

    ! Finding (0,1) and (0,2) excitations 
    k_b = psi_bilinear_matrix_order_transp_reverse(k_a)
    krow = psi_bilinear_matrix_transp_rows(k_b) !this is unnecessary, technically
    kcol = psi_bilinear_matrix_transp_columns(k_b)
    
    ! loop over same alpha different beta
    n_buffer = 0
    do i = psi_bilinear_matrix_transp_rows_loc(krow), psi_bilinear_matrix_transp_rows_loc(krow+1)-1
        tidx = psi_bilinear_matrix_transp_order(i)
        lidx = psi_bilinear_matrix_order(tidx)

        ! check if determinant is in upper half of reduced Hamiltonian matrix
        if ((lidx <= kidx) .or. (lidx > N_det_l)) cycle
        
        lcol = psi_bilinear_matrix_transp_columns(i)
        lrow = psi_bilinear_matrix_transp_rows(i)
        
        tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
        tmp_det(:,2) = psi_det_beta_unique (:, lcol)

        ! add determinant to buffer
        ! buffer is list of beta spin determinants
        n_buffer += 1
        buffer(:,n_buffer) = tmp_det(:,2)
        idx(n_buffer) = lidx
    end do
    
    allocate(singles_b(n_buffer), doubles_bb(n_buffer))

    ! grab indices of all determinants in buffer related to ref_det by (0,1) or (0,2) excitations 
    sdet_b = ref_det(:,2)
    call get_all_spin_singles_and_doubles(buffer, idx, sdet_b, &
                        N_int, n_buffer, singles_b, doubles_bb, n_singles_b, n_doubles_bb)
                        
    deallocate(buffer, idx)

    ! Finding (1,1) excitations
    allocate(buffer(N_int, N_det), idx(N_det))
    allocate(doubles_ab(N_det))

    n_buffer = 0
    n_buffer_old = 0
    n_doubles_ab_tot = 0

    ! starting from list of beta singles does not give all the (1,1) excitations
    ! so we need to search over either all beta or all alpha at some point
    ! start from (0,1) to excite to (1,1)
    do j = 1, n_det_beta_unique

        kcol = psi_bilinear_matrix_columns(k)

        tmp_det2(:,2) = psi_det_beta_unique (:, j)

        ! check if a single excitation different
        tdegree_beta = 0
        do i = 1, N_int
            tdegree_beta += popcnt(ieor(tmp_det2(i,2), ref_det(i,2)))
        end do

        if (tdegree_beta /= 2) then
            cycle
        end if

        ! loop over same beta different alpha
        do i = psi_bilinear_matrix_columns_loc(j), psi_bilinear_matrix_columns_loc(j+1)-1
            lidx = psi_bilinear_matrix_order(i)

            ! check if determinant is in upper half of reduced Hamiltonian matrix
            if ((lidx <= kidx) .or. (lidx > N_det_l)) cycle

            lcol = psi_bilinear_matrix_columns(i)
            lrow = psi_bilinear_matrix_rows(i)

            tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
            tmp_det(:,2) = psi_det_beta_unique (:, lcol)

            ! add determinant to buffer
            ! buffer is list of alpha spin determinants
            n_buffer += 1
            buffer(:,n_buffer) = tmp_det(:,1)
            idx(n_buffer) = lidx
        end do

        sdet_a = ref_det(:,1)

        ! all determinants are (X,1) excitations from ref_det
        ! so we just need to check alpha channel now
        ! grab indices of all determinants in buffer related to ref_det by (1,1) excitations 
        call get_all_spin_singles(buffer(:,n_buffer_old+1:n_buffer), idx(n_buffer_old+1:n_buffer),&
                                sdet_a, N_int, n_buffer-n_buffer_old,&
                                doubles_ab(n_doubles_ab_tot+1:n_doubles_ab_tot+n_buffer-n_buffer_old),&
                                n_doubles_ab)


        n_buffer_old = n_buffer
        n_doubles_ab_tot += n_doubles_ab
    end do


    !!! Finding triple excitations of form (3,0), (2,1), and (1,2)
    ! start from (0,0), (0,1), (0,2) to excite to (3,0), (2,1), (1,2)
    n_triples_tot = 0
    allocate(triples(N_det_l), triples_case(N_det_l))

    do j = 1, n_det_beta_unique

        kcol = psi_bilinear_matrix_columns(k)

        tmp_det2(:,2) = psi_det_beta_unique (:, j)

        ! check if a single excitation different
        tdegree_beta = 0
        do i = 1, N_int
            tdegree_beta += popcnt(ieor(tmp_det2(i,2), ref_det(i,2)))
        end do

        ! exclude anything with beta excitation > 2
        if (tdegree_beta > 4) then
            cycle
        end if

        ! loop over same beta different alpha
        do i = psi_bilinear_matrix_columns_loc(j), psi_bilinear_matrix_columns_loc(j+1)-1
            lidx = psi_bilinear_matrix_order(i)

            ! check if determinant is in upper half of reduced Hamiltonian matrix
            if ((lidx <= kidx) .or. (lidx > N_det_l)) cycle

            lcol = psi_bilinear_matrix_columns(i)
            lrow = psi_bilinear_matrix_rows(i)

            tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
            tmp_det(:,2) = psi_det_beta_unique (:, lcol)

            ! tmp_det(:,2) and tmp_det2(:,2) should be the same?
            tdegree_alpha = 0
            do i = 1, N_int
                tdegree_alpha += popcnt(ieor(tmp_det(i,1), ref_det(i,1)))
            end do

            ! add to triples buffer 
            if (tdegree_alpha + tdegree_beta == 6) then
                n_triples_tot += 1
                triples(n_triples_tot) = lidx
                select case (tdegree_beta)
                    case (0)
                        triples_case(n_triples_tot) = 3
                    case (2)
                        triples_case(n_triples_tot) = 8
                    case (4)
                        triples_case(n_triples_tot) = 7
                    case default
                        write(*, ("(A40, I16, I16, I16)")),"Error in triple selection on krow/beta det/alpha det: " k_a, j, i
                        exit
                end select
            end if

        end do
    end do


    !!! Create final buffers
    ! add all indices into list and sort
    ! number of off-diagonal terms needed to caclulate to fill this row in sparse Hamlitonian
    n_off_diagonal = n_singles_a + n_singles_b + n_doubles_aa + n_doubles_ab_tot + n_doubles_bb + n_triples_tot

    allocate(all_degree(n_off_diagonal+1), all_idx(n_off_diagonal+1), srt_idx(n_off_diagonal+1))

    do i = 1, n_off_diagonal+1
        srt_idx(i) = i
    end do
  
    ! a - 1; aa - 2; aaa -3; b - 4; bb - 5; ab - 6; abb - 7; aab - 8
    n_offset = 0
    all_idx(n_offset+1:n_offset+1)               = kidx
    all_degree(n_offset+1:n_offset+1)            = 0
                            n_offset            += 1

    all_idx(n_offset+1:n_offset+n_singles_a)     = singles_a(:n_singles_a)
    all_degree(n_offset+1:n_offset+n_singles_a)  = 1
                            n_offset            += n_singles_a

    all_idx(n_offset+1:n_offset+n_singles_b)     = singles_b(:n_singles_b)
    all_degree(n_offset+1:n_offset+n_singles_b)  = 3
                            n_offset            += n_singles_b

    all_idx(n_offset+1:n_offset+n_doubles_aa)    = doubles_aa(:n_doubles_aa)
    all_degree(n_offset+1:n_offset+n_doubles_aa) = 2
                            n_offset            += n_doubles_aa

    all_idx(n_offset+1:n_offset+n_doubles_ab_tot)    = doubles_ab(:n_doubles_ab_tot)
    all_degree(n_offset+1:n_offset+n_doubles_ab_tot) = 6
                            n_offset                += n_doubles_ab_tot

    all_idx(n_offset+1:n_offset+n_doubles_bb)     = doubles_bb(:n_doubles_bb)
    all_degree(n_offset+1:n_offset+n_doubles_bb)  = 5
                            n_offset             += n_doubles_bb

    all_idx(n_offset+1:n_offset+n_triples_tot)    = triples(:n_triples_tot)
    all_degree(n_offset+1:n_offset+n_triples_tot) = triples_case(:n_triples_tot)

    call insertion_isort(all_idx, srt_idx, n_off_diagonal+1)

    columns(:n_off_diagonal+1) = all_idx

    do i = 1, n_off_diagonal+1
        exc_degree(i) = all_degree(srt_idx(i))
    end do

    nnz = n_off_diagonal + 1
    row = kidx
    
    deallocate(buffer, idx, singles_a, doubles_aa,&
                singles_b, doubles_bb, doubles_ab,&
                triples, triples_case,&
                srt_idx, all_idx)
end

subroutine calc_sparse_dH(H_p, H_c, H_v, sze, nnz, dets)
    BEGIN_DOC
    ! Calculate the entries to the sparse Hamiltonian, given the structure
    ! in CSR format
    ! H_p are the row pointers
    ! H_c are the column indices
    ! H_v are the values of the Hamiltonian
    ! sze is the square side length of the Hamiltonian
    ! nnz is the total number of non zeros
    ! dets is the set of (excited) determinants used to calculate the matrix entries
    END_DOC

    integer :: i, j
    integer, intent(in) :: sze, nnz, H_p(sze+1), H_c(nnz)
    integer(bit_kind), intent(in) :: dets(N_int, 2, sze)
    double precision, intent(out) :: H_v(nnz)
    double precision :: hij

    H_v = 0.d0

    ! force provide early so that threads don't each try to provide
    call i_H_j(dets(:,:,1), dets(:,:,1), N_int, hij) 

    !$OMP PARALLEL PRIVATE(i, j, hij) SHARED(H_c, H_p, H_v, sze, dets)
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, sze                
        ! loop over columns
        do j = H_p(i), H_p(i+1)-1
            call i_H_j(dets(:,:,i),&
                       dets(:,:,H_c(j)), N_int, hij)
            H_v(j) = hij
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
end

subroutine calc_sparse_zH(H_p, H_c, H_v, sze, nnz, dets)
    BEGIN_DOC
    ! Calculate the entries to the sparse Hamiltonian, given the structure
    ! in CSR format
    ! H_p are the row pointers
    ! H_c are the column indices
    ! H_v are the values of the Hamiltonian
    ! sze is the square side length of the Hamiltonian
    ! nnz is the total number of non zeros
    ! dets is the set of (excited) determinants used to calculate the matrix entries
    END_DOC

    integer :: i, j
    integer, intent(in) :: sze, nnz, H_p(sze+1), H_c(nnz)
    integer(bit_kind), intent(in) :: dets(N_int, 2, sze)
    complex*16, intent(out) :: H_v(nnz)
    complex*16 :: hij

    H_v = 0.d0

    ! force provide early so that threads don't each try to provide
    call i_h_j_complex(dets(:,:,1), dets(:,:,1), N_int, hij) 

    !$OMP PARALLEL PRIVATE(i, j, hij) SHARED(H_c, H_p, H_v, sze, dets)
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, sze                
        ! loop over columns
        do j = H_p(i), H_p(i+1)-1
            call i_h_j_complex(dets(:,:,i),&
                       dets(:,:,H_c(j)), N_int, hij)
            H_v(j) = hij
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
end