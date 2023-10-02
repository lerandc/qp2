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

    integer,          intent(in)  :: sze
    integer(kind=8),  intent(in)  :: nnz, A_p(sze+1)
    integer,          intent(in)  :: A_c(nnz)
    double precision, intent(in)  :: A_v(nnz), x(sze)
    double precision, intent(out) :: y(sze)
    double precision, allocatable :: y_t(:)
    integer(kind=8)               :: i, j

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

subroutine sparse_csr_dmv_row_part(A_v, A_c, A_p, x, y, sze, nnz, row_starts, n_threads)
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

    integer,          intent(in)  :: sze, n_threads, row_starts(n_threads)
    integer(kind=8),  intent(in)  :: nnz, A_p(sze+1)
    integer,          intent(in)  :: A_c(nnz)
    double precision, intent(in)  :: A_v(nnz), x(sze)
    double precision, intent(out) :: y(sze)
    double precision, allocatable :: y_t(:)
    integer(kind=8)               :: i, j 
    integer :: OMP_get_thread_num, ID


    !$OMP PARALLEL PRIVATE(i, j, y_t, ID) SHARED(y, x, A_c, A_p, A_v, sze, row_starts)
    allocate(y_t(sze))
    y_t = 0.d0
    !$ ID = OMP_get_thread_num() + 1

    ! loop over rows
    do i = row_starts(ID), row_starts(ID+1) - 1
        
        
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
    integer, allocatable:: coo_s(:)
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

    integer,         intent(in) :: sze 
    integer(kind=8), intent(in) :: nnz_in, A_p(sze+1)
    integer, intent(in) :: A_c(nnz_in), I_k(sze)
    integer, intent(out):: B_c(nnz_in)
    integer(kind=8), intent(out) :: B_p(sze+1)
    integer :: lcol
    integer(kind=8) :: i, j, k, ii, kk, nnz_cnt, old_row, nnz_tot, scn_a
    integer              :: n_vals, n_threads, ID, nnz
    integer :: OMP_get_num_threads, OMP_get_thread_num
    integer, allocatable :: nnz_arr(:), coo_r(:), coo_c(:)
    integer, allocatable :: coo_r_all(:), coo_c_all(:), coo_r_t(:), coo_c_t(:)
    integer(kind=8), allocatable :: coo_s(:)
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
    !$OMP BARRIER

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
    !$OMP BARRIER
    
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
        B_c(B_p(i):B_p(i+1)-1) = coo_c_all(coo_s(i):coo_s(i) + (B_p(i+1) - B_p(i) -1) ) 
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
    integer(kind=8), intent(in)   :: sze
    integer, intent(out)          :: csr_c(sze)
    integer(kind=8), intent(out)   :: csr_s(N_det_l+1)

    integer              :: l_row
    integer              :: nnz, nnz_cnt
    integer              :: n_vals, n_vals_row, n_threads, ID
    integer(kind=8)      :: nnz_tot, i, j, k, ii, kk, scn_a, old_row
    integer :: OMP_get_num_threads, OMP_get_thread_num
    integer, allocatable :: nnz_arr(:), coo_r(:), coo_c(:), l_cols(:), coo_n(:)
    integer, allocatable :: coo_r_all(:), coo_c_all(:), coo_r_t(:), coo_c_t(:)
    integer(kind=8), allocatable:: coo_s(:) 
    double precision     :: frac
    

    !$OMP PARALLEL SHARED(nnz_tot, nnz_arr, n_threads, N_det, N_det_l, nnz_max_per_row, n_vals_row,&
    !$OMP                 coo_r_all, coo_c_all, csr_s, csr_c, coo_s, coo_n)& 
    !$OMP PRIVATE(i, j , old_row, k, ii, kk, scn_a, ID, nnz, nnz_cnt, coo_r, coo_c, coo_r_t, coo_c_t, l_cols, l_row) 

    !$ ID = OMP_get_thread_num() + 1
    
    !$OMP SINGLE
    !$ n_threads = OMP_get_num_threads()
    allocate(coo_s(N_det_l), coo_n(N_det_l))
    allocate(nnz_arr(n_threads))
    nnz_tot = 0
    coo_n = 0
    !$OMP END SINGLE
    !$OMP BARRIER

    ! initial allocation sizes for vectors
    frac = 0.2
    n_vals = max(nint(N_det_l*N_det_l*frac/n_threads), 128)
    n_vals_row = min(N_det_l, nnz_max_per_row)

    allocate(coo_r(n_vals))
    allocate(coo_c(n_vals))
    allocate(l_cols(n_vals_row))
        
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
        if ((nnz_cnt + nnz) > size(coo_r, 1)) then
            allocate(coo_r_t(nnz_cnt + 10*nnz))
            allocate(coo_c_t(nnz_cnt + 10*nnz))
            
            coo_r_t(:nnz_cnt) = coo_r(:nnz_cnt)
            coo_c_t(:nnz_cnt) = coo_c(:nnz_cnt)
            
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
    coo_r_all(k+1:k+nnz_arr(ID)) = coo_r(:nnz_cnt)
    coo_c_all(k+1:k+nnz_arr(ID)) = coo_c(:nnz_cnt)
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
        csr_c(csr_s(i):csr_s(i+1)-1) = coo_c_all(coo_s(i):coo_s(i) + (csr_s(i+1) - csr_s(i) - 1))
    end do
    !$OMP END DO

    deallocate(coo_r, coo_c, l_cols)

    !$OMP SINGLE
    deallocate(coo_s, coo_n, nnz_arr, coo_r_all, coo_c_all)
    !$OMP END SINGLE
    !$OMP END PARALLEL

end

subroutine write_degree(col, exc_degree)
    implicit none
    integer, intent(in) :: col, exc_degree
    integer degree_a, degree_b
    select case (exc_degree)
    case (0)
        degree_a = 0
        degree_b = 0
    case (1)
        degree_a = 1
        degree_b = 0
    case (2)
        degree_a = 2
        degree_b = 0
    case (3)
        degree_a = 3
        degree_b = 0
    case (4)
        degree_a = 0
        degree_b = 1
    case (5)
        degree_a = 0
        degree_b = 2
    case (6)
        degree_a = 1
        degree_b = 1
    case (7)
        degree_a = 1
        degree_b = 2
    case (8)
        degree_a = 2
        degree_b = 1
    end select
    write(*, "(I8, I8, I8)"), col, degree_a, degree_b
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
    integer(kind=8), intent(in)   :: sze
    integer, intent(out)          :: csr_c(sze), csr_e(sze)
    integer(kind=8), intent(out)   :: csr_s(N_det_l+1)

    integer              :: l_row
    integer              :: nnz, nnz_cnt
    integer              :: n_vals, n_vals_row, n_threads, ID
    integer              :: stat
    integer(kind=8)      :: nnz_tot, i, j, k, ii, kk, scn_a, old_row
    integer :: OMP_get_num_threads, OMP_get_thread_num
    integer, allocatable :: nnz_arr(:), coo_r(:), coo_c(:), coo_e(:), l_cols(:), l_exc(:), coo_n(:)
    integer, allocatable :: coo_r_all(:), coo_c_all(:), coo_e_all(:), coo_r_t(:), coo_c_t(:), coo_e_t(:)
    integer(kind=8), allocatable :: coo_s(:)
    double precision     :: frac
    
    !$OMP PARALLEL SHARED(nnz_tot, nnz_arr, n_threads, N_det, N_det_l, n_vals_row,&
    !$OMP                 coo_r_all, coo_c_all, coo_e_all, csr_s, csr_c, csr_e, coo_s, coo_n, sze)& 
    !$OMP PRIVATE(stat, i,j,k,ii,kk,old_row, scn_a, ID, nnz, nnz_cnt, coo_r, coo_c, coo_e, coo_r_t, coo_c_t, coo_e_t, l_cols, l_exc, l_row, frac) 

    !$ ID = OMP_get_thread_num() + 1
    
    !$OMP SINGLE
    !$ n_threads = OMP_get_num_threads()
    allocate(nnz_arr(n_threads))
    allocate(coo_s(N_det_l))
    allocate(coo_n(N_det_l))
    coo_n = 0
    nnz_tot = 0
    !$OMP END SINGLE
    
    ! ! initial allocation sizes for vectors
    frac = 0.2
    n_vals = max(nint(N_det_l*N_det_l*frac/n_threads), 128)
    n_vals_row = N_det_l !min(N_det_l, nnz_max_per_row)

    allocate(coo_r(n_vals), coo_c(n_vals), coo_e(n_vals), l_cols(n_vals_row), l_exc(n_vals_row))
    nnz_cnt = 0
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det ! this loop needs to go over all the determinants, since this loop is not in determinant order but rather k_a order
        nnz = 0
        l_cols = 0
        l_exc = 0

        ! check if row in N+1/N-1 determinant space
        ! if so, grab all the indices of all nonzero columns for that row in the upper half of H
        call get_all_sparse_columns_with_triples(i, l_cols, l_exc, l_row, nnz, n_vals_row, N_det_l)

        if (nnz == 0) cycle

        ! reallocate arrays if necessary
        if ((nnz_cnt + nnz) > size(coo_r, 1)) then
            allocate(coo_r_t(nnz_cnt + 10*nnz))
            allocate(coo_c_t(nnz_cnt + 10*nnz))
            allocate(coo_e_t(nnz_cnt + 10*nnz))
            
            coo_r_t(:nnz_cnt) = coo_r(:nnz_cnt)
            coo_c_t(:nnz_cnt) = coo_c(:nnz_cnt)
            coo_e_t(:nnz_cnt) = coo_e(:nnz_cnt)
            
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
    coo_r_all(k+1:k+nnz_arr(ID)) = coo_r(:nnz_cnt)
    coo_c_all(k+1:k+nnz_arr(ID)) = coo_c(:nnz_cnt)
    coo_e_all(k+1:k+nnz_arr(ID)) = coo_e(:nnz_cnt)
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
        csr_c(csr_s(i):csr_s(i+1)-1) = coo_c_all(coo_s(i):coo_s(i) + (csr_s(i+1) - csr_s(i) - 1))
        csr_e(csr_s(i):csr_s(i+1)-1) = coo_e_all(coo_s(i):coo_s(i) + (csr_s(i+1) - csr_s(i) - 1))
    end do
    !$OMP END DO

    deallocate(coo_r, coo_c, coo_e, l_cols, l_exc)

    !$OMP SINGLE
    deallocate(coo_s, coo_n, nnz_arr, coo_r_all, coo_c_all, coo_e_all)
    !$OMP END SINGLE
    !$OMP END PARALLEL

end

subroutine uH_structure_from_gH(H_p, H_c, H_e, uH_p, uH_c, N_det_g, N_det_u,nnz_g,nnz_max,N_exc,&
                                hash_alpha, hash_beta, hash_vals, hash_prime, hash_table_size,&
                                I_cut_k, I_det_k, I_det_ind)
    BEGIN_DOC
    ! Calculate the N+1/N-1 Hamiltonian from the ground state Hamiltonian with triples in CSR format
    ! N+1/N-1 hamiltonian constructed on the union of the excitation spaces generated by a set of 
    ! creation/annhiliation operators
    END_DOC

    implicit none
    
    !! input data
    integer(kind=8), intent(in) :: nnz_g, nnz_max
    integer, intent(in) :: N_det_g, N_det_u, N_exc
    integer, intent(in) :: H_c(nnz_g), H_e(nnz_g)
    integer(kind=8), intent(in) :: H_p(N_det_g+1)
    
    integer, intent(in) :: hash_table_size, hash_prime
    integer, intent(in) :: hash_alpha(hash_table_size), hash_beta(hash_table_size), hash_vals(hash_table_size)
    logical             :: hash_success
    
    integer, intent(in) :: I_cut_k(N_det_g, N_exc), I_det_ind(N_det_g, N_exc)
    integer(bit_kind), intent(in)   :: I_det_k(N_int, 2, N_det_g, N_exc)
    
    !! output data
    integer, intent(out) :: uH_c(nnz_max)
    integer(kind=8), intent(out) :: uH_p(N_det_u+1)

    !! hash functions
    integer :: hash_value

    !! OMP calls
    integer :: OMP_get_num_threads, OMP_get_thread_num

    !! routine variables
    integer(kind=8)     :: i, ii, j, jj, k, kk, scn_a, pidx, block_start, block_end, umax_block_size
    integer             :: i_exc, j_exc, row, col, target_row, target_col, cur_exc_degree
    integer*1           :: bit_degree, criterion_a_mask, criterion_b_mask, criterion_c_mask, bit_map(8)
    integer             :: target_N, n_threads, ID
    integer(bit_kind)   :: target_row_det(N_int, 2), target_col_det(N_int, 2)
    logical             :: block_a, block_b, criterion_a, criterion_b, criterion_c, pass_a, pass_b, pass_c, add_pair

    integer              :: buffer_count, nts_buffer_count, max_new_entries, buffer_total, ubuffer_total
    integer              :: nts_block_start, nts_block_end, block_total, max_block_size

    integer, allocatable :: coo_r(:), coo_c(:), coo_n(:), coo_r_nts(:), coo_c_nts(:), coo_n_nts(:)
    integer, allocatable :: u_coo_r(:), u_coo_c(:), u_coo_n(:)
    integer, allocatable :: coo_c_all(:), t_coo_r(:), t_coo_c(:), sort_idx(:)
    integer(kind=8), allocatable :: coo_n_all(:), nnz_arr(:,:), u_coo_n_all(:),  u_block_rows(:,:)
    integer, allocatable :: row_starts(:), pointer_blocks(:)
    integer              :: row_start, row_end, old_row

    !! profiling
    double precision :: t0, t1, t_tot, t_all, t_main_loop

    criterion_a_mask = 144
    criterion_b_mask = 231
    criterion_c_mask = 220
    bit_map = (/128, 64, 32, 16, 8, 4, 2, 1/)

    ! coo_x stores entries blocked by row; coo_x_nts stores entries in unsorted order
    !$OMP PARALLEL SHARED(N_det_g, N_det_u, nnz_g, nnz_max, N_exc, H_p, H_c, H_e,&
    !$OMP                 uH_p, uH_c, pointer_blocks, row_starts, target_N, n_threads,&
    !$OMP                 hash_table_size, hash_prime, hash_alpha, hash_beta, hash_vals,&
    !$OMP                 I_cut_k, I_det_k, I_det_ind,&
    !$OMP                 coo_c_all, coo_n_all, u_coo_n_all, nnz_arr, umax_block_size, t_all,&
    !$OMP                 criterion_a_mask, criterion_b_mask, criterion_c_mask, bit_map)&
    !$OMP PRIVATE(i, ii, j, jj, k, kk, i_exc, j_exc, row, pidx, col, target_row, target_col, cur_exc_degree,&
    !$OMP         ID, scn_a, target_row_det, target_col_det, block_a, block_b, criterion_a, criterion_b, criterion_c, add_pair,&
    !$OMP         buffer_count, nts_buffer_count, max_new_entries, buffer_total, ubuffer_total,&
    !$OMP         block_start, block_end, nts_block_start, nts_block_end, block_total, max_block_size,&
    !$OMP         coo_r, coo_c, coo_n, coo_r_nts, coo_c_nts, coo_n_nts, u_coo_r, u_coo_c, u_coo_n, u_block_rows,&
    !$OMP         t_coo_r, t_coo_c, sort_idx, row_start, row_end, old_row, t0, t1, t_tot, t_main_loop,&
    !$OMP         bit_degree, pass_a, pass_b, pass_c)
    
    !$ ID = OMP_get_thread_num() + 1
    
    !! calculate row partitions
    !$OMP SINGLE
    !$ n_threads = OMP_get_num_threads()
    allocate(pointer_blocks(N_det_g+1), row_starts(n_threads+1), coo_n_all(N_det_u+1))
    allocate(nnz_arr(N_det_u, n_threads+1))
    nnz_arr = 0
    nnz_arr(:, 1) = 1
    coo_n_all = 0
    pointer_blocks = 0
    row_starts = 0
    target_N = nnz_g / n_threads
    t_all = 0.0
    !$OMP END SINGLE
    !$OMP BARRIER

    t_tot = 0.0

    ! split the work so that different groups of contiguous rows have roughly equal entries
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det_g+1
        ! if ((H_p(i)/target_N + 1) < 0)
        pointer_blocks(i) = H_p(i) / target_N + 1
    end do
    !$OMP END DO
    !$OMP BARRIER

    ii = (N_det_g / n_threads)*(ID-1) + 1 ! start of range
    kk = (N_det_g / n_threads) * ID + 3 ! end of range, slight overlap to catch boundary pointers
    if (ID == n_threads) kk = N_det_g
    old_row = pointer_blocks(ii)
    
    do i = ii, kk 
        if(pointer_blocks(i) .ne. old_row) then 
            row_starts(pointer_blocks(i)) = i
            old_row = pointer_blocks(i)
        end if 
    end do
    
    !$OMP SINGLE
    row_starts(1) = 1
    row_starts(n_threads+1) = N_det_g+1
    !$OMP END SINGLE
    !$OMP BARRIER

    row_start = row_starts(ID)
    row_end = row_starts(ID+1) - 1
    if (ID == n_threads) row_end = N_det_g
    
    allocate(coo_r(8192), coo_c(8192), coo_n(N_det_u+1), coo_r_nts(8192), coo_c_nts(8192), coo_n_nts(N_det_u+1)) ! private
    buffer_count = 0
    nts_buffer_count = 0
    coo_n = 0
    coo_n_nts = 0

    !! this loop accomplishes the equivalent task as the loop over determinants in get_sparsity_structure()
    !! loop over excitation pairs

    call wall_time(t0)
    do i_exc = 1, N_exc
        do j_exc = 1, N_exc

            !! loop over matrix entries
            !! since this loop is done multiple times, can be sped up with row partitioning
            do row = row_start, row_end

                ! row cannot accept current i excitation
                if (I_cut_k(row, i_exc) == 0) then
                    cycle
                end if
                
                max_new_entries = H_p(row+1) - H_p(row)
                !! reallocate buffers if necessary
                if (max_new_entries > (size(coo_r,1) - buffer_count - 100)) then ! reallocate sorted buffer
                    allocate(t_coo_r(buffer_count + 4 *max_new_entries), t_coo_c(buffer_count + 4*max_new_entries))

                    t_coo_r(:buffer_count) = coo_r(:buffer_count)
                    t_coo_c(:buffer_count) = coo_c(:buffer_count)

                    call move_alloc(t_coo_r, coo_r)
                    call move_alloc(t_coo_c, coo_c)
                end if
                
                if (max_new_entries > (size(coo_r_nts,1) - nts_buffer_count - 100)) then ! reallocate unsorted buffer
                    allocate(t_coo_r(nts_buffer_count + 4 *max_new_entries), t_coo_c(nts_buffer_count + 4*max_new_entries))
    
                    t_coo_r(:nts_buffer_count) = coo_r_nts(:nts_buffer_count)
                    t_coo_c(:nts_buffer_count) = coo_c_nts(:nts_buffer_count)
    
                    call move_alloc(t_coo_r, coo_r_nts)
                    call move_alloc(t_coo_c, coo_c_nts)
                end if

                ! target_row_det = I_det_k(:, :, row, i_exc)
                ! target_row = hash_value(hash_alpha, hash_beta, hash_vals, hash_prime, hash_table_size, N_exc,&
                !                         target_row_det(:,1), target_row_det(:,2))

                target_row = I_det_ind(row, i_exc)

                ! if (target_row /= I_det_ind(row, i_exc)) then
                !     print *, i_exc, row, target_row, I_det_ind(row, i_exc)
                ! end if

                !! handle diagonal excitations out of loop
                if (i_exc == j_exc) then
                    buffer_count += 1
                    coo_r(buffer_count) = target_row
                    coo_c(buffer_count) = target_row
                else if (( j_exc > i_exc ) .and. ( I_cut_k(row, j_exc) == 1 )) then
                    
                    ! target_col_det = I_det_k(:, :, row, j_exc)
                    ! target_col = hash_value(hash_alpha, hash_beta, hash_vals, hash_prime, hash_table_size, N_exc,&
                    !                         target_col_det(:,1), target_col_det(:,2))
                    
                    target_col = I_det_ind(row, j_exc)

                    ! if (target_col /= I_det_ind(row, j_exc)) then
                    !     print *, j_exc, row, target_col, I_det_ind(row, j_exc)
                    ! end if
                    
                    if (target_col < target_row) then
                        nts_buffer_count += 1
                        coo_r_nts(nts_buffer_count) = target_col
                        coo_c_nts(nts_buffer_count) = target_row
                    else if (target_col > target_row) then
                        buffer_count += 1
                        coo_r(buffer_count) = target_row
                        coo_c(buffer_count) = target_col
                    end if

                end if
                
                
                !! handle off diagonal excitations
                do pidx = H_p(row) + 1, H_p(row+1) - 1
                    
                    
                    !! add in new entries from determinant pair
                    col = H_c(pidx)

                    
                    if (I_cut_k(col, j_exc) == 1) then
                        ! target_col_det = I_det_k(:, :, col, j_exc)
                        ! target_col = hash_value(hash_alpha, hash_beta, hash_vals, hash_prime, hash_table_size, N_exc,&
                        !                         target_col_det(:,1), target_col_det(:,2))

                        ! if (target_col /= I_det_ind(col, j_exc)) then
                        !     print *, j_exc, col, target_col, I_det_ind(col, j_exc)
                        ! end if
                        target_col = I_det_ind(col, j_exc)

                        ! if (target_col == -1) then
                        !     print *, "off_diag", ID, i_exc, j_exc, row, col
                        ! end if
                        
                        cur_exc_degree = H_e(pidx)
                        bit_degree = bit_map(H_e(pidx))

                        if (i_exc == j_exc) then
                            ! add in edges from intrinsic space
                            ! avoid any triple excitation pairs
                            if ((cur_exc_degree < 7) .and. (cur_exc_degree /= 3 )) then
                                if (target_col < target_row) then
                                    nts_buffer_count += 1
                                    coo_r_nts(nts_buffer_count) = target_col
                                    coo_c_nts(nts_buffer_count) = target_row
                                else if (target_col > target_row) then
                                    buffer_count += 1
                                    coo_r(buffer_count) = target_row
                                    coo_c(buffer_count) = target_col
                                end if
                            end if

                        else
                            !! following schema for determining how many determinants to add relies on fact that
                            !! applied excitations are within the same spin channel

                            ! ! both excitations are accepted in both source determinants, and alpha degree increases by 1
                            ! ! i, j are common alpha-occupied shells
                            ! criterion_a = (I_cut_k(row, i_exc) == I_cut_k(col, i_exc)) .and. (I_cut_k(row, j_exc) == I_cut_k(col, j_exc))

                            ! ! source determinants accept exclusive excitations, and alpha degree reduces by 1
                            ! ! i, j are exclusive alpha-occupied shells 
                            ! criterion_b = (I_cut_k(row, i_exc) /= I_cut_k(col, i_exc)) .and. (I_cut_k(row, j_exc) /= I_cut_k(col, j_exc))

                            
                            ! ! one source determinant can accept both excitations, while other can accept only one; overall degree stays same
                            ! ! i/j are common alpha-occupied alpha; j/i are exclusive alpha-occupied
                            ! criterion_c = ( (I_cut_k(row, i_exc) /= I_cut_k(col, i_exc)) /= (I_cut_k(row, j_exc) /= I_cut_k(col, j_exc)) )&
                            !         .and. ( (I_cut_k(row, i_exc) == I_cut_k(col, i_exc)) /= (I_cut_k(row, j_exc) == I_cut_k(col, j_exc)) )
                                    
                                    
                            block_a = I_cut_k(row, i_exc) == I_cut_k(col, i_exc)
                            block_b = I_cut_k(row, j_exc) == I_cut_k(col, j_exc)
                            
                            ! both excitations are accepted in both source determinants, and alpha degree increases by 1
                            ! i, j are common alpha-occupied shells
                            criterion_a = block_a .and. block_b


                            ! source determinants accept exclusive excitations, and alpha degree reduces by 1
                            ! i, j are exclusive alpha-occupied shells 
                            criterion_b = (.not. block_a) .and. (.not. block_b)

                            
                            ! one source determinant can accept both excitations, while other can accept only one; overall degree stays same
                            ! i/j are common alpha-occupied alpha; j/i are exclusive alpha-occupied
                            criterion_c = block_a /= block_b

                            pass_a = iand(bit_degree, criterion_a*criterion_a_mask)
                            pass_b = iand(bit_degree, criterion_b*criterion_b_mask)
                            pass_c = iand(bit_degree, criterion_c*criterion_c_mask)
                            add_pair = ior(ior(pass_a, pass_b), pass_c)
                            ! add_pair = .false.
                            ! select case (cur_exc_degree)
                            !     case (1) ! (1,0) -> (1,0), (2,0)
                            !         !add_pair = criterion_c .or. ~criterion_c ! if c then (1,0), else (2,0)
                            !         add_pair = .true. 
                            !     case (2) ! (2,0) -> (1,0), (2,0)
                            !         add_pair = criterion_b .or. criterion_c 
                            !     case (3) ! (3,0) -> (2,0)
                            !         add_pair = criterion_b
                            !     case (4) ! (0,1) -> (0,1), (1,1) 
                            !         add_pair = criterion_a .or. criterion_c
                            !     case (5) ! (0,2) -> (0,2)
                            !         add_pair = criterion_c
                            !     case (6) ! (1,1) -> (0,1), (1,1)
                            !         add_pair = criterion_b .or. criterion_c
                            !     case (7) ! (1,2) -> (0,2)
                            !         add_pair = criterion_b
                            !     case (8) ! (2,1) -> (1,1)
                            !         add_pair = criterion_b
                            ! end select

                            if ( add_pair ) then
                                if (target_col < target_row) then
                                    nts_buffer_count += 1
                                    coo_r_nts(nts_buffer_count) = target_col
                                    coo_c_nts(nts_buffer_count) = target_row
                                else if (target_col > target_row) then
                                    buffer_count += 1
                                    coo_r(buffer_count) = target_row
                                    coo_c(buffer_count) = target_col
                                end if
                            end if

                        end if 

                    end if

                end do

            end do

        end do
    end do
    call wall_time(t1)
    t_main_loop = t1-t0
    !! sort buffers per task; reduce combined buffers to unique buffers, blocked by per row

    ! standard buffer is blocked by row but blocks aren't sorted
    ! NTS is blocked by col, rows aren't sorted
    call wall_time(t0)
    call double_sort(coo_r, coo_c, buffer_count, coo_n, N_det_u)
    call double_sort(coo_r_nts, coo_c_nts, nts_buffer_count, coo_n_nts, N_det_u)
    call wall_time(t1)
    t_tot += t1-t0

    ! combine buffers by row block
    buffer_total = nts_buffer_count + buffer_count
    allocate(u_coo_c(buffer_total), u_coo_n(N_det_u+1))
    u_coo_c = 0
    u_coo_n = 0
    u_coo_n(1) = 1

    ! scan block sizes to get a maximum block allocation
    max_block_size = 1
    do i = 1, N_det_u
        block_total = (coo_n(i+1) - coo_n(i)) + (coo_n_nts(i+1) - coo_n_nts(i))
        max_block_size = merge(block_total, max_block_size, block_total > max_block_size)
    end do
    
    ! reduce to unique entires by block and add to ubuffers
    allocate(t_coo_c(max_block_size), sort_idx(max_block_size))
    do i = 1, N_det_u
        block_start = coo_n(i)
        block_end = coo_n(i+1)-1
        nts_block_start = coo_n_nts(i)
        nts_block_end = coo_n_nts(i+1)-1

        ! print *, i, block_start, block_end, nts_block_start, nts_block_end
        block_total = (coo_n(i+1) - coo_n(i)) + (coo_n_nts(i+1) - coo_n_nts(i))

        if (block_total > 0) then
            ! combine blocks into buffer
            t_coo_c(1 : coo_n(i+1) - coo_n(i)) = coo_c(block_start:block_end)
            t_coo_c(coo_n(i+1) - coo_n(i) + 1 : block_total) = coo_c_nts(nts_block_start:nts_block_end)

            ! sort, then get unique entries
            do j = 1, block_total
                sort_idx(j) = j
            end do

            call wall_time(t0)
            call quick_isort(t_coo_c, sort_idx, block_total)
            call wall_time(t1)
            call unique_from_sorted_buffer(t_coo_c, block_total, ubuffer_total)
            t_tot += t1-t0

            ! add into compact unique buffers
            u_coo_n(i+1) = u_coo_n(i) + ubuffer_total
            u_coo_c(u_coo_n(i):u_coo_n(i+1)-1) = t_coo_c(:ubuffer_total)

            nnz_arr(i, ID+1) = ubuffer_total
        else
            u_coo_n(i+1) = u_coo_n(i)
            nnz_arr(i, ID+1) = 0
        end if

    end do
    deallocate(t_coo_c, sort_idx)

    !$OMP SINGLE
    coo_n_all = 0
    coo_n_all(1) = 1
    !$OMP END SINGLE
    allocate(u_block_rows(N_det_u, 2))
    u_block_rows = 0
    !$OMP BARRIER

    !! combine buffers into shared memory from unique components by task
    !$OMP SINGLE
    print *, "reducing unique buffers into shared memory"
    !$OMP END SINGLE
    !$OMP BARRIER

    ! reduce counts over threads and set pointer ranges
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det_u
        ! get total number in block across threads
        scn_a = 0
        !$OMP SIMD REDUCTION(inscan, +:scn_a)
        do j = 1, n_threads + 1
            scn_a = scn_a + nnz_arr(i, j)
            !$OMP SCAN INCLUSIVE(scn_a)
            nnz_arr(i, j) = scn_a
        end do
        coo_n_all(i+1) = scn_a - 1
    end do
    !$OMP END DO
    !$OMP BARRIER

    !$OMP SINGLE
    scn_a = 0
    !$OMP SIMD REDUCTION(inscan, +:scn_a)
    do i = 1, N_det_u+1
        scn_a = scn_a + coo_n_all(i)
        !$OMP SCAN INCLUSIVE(scn_a)
        coo_n_all(i) = scn_a 

    end do
    !$OMP END SINGLE
    !$OMP BARRIER

    ! TODO: better checking of the coo_n blocks, something might not be set correctly, but only seems 
    ! to affect the end block

    do i = 1, N_det_u
        u_block_rows(i, 1) = coo_n_all(i) + nnz_arr(i, ID) - 1
        if (ID == n_threads) then
            u_block_rows(i, 2) = coo_n_all(i+1)
        else
            u_block_rows(i, 2) = coo_n_all(i) + nnz_arr(i, ID + 1) - 1
        end if
    end do

    ! while work is finishing, allocate some arrays and get a buffer size for temp array allocation later
    !$OMP SINGLE
    allocate(u_coo_n_all(N_det_u), coo_c_all(coo_n_all(N_det_u+1)))
    !$OMP END SINGLE

    !$OMP SINGLE
    umax_block_size = 1
    do i = 1, N_det_u 
        kk = coo_n_all(i+1) - coo_n_all(i)
        umax_block_size = merge(kk, umax_block_size, kk > umax_block_size)
    end do
    print *, "umax block size", umax_block_size
    !$OMP END SINGLE

    !$OMP BARRIER
    
    ! at this point, nnz_arr contains pointers for where each block owned by each thread begins
    ! assemble into shared memory, then reduce to unique entries
    ! in future implementation, this could be done with RMA or other global address space schemes
    do i = 1, N_det_u - 1
        block_start = u_block_rows(i, 1)
        block_end = u_block_rows(i, 2) - 1
        
        ! write(*, '(I4, I6, I8, I8, I8, I8, I8, I8)'), ID, i, block_start, block_end, block_end-block_start + 1, u_coo_n(i), u_coo_n(i+1), u_coo_n(i+1) - u_coo_n(i)
        coo_c_all(block_start:block_end) = u_coo_c(u_coo_n(i):u_coo_n(i+1)-1)
    end do

    !$OMP SINGLE
    ! since we're storing only the upper triangle, we can skip this is in the above loop
    coo_c_all(coo_n_all(N_det_u)) = N_det_u
    !$OMP END SINGLE
    !$OMP BARRIER

    !$OMP SINGLE
    print *, "Getting unique entries per shared row block"
    !$OMP END SINGLE
    !$OMP BARRIER

    ! get unique entries per shared row block
    allocate(sort_idx(umax_block_size))
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det_u
        do j = 1, umax_block_size
            sort_idx(j) = j
        end do

        block_start = coo_n_all(i)
        block_end = coo_n_all(i+1) - 1
        if ((block_end - block_start + 1) > 0) then ! this shouldn't be needed
            call wall_time(t0)
            call quick_isort(coo_c_all(block_start:block_end), sort_idx, block_end-block_start+1)
            call wall_time(t1)
            t_tot += t1-t0
            call unique_from_sorted_buffer(coo_c_all(block_start:block_end), block_end-block_start+1, ubuffer_total)
            u_coo_n_all(i) = ubuffer_total
        else
            u_coo_n_all(i) = 0
        end if
    end do
    !$OMP END DO
    deallocate(sort_idx)
    !$OMP BARRIER

    !$OMP SINGLE
    print *, "Assembling into CSR format"
    !$OMP END SINGLE
    !$OMP BARRIER

    !! assemble into CSR format

    ! set up row pointers
    !$OMP SINGLE
    uH_p(1) = 1
    uH_p(2:N_det_u+1) = u_coo_n_all
    scn_a = 0
    !$OMP SIMD REDUCTION(inscan, +:scn_a)
    do i = 1, N_det_u+1
        scn_a = scn_a + uH_p(i)
        !$OMP SCAN INCLUSIVE(scn_a)
        uH_p(i) = scn_a
    end do

    print *, "Total nonzero entries in union Hamiltonian: ", uH_p(N_det_u+1)-1, " max size: ", nnz_max
    !$OMP END SINGLE
    !$OMP BARRIER

    ! fill in matrix
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det_u
        block_start = coo_n_all(i)
        block_end = block_start + u_coo_n_all(i) - 1
        uH_c(uH_p(i):uH_p(i+1)-1) = coo_c_all(block_start:block_end)
    end do
    !$OMP END DO
    !$OMP BARRIER

    !! remaining clean up
    !$OMP SINGLE
    ! deallocate(pointer_blocks, row_starts, coo_n_all, coo_c_all, u_coo_n_all, nnz_arr)
    deallocate(pointer_blocks, row_starts, coo_n_all, nnz_arr)
    !$OMP END SINGLE

    deallocate(coo_r, coo_c, coo_n, coo_r_nts, coo_c_nts, coo_n_nts, u_coo_c, u_coo_n, u_block_rows)
    !$OMP BARRIER

    !$OMP CRITICAL
    t_all += t_tot
    ! print *, ID, t_tot
    !$OMP END CRITICAL
    !$OMP BARRIER
    !$OMP SINGLE
    print *, "Average time per thread spent sorting: ", t_all/n_threads
    t_all = 0
    !$OMP END SINGLE
    !$OMP BARRIER
    !$OMP CRITICAL
    t_all += t_main_loop
    ! print *, ID, t_main_loop
    !$OMP END CRITICAL
    !$OMP BARRIER
    !$OMP SINGLE
    print *, "Average time per thread spent in main loop: ", t_all/n_threads
    !$OMP END SINGLE
    !$OMP END PARALLEL
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

        ! kcol = psi_bilinear_matrix_columns(k)

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


    call quick_isort(all_idx, srt_idx, n_off_diagonal+1)

    columns(:n_off_diagonal+1) = all_idx
    nnz = n_off_diagonal + 1
    row = kidx
    
    deallocate(buffer, idx, singles_a, doubles_aa,&
                singles_b, doubles_bb, doubles_ab,&
                srt_idx, all_idx)
end

subroutine get_all_sparse_columns_with_triples(k_a, columns, exc_degree, row, nnz, nnz_max, N_det_l)
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
    ! allocate(buffer(N_int, N_det), idx(N_det))
    allocate(doubles_ab(N_det))

    n_buffer = 0
    n_buffer_old = 0
    n_doubles_ab_tot = 0

    ! starting from list of beta singles does not give all the (1,1) excitations
    ! so we need to search over either all beta or all alpha at some point
    ! start from (0,1) to excite to (1,1)
    ! do j = 1, n_det_beta_unique

    !     ! kcol = psi_bilinear_matrix_columns(j)

    !     tmp_det2(:,2) = psi_det_beta_unique (:, j)

    !     ! check if a single excitation different
    !     tdegree_beta = 0
    !     do i = 1, N_int
    !         tdegree_beta += popcnt(ieor(tmp_det2(i,2), ref_det(i,2)))
    !     end do

    !     if (tdegree_beta /= 2) then
    !         cycle
    !     end if

    !     ! loop over same beta different alpha
    !     do i = psi_bilinear_matrix_columns_loc(j), psi_bilinear_matrix_columns_loc(j+1)-1
    !         lidx = psi_bilinear_matrix_order(i)

    !         ! check if determinant is in upper half of reduced Hamiltonian matrix
    !         if ((lidx <= kidx) .or. (lidx > N_det_l)) cycle

    !         lcol = psi_bilinear_matrix_columns(i)
    !         lrow = psi_bilinear_matrix_rows(i)

    !         tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
    !         tmp_det(:,2) = psi_det_beta_unique (:, lcol)

    !         ! add determinant to buffer
    !         ! buffer is list of alpha spin determinants
    !         n_buffer += 1
    !         buffer(:,n_buffer) = tmp_det(:,1)
    !         idx(n_buffer) = lidx
    !     end do

    !     sdet_a = ref_det(:,1)

    !     ! all determinants are (X,1) excitations from ref_det
    !     ! so we just need to check alpha channel now
    !     ! grab indices of all determinants in buffer related to ref_det by (1,1) excitations 
    !     call get_all_spin_singles(buffer(:,n_buffer_old+1:n_buffer), idx(n_buffer_old+1:n_buffer),&
    !                             sdet_a, N_int, n_buffer-n_buffer_old,&
    !                             doubles_ab(n_doubles_ab_tot+1:n_doubles_ab_tot+n_buffer-n_buffer_old),&
    !                             n_doubles_ab)


    !     n_buffer_old = n_buffer
    !     n_doubles_ab_tot += n_doubles_ab
    ! end do


    !!! Finding triple excitations of form (3,0), (2,1), and (1,2)
    ! also, find (1,1) at the same time
    ! start from (0,0), (0,1), (0,2) to excite to (3,0), (2,1), (1,2)
    n_triples_tot = 0
    allocate(triples(N_det_l), triples_case(N_det_l))

    integer :: exc_map(0:4)

    exc_map = (/3, 0, 8, 0, 7/)
    do j = 1, n_det_beta_unique

        ! kcol = psi_bilinear_matrix_columns(k)

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
            do k = 1, N_int
                tdegree_alpha += popcnt(ieor(tmp_det(k,1), ref_det(k,1)))
            end do

            ! add to triples buffer 
            if (tdegree_alpha + tdegree_beta == 6) then
                n_triples_tot += 1
                triples(n_triples_tot) = lidx
                triples_case(n_triples_tot) = exc_map(tdegree_beta)
                ! select case (tdegree_beta)
                !     case (0)
                !         triples_case(n_triples_tot) = 3
                !     case (2)
                !         triples_case(n_triples_tot) = 8
                !     case (4)
                !         triples_case(n_triples_tot) = 7
                !     case default
                !         write(*, ("(A40, I16, I16, I16)")),"Error in triple selection on krow/beta det/alpha det: ", k_a, j, i
                !         exit
                ! end select
            else if ((tdegree_alpha == 2) .and. (tdegree_beta == 2)) then ! add (1,1) excitations in same loop
                n_doubles_ab_tot += 1
                doubles_ab(n_doubles_ab_tot) = lidx
            end if

        end do
    end do


    !!! Create final buffers
    ! add all indices into list and sort
    ! number of off-diagonal terms needed to caclulate to fill this row in sparse Hamlitonian
    n_off_diagonal = n_singles_a + n_singles_b + n_doubles_aa + n_doubles_ab_tot + n_doubles_bb + n_triples_tot

    ! allocate(all_degree(n_off_diagonal+1), all_idx(n_off_diagonal+1), srt_idx(n_off_diagonal+1))

    ! do i = 1, n_off_diagonal+1
    !     srt_idx(i) = i
    ! end do
  
    ! a - 1; aa - 2; aaa -3; b - 4; bb - 5; ab - 6; abb - 7; aab - 8
    n_offset = 0
    columns(n_offset+1:n_offset+1)               = kidx
    exc_degree(n_offset+1:n_offset+1)            = 0
                            n_offset            += 1

    columns(n_offset+1:n_offset+n_singles_a)     = singles_a(:n_singles_a)
    exc_degree(n_offset+1:n_offset+n_singles_a)  = 1
                            n_offset            += n_singles_a

    columns(n_offset+1:n_offset+n_singles_b)     = singles_b(:n_singles_b)
    exc_degree(n_offset+1:n_offset+n_singles_b)  = 3
                            n_offset            += n_singles_b

    columns(n_offset+1:n_offset+n_doubles_aa)    = doubles_aa(:n_doubles_aa)
    exc_degree(n_offset+1:n_offset+n_doubles_aa) = 2
                            n_offset            += n_doubles_aa

    columns(n_offset+1:n_offset+n_doubles_ab_tot)    = doubles_ab(:n_doubles_ab_tot)
    exc_degree(n_offset+1:n_offset+n_doubles_ab_tot) = 6
                            n_offset                += n_doubles_ab_tot

    columns(n_offset+1:n_offset+n_doubles_bb)     = doubles_bb(:n_doubles_bb)
    exc_degree(n_offset+1:n_offset+n_doubles_bb)  = 5
                            n_offset             += n_doubles_bb

    columns(n_offset+1:n_offset+n_triples_tot)    = triples(:n_triples_tot)
    exc_degree(n_offset+1:n_offset+n_triples_tot) = triples_case(:n_triples_tot)

    ! call quick_isort(all_idx, srt_idx, n_off_diagonal+1)

    ! columns(:n_off_diagonal+1) = all_idx

    ! do i = 1, n_off_diagonal+1
        ! exc_degree(i) = all_degree(srt_idx(i))
    ! end do

    nnz = n_off_diagonal + 1
    row = kidx
    
    ! deallocate(buffer, idx, singles_a, doubles_aa,&
    deallocate(singles_a, doubles_aa,&
                singles_b, doubles_bb, doubles_ab,&
                triples, triples_case)!,&
                ! srt_idx, all_idx, all_degree)
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

    integer(kind=8) :: i, j
    integer, intent(in) :: sze
    integer(kind=8), intent(in) :: nnz, H_p(sze+1)
    integer, intent(in) :: H_c(nnz)
    integer(bit_kind), intent(in) :: dets(N_int, 2, sze)
    double precision, intent(out) :: H_v(nnz)
    double precision :: hij

    !$OMP PARALLEL PRIVATE(i, j, hij) SHARED(H_c, H_p, H_v, sze, dets)
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, sze                
        ! loop over columns
        do j = H_p(i), H_p(i+1)-1
            call i_H_j(dets(:,:,i), dets(:,:,H_c(j)), N_int, hij)
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

subroutine unique_from_sorted_buffer(arr, n_max, n_out)
    BEGIN_DOC
    ! Pack unique values into front of array, in-place
    ! Assumes array is sorted
    END_DOC
    implicit none

    integer, intent(in)   :: n_max
    integer, intent(inout) :: arr(n_max)
    integer, intent(out)  :: n_out
    integer :: i


    n_out = 1
    do i = 2, n_max
        if (arr(i) == arr(i-1)) then
            cycle
        end if
        n_out += 1
        arr(n_out) = arr(i)
    end do

end

subroutine double_sort(rows, cols, n_max, row_starts, N_det_u)
    implicit none

    integer, intent(in)    :: n_max, N_det_u
    integer, intent(inout) :: rows(n_max), cols(n_max), row_starts(N_det_u+1)

    integer :: i, j, cur_row, prev_row, n_cols
    integer, allocatable :: tcols(:), sort_idx(:)

    if (n_max == 0) return

    allocate(sort_idx(n_max), tcols(n_max))
    do i = 1, n_max
        sort_idx(i) = i
    end do

    call quick_isort(rows, sort_idx, n_max) ! rows are sorted
  
    ! sort cols by row, and calculate all row ptrs
    prev_row = rows(1)
    row_starts(1) = 1
    if (prev_row > 1) row_starts(:prev_row) = 1

    tcols = cols
    do i = 1, n_max 
        cols(i) = tcols(sort_idx(i))

        cur_row = rows(i)
        if (cur_row /= prev_row) then
            row_starts(prev_row+1:cur_row) = i
            prev_row = cur_row
        end if
    end do
    row_starts(cur_row+1:) = n_max + 1

    ! sort cols within row block
    do i = 1, N_det_u
        n_cols = row_starts(i+1) - row_starts(i)
        do j = 1, n_cols
            sort_idx(j) = j
        end do

        if (n_cols > 0) then
            call quick_isort(cols(row_starts(i):row_starts(i+1)-1), sort_idx(:n_cols), n_cols)
        end if
    end do


    deallocate(sort_idx, tcols)
end