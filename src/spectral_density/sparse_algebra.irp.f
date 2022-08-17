BEGIN_PROVIDER [integer, nnz_max_per_row]
&BEGIN_PROVIDER [integer, nnz_max_tot]

    BEGIN_DOC
    ! Total possible number of entries in each column/row in full determinant space
    ! Vastly overestimates actual number needed
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
    END_DOC

    integer,          intent(in)  :: sze, nnz, A_c(nnz), A_p(sze+1)
    double precision, intent(in)  :: A_v(nnz), x(sze)
    double precision, intent(out) :: y(sze)
    double precision, allocatable :: y_t(:)
    integer                       :: i, j

    ! loop over rows
    !$OMP PARALLEL PRIVATE(i, j, y_t) SHARED(y, x, A_c, A_p, A_v, sze)
    allocate(y_t(sze))
    y_t = 0.d0
    !$OMP BARRIER
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, sze
        ! loop over columns

        
        ! make sure row actually is in reduced determinant space
        if (A_p(i+1) - A_p(i) > 0) then
            ! calculate diagonal separately to avoid double counting
            ! first entry per column is guaranteed to be diagonal since all diagonal
            ! elements of H are nonzero
            j = A_p(i)
            y_t(i) = y_t(i) + A_v(j) * x(A_c(j))
        end if

        do j = A_p(i)+1, A_p(i+1)-1
            ! calculate element of owned row
            y_t(i) = y_t(i) + A_v(j) * x(A_c(j))

            ! calculate element of owned column
            y_t(A_c(j)) = y_t(A_c(j)) + A_v(j) * x(i)
        end do
    end do
    !$OMP END DO

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

    ! loop over rows
    !$OMP PARALLEL PRIVATE(i, j, y_t) SHARED(y, x, A_c, A_p, A_v, sze)

    allocate(y_t(sze))
    y_t = (0.d0, 0.d0)
    !$OMP BARRIER
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, sze
        ! loop over columns


            ! make sure row actually is in reduced determinant space
        if (A_p(i+1) - A_p(i) > 0) then
            ! calculate diagonal separately to avoid double counting
            ! first entry per column is guaranteed to be diagonal since all diagonal
            ! elements of H are nonzero
            j = A_p(i)
            y_t(i) = y_t(i) + A_v(j) * x(A_c(j))
        end if

        do j = A_p(i)+1, A_p(i+1)-1
            ! calculate element of owned row
            y_t(i) = y_t(i) + A_v(j) * x(A_c(j))

            ! calculate element of owned column
            y_t(A_c(j)) = y_t(A_c(j)) + conjg(A_v(j)) * x(i)
        end do
    end do
    !$OMP END DO

    !$OMP CRITICAL
    y = y + y_t
    !$OMP END CRITICAL

    deallocate(y_t)
    !$OMP END PARALLEL
end

subroutine form_sparse_dH(csr_s, csr_c, csr_v, sze, dets, iorb, ispin, ac_type, N_det_l)
    ! use MKL_SPBLAS

    implicit none
    BEGIN_DOC
    ! Form a compressed sparse row matrix representation of the Hamiltonian
    ! in the space of the determinants
    END_DOC

    integer, intent(in)           :: iorb, ispin, N_det_l
    integer(kind=8), intent(in)        :: sze
    integer(bit_kind), intent(in) :: dets(N_int, 2, N_det_l)
    logical, intent(in)           :: ac_type
    integer, intent(out)          :: csr_s(N_det_l+1), csr_c(sze)
    double precision, intent(out) :: csr_v(sze)

    integer              :: n, i, j, k, l_row, old_row
    integer              :: nnz, nnz_cnt, nnz_tot
    integer              :: n_vals, n_vals_row, n_threads, ID
    integer              :: nnz_csr, ii, scn_a, kk
    integer :: OMP_get_num_threads, OMP_get_thread_num
    integer, allocatable :: nnz_arr(:), coo_r(:), coo_c(:), l_cols(:)
    integer, allocatable :: coo_r_all(:), coo_c_all(:), coo_r_t(:), coo_c_t(:)
    ! integer             :: coo_s(N_det_l), coo_n(N_det_l), coo_c_n(N_det_l), coo_c_n_all(N_det_l)
    integer, allocatable:: coo_s(:), coo_n(:)
    double precision     :: hij, frac
    double precision, allocatable :: coo_v(:), coo_v_t(:), coo_v_all(:)
    
    ! force provide early so that threads don't each try to provide
    call i_H_j(dets(:,:,1), dets(:,:,1), N_int, hij) 

    !$OMP PARALLEL SHARED(nnz_tot, nnz_arr, nnz_csr, n_threads, dets, psi_det, N_det, N_det_l, N_int, nnz_max_per_row, n_vals_row,&
    !$OMP                 coo_r_all, coo_c_all, coo_v_all, csr_s, csr_c, csr_v, coo_s, coo_n)& 
    !$OMP PRIVATE(i,j,old_row, k,ii,kk, scn_a, ID, hij, nnz, nnz_cnt, coo_r, coo_c, coo_v, coo_r_t, coo_c_t, coo_v_t, l_cols, l_row) 

    !$ n_threads = OMP_get_num_threads()
    !$ ID = OMP_get_thread_num() + 1

    !$OMP SINGLE
    allocate(coo_s(N_det_l), coo_n(N_det_l))
    coo_n = 0
    !$OMP END SINGLE
    !$OMP BARRIER

    frac = 0.2
    n_vals = max(nint(N_det_l*N_det_l*frac/n_threads), 128)
    n_vals_row = nnz_max_per_row

    allocate(coo_r(n_vals))
    allocate(coo_c(n_vals))
    allocate(coo_v(n_vals))
    allocate(l_cols(n_vals_row))
    
    !$OMP SINGLE
    allocate(nnz_arr(n_threads))
    nnz_tot = 0
    !$OMP END SINGLE
    
    ! !$OMP SINGLE
    ! print *, "## Calculating nonzero entries"
    ! !$OMP END SINGLE
    
    nnz_cnt = 0
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det ! this loop needs to go over all the determinants, since this loop is not in determiant order but rather k_a order
        nnz = 0
        l_cols = 0
        call get_sparse_columns(i, l_cols, l_row, nnz, n_vals_row,&
                                 iorb, ispin, ac_type, N_det_l)

        if (nnz == 0) cycle

        ! reallocate arrays if necessary
        if (nnz_cnt + nnz > size(coo_r, 1)) then
            allocate(coo_r_t(nnz_cnt + 10*nnz))
            allocate(coo_c_t(nnz_cnt + 10*nnz))
            allocate(coo_v_t(nnz_cnt + 10*nnz))
            
            coo_r_t(:nnz_cnt) = coo_r
            coo_c_t(:nnz_cnt) = coo_c
            coo_v_t(:nnz_cnt) = coo_v
            
            call move_alloc(coo_r_t, coo_r)
            call move_alloc(coo_c_t, coo_c)
            call move_alloc(coo_v_t, coo_v)
        end if
        
        
        coo_n(l_row) = nnz ! store for later
        do j = 1, nnz
            nnz_cnt += 1

            call i_H_j(dets(:,:,l_row),&
                       dets(:,:,l_cols(j)), N_int, hij)
            
            coo_r(nnz_cnt) = l_row
            coo_c(nnz_cnt) = l_cols(j)
            coo_v(nnz_cnt) = hij
        end do
    end do
    !$OMP END DO

    nnz_arr(ID) = nnz_cnt

    !$OMP CRITICAl
    nnz_tot = nnz_tot + nnz_cnt
    !$OMP END CRITICAl
    !$OMP BARRIER
    
    !$OMP SINGLE
    print *, "Total non-zero entries in Hamiltonian: ", nnz_tot, " max size:", sze
    ! print *, "## Constructing pointer arrays"
    !$OMP END SINGLE
    
    !$OMP SINGLE
    allocate(coo_r_all(nnz_tot))
    allocate(coo_c_all(nnz_tot))
    allocate(coo_v_all(nnz_tot))
    !$OMP END SINGLE
    
    k = 0
    do i = 1, ID-1
        k += nnz_arr(i)
    end do

    coo_r_all(k+1:k+nnz_arr(ID)) = coo_r
    coo_c_all(k+1:k+nnz_arr(ID)) = coo_c
    coo_v_all(k+1:k+nnz_arr(ID)) = coo_v
    !$OMP BARRIER

    ! calculate the starting index of each row in COO, since they were not sorted
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

    !$OMP SINGLE
    coo_s(coo_r_all(1)) = 1
    !$OMP END SINGLE

    ! calculate CSR pointer ranges
    ! row i data goes from csr_s(i) to csr_s(i+1) - 1
    ! first, count all entries in parallel, then perform scan to set pointer ranges
    ! then reduce with inclsuive scan

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

    ! !$OMP SINGLE
    ! print *, "## Constructing CSR arrays"
    ! !$OMP END SINGLE

    ! loop through rows and construct CSR matrix
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det_l
        csr_v(csr_s(i):csr_s(i+1)-1) = coo_v_all(coo_s(i):coo_s(i+1)-1)
        csr_c(csr_s(i):csr_s(i+1)-1) = coo_c_all(coo_s(i):coo_s(i+1)-1)
    end do
    !$OMP END DO

    deallocate(coo_r, coo_c, coo_v, l_cols)

    !$OMP SINGLE
    deallocate(coo_s, coo_n, nnz_arr, coo_r_all, coo_c_all, coo_v_all)
    !$OMP END SINGLE
    !$OMP END PARALLEL

end

subroutine form_sparse_zH(csr_s, csr_c, csr_v, sze, dets, iorb, ispin, ac_type, N_det_l)
    ! use MKL_SPBLAS

    implicit none
    BEGIN_DOC
    ! Form a compressed sparse row matrix representation of the Hamiltonian
    ! in the space of the determinants
    END_DOC

    integer, intent(in)           :: iorb, ispin, N_det_l
    integer(kind=8), intent(in)        :: sze
    integer(bit_kind), intent(in) :: dets(N_int, 2, N_det_l)
    logical, intent(in)           :: ac_type
    integer, intent(out)          :: csr_s(N_det_l+1), csr_c(sze)
    complex*16, intent(out) :: csr_v(sze)

    integer              :: n, i, j, k, l_row, old_row
    integer              :: nnz, nnz_cnt, nnz_tot
    integer              :: n_vals, n_vals_row, n_threads, ID
    integer              :: nnz_csr, ii, scn_a, kk
    integer :: OMP_get_num_threads, OMP_get_thread_num
    integer, allocatable :: nnz_arr(:), coo_r(:), coo_c(:), l_cols(:)
    integer, allocatable :: coo_r_all(:), coo_c_all(:), coo_r_t(:), coo_c_t(:)
    ! integer             :: coo_s(N_det), coo_n(N_det), coo_c_n(N_det), coo_c_n_all(N_det)
    integer, allocatable:: coo_s(:), coo_n(:)
    double precision     :: frac
    complex*16           :: hij
    complex*16, allocatable :: coo_v(:), coo_v_t(:), coo_v_all(:)
    

    ! force provide early so that threads don't each try to provide
    call i_h_j_complex(dets(:,:,1), dets(:,:,1), N_int, hij) 

    !$OMP PARALLEL SHARED(nnz_tot, nnz_arr, nnz_csr, n_threads, dets, psi_det, N_det, N_det_l, N_int, nnz_max_per_row, n_vals_row,&
    !$OMP                 coo_r_all, coo_c_all, coo_v_all, csr_s, csr_c, csr_v, coo_s, coo_n)& 
    !$OMP PRIVATE(i,j,old_row, k,ii,kk, scn_a, ID, hij, nnz, nnz_cnt, coo_r, coo_c, coo_v, coo_r_t, coo_c_t, coo_v_t, l_cols, l_row) 

    !$ n_threads = OMP_get_num_threads()
    !$ ID = OMP_get_thread_num() + 1

    !$OMP SINGLE
    allocate(coo_s(N_det_l), coo_n(N_det_l))
    coo_n = 0
    !$OMP END SINGLE
    !$OMP BARRIER

    frac = 0.2
    n_vals = max(nint(N_det_l*N_det_l*frac/n_threads), 128)
    n_vals_row = nnz_max_per_row

    allocate(coo_r(n_vals))
    allocate(coo_c(n_vals))
    allocate(coo_v(n_vals))
    allocate(l_cols(n_vals_row))
    
    !$OMP SINGLE
    allocate(nnz_arr(n_threads))
    nnz_tot = 0
    !$OMP END SINGLE
    
    ! !$OMP SINGLE
    ! print *, "## Calculating nonzero entries"
    ! !$OMP END SINGLE
    
    nnz_cnt = 0
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det ! this loop needs to go over all the determinants, since this loop is not in determiant order but rather k_a order
        nnz = 0
        l_cols = 0
        call get_sparse_columns(i, l_cols, l_row, nnz, n_vals_row,&
                                 iorb, ispin, ac_type, N_det_l)

        if (nnz == 0) cycle

        ! reallocate arrays if necessary
        if (nnz_cnt + nnz > size(coo_r, 1)) then
            allocate(coo_r_t(nnz_cnt + 10*nnz))
            allocate(coo_c_t(nnz_cnt + 10*nnz))
            allocate(coo_v_t(nnz_cnt + 10*nnz))
            
            coo_r_t(:nnz_cnt) = coo_r
            coo_c_t(:nnz_cnt) = coo_c
            coo_v_t(:nnz_cnt) = coo_v
            
            call move_alloc(coo_r_t, coo_r)
            call move_alloc(coo_c_t, coo_c)
            call move_alloc(coo_v_t, coo_v)
        end if
        
        
        coo_n(l_row) = nnz ! store for later
        do j = 1, nnz
            nnz_cnt += 1

            call i_h_j_complex(dets(:,:,l_row),&
                               dets(:,:,l_cols(j)), N_int, hij)
            
            coo_r(nnz_cnt) = l_row
            coo_c(nnz_cnt) = l_cols(j)
            coo_v(nnz_cnt) = hij
        end do
    end do
    !$OMP END DO

    nnz_arr(ID) = nnz_cnt

    !$OMP CRITICAl
    nnz_tot = nnz_tot + nnz_cnt
    !$OMP END CRITICAl
    !$OMP BARRIER
    
    !$OMP SINGLE
    print *, "Total non-zero entries: ", nnz_tot, " max size:", sze
    ! print *, "## Constructing pointer arrays"
    !$OMP END SINGLE
    
    !$OMP SINGLE
    allocate(coo_r_all(nnz_tot))
    allocate(coo_c_all(nnz_tot))
    allocate(coo_v_all(nnz_tot))
    !$OMP END SINGLE
    
    k = 0
    do i = 1, ID-1
        k += nnz_arr(i)
    end do

    coo_r_all(k+1:k+nnz_arr(ID)) = coo_r
    coo_c_all(k+1:k+nnz_arr(ID)) = coo_c
    coo_v_all(k+1:k+nnz_arr(ID)) = coo_v
    !$OMP BARRIER

    ! calculate the starting index of each row in COO, since they were not sorted
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

    !$OMP SINGLE
    coo_s(coo_r_all(1)) = 1
    !$OMP END SINGLE

    ! calculate CSR pointer ranges
    ! row i data goes from csr_s(i) to csr_s(i+1) - 1
    ! first, count all entries in parallel, then perform scan to set pointer ranges
    ! then reduce with inclsuive scan

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

    ! !$OMP SINGLE
    ! print *, "## Constructing CSR arrays"
    ! !$OMP END SINGLE

    ! loop through rows and construct CSR matrix
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, N_det_l
        csr_v(csr_s(i):csr_s(i+1)-1) = coo_v_all(coo_s(i):coo_s(i+1)-1)
        csr_c(csr_s(i):csr_s(i+1)-1) = coo_c_all(coo_s(i):coo_s(i+1)-1)
    end do
    !$OMP END DO

    deallocate(coo_r, coo_c, coo_v, l_cols)

    !$OMP SINGLE
    deallocate(coo_s, coo_n, nnz_arr, coo_r_all, coo_c_all, coo_v_all)
    !$OMP END SINGLE
    !$OMP END PARALLEL

end

subroutine get_sparse_columns(k_a, columns, row, nnz, nnz_max, iorb, ispin, ac_type, N_det_l)
    ! this whole function should be within a parallel loop over N_det
    ! the final output is the list of indices of off diagonal terms, sorted by column
    ! for the thread in the mainloop to calculate the (non-zero) matrix elements H_i, j>i
    ! for the hamiltonian in the space of the set of determinants
    ! ac_type == F if adding electron, T if removing
    implicit none

    integer, intent(in)      :: k_a, nnz_max, iorb, ispin, N_det_l
    logical, intent(in)      :: ac_type
    integer, intent(out)     :: nnz, columns(nnz_max), row
    integer :: i, j, k, k_b
    integer :: krow, kcol, lrow, lcol
    integer :: lidx, kidx, tidx, n_buffer
    integer :: n_singles_a, n_singles_b, n_doubles_aa, n_doubles_bb, n_doubles_ab
    integer :: n_buffer_old, n_doubles_ab_tot
    integer :: n_off_diagonal, n_offset, tdegree
    logical :: is_filled

    integer(bit_kind) :: ref_det(N_int, 2), tmp_det(N_int, 2), tmp_det2(N_int, 2), sdet_a(N_int), sdet_b(N_int)
    integer(bit_kind), allocatable    :: buffer(:,:)
    integer, allocatable              :: singles_a(:), singles_b(:)
    integer, allocatable              :: doubles_aa(:), doubles_bb(:), doubles_ab(:)
    integer, allocatable              :: idx(:), all_idx(:), srt_idx(:)
    
    allocate(buffer(N_int, N_det_l), idx(N_det_l))

    n_singles_a = 0
    n_singles_b = 0
    n_doubles_aa = 0
    n_doubles_ab = 0
    n_doubles_bb = 0
    
    kidx = psi_bilinear_matrix_order(k_a)

    if (kidx > N_det_l) return ! determinant not included in this subset

    krow = psi_bilinear_matrix_rows(k_a)
    kcol = psi_bilinear_matrix_columns(k_a)
    ref_det(:,1) = psi_det_alpha_unique(:, krow)
    ref_det(:,2) = psi_det_beta_unique (:, kcol)

    call orb_is_filled(ref_det, iorb, ispin, N_int, is_filled)
    if (is_filled .neqv. ac_type) return 
    
    ! Finding (1,0) and (2,0) excitations
    ! loop over same beta different alpha
    n_buffer = 0
    do i = psi_bilinear_matrix_columns_loc(kcol), psi_bilinear_matrix_columns_loc(kcol+1)-1
        lidx = psi_bilinear_matrix_order(i)
        
        if ((lidx <= kidx) .or. (lidx > N_det_l)) cycle
        ! if (lidx > N_det_l) cycle

        lcol = psi_bilinear_matrix_columns(i)
        lrow = psi_bilinear_matrix_rows(i)
        
        tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
        tmp_det(:,2) = psi_det_beta_unique (:, lcol)

        call orb_is_filled(tmp_det, iorb, ispin, N_int, is_filled)

        if (is_filled .neqv. ac_type) cycle 
        
        ! add determinant to buffer
        ! buffer is list of alpha spin determinants
        n_buffer += 1
        buffer(:,n_buffer) = tmp_det(:,1)
        idx(n_buffer) = lidx
    end do

    allocate(singles_a(n_buffer), doubles_aa(n_buffer))

    sdet_a = ref_det(:,1)
    call get_all_spin_singles_and_doubles(buffer, idx, sdet_a, &
                        N_int, n_buffer, singles_a, doubles_aa, n_singles_a, n_doubles_aa)

    deallocate(buffer, idx)
    allocate(buffer(N_int, N_det), idx(N_det))

    ! print *, "Getting beta singles/doubles"
    ! Finding (0,1) and (0,2) excitations 
    k_b = psi_bilinear_matrix_order_transp_reverse(k_a)
    krow = psi_bilinear_matrix_transp_rows(k_b) !this is unnecessary, technically
    kcol = psi_bilinear_matrix_transp_columns(k_b)
    
    ! loop over same alpha different beta
    n_buffer = 0
    do i = psi_bilinear_matrix_transp_rows_loc(krow), psi_bilinear_matrix_transp_rows_loc(krow+1)-1
        tidx = psi_bilinear_matrix_transp_order(i)
        lidx = psi_bilinear_matrix_order(tidx)

        if ((lidx <= kidx) .or. (lidx > N_det_l)) cycle
        ! if (lidx > N_det_l) cycle
        
        lcol = psi_bilinear_matrix_transp_columns(i)
        lrow = psi_bilinear_matrix_transp_rows(i)
        
        tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
        tmp_det(:,2) = psi_det_beta_unique (:, lcol)

        call orb_is_filled(tmp_det, iorb, ispin, N_int, is_filled)

        if (is_filled .neqv. ac_type) cycle 
                
        ! add determinant to buffer
        ! buffer is list of beta spin determinants
        n_buffer += 1
        buffer(:,n_buffer) = tmp_det(:,2)
        idx(n_buffer) = lidx
    end do
    
    allocate(singles_b(n_buffer), doubles_bb(n_buffer))
    sdet_b = ref_det(:,2)
    call get_all_spin_singles_and_doubles(buffer, idx, sdet_b, &
                        N_int, n_buffer, singles_b, doubles_bb, n_singles_b, n_doubles_bb)
                        
    deallocate(buffer, idx)

    ! print *, "Getting alpha beta doubles"
    ! Finding (1,1) excitations
    allocate(buffer(N_int, N_det), idx(N_det))
    allocate(doubles_ab(N_det))

    n_buffer = 0
    n_buffer_old = 0
    n_doubles_ab_tot = 0

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

            if ((lidx <= kidx) .or. (lidx > N_det_l)) cycle
            ! if (lidx > N_det_l) cycle

            lcol = psi_bilinear_matrix_columns(i)
            lrow = psi_bilinear_matrix_rows(i)

            tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
            tmp_det(:,2) = psi_det_beta_unique (:, lcol)

            call orb_is_filled(tmp_det, iorb, ispin, N_int, is_filled)

            if (is_filled .neqv. ac_type) cycle 

            ! add determinant to buffer
            ! buffer is list of alpha spin determinants
            n_buffer += 1
            buffer(:,n_buffer) = tmp_det(:,1)
            idx(n_buffer) = lidx
        end do

        sdet_a = ref_det(:,1)
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

subroutine orb_is_filled(key_ref,iorb,ispin,Nint,is_filled)
    use bitmasks
    implicit none
    BEGIN_DOC
    ! determine whether iorb, ispin is filled in key_ref
    ! key_ref has alpha and beta parts
    END_DOC
    integer, intent(in)            :: iorb, ispin, Nint
    integer(bit_kind), intent(in) :: key_ref(Nint,2)
    logical, intent(out) :: is_filled
    
    integer                        :: k,l
    
    ASSERT (iorb > 0)
    ASSERT (ispin > 0)
    ASSERT (ispin < 3)
    ASSERT (Nint > 0)
    
    ! k is index of the int where iorb is found
    ! l is index of the bit where iorb is found
    k = ishft(iorb-1,-bit_kind_shift)+1
    ASSERT (k >0)
    l = iorb - ishft(k-1,bit_kind_shift)-1
    ASSERT (l >= 0)
    is_filled = btest(key_ref(k,ispin),l)  
  end