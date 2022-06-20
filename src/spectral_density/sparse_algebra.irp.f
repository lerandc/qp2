! faster loop for getting nonzero entries in H matrix

! use sorted determinants
! form maps between unique alpha/beta to full determinant space
! have it be CSR pointer style: ranges are separated by unique alpha
! thus map is length N_unique_alpha_dets?
! how to perform a look up for where beta fall?

! get reference determinant from bilinear matrix indexing

! check if it can accept homo/lumo hole/particle

! if so
! calculate diagonal component


! generate list of dets A related by single alpha excitation
! generate list of dets B related by single beta excitation
! generate list of dets AA related by double alpha excitation
! generate list of dets BB related by double beta excitation
! from A, B, generate list of dets AB and BA related by double alpha/beta excitation

! when generating list, only add dets that:
! 1) are greater in index in sorting order
! 2) can accept hole/particle

! iterate over lists and calculate H_ij calculations for determinant pairs
! using maps, populate COO sparse format at all coordinates where det pairs correspond to
! populate in row -> column order (this should be implicit in the sorted ordering)

BEGIN_PROVIDER [integer, nnz_max_per_row]
&BEGIN_PROVIDER [integer, nnz_max_tot]

    BEGIN_DOC
    ! Total possible number of entries in each column/row in full determinant space
    ! Vastly overestimates actual number needed
    END_DOC
    implicit none

    integer     :: s_a, s_b, d_a, d_b

    s_a = elec_alpha_num * (mo_num - elec_alpha_num)
    s_b = elec_beta_num * (mo_num - elec_beta_num)
    
    d_a = ( elec_alpha_num * (elec_alpha_num - 1) / 2) * &
    ( (mo_num - elec_alpha_num) * (mo_num - elec_alpha_num - 1) / 2)
    
    d_b = ( elec_beta_num * (elec_beta_num - 1) / 2) * &
    ( (mo_num - elec_beta_num) * (mo_num - elec_beta_num - 1) / 2)
    
    write(*, '(I10, I10, I10, I10, I10)'),&
            elec_alpha_num, elec_beta_num, mo_num, s_a, s_b

    nnz_max_per_row = 1 + s_a + s_b + s_a*s_b + d_a + d_b
    nnz_max_tot = N_det * nnz_max_per_row

END_PROVIDER

! to convert to provider, provide all 4 array for intel MKL CSR representation
subroutine form_sparse_dH(finished)
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
    double precision, allocatable :: coo_v(:), coo_v_t(:), coo_v_all(:)

    ! force provide early so that threads don't each try to provide
    call i_H_j(psi_det(1,1,1), psi_det(1,1,1), N_int, hij) 
    
    !$OMP PARALLEL SHARED(n, nnz_tot, nnz_arr, n_threads, psi_det, N_int) &
    !$OMP PRIVATE(i,j,k, hij, nnz, ID, coo_r, coo_c, coo_v, coo_r_t, coo_c_t, coo_v_t) 
    !$ n_threads = OMP_get_num_threads()
    !$ ID = OMP_get_thread_num() + 1
    n = N_det
    frac = 0.2
    n_vals = max(nint(n*n*frac/n_threads), 128)

    allocate(coo_r(n_vals))
    allocate(coo_c(n_vals))
    allocate(coo_v(n_vals))

    !$OMP SINGLE
    allocate(nnz_arr(n_threads))
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

    ! ! convert to sparse matrix in CSR format
    ! type(SPARSE_MATRIX_T) :: sH
    ! integer               :: status, mkl_d_create_coo

    ! status = mkl_d_create_coo(sH, SPARSE_INDEX_BASE_ONE, &
    !                           n, n, nnz_tot,&
    !                           coo_r_all, coo_c_all, coo_v_all)

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
    
    !$OMP DO SCHEDULE(GUIDED) 
    do i = 2, n+1
        k = 0 ! count all entries with column i-1
        do j = 1, nnz_tot
            k = k + merge(1,0,coo_c_all(j)==i-1)
        end do
        csr_s(i) = coo_n(i-1) + k - 1 ! account for double counting
    end do
    !$OMP END DO

    ! need to strongly enforce synchronization here
    !$OMP BARRIER
    !$OMP SINGLE
    scn_a = 0
    !$OMP SIMD REDUCTION(inscan, +:scn_a)
    do i = 1, n+1
        scn_a = scn_a + csr_s(i)
        !$OMP SCAN INCLUSIVE(scn_a)
        csr_s(i) = scn_a
    end do
    !$OMP END SINGLE
    !$OMP BARRIER
    

    ! loop through rows and construct CSR matrix
    !$OMP DO SCHEDULE(GUIDED)
    do i = 1, n
        k = 0
        do j = 1, i-1 !TODO: possible to only loop over number of entries?
            ! get pre-diagonal entries
            ! search jth row for column i
            do ii = coo_s(j)+1, coo_s(j) + coo_n(j) -1
                if (coo_c_all(ii) == i) then
                    csr_v(csr_s(i) + k) = coo_v_all(ii)
                    csr_c(csr_s(i) + k) = j !coo_c_all(ii)
                    k += 1
                    exit
                end if
            end do
        ! get diagonal+ entries
        end do
        
        do ii = coo_s(i), coo_s(i) + coo_n(i) - 1
            csr_v(csr_s(i) + k) = coo_v_all(ii)
            csr_c(csr_s(i) + k) = coo_c_all(ii)
            k += 1
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    print *, 'CSR representation'
    print *, csr_v
    print *, csr_c
    print *, csr_s


    finished = .TRUE.
end

subroutine sparse_csr_dmv(A_v, A_c, A_p, x, y, sze, nnz)
    implicit none
    BEGIN_DOC
    ! Compute the matrix vector product y = Ax with sparse A in CSR format
    ! Not general. Assumes A is sze x sze; x, y are sze
    END_DOC

    integer,          intent(in)  :: sze, nnz, A_c(nnz), A_p(sze+1)
    double precision, intent(in)  :: A_v(nnz), x(sze)
    double precision, intent(out) :: y(sze)
    integer                      :: i, j

    ! loop over rows
    !$OMP PARALLEL DO PRIVATE(i, j, sze) SHARED(y, x, A_c, A_p, A_v)&
    !$OMP SCHEDULE(GUIDED)
    do i = 1, sze
        ! loop over columns
        do j = A_p(i), A_p(i+1)-1
            y(i) = y(i) + A_v(j) * x(A_c(j))
        end do
    end do
    !$OMP END PARALLEL DO
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
    integer                      :: i, j

    ! loop over rows
    !$OMP PARALLEL DO PRIVATE(i, j, sze) SHARED(y, x, A_c, A_p, A_v)&
    !$OMP SCHEDULE(GUIDED)
    do i = 1, sze
        ! loop over columns
        do j = A_p(i), A_p(i+1)-1
            y(i) = y(i) + A_v(j) * x(A_c(j))
        end do
    end do
    !$OMP END PARALLEL DO
end


subroutine get_sparse_columns(k_a, columns, nnz, nnz_max)
    ! needs to be a function, not a subroutine, if output size is unknown

    ! this whole function should be within a parallel loop over N_det
    ! the final output is the list of indices of off diagonal terms, sorted by column
    ! for the thread in the mainloop to calculate the (non-zero) matrix elements H_i, j>i
    ! for the hamiltonian in the space of the set of determinants
    implicit none

    integer, intent(in)      :: k_a, nnz_max
    integer, intent(out)     :: nnz, columns(nnz_max)
    integer :: i, j, k, k_b
    integer :: krow, kcol, lrow, lcol
    integer :: lidx, kidx, tidx, n_buffer
    integer :: n_singles_a, n_singles_b, n_doubles_aa, n_doubles_bb, n_doubles_ab
    integer :: n_buffer_old, n_doubles_ab_tot, n_doubles_ab_max, n_off_diagonal, n_offset

    integer(bit_kind) :: ref_det(N_int, 2), tmp_det(N_int, 2), tmp_det2(N_int, 2), sdet_a(N_int), sdet_b(N_int)
    integer(bit_kind), allocatable    :: buffer(:,:)
    integer, allocatable              :: singles_a(:), singles_b(:), doubles_aa(:), doubles_bb(:), doubles_ab(:)
    integer, allocatable              :: idx(:), all_idx(:), srt_idx(:)
    
    allocate(buffer(N_int, N_det), idx(N_det))

    n_singles_a = 0
    n_singles_b = 0
    n_doubles_aa = 0
    n_doubles_ab = 0
    n_doubles_bb = 0
    
    kidx = psi_bilinear_matrix_order(k_a)
    krow = psi_bilinear_matrix_rows(k_a)
    kcol = psi_bilinear_matrix_columns(k_a)
    ref_det(:,1) = psi_det_alpha_unique(:, krow)
    ref_det(:,2) = psi_det_beta_unique (:, kcol)
    
    ! check if determinant can accept electron
    ! iorb = elec_alpha_num + 1
    ! ispin = 1
    ! call orb_is_filled(tmp_det, iorb, ispin, N_int, filled)

    ! Finding (1,0) and (2,0) excitations
    ! loop over same beta different alpha
    n_buffer = 0
    do i = psi_bilinear_matrix_columns_loc(kcol), psi_bilinear_matrix_columns_loc(kcol+1)-1
        lidx = psi_bilinear_matrix_order(i)

        if (lidx <= kidx) cycle

        lcol = psi_bilinear_matrix_columns(i)
        lrow = psi_bilinear_matrix_rows(i)

        tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
        tmp_det(:,2) = psi_det_beta_unique (:, lcol)

        ! call orb_is_filled(tmp_det, iorb, ispin, N_int, filled)
        
        ! if (filled) cycle

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


    ! Finding (0,1) and (0,2) excitations 
    k_b = psi_bilinear_matrix_order_transp_reverse(k_a)
    krow = psi_bilinear_matrix_transp_rows(k_b) !this is unnecessary, technically
    kcol = psi_bilinear_matrix_transp_columns(k_b)
    
    ! check if determinant can accept electron
    ! iorb = elec_alpha_num + 1
    ! ispin = 1
    ! call orb_is_filled(tmp_det, iorb, ispin, N_int, filled)
    ! call print_det(tmp_det, N_int)

    ! loop over same alpha different beta
    n_buffer = 0
    do i = psi_bilinear_matrix_transp_rows_loc(krow), psi_bilinear_matrix_transp_rows_loc(krow+1)-1
        tidx = psi_bilinear_matrix_order_transp_reverse(i)
        lidx = psi_bilinear_matrix_order(tidx)

        if (lidx <= kidx) cycle

        lcol = psi_bilinear_matrix_transp_columns(i)
        lrow = psi_bilinear_matrix_transp_rows(i)
        
        tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
        tmp_det(:,2) = psi_det_beta_unique (:, lcol)
        
        ! call orb_is_filled(tmp_det, iorb, ispin, N_int, filled)
        
        ! if (filled) cycle

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


    ! Finding (1,1) excitations
    deallocate(buffer, idx)
    allocate(buffer(N_int, N_det), idx(N_det))

    ! need to iterate over single excitations over a small set of determinants
    ! and need to add new (1,1) determinants, but without transferring allocation
    ! after each new set

    ! create long term storage buffers that upper bound total ab excitations
    n_doubles_ab_max = 0
    do j = 1, n_singles_b
        k = singles_b(j)
        kcol = psi_bilinear_matrix_columns(k)
        n_doubles_ab_max += psi_bilinear_matrix_columns_loc(kcol+1)&
                        - psi_bilinear_matrix_columns_loc(kcol)
    end do

    n_doubles_ab_max -= n_singles_b ! these won't get included in buffers but are in pointer ranges
    
    ! n_doubles_ab_max will always be larger than the sum of n_buffer across all iterations
    ! n_doubles_ab_tot will always be smaller than n_buffer 
    ! in the end, only first n_doubles_ab_tot values should be needed

    allocate(doubles_ab(n_doubles_ab_max))

    n_buffer = 0
    n_buffer_old = 0
    n_doubles_ab_tot = 0

    ! start from (0,1) to excite to (1,1)
    do j = 1, n_singles_b

        k = singles_b(j)
        krow = psi_bilinear_matrix_rows(k)
        kcol = psi_bilinear_matrix_columns(k)

        tmp_det2(:,1) = psi_det_alpha_unique(:, krow)
        tmp_det2(:,2) = psi_det_beta_unique (:, kcol)

        ! loop over same beta different alpha
        do i = psi_bilinear_matrix_columns_loc(kcol), psi_bilinear_matrix_columns_loc(kcol+1)-1
            lidx = psi_bilinear_matrix_order(i)

            if (lidx <= kidx) cycle

            lcol = psi_bilinear_matrix_columns(i)
            lrow = psi_bilinear_matrix_rows(i)

            tmp_det(:,1) = psi_det_alpha_unique(:, lrow)
            tmp_det(:,2) = psi_det_beta_unique (:, lcol)

            ! call orb_is_filled(tmp_det, iorb, ispin, N_int, filled)
            
            ! if (filled) cycle

            ! add determinant to buffer
            ! buffer is list of alpha spin determinants
            n_buffer += 1
            buffer(:,n_buffer) = tmp_det(:,1)
            idx(n_buffer) = lidx
        end do

        sdet_a = tmp_det2(:,1)
        call get_all_spin_singles(buffer(:,n_buffer_old+1:n_buffer), idx(n_buffer_old+1:n_buffer),&
                                sdet_a, N_int, n_buffer-n_buffer_old,&
                                doubles_ab(n_doubles_ab_tot+1:n_doubles_ab_tot+n_buffer-n_buffer_old),&
                                n_doubles_ab)


        n_buffer_old = n_buffer
        n_doubles_ab_tot += n_doubles_ab
        
    end do

    ! add all indices into list and sort
    ! number of off-diagonal terms needed to caclulate to fill this row in sparse Hamlitonian
    n_off_diagonal = n_singles_a + n_singles_b + n_doubles_aa + n_doubles_ab + n_doubles_bb

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
    all_idx(n_offset+1:n_offset+n_doubles_ab)   = doubles_ab(:n_doubles_ab)
                            n_offset           += n_doubles_ab
    all_idx(n_offset+1:n_offset+n_doubles_bb)   = doubles_bb(:n_doubles_bb)
                            n_offset           += n_doubles_bb

    call insertion_isort(all_idx, srt_idx, n_off_diagonal)

    deallocate(buffer, idx, singles_a, doubles_aa,&
                singles_b, doubles_bb, doubles_ab,&
                srt_idx)

    columns = all_idx
    nnz = n_off_diagonal

end

! Utils
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