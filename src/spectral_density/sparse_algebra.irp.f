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
    
    !$OMP PARALLEL SHARED(n, nnz_tot, nnz_arr, n_threads) &
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

    !$OMP SINGLE
    !$OMP SIMD REDUCTION(inscan, +:scn_a)
    do i = 1, n+1
        scn_a = scn_a + csr_s(i)
        !$OMP SCAN INCLUSIVE(scn_a)
        csr_s(i) = scn_a
    end do
    !$OMP END SINGLE
    
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


subroutine test_unique_looping(finished)
    implicit none

    logical :: finished
    integer :: i, j, k, k_a, k_b, l_a, l_b, maxab
    integer :: krow, kcol, lrow, lcol
    integer(bit_kind) :: tmp_det(N_int, 2), tmp_det2(N_int, 2), tmp_det3(N_int, 2), sdet_a(N_int), sdet_b(N_int)
    ! only need two tmp dets, change first to ref_det

    ! print *, k_a, krow, kcol

    ! call print_det(tmp_det, N_int)

    ! do i = 1, 20
    !     tmp_det(1:N_int, :) = psi_det_sorted_bit(:,:,i)
    !     print *, tmp_det(:,1), tmp_det(:,2)
    ! end do

    ! print *, "----------------------"

    ! do i = 1, 20
    !     sdet_a = psi_det_alpha_unique(:,i)
    !     sdet_b = psi_det_beta_unique(:,i)
    !     print *, sdet_a, sdet_b
    ! end do

    ! print *, "----------------------"

    ! do i = 1, 20
    !     print *, psi_bilinear_matrix_rows(i), psi_bilinear_matrix_columns(i), psi_bilinear_matrix_order_reverse(i), psi_bilinear_matrix_order(i)
    ! end do

    ! print *, "----------------------"
    ! do i = 1, 20
    !     krow = psi_bilinear_matrix_rows(i)
    !     kcol = psi_bilinear_matrix_columns(i)
    !     tmp_det(:,1) = psi_det_alpha_unique(1:N_int, krow)
    !     tmp_det(:,2) = psi_det_beta_unique (1:N_int, kcol)
    !     tmp_det2 = psi_det(:, :, psi_bilinear_matrix_order(i))
    !     print *, i, krow, kcol
    !     print *, tmp_det(:,1), tmp_det2(:,1), tmp_det(:,2), tmp_det2(:,2)
    ! end do

    ! print *, "----------------------"

    ! do i = 1, 20
    !     print *, psi_bilinear_matrix_transp_rows(i), psi_bilinear_matrix_transp_columns(i), psi_bilinear_matrix_order_transp_reverse(i), psi_bilinear_matrix_transp_order(i)
    ! end do

    ! print *, "----------------------"
    ! do i = 1, 20
    !     k_b = psi_bilinear_matrix_order_transp_reverse(i)
    !     krow = psi_bilinear_matrix_transp_rows(k_b)
    !     kcol = psi_bilinear_matrix_transp_columns(k_b)
    !     tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int, krow)
    !     tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int, kcol)
    !     print *, i
    !     call print_det(tmp_det, N_int)
    !     ! tmp_det2 = psi_det(:, :, psi_bilinear_matrix_order(i))
    !     ! write(*, '(I8, I8, I8)'), i, k_b, psi_bilinear_matrix_order_transp_reverse(k_b)
    !     ! write(*, '(I12, I12, I12, I12)'), tmp_det(:,1), tmp_det2(:,1), tmp_det(:,2), tmp_det2(:,2)
    ! end do

    ! print *, "----------------------"
    ! do i = 1, maxab
    !     print *, psi_bilinear_matrix_columns_loc(i), psi_bilinear_matrix_transp_rows_loc(i)
    ! end do

    ! search in tranp (alpha major) for beta single/double excitations
    ! only add dets to buffer if index in determinant space > index of reference determinant

    integer :: lidx, kidx, tidx, ispin, iorb, n_buffer
    integer :: n_singles_a, n_singles_b, n_doubles_aa, n_doubles_bb, n_doubles_ab
    integer(bit_kind), allocatable  :: buffer(:,:)
    integer,   allocatable          :: idx(:), singles_a(:), singles_b(:), doubles_aa(:), doubles_bb(:), doubles_ab(:)
    maxab = max(N_det_alpha_unique, N_det_beta_unique) + 1
    
    allocate(buffer(N_int, N_det), idx(N_det))
    
    logical :: filled
    
    k_a = 5
    kidx = psi_bilinear_matrix_order(k_a)

    krow = psi_bilinear_matrix_rows(k_a)
    kcol = psi_bilinear_matrix_columns(k_a)
    tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int, krow)
    tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int, kcol)
    
    ! check if determinant can accept electron
    iorb = elec_alpha_num + 1
    ispin = 1
    call orb_is_filled(tmp_det, iorb, ispin, N_int, filled)
    print *, "kidx: ", kidx
    print *, "Filled?", filled
    call print_det(tmp_det, N_int)


    
    print *, "--------- Finding (1,0) and (2,0) excitations ----------"
    ! loop over same beta different alpha
    n_buffer = 0
    do i = psi_bilinear_matrix_columns_loc(kcol), psi_bilinear_matrix_columns_loc(kcol+1)-1
        lidx = psi_bilinear_matrix_order(i)

        print *, i - psi_bilinear_matrix_columns_loc(kcol)

        if (lidx <= kidx) then
            print *, i - psi_bilinear_matrix_columns_loc(kcol), "idx bound"
            cycle ! only work on upper triangle
        end if

        lcol = psi_bilinear_matrix_columns(i)
        lrow = psi_bilinear_matrix_rows(i)

        tmp_det2(1:N_int,1) = psi_det_alpha_unique(1:N_int, lrow)
        tmp_det2(1:N_int,2) = psi_det_beta_unique (1:N_int, lcol)

        call orb_is_filled(tmp_det2, iorb, ispin, N_int, filled)
        
        if (filled) then
            print *, i - psi_bilinear_matrix_columns_loc(kcol), "can't accept electron"
            cycle ! new determinant cannot accept electron
        end if

        ! add determinant to buffer
        ! buffer is list of alpha spin determinants
        n_buffer += 1
        call print_det(tmp_det2, N_int)
        buffer(:,n_buffer) = tmp_det2(:,1)
        idx(n_buffer) = lidx

    end do

    print *, "n_buffer", n_buffer
    allocate(singles_a(n_buffer), doubles_aa(n_buffer))
    sdet_a = tmp_det(:,1)
    call get_all_spin_singles_and_doubles(buffer, idx, sdet_a, &
                        N_int, n_buffer, singles_a, doubles_aa, n_singles_a, n_doubles_aa)

    print *, n_singles_a, n_doubles_aa
    print *, "---------------"
    do i = 1, n_buffer
        write(*, '(I10, I10, I10, I10)'), i, idx(i), singles_a(i), doubles_aa(i)
    end do

    deallocate(buffer, idx)

    allocate(buffer(N_int, N_det), idx(N_det))


    print *, "--------- Finding (0,1) and (0,2) excitations ----------"
    k_b = psi_bilinear_matrix_order_transp_reverse(k_a)
    krow = psi_bilinear_matrix_transp_rows(k_b)
    kcol = psi_bilinear_matrix_transp_columns(k_b)
    tmp_det(1:N_int,1) = psi_det_alpha_unique(1:N_int, krow)
    tmp_det(1:N_int,2) = psi_det_beta_unique (1:N_int, kcol)
    
    ! check if determinant can accept electron
    iorb = elec_alpha_num + 1
    ispin = 1
    call orb_is_filled(tmp_det, iorb, ispin, N_int, filled)
    print *, "Filled?", filled
    call print_det(tmp_det, N_int)

    ! loop over same alpha different beta
    n_buffer = 0
    do i = psi_bilinear_matrix_transp_rows_loc(krow), psi_bilinear_matrix_transp_rows_loc(krow+1)-1
        tidx = psi_bilinear_matrix_order_transp_reverse(i)
        lidx = psi_bilinear_matrix_order(tidx)

        print *, i - psi_bilinear_matrix_transp_rows_loc(krow)
        
        if (lidx <= kidx) then
            print *, i - psi_bilinear_matrix_transp_rows_loc(krow), "idx bound"
            cycle ! only work on upper triangle
        end if
        lcol = psi_bilinear_matrix_transp_columns(i)
        lrow = psi_bilinear_matrix_transp_rows(i)
        
        tmp_det2(1:N_int,1) = psi_det_alpha_unique(1:N_int, lrow)
        tmp_det2(1:N_int,2) = psi_det_beta_unique (1:N_int, lcol)
        
        call orb_is_filled(tmp_det2, iorb, ispin, N_int, filled)
        
        if (filled) then
            print *, i - psi_bilinear_matrix_transp_rows_loc(krow), "can't accept electron"
            cycle ! new determinant cannot accept electron
        end if
        ! add determinant to buffer
        ! buffer is list of beta spin determinants
        n_buffer += 1
        call print_det(tmp_det2, N_int)
        buffer(:,n_buffer) = tmp_det2(:,2)
        idx(n_buffer) = lidx
    end do
    
    print *, "n_buffer", n_buffer
    allocate(singles_b(n_buffer), doubles_bb(n_buffer))
    sdet_b = tmp_det(:,2)
    call get_all_spin_singles_and_doubles(buffer, idx, sdet_b, &
                        N_int, n_buffer, singles_b, doubles_bb, n_singles_b, n_doubles_bb)

    print *, n_singles_b, n_doubles_bb
    print *, "---------------"
    do i = 1, n_buffer
        write(*, '(I10, I10, I10, I10)'), i, idx(i), singles_b(i), doubles_bb(i)
    end do

    ! by now, have acquired (0,1), (0,2), (1,0), (2,0) excitations
    ! need just (1,1)
    ! can loop either over single a or single b to get doubles ab
    ! do I need another index or do I have enough info already?

    print *, "--------- Finding (1,1) excitations ----------"
    do j = 1, n_singles_b
        deallocate(buffer, idx)
        allocate(buffer(N_int, N_det), idx(N_det))

        k = singles_b(j)
        krow = psi_bilinear_matrix_rows(k)
        kcol = psi_bilinear_matrix_columns(k)

        tmp_det3(1:N_int,1) = psi_det_alpha_unique(1:N_int, krow)
        tmp_det3(1:N_int,2) = psi_det_beta_unique (1:N_int, kcol)

        print *, "--- Local ref det"

        call print_det(tmp_det3, N_int)
        ! loop over same beta different alpha
        n_buffer = 0
        do i = psi_bilinear_matrix_columns_loc(kcol), psi_bilinear_matrix_columns_loc(kcol+1)-1
            lidx = psi_bilinear_matrix_order(i)

            print *, i - psi_bilinear_matrix_columns_loc(kcol)

            if (lidx <= kidx) then
                print *, i - psi_bilinear_matrix_columns_loc(kcol), "idx bound"
                cycle ! only work on upper triangle
            end if

            lcol = psi_bilinear_matrix_columns(i)
            lrow = psi_bilinear_matrix_rows(i)

            tmp_det2(1:N_int,1) = psi_det_alpha_unique(1:N_int, lrow)
            tmp_det2(1:N_int,2) = psi_det_beta_unique (1:N_int, lcol)

            call orb_is_filled(tmp_det2, iorb, ispin, N_int, filled)
            
            if (filled) then
                print *, i - psi_bilinear_matrix_columns_loc(kcol), "can't accept electron"
                cycle ! new determinant cannot accept electron
            end if

            ! add determinant to buffer
            ! buffer is list of alpha spin determinants
            n_buffer += 1
            call print_det(tmp_det2, N_int)
            buffer(:,n_buffer) = tmp_det2(:,1)
            idx(n_buffer) = lidx

        end do

        
        print *, "n_buffer", n_buffer
        allocate(doubles_ab(n_buffer))
        sdet_a = tmp_det3(:,1)
        call get_all_spin_singleS(buffer, idx, sdet_a, &
                            N_int, n_buffer, doubles_ab, n_doubles_ab)
    
        print *, n_doubles_ab
        print *, "---------------"
        do i = 1, n_buffer
            write(*, '(I10, I10, I10, I10)'), i, idx(i), doubles_ab(i)
        end do

        deallocate(doubles_ab)

    end do

    ! now, need to just add all indices into list; sort; remove duplicates
    ! though if number of duplicates is small, then might be better to just
    ! check in calculation loop to see if previous j equal current j and cycle

    finished = .TRUE.
end

! KEEP IN MIND: how to sort by j? is that necessary or worth it? probably better locality in sparse MV to have it ordered in columns
! might need to manually sort it

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