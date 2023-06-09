logical function det_equal(alpha_a, beta_a, alpha_b, beta_b)
    implicit none

    integer(bit_kind), intent(in) :: alpha_a(N_int), beta_a(N_int), alpha_b(N_int), beta_b(N_int)
    integer :: i

    do i = 1, N_int
        if ((alpha_a(i) /= alpha_b(i)) .or. (beta_a(i) /= beta_b(i))) then
            det_equal = .false.
            return
        end if
    end do

    det_equal = .true.
    return

end function det_equal

integer function hash_index(det_alpha, det_beta, hash_prime, ht_size, n_orbs)
    implicit none
    BEGIN_DOC
    ! Hash function for table of excited determinants
    END_DOC

    integer, intent(in) :: hash_prime, ht_size, n_orbs
    integer(bit_kind), intent(in) :: det_alpha(N_int), det_beta(N_int)
    
    !! need to adapt modulo to work with integers represented as multiple ints,  ie., bit_kind > 1
    hash_index = modulo(det_alpha(1), hash_prime) + (det_beta(1) / (ht_size / n_orbs))
    hash_index = modulo(hash_index, ht_size)
    return

end function hash_index

integer function hash_value( hash_alpha, hash_beta, hash_vals, hash_prime, ht_size, n_orbs, det_alpha, det_beta)
    implicit none
    BEGIN_DOC
    ! Hash function for table of excited determinants
    END_DOC
    
    integer, intent(in) :: det_alpha, det_beta, hash_prime, ht_size, n_orbs
    integer, intent(inout) :: hash_vals(ht_size)
    integer(bit_kind), intent(inout) :: hash_alpha(N_int, ht_size), hash_beta(N_int, ht_size)
    integer :: idx, hash_index
    logical :: det_equal

    idx = hash_index(det_alpha, det_beta, hash_prime, ht_size, n_orbs)
    
    ! search for value with linear probing
    ! rely on table not being full to avoid infinite loop when key doesn't exist
    ! TODO: hash_alpha(1, idx) can be zero if the none of the electrons are in the first integer
    do while ( (hash_alpha(1, idx) > 0))
        if ( det_equal(hash_alpha(:,idx), hash_beta(:,idx), det_alpha, det_beta) ) then
            hash_value = hash_vals(idx)
            return
        end if
        idx += 1
        idx = modulo(idx, ht_size)
    end do
    
    hash_value = -1
    return
    
end function hash_value

subroutine insert_key_val_pair(hash_alpha, hash_beta, hash_vals, hash_prime, ht_size, n_orbs, det_alpha, det_beta, target_row, success)
    implicit none
    BEGIN_DOC
    ! Hash function for table of excited determinants
    END_DOC
    
    integer, intent(in)    :: ht_size, target_row, hash_prime, n_orbs
    logical, intent(out)   :: success
    integer, intent(inout) :: hash_vals(ht_size)
    integer(bit_kind), intent(in) :: det_alpha(N_int), det_beta(N_int)
    integer(bit_kind), intent(inout) :: hash_alpha(N_int, ht_size), hash_beta(N_int, ht_size)
    integer :: idx, hash_index
    logical :: det_equal
    
    idx = hash_index(det_alpha, det_beta, hash_prime, ht_size, n_orbs)
    
    ! search for value with linear probing
    ! rely on table not being full to avoid infinite loop when key doesn't exist
    do while ((hash_alpha(1, idx) > 0))
        if ( det_equal(hash_alpha(:,idx), hash_beta(:,idx), det_alpha, det_beta) ) then
            ! determinant exists already in the table
            success = .false.
            return
        end if
        idx += 1
        idx = modulo(idx, ht_size)
    end do

    hash_alpha(:, idx) = det_alpha
    hash_beta(:, idx) = det_beta
    hash_vals(idx) = target_row

    success = .true.
end

subroutine build_hash_table(hash_alpha, hash_beta, hash_vals, hash_prime, ht_size, I_cut, I_det, n_orb, n_det_l, det_basis, n_det_out)
    implicit none
    BEGIN_DOC
    ! Construct the hash table
    END_DOC

    integer, intent(in)     :: ht_size, hash_prime, n_orb, n_det_l, I_cut(n_det_l, n_orb)
    integer(bit_kind), intent(in) :: I_det(N_int, 2, n_det_l, n_orb)
    integer, intent(inout)  :: n_det_out
    integer, intent(out)    :: hash_vals(ht_size)
    integer(bit_kind), intent(out) :: hash_alpha(N_int, ht_size), hash_beta(N_int, ht_size)
    integer(bit_kind), intent(out) :: det_basis(N_int, 2, ht_size)
    logical                 :: hash_success

    integer :: i, iorb

    n_det_out = 0

    ! iterate over each set of determinants in excitation order
    do iorb = 1, n_orb
        do i = 1, n_det_l
            if (I_cut(i, iorb) == 1) then

                call insert_key_val_pair(hash_alpha, hash_beta, hash_vals, hash_prime, ht_size, n_orb, &
                                         I_det(:,1,i,iorb), I_det(:,2,i,iorb), n_det_out, hash_success)

                if (hash_success) then
                    n_det_out += 1
                    det_basis(:, :, n_det_out) = I_det(:,:,i,iorb)
                end if

            end if
        end do
    end do


end