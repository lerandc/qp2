integer function hash_index(det_alpha, det_beta, hash_prime, ht_size, n_orbs)
    implicit none
    BEGIN_DOC
    ! Hash function for table of excited determinants
    END_DOC

    integer, intent(in) :: det_alpha, det_beta, hash_prime, ht_size, n_orbs
    
    !! need to adapt modulo to work with integers represented as multiple ints,  ie., bit_kind > 1
    !! need to also make sure it works for bit kind integers
    hash_index = modulo(det_alpha, hash_prime) + (det_beta / (ht_size / n_orbs))
    hash_index = modulo(hash_index, ht_size)
    return

end function hash_index

integer function hash_value( hash_alpha, hash_beta, hash_vals, hash_prime, ht_size, det_alpha, det_beta)
    implicit none
    BEGIN_DOC
    ! Hash function for table of excited determinants
    END_DOC
    
    integer, intent(in) :: det_alpha, det_beta, hash_prime, ht_size
    integer, intent(inout) :: hash_vals(ht_size)
    integer(bit_kind), intent(inout) :: hash_alpha(N_int, ht_size), hash_beta(N_int, ht_size)
    integer:: idx, hash_index

    idx = hash_index( det_alpha, det_beta, hash_prime)
    
    ! search for value with linear probing
    ! rely on table not being full to avoid infinite loop when key doesn't exist
    do while ( (hash_alpha(idx) > 0))
        if ((hash_alpha(idx) == det_alpha) .and. (hash_beta(idx) == det_beta)) then
            hash_value = hash_vals(idx)
            return
        end if
        idx += 1
        idx = modulo(idx, ht_size)
    end do
    
    hash_value = -1
    return
    
end function hash_value

subroutine insert_key_val_pair(hash_alpha ,hash_beta,hash_vals,hash_prime,ht_size,det_alpha,det_beta,target_row,success)
    implicit none
    BEGIN_DOC
    ! Hash function for table of excited determinants
    END_DOC
    
    integer, intent(in)    :: ht_size, det_alpha, det_beta, target_row, hash_prime
    logical, intent(out)   :: success
    integer, intent(inout) :: hash_vals(ht_size)
    integer(bit_kind), intent(inout) :: hash_alpha(N_int, ht_size), hash_beta(N_int, ht_size)
    integer                :: idx, hash_index
    
    idx = hash_index(det_alpha, det_beta, hash_prime)
    
    ! search for value with linear probing
    ! rely on table not being full to avoid infinite loop when key doesn't exist
    do while ((hash_alpha(idx) > 0))
        if ((hash_alpha(idx) == det_alpha) .and. (hash_beta(idx) == det_beta)) then
            ! determinant exists already in the table
            success = .false.
            return
        end if
        idx += 1
        idx = modulo(idx, ht_size)
    end do

    hash_alpha(idx) = det_alpha
    hash_beta(idx) = det_beta
    hash_vals(idx) = target_row

    success = .true.
end

subroutine build_hash_table( hash_alpha, hash_beta, hash_vals, hash_prime, ht_size, I_cut, I_det, n_orb, n_det_l, det_basis, n_det_out)
    implicit none
    BEGIN_DOC
    ! Construct the hash table
    END_DOC

    integer, intent(in)     :: ht_size, hash_prime, n_orb, n_det_l, I_cut(n_det_l, n_orb),
    integer(bit_kind), intent(in) :: I_det(N_int, 2, n_det_l, n_orb)
    integer, intent(out)    :: n_det_out, hash_vals(ht_size)
    integer(bit_kind), intent(out) :: hash_alpha(N_int, ht_size), hash_beta(N_int, ht_size)
    integer(bit_kind), intent(out) :: det_basis(N_int, 2, n_det_l*n_orb)
    logical                 :: hash_success


    n_det_out = 0

    ! iterate over each set of determinants in excitation order
    do iorb = 1, n_orb
        do i = 1, n_det_l
            if (I_cut(i, iorb) == 1) then
                
                call insert_key_val_pair(hash_alpha, hash_beta, hash_vals, hash_prime, ht_size, &
                                         I_det(:,1,i,iorb), I_det_k(:,2,i,iorb), n_det_out, hash_success)

                if (hash_success) then
                    n_det_out += 1
                    det_basis(:, :, n_det_out) = I_det(:,:,i,iorb)
                end if

            end if
    end do


end