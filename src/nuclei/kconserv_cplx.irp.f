BEGIN_PROVIDER [integer, kpt_pair_num]
  implicit none
  kpt_pair_num = shiftr(kpt_num*kpt_num+kpt_num,1)
END_PROVIDER

BEGIN_PROVIDER [integer, kconserv, (kpt_num,kpt_num,kpt_num)]
  implicit none
  BEGIN_DOC
  ! Information about k-point symmetry
  !
  ! for k-points I,J,K: kconserv(I,J,K) gives L such that
  ! k_I + k_J = k_K + k_L
  ! two-electron integrals of the form <ij|kx>
  ! (where i,j,k have momentum k_I, k_J, k_K)
  ! will only be nonzero if x has momentum k_L (as described above)
  !
  END_DOC
  integer                        :: i,j,k,l

  if (read_kconserv) then
    call ezfio_get_nuclei_kconserv(kconserv)
    print *,  'kconserv read from disk'
  else
    call set_kconserv(kconserv)
    !print*,'kconserv must be provided'
    !stop -1
  endif
  if (write_kconserv) then
    call ezfio_set_nuclei_kconserv(kconserv)
    print *,  'kconserv written to disk'
  endif
END_PROVIDER

BEGIN_PROVIDER [integer, kpt_pair_map, (kpt_num,kpt_num)]
  implicit none
  BEGIN_DOC
  ! Information about k-point symmetry
  !
  ! for k-points I,K: kpt_pair_map(I,K) = \alpha
  ! where Q_{\alpha} = k_I - k_K
  ! 
  END_DOC

  if (read_kpt_symm) then
    call ezfio_get_nuclei_kpt_pair_map(kpt_pair_map)
    print *,  'kpt_pair_map read from disk'
  else
    print*,'kpt_pair_map must be provided'
    stop -1
  endif
  if (write_kpt_symm) then
    call ezfio_set_nuclei_kpt_pair_map(kpt_pair_map)
    print *,  'kpt_pair_map written to disk'
  endif
END_PROVIDER

BEGIN_PROVIDER [integer, kpt_inv, (kpt_num)]
  implicit none
  BEGIN_DOC
  ! Information about k-point symmetry
  !
  ! for k-point I: kpt_inv(I) = K
  ! where k_I + k_K = 0 (mod G)
  ! 
  END_DOC

  if (read_kpt_symm) then
    call ezfio_get_nuclei_kpt_inv(kpt_inv)
    print *,  'kpt_inv read from disk'
  else
    print*,'kpt_inv must be provided'
    stop -1
  endif
  if (write_kpt_symm) then
    call ezfio_set_nuclei_kpt_inv(kpt_inv)
    print *,  'kpt_inv written to disk'
  endif
END_PROVIDER

BEGIN_PROVIDER [integer, kpt_sparse_map, (kpt_num)]
  implicit none
  BEGIN_DOC
  ! Information about k-point symmetry
  !
  ! for k-point I: if kpt_sparse_map(I) = j
  ! if j>0: data for k_I is stored at index j in chol_ints
  ! if j<0: data for k_I is conj. transp. of data at index j in chol_{ao,mo}_integrals_complex
  ! 
  END_DOC

  if (read_kpt_symm) then
    call ezfio_get_nuclei_kpt_sparse_map(kpt_sparse_map)
    print *,  'kpt_sparse_map read from disk'
  else
    print*,'kpt_sparse_map must be provided'
    stop -1
  endif
  if (write_kpt_symm) then
    call ezfio_set_nuclei_kpt_sparse_map(kpt_sparse_map)
    print *,  'kpt_sparse_map written to disk'
  endif
END_PROVIDER

subroutine double_allowed_kpts(kh1,kh2,kp1,kp2,is_allowed)
  implicit none
  integer, intent(in) :: kh1,kh2,kp1,kp2
  logical, intent(out) :: is_allowed

  is_allowed = (kconserv(kh1,kh2,kp1) == kp2)
end subroutine

subroutine set_kconserv(kcon)
  implicit none
  integer, intent(out) :: kcon(kpt_num,kpt_num,kpt_num)
  integer :: i,j,k,qij

  do i=1,kpt_num
    do k=1,kpt_num
      ! Q = k_I - k_K
      qik = kpt_pair_map(i,k)
      do j=1,kpt_num
        ! k_L = k_J - (-(k_I - k_K))
        kcon(i,j,k) = kpt_pair_map(j,kpt_inv(qik))
      enddo
    enddo
  enddo
end subroutine
