BEGIN_PROVIDER [complex*16, ao_ortho_lowdin_coef_complex, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
! matrix of the coefficients of the mos generated by the
! orthonormalization by the S^{-1/2} canonical transformation of the aos
! ao_ortho_lowdin_coef(i,j) = coefficient of the ith ao on the jth ao_ortho_lowdin orbital
  END_DOC
  integer                        :: i,j,k,l
  complex*16, allocatable  :: tmp_matrix(:,:)
  allocate (tmp_matrix(ao_num,ao_num))
  tmp_matrix(:,:) = (0.d0,0.d0)
  do j=1, ao_num
    tmp_matrix(j,j) = (1.d0,0.d0)
  enddo
  call ortho_lowdin_complex(ao_overlap_complex,ao_num,ao_num,tmp_matrix,ao_num,ao_num)
  do i=1, ao_num
    do j=1, ao_num
      ao_ortho_lowdin_coef_complex(j,i) = tmp_matrix(i,j)
    enddo
  enddo
  deallocate(tmp_matrix)
END_PROVIDER

BEGIN_PROVIDER [complex*16, ao_ortho_lowdin_overlap_complex, (ao_num,ao_num)]
  implicit none
  BEGIN_DOC
! overlap matrix of the ao_ortho_lowdin
! supposed to be the Identity
  END_DOC
  integer                        :: i,j,k,l
  complex*16               :: c
  do j=1, ao_num
    do i=1, ao_num
      ao_ortho_lowdin_overlap_complex(i,j) = (0.d0,0.d0)
    enddo
  enddo
  do k=1, ao_num
    do j=1, ao_num
      c = (0.d0,0.d0)
      do l=1, ao_num
        c +=  dconjg(ao_ortho_lowdin_coef_complex(j,l)) * ao_overlap_complex(k,l)
      enddo
      do i=1, ao_num
        ao_ortho_lowdin_overlap_complex(i,j) += ao_ortho_lowdin_coef_complex(i,k) * c
      enddo
    enddo
  enddo
END_PROVIDER
