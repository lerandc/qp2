
BEGIN_PROVIDER [complex*16, FPS_SPF_Matrix_AO_complex, (AO_num, AO_num)]
  implicit none
  BEGIN_DOC
  !   Commutator FPS - SPF
  END_DOC
  complex*16, allocatable  :: scratch(:,:)
  allocate(                                                          &
      scratch(AO_num, AO_num)                                  &
      )

  ! Compute FP

  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
      (1.d0,0.d0),                                                          &
      Fock_Matrix_AO_complex,Size(Fock_Matrix_AO_complex,1),                &
      SCF_Density_Matrix_AO_complex,Size(SCF_Density_Matrix_AO_complex,1),  &
      (0.d0,0.d0),                                                          &
      scratch,Size(scratch,1))

  ! Compute FPS

  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
      (1.d0,0.d0),                                                   &
      scratch,Size(scratch,1),                                       &
      AO_Overlap_complex,Size(AO_Overlap_complex,1),                 &
      (0.d0,0.d0),                                                  &
      FPS_SPF_Matrix_AO_complex,Size(FPS_SPF_Matrix_AO_complex,1))

  ! Compute SP

  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
      (1.d0,0.d0),                                                          &
      AO_Overlap_complex,Size(AO_Overlap_complex,1),                        &
      SCF_Density_Matrix_AO_complex,Size(SCF_Density_Matrix_AO_complex,1),  &
      (0.d0,0.d0),                                                          &
      scratch,Size(scratch,1))

  ! Compute FPS - SPF

  call zgemm('N','N',AO_num,AO_num,AO_num,                           &
      (-1.d0,0.d0),                                                  &
      scratch,Size(scratch,1),                                       &
      Fock_Matrix_AO_complex,Size(Fock_Matrix_AO_complex,1),         &
      (1.d0,0.d0),                                                   &
      FPS_SPF_Matrix_AO_complex,Size(FPS_SPF_Matrix_AO_complex,1))

END_PROVIDER

BEGIN_PROVIDER [complex*16, FPS_SPF_Matrix_MO, (mo_num, mo_num)]
  implicit none
  begin_doc
!   Commutator FPS - SPF in MO basis
  end_doc
  call ao_to_mo_complex(FPS_SPF_Matrix_AO_complex, size(FPS_SPF_Matrix_AO_complex,1), &
     FPS_SPF_Matrix_MO_complex, size(FPS_SPF_Matrix_MO_complex,1))
END_PROVIDER


 BEGIN_PROVIDER [ double precision, eigenvalues_Fock_matrix_AO_complex, (AO_num) ]
&BEGIN_PROVIDER [ complex*16, eigenvectors_Fock_matrix_AO_complex, (AO_num,AO_num) ]
   !TODO: finish this provider; write provider for S_half_inv_complex
   BEGIN_DOC
   ! Eigenvalues and eigenvectors of the Fock matrix over the AO basis
   END_DOC

   implicit none

   double precision, allocatable  :: scratch(:,:),work(:),Xt(:,:)
   integer                        :: lwork,info
   integer                        :: i,j

   lwork = 3*AO_num - 1
   allocate(                                                         &
       scratch(AO_num,AO_num),                                 &
       work(lwork),                                                  &
       Xt(AO_num,AO_num)                                             &
       )

! Calculate Xt

  do i=1,AO_num
    do j=1,AO_num
      Xt(i,j) = S_half_inv(j,i)
    enddo
  enddo

! Calculate Fock matrix in orthogonal basis: F' = Xt.F.X

  call dgemm('N','N',AO_num,AO_num,AO_num,     &
       1.d0,                                   &
       Fock_matrix_AO,size(Fock_matrix_AO,1),  &
       S_half_inv,size(S_half_inv,1),        &
       0.d0,                                   &
       eigenvectors_Fock_matrix_AO,size(eigenvectors_Fock_matrix_AO,1))

  call dgemm('N','N',AO_num,AO_num,AO_num,                              &
       1.d0,                                                            &
       Xt,size(Xt,1),                                                   &
       eigenvectors_Fock_matrix_AO,size(eigenvectors_Fock_matrix_AO,1), &
       0.d0,                                                            &
       scratch,size(scratch,1))

! Diagonalize F' to obtain eigenvectors in orthogonal basis C' and eigenvalues

   call dsyev('V','U',AO_num,       &
        scratch,size(scratch,1),    &
        eigenvalues_Fock_matrix_AO, &
        work,lwork,info)

   if(info /= 0) then
     print *,  irp_here//' failed : ', info
     stop 1
   endif

! Back-transform eigenvectors: C =X.C'

  call dgemm('N','N',AO_num,AO_num,AO_num,     &
       1.d0,                                   &
       S_half_inv,size(S_half_inv,1),        &
       scratch,size(scratch,1),                &
       0.d0,                                   &
       eigenvectors_Fock_matrix_AO,size(eigenvectors_Fock_matrix_AO,1))

END_PROVIDER

