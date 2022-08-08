! modules for calculating relevant finite size corrections for 
! excited state properties calculated for solids



! following along Yang, et. al, Phys. Rev. B 101, 085115 (2020).
! TODOs
! 1) calculate the dielectric constant
! needs:
!   a) fourier transforms of electron density -> through one_e_dm_mo_kpts ?
!   b) calculation of structure factors -> straight forward
!   c) extraploation of structure factor bound to low k -> fit low k with QR factorization
!
! 2) calculate the Madelung constant
! needs:
!   a) set up of k-space grid
!   b) integration of columb singularities over k-space grid
! maybe this could be imported from pyscf?


! BEGIN_PROVIDER [complex*16, rho_k, (kpt_num)]
!     implicit none

!     integer :: i, j, k
    
!     ! I think this should be an upper triangle sum of the one_e_dm_mo_kpts matrices
!     rho_k = (0.d0, 0.d0)
!     do k = 1, kpt_num
!         do i = 1, mo_num_per_kpt
!             do j = i, mo_num_per_kpt
!                 rho_k(k) += one_e_dm_mo_kpts(i,j,k)
!             end do
!         end do
!     end do

! END_PROVIDER


! BEGIN_PROVIDER [double precision, structure_factors, (kpt_num)]
!     implicit none
!     ! electron density is a real valued function
!     ! so its fourier transform has conjugate inversion symmetry
!     integer :: e_n 
!     e_n = (elec_alpha_num + elec_beta_num)/kpt_num
!     structure_factors = real(rho_k * dconjg(rho_k), kind=16) / (e_n)
! END_PROVIDER


! BEGIN_PROVIDER [double precision, dielectric_eps]
!     ! hopefully we can use MKL routines here
!     implicit none
!     print *, "building dielectric"

!     ! number of kpts probably on the lower end?
!     ! should also adjust this to only take into account the smallest N kpts
!     double precision :: pi, n_e, k2(kpt_num-1), gamma_k(kpt_num-1)
!     double precision :: b(kpt_num-1,1), A(kpt_num-1, 2), work(16)
!     integer :: i, e_n, info

!     pi = acos(-1.d0)
!     e_n = (elec_alpha_num + elec_beta_num)/kpt_num

!     ! set up gamma array
!     ! don't include gamma point
!     do i = 2, kpt_num
!         print *, i
!         k2(i) = kpts(i,1)**2.0 + kpts(i,2)**2.0 + kpts(i,3)**2.0
!         gamma_k(i) = 8.d0 * pi* e_n * structure_factors(i) / k2(i)
!     end do

!     ! set up vector b of 1 - gamma_k ^2 values
!     b(:,1) = 1 - gamma_k **2.0 

!     print *, structure_factors
!     print *, "---------------"
!     print *, b

!     ! set up matrix m x n matrix  A of coefficients for linear regression
!     A(:,1) = k2
!     A(:,2) = 1.d0 

!     ! calculate x = R^-1 Q^T b
!     call dgels('N',kpt_num-1, 2, 1,&
!                 A, kpt_num-1, b, kpt_num-1, &
!                 work, 16, info) ! since this is such a small least squares problem LWORK can be small

!     if(info == 0) then
!         dielectric_eps = b(2,1) !  m >= n, solution to linear system overwritten in first n rows of column
!     else
!         print *, info, " DGELS has failed"
!         dielectric_eps = 1
!     end if

! END_PROVIDER

subroutine calc_dielectric_eps(d_eps)

    implicit none
    print *, "building dielectric"

    ! number of kpts probably on the lower end?
    ! should also adjust this to only take into account the smallest N kpts
    double precision :: pi, n_e, k2(kpt_num-1), gamma_k(kpt_num-1)
    double precision :: b(kpt_num-1,1), A(kpt_num-1, 2), work(16)
    double precision, intent(out) :: d_eps
    complex*16  :: rho_k(kpt_num)
    double precision :: structure_factors(kpt_num)
    integer :: i,j,k, e_n, info

    
    e_n = (elec_alpha_num + elec_beta_num)/kpt_num

    rho_k = (0.d0, 0.d0)
    do k = 1, kpt_num
        do i = 1, mo_num_per_kpt
            do j = i, mo_num_per_kpt
                rho_k(k) += one_e_dm_mo_kpts(i,j,k)
            end do
        end do
    end do

    e_n = (elec_alpha_num + elec_beta_num)/kpt_num
    structure_factors = real(rho_k * dconjg(rho_k), kind=16) / (e_n)

    pi = acos(-1.d0)
    
    ! set up gamma array
    ! don't include gamma point
    do i = 2, kpt_num
        k2(i) = kpts(i,1)**2.0 + kpts(i,2)**2.0 + kpts(i,3)**2.0
        gamma_k(i) = 8.d0 * pi* e_n * structure_factors(i) / k2(i)
    end do

    ! set up vector b of 1 - gamma_k ^2 values
    b(:,1) = 1 - gamma_k **2.0 

    print *, rho_k(2:)
    print *, "---------------"
    print *, structure_factors(2:)
    print *, "---------------"
    print *, gamma_k
    print *, "---------------"
    print *, b

    ! set up matrix m x n matrix  A of coefficients for linear regression
    A(:,1) = k2
    A(:,2) = 1.d0 

    ! calculate x = R^-1 Q^T b
    call dgels('N',kpt_num-1, 2, 1,&
                A, kpt_num-1, b, kpt_num-1, &
                work, 16, info) ! since this is such a small least squares problem LWORK can be small

    if(info == 0) then
        d_eps = b(2,1) !  m >= n, solution to linear system overwritten in first n rows of column
    else
        print *, info, " DGELS has failed"
        d_eps = 1
    end if

end

