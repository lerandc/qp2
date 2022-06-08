!
!
! Complex implementations
!
!

subroutine lanczos_tridiag_c(H, u0, alpha, beta, k, sze)
    implicit none
    BEGIN_DOC
    ! Function that performs lanczos triadiaglonization of a hermitian matrix
    ! Takes H matrix and u initial vector as inputs
    ! Outputs triadiagonlized Hermitian matrix as alpha and beta vectors
    ! Outputs up to k elements of tridiagonlization
    END_DOC

    integer, intent(in)             :: k, sze
    integer                         :: i, incx, incy
    complex*16, intent(in)          :: H(sze, sze), u0(sze)
    complex*16, allocatable         :: z(:), u(:), v(:) 
    complex*16, intent(out)         :: alpha(k)
    complex*16                      :: zdotc
    double precision, intent(out)   :: beta(k)
    double precision                :: dznrm2
    
    incx = 1 
    incy = 1

    allocate (z(sze), u(sze), v(sze))

    ! initialize vectors
    u = u0
    v = (0.d0, 0.d0)
    alpha = (0.d0, 0.d0)
    beta = 0.d0
    
    do i = 1, k
        z = 0
        call zhemv('U', sze, (1.d0,0.d0), H, sze, u, incx, (1.d0,0.d0), z, incy)
        alpha(i) = zdotc(sze, z, incx, u, incy)
        
        if (i == k) then
            exit
        end if

        z = z - alpha(i) * u - beta(i)*v
        
        beta(i+1) = dznrm2(sze, z, incx)
        if (beta(i+1) < 1e-16) then ! some small number
        ! add some type of escape or warning for beta(i) = 0
            exit
        end if

        v = u
        u = z / beta(i+1)
    enddo

    deallocate(z, u, v)

end

subroutine lanczos_tridiag_reortho_c(H, u0, uu, alpha, beta, k, sze)
    implicit none
    BEGIN_DOC
    ! Function that performs lanczos triadiaglonization of a hermitian matrix
    ! Performs reorthogonalization of intermediate vectors
    ! Outputs triadiagonlized Hermitian matrix as alpha and beta vectors
    ! Outputs up to k elements of tridiagonlization
    END_DOC

    integer, intent(in)             :: k, sze
    integer                         :: i, ii, j, incx, incy
    complex*16, intent(in)          :: H(sze, sze), u0(sze)
    complex*16                      :: z(sze)
    complex*16, intent(out)         :: alpha(k), uu(sze,k)
    complex*16                      :: zdotc
    double precision, intent(out)   :: beta(k)
    double precision                :: dznrm2, coef
    
    incx = 1 
    incy = 1

    ! allocate (z(sze), uu(sze, k))

    ! initialize vectors
    uu(:,1) = u0
    alpha = (0.d0, 0.d0)
    beta = 0.d0
    
    do i = 1, k
        z = (0.d0, 0.d0)
        call zhemv('U', sze, (1.d0,0.d0), H, sze, uu(:,i), incx, (0.d0,0.d0), z, incy)
        alpha(i) = zdotc(sze, uu(:,i), incx, z, incy)

        if (i == k) then
            exit
        end if

        do ii = 1, 2 ! repeat process twice
            do j = 1, i
                ! TODO: check to make sure that this loop can be done in this way, i.e., that the sum need not be completed before subtracting from z
                ! in exact arithmetic, it is equivalent; does order matter for roundoff?
                coef = zdotc(sze, z, incx, uu(:,j), incy)
                z = z - coef * uu(:,j)
            enddo
        enddo 

        beta(i+1) = dznrm2(sze, z, incx)
        if (beta(i+1) < 1e-16) then ! some small number
        ! add some type of escape or warning for beta(i) = 0
            exit
        end if

        uu(:,i+1) = z / beta(i+1)
    enddo

    ! deallocate(z, uu)

end

!
!
! Real implementations
!
!


subroutine lanczos_tridiag_r(H, u0, uu, alpha, beta, k, sze)
    implicit none
    BEGIN_DOC
    ! Function that performs lanczos triadiaglonization of a hermitian matrix
    ! Takes H matrix and u initial vector as inputs
    ! Outputs triadiagonlized Hermitian matrix as alpha and beta vectors
    ! Outputs up to k elements of tridiagonlization
    END_DOC

    integer, intent(in)             :: k, sze
    integer                         :: i, incx, incy
    double precision, intent(in)    :: H(sze, sze), u0(sze)
    double precision, allocatable   :: z(:), u(:), v(:) 
    double precision, intent(out)   :: alpha(k), beta(k), uu(sze,k)
    double precision                :: ddot
    double precision                :: dnrm2
    
    incx = 1 
    incy = 1

    allocate (z(sze), u(sze), v(sze))

    ! initialize vectors
    u = u0
    uu(:,1) = u0
    v = 0.d0
    alpha = 0.d0
    beta = 0.d0
    
    do i = 1, k
        z = 0
        call dsymv('U', sze, 1.d0, H, sze, u, incx, 0.d0, z, incy)
        alpha(i) = ddot(sze, z, incx, u, incy)
        
        if (i == k) then
            exit
        end if

        z = z - alpha(i) * u - beta(i)*v
        
        beta(i+1) = dnrm2(sze, z, incx)
        if (beta(i+1) < 1e-16) then ! some small number
        ! add some type of escape or warning for beta(i) = 0
            exit
        end if

        v = u
        u = z / beta(i+1)
        uu(:,i+1) = u
    enddo

    deallocate(z, u, v)

end

subroutine lanczos_tridiag_reortho_r(H, u0, uu, alpha, beta, k, sze)
    implicit none
    BEGIN_DOC
    ! Function that performs lanczos triadiaglonization of a hermitian matrix
    ! Performs reorthogonalization of intermediate vectors
    ! Outputs triadiagonlized Hermitian matrix as alpha and beta vectors
    ! Outputs up to k elements of tridiagonlization
    END_DOC

    integer, intent(in)             :: k, sze
    integer                         :: i, ii, j, incx, incy
    double precision, intent(in)    :: H(sze, sze), u0(sze)
    double precision                :: z(sze)
    double precision, intent(out)   :: alpha(k), beta(k), uu(sze,k)
    double precision                :: ddot, coef
    double precision                :: dnrm2
    
    incx = 1 
    incy = 1

    ! allocate (z(sze))

    ! initialize vectors
    uu(:,1) = u0
    alpha = 0.d0
    beta = 0.d0
    
    do i = 1, k
        z = 0
        call dsymv('U', sze, 1.d0, H, sze, uu(:,i), incx, 0.d0, z, incy)
        alpha(i) = ddot(sze, z, incx, uu(:,i), incy)

        if (i == k) then
            exit
        end if

        do ii = 1, 2 ! repeat process twice
            do j = 1, i
                ! TODO: check to make sure that this loop can be done in this way, i.e., that the sum need not be completed before subtracting from z
                ! in exact arithmetic, it is equivalent; does order matter for roundoff?
                coef = ddot(sze, z, incx, uu(:,j), incy)
                z = z - coef * uu(:,j)
            enddo
        enddo 
        
        beta(i+1) = dnrm2(sze, z, incx)
        if (beta(i+1) < 1e-16) then ! some small number
        ! add some type of escape or warning for beta(i) = 0
            print *, "escaping early"
            exit
        end if

        uu(:,i+1) = z / beta(i+1)
    enddo

    ! deallocate(z)

end


! subroutine reorthogonalize_vector_c(z, q, sze, j)
!     implicit none

!     integer, intent(in)         :: j, sze
!     integer                     :: i, incx, incy
!     complex*16, intent(in)      :: q
!     complex*16, intent(inout)   :: z
!     complex*16                  :: zdotu
    
!     incx = 1 
!     incy = 1

!     do ii = 1, 2 ! repeat process twice
!         do i = 2, j
!             ! TODO: check to make sure that this loop can be done in this way, i.e., that the sum need not be completed before subtracting from z
!             ! in exact arithmetic, it is equivalent; does order matter for roundoff?
!             coef = zdotu(sze, z, incx, q(:, i), incy)
!             z = z - coef * q(:,i)
!         enddo
!     enddo 

! end

! subroutine reorthogonalize_vector_r(z, q, sze, j)
!     implicit none

!     integer, intent(in)               :: j, sze
!     integer                           :: i, incx, incy
!     double precision, intent(in)      :: q
!     double precision, intent(inout)   :: z
!     double precision,                 :: ddot
    
!     incx = 1 
!     incy = 1

!     do ii = 1, 2 ! repeat process twice
!         do i = 2, j
!             ! TODO: check to make sure that this loop can be done in this way, i.e., that the sum need not be completed before subtracting from z
!             ! in exact arithmetic, it is equivalent; does order matter for roundoff?
!             coef = ddot(sze, z, incx, q(:, i), incy)
!             z = z - coef * q(:,i)
!         enddo
!     enddo 

! end