!
!
! Complex implementations
!
!

subroutine lanczos_tridiag_reortho_cb(H, u0, uu, alpha, beta, k, sze)
    implicit none
    BEGIN_DOC
    ! Function that performs lanczos triadiaglonization of a hermitian matrix
    ! Performs reorthogonalization of intermediate vectors
    ! Outputs triadiagonlized Hermitian matrix as alpha and beta vectors
    ! Outputs up to k elements of tridiagonlization and orthogonal basis (testing only)
    END_DOC

    integer, intent(in)             :: k, sze
    integer                         :: i, ii, j, incx, incy
    complex*16, intent(in)          :: H(sze, sze), u0(sze)
    complex*16                      :: z(sze)
    complex*16, intent(out)         :: uu(sze,k)
    complex*16                      :: zdotc
    double precision, intent(out)   :: beta(k), alpha(k)
    double precision                :: dznrm2, coef
    
    incx = 1 
    incy = 1

    ! initialize vectors
    uu(:,1) = u0
    alpha = (0.d0, 0.d0)
    beta = 0.d0
    
    do i = 1, k
        z = (0.d0, 0.d0)
        call zhemv('U', sze, (1.d0,0.d0), H, sze, uu(:,i), incx, (0.d0,0.d0), z, incy)
        alpha(i) = real(zdotc(sze, uu(:,i), incx, z, incy), kind=16)

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


end

subroutine lanczos_tridiag_reortho_c(H, u0, alpha, beta, k, sze)
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
    complex*16                      :: z(sze), uu(sze,k)
    complex*16                      :: zdotc
    double precision, intent(out)   :: beta(k), alpha(k)
    double precision                :: dznrm2, coef
    
    incx = 1 
    incy = 1

    ! initialize vectors
    uu(:,1) = u0
    alpha = (0.d0, 0.d0)
    beta = 0.d0
    
    do i = 1, k
        z = (0.d0, 0.d0)
        call zhemv('U', sze, (1.d0,0.d0), H, sze, uu(:,i), incx, (0.d0,0.d0), z, incy)
        alpha(i) = real(zdotc(sze, uu(:,i), incx, z, incy), kind=16)

        if (i == k) then
            exit
        end if

        do ii = 1, 2 ! repeat process twice
            !$OMP PARALLEL DO PRIVATE(j, coef) SHARED(sze, incx, incy, uu) &
            !$OMP REDUCTION(+:z)
            do j = 1, i
                coef = zdotc(sze, z, incx, uu(:,j), incy)
                z = z - coef * uu(:,j)
            enddo
            !$OMP END PARALLEL DO
        enddo 

        beta(i+1) = dznrm2(sze, z, incx)
        if (beta(i+1) < 1e-16) then ! some small number
        ! add some type of escape or warning for beta(i) = 0
            exit
        end if

        uu(:,i+1) = z / beta(i+1)
    enddo


end

subroutine lanczos_tridiag_sparse_reortho_c(H_v, H_c, H_p, u0, alpha, beta, k, nnz, sze)
    implicit none
    BEGIN_DOC
    ! Function that performs lanczos triadiaglonization of a hermitian matrix
    ! Performs reorthogonalization of intermediate vectors
    ! Outputs triadiagonlized Hermitian matrix as alpha and beta vectors
    ! Outputs up to k elements of tridiagonlization
    END_DOC

    integer, intent(in)             :: k, sze, nnz, H_c(nnz), H_p(sze+1)
    integer                         :: i, ii, j, incx, incy
    complex*16, intent(in)          :: H_v(nnz), u0(sze)
    complex*16                      :: z(sze), uu(sze,k)
    complex*16                      :: zdotc
    double precision, intent(out)   :: beta(k), alpha(k)
    double precision                :: dznrm2, coef
    
    incx = 1 
    incy = 1

    ! initialize vectors
    uu(:,1) = u0
    alpha = (0.d0, 0.d0)
    beta = 0.d0
    
    do i = 1, k
        z = (0.d0, 0.d0)
        call sparse_csr_zmv(H_v, H_c, H_p, uu(:,i), z, sze, nnz)
        alpha(i) = real(zdotc(sze, uu(:,i), incx, z, incy), kind=16)

        if (i == k) then
            exit
        end if

        do ii = 1, 2 ! repeat process twice
            !$OMP PARALLEL DO PRIVATE(j, coef) SHARED(sze, incx, incy, uu) &
            !$OMP REDUCTION(+:z)
            do j = 1, i
                coef = zdotc(sze, z, incx, uu(:,j), incy)
                z = z - coef * uu(:,j)
            enddo
            !$OMP END PARALLEL DO
        enddo 

        beta(i+1) = dznrm2(sze, z, incx)
        if (beta(i+1) < 1e-16) then ! some small number
        ! add some type of escape or warning for beta(i) = 0
            exit
        end if

        uu(:,i+1) = z / beta(i+1)
    enddo


end

!
!
! Real implementations
!
!

subroutine lanczos_tridiag_reortho_rb(H, u0, uu, alpha, beta, k, sze)
    implicit none
    BEGIN_DOC
    ! Function that performs lanczos triadiaglonization of a hermitian matrix
    ! Performs reorthogonalization of intermediate vectors
    ! Outputs triadiagonlized Hermitian matrix as alpha and beta vectors
    ! Outputs up to k elements of tridiagonlization and orthogonal basis (testing only)
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

end

subroutine lanczos_tridiag_reortho_r(H, u0, alpha, beta, k, sze)
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
    double precision                :: z(sze), uu(sze,k)
    double precision, intent(out)   :: alpha(k), beta(k)
    double precision                :: ddot, coef
    double precision                :: dnrm2
    
    incx = 1 
    incy = 1

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
            !$OMP PARALLEL DO PRIVATE(j, coef) SHARED(sze, incx, incy, uu) &
            !$OMP REDUCTION(+:z)
            do j = 1, i
                coef = ddot(sze, z, incx, uu(:,j), incy)
                z = z - coef * uu(:,j)
            enddo
            !$OMP END PARALLEL DO
        enddo 
        
        beta(i+1) = dnrm2(sze, z, incx)
        if (beta(i+1) < 1e-16) then ! some small number
        ! add some type of escape or warning for beta(i) = 0
            print *, "escaping early"
            exit
        end if

        uu(:,i+1) = z / beta(i+1)
    enddo

end

subroutine lanczos_tridiag_sparse_reortho_r(H_v, H_c, H_p, u0, alpha, beta, k, nnz, sze)
    implicit none
    BEGIN_DOC
    ! Function that performs lanczos triadiaglonization of a hermitian matrix
    ! Uses CSR sparse H 
    ! Performs reorthogonalization of intermediate vectors
    ! Outputs triadiagonlized Hermitian matrix as alpha and beta vectors
    ! Outputs up to k elements of tridiagonlization
    END_DOC

    integer, intent(in)             :: k, sze, nnz, H_c(nnz), H_p(sze+1)
    integer                         :: i, ii, j, incx, incy
    double precision, intent(in)    :: H_v(nnz), u0(sze)
    double precision                :: z(sze), uu(sze,k)
    double precision, intent(out)   :: alpha(k), beta(k)
    double precision                :: ddot, coef
    double precision                :: dnrm2
    
    incx = 1 
    incy = 1

    ! initialize vectors
    uu(:,1) = u0
    alpha = 0.d0
    beta = 0.d0
    
    do i = 1, k
        z = 0
        call sparse_csr_dmv(H_v, H_c, H_p, uu(:,i), z, sze, nnz)
        alpha(i) = ddot(sze, z, incx, uu(:,i), incy)

        if (i == k) then
            exit
        end if

        do ii = 1, 2 ! repeat process twice
            !$OMP PARALLEL DO PRIVATE(j, coef) SHARED(sze, incx, incy, uu) &
            !$OMP REDUCTION(+:z)
            do j = 1, i
                coef = ddot(sze, z, incx, uu(:,j), incy)
                z = z - coef * uu(:,j)
            enddo
            !$OMP END PARALLEL DO
        enddo 
        
        beta(i+1) = dnrm2(sze, z, incx)
        if (beta(i+1) < 1e-16) then ! some small number
        ! add some type of escape or warning for beta(i) = 0
            print *, "escaping early"
            exit
        end if

        uu(:,i+1) = z / beta(i+1)
    enddo

end