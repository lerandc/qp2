program spectral_density
    implicit none
    BEGIN_DOC
    ! Program that calculates the spectral density.
    END_DOC
    !print *, cfraction_test
    print *, lanczos_test
    print *, lanczos_alpha
    print *, lanczos_beta
end

BEGIN_PROVIDER [integer, cfraction_test]
    implicit none
    use cfraction

    !!!! Continued fraction tests
    ! Real tests
    double precision, allocatable    :: a_r(:), b_r(:)
    double precision                 :: a0_r, phi, phi_ref, pi, pi_ref
    integer                          :: i, N

    N = 97
    allocate(a_r(N), b_r(N))

    a_r = 1
    b_r = 1
    a0_r = 1

    print *, "Calculating phi directly"
    phi_ref = (1.d0+sqrt(5.d0))/2.d0
    print *, phi_ref

    print *, "Calculating phi with continued fraction expansion"
    phi = cfraction_r(a0_r, a_r, b_r, N)
    print *, phi

    print *, "Calculating pi directly"
    pi_ref = acos(-1.d0)
    print *, pi_ref

    print *, "Calculating pi with continued fraction expansion"
    a0_r = 3
    b_r = (/7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2,&
            2, 1, 84, 2, 1, 1, 15, 3, 13, 1, 4, 2, 6, 6, 99, 1, 2, 2, 6,&
            3, 5, 1, 1, 6, 8, 1, 7, 1, 2, 3, 7, 1, 2, 1, 1, 12, 1, 1, 1,&
            3, 1, 1, 8, 1, 1, 2, 1, 6, 1, 1, 5, 2, 2, 3, 1, 2, 4, 4, 16,&
            1, 161, 45, 1, 22, 1, 2, 2, 1, 4, 1, 2, 24, 1, 2, 1, 3, 1,&
            2, 1/)

    pi = cfraction_r(a0_r, a_r, b_r, N)
    print *, pi

    deallocate(a_r, b_r)

    ! Complex tests
    complex*16, allocatable     :: a_c(:), b_c(:)
    complex*16                  :: a0_c, ref_val, test_val, x

    N = 100

    allocate(a_c(N), b_c(N))

    x = (2.3, -1.4) ! random value
    print *, x
    a0_c = (1.d0, 0.d0)
    a_c(1) =  x
    b_c(1) = 1

    do i = 2, N
        a_c(i) = -(i-1)*x
        b_c(i) = i + x
    enddo

    print *, "Calculating reference exp(x)"
    ref_val = exp(x)
    print *, ref_val

    print *, "Calculating test exp(x) with complex continued fraction"
    test_val = cfraction_c(a0_c, a_c, b_c, N)
    print *, test_val
    deallocate(a_c, b_c)

    cfraction_test = 1
END_PROVIDER

BEGIN_PROVIDER [integer, lanczos_test]
&BEGIN_PROVIDER [double precision, lanczos_alpha, (lanczos_N)]
&BEGIN_PROVIDER [double precision, lanczos_beta, (lanczos_N)]
    implicit none
    lanczos_test = 23

    !!! Real tests
    integer                       :: k, sze, i
    double precision, allocatable :: H(:,:), uu(:,:), uu_h(:,:), Q(:,:), u(:), alpha(:), beta(:)
    double precision                :: dnrm2

    sze = 100
    k = 100
    allocate(H(sze,sze), uu(sze,sze), uu_h(sze,sze), Q(sze, sze), u(sze), alpha(sze), beta(sze))

    H = 0
    u = 1.0
    u(1) = 2.0
    u = u / dnrm2(sze, u, 1)

    print *, 0.d0

    ! test evals spans -1000,1000
    ! logspaced

    ! set the exponents first
    do i = 1, sze-1
        H(i,i) = sze - i + 1
        print *, i, H(i,i)
    end do

    ! print *, H

    ! call lanczos_tridiag_r(H, u, alpha, beta, k, sze)

    ! print *, "-----------"
    ! do i = 1, sze
    !     print *, alpha(i), beta(i)
    ! end do 


    call lanczos_tridiag_reortho_r(H, u, uu, alpha, beta, k, sze)

    ! print *, "-----------"
    ! do i = 1, sze
    !     print *, alpha(i), beta(i)
    ! end do 

    lanczos_alpha = alpha
    lanczos_beta = beta

    !!! Complex tests
END_PROVIDER