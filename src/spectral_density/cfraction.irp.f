
complex*16 function cfraction_c(a0, a, b, N)
    BEGIN_DOC
    ! Compute continued fraction expansions with Steed's Algorithm
    ! See https://doi.org/10.1016/0010-4655(85)90025-6 or
    ! https://dlmf.nist.gov/3.10
    ! N is number of convergents
    ! a0 is constant out front
    ! arrays a,b contain coefficients of continued fraction itself
    END_DOC
    implicit none
        
    integer, intent(in)     :: N
    complex*16, intent(in)  :: a0, a(N), b(N)
    complex*16              :: Cn, dC, Dn
    integer                 :: i

    Cn = a0
    Dn = (1.d0, 0.d0)/b(1)
    dC = a(1)*Dn
    Cn = Cn + dC

    do i = 2, N
        Dn = (1.d0, 0.d0)/ (Dn*a(i) + b(i))
        dC = (b(i) * Dn - (1.d0, 0.d0)) * dC
        Cn = Cn + dC 
    enddo
        
        cfraction_c = Cn
        return
    
end function cfraction_c


double precision function cfraction_r(a0, a, b, N)
    BEGIN_DOC
    ! Compute continued fraction expansions with Steed's Algorithm
    ! See https://doi.org/10.1016/0010-4655(85)90025-6 or
    ! https://dlmf.nist.gov/3.10
    ! N is number of convergents
    ! a0 is constant out front
    ! arrays a,b contain coefficients of continued fraction itself
    END_DOC
    implicit none

    integer, intent(in)             :: N
    double precision, intent(in)    :: a0, a(N), b(N)
    double precision                :: Cn, dC, Dn
    integer                         :: i


    Cn = a0
    Dn = 1.d0/b(1)
    dC = a(1)*Dn
    Cn = Cn + dC

    do i = 2, N
        Dn = 1.d0/ (Dn*a(i) + b(i))
        dC = (b(i) * Dn - 1.d0) * dC
        Cn = Cn + dC 
    enddo

    cfraction_r = Cn
    return

end function cfraction_r