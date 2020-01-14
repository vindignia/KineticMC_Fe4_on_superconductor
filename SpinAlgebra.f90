Module  SpinAlgebra
Contains

    pure function spinModulus(matrix_size)
        implicit none
        INTEGER (Kind=4),intent(in) ::  matrix_size
        REAL (Kind=8) :: spinModulus
        spinModulus = Real((matrix_size-1), Kind=8)/2.d0
    end function spinModulus

    !---------------------------------------------------------

    pure function dimMatrix(S)
        implicit none
        REAL (Kind=8),intent(in) :: S
        INTEGER (Kind=4) ::  dimMatrix
        dimMatrix = Int((2*S), Kind=8) + 1
    end function dimMatrix

    !---------------------------------------------------------

    pure function S_z_coeff(S, m)
        implicit none
        REAL (Kind=8), intent(in)   :: S
        INTEGER (Kind=4),intent(in) :: m
        REAL (Kind=8) 	            :: S_z_coeff
        S_z_coeff = S + 1.d0 - real(m, kind=8)
    end function S_z_coeff

    !---------------------------------------------------------

    pure function S_plus_coeff(S, m)
        implicit none
        REAL (Kind=8), intent(in)   :: S
        INTEGER (Kind=4),intent(in) :: m
        REAL (Kind=8) 	            :: S_plus_coeff,arg,Ms
        Ms =  S_z_coeff(S, m)
        if (abs(Ms) .gt. S)then
            S_plus_coeff =0.d0
        else
            arg = S*(S+1.d0) - Ms*(Ms+1.d0)
            S_plus_coeff = dsqrt(arg)
        end if
    end function S_plus_coeff

    !---------------------------------------------------------

    pure function S_minus_coeff(S, m)
        implicit none
        REAL (Kind=8), intent(in)   :: S
        INTEGER (Kind=4),intent(in) :: m
        REAL (Kind=8) 	            :: S_minus_coeff,arg,Ms
        Ms =  S_z_coeff(S, m)
        if (abs(Ms) .gt. S)then
            S_minus_coeff =0.d0
        else
            arg = S*(S+1.d0) - Ms*(Ms-1.d0)
            S_minus_coeff = dsqrt(arg)
        end if
    end function S_minus_coeff


    !---------------------------------------------------------

    pure function S_z(Nd,eigenVect,p,q) ! <p| ... |q>
        implicit none
        INTEGER (Kind=4),intent(in) :: Nd, p,q
        COMPLEX (Kind=8),intent(in) :: eigenVect(Nd,Nd)
        COMPLEX (Kind=8)            :: S_z ! output
        !----
        INTEGER (Kind=4) :: m
        REAL (Kind=8) :: S, Ms
        COMPLEX (Kind=8) :: sum, arg
        REAL (Kind=8) 	            :: coeff

        S =spinModulus(Nd)
        sum = (0.d0, 0.d0)
        if(p.eq.q) then
            do m = 1,Nd
                arg = S_z_coeff(S, m)*dconjg(eigenVect(m,p))*eigenVect(m,q)
                sum = sum + arg
            end do
        endif
        S_z = sum
    end function S_z


    !---------------------------------------------------------

    pure function S_plus(Nd,eigenVect,p,q) ! <p| ... |q>
        implicit none
        INTEGER (Kind=4),intent(in) :: Nd, p,q
        COMPLEX (Kind=8),intent(in) :: eigenVect(Nd,Nd)
        COMPLEX (Kind=8)            :: S_plus ! output
        !----
        INTEGER (Kind=4) :: m
        REAL (Kind=8) :: S, Ms
        COMPLEX (Kind=8) :: sum, arg
        REAL (Kind=8) 	            :: coeff

        S =spinModulus(Nd)
        sum = (0.d0, 0.d0)
        do m = 2,Nd     ! ACHTUNG |M_s + 1> = | m-1>
            arg = S_plus_coeff(S, m)*dconjg(eigenVect(m-1,p))*eigenVect(m,q)
            sum = sum + arg
        end do
        S_plus = sum
    end function S_plus


    !---------------------------------------------------------

    pure function S_minus(Nd,eigenVect,p,q) ! <p| ... |q>
        implicit none
        INTEGER (Kind=4),intent(in) :: Nd, p,q
        COMPLEX (Kind=8),intent(in) :: eigenVect(Nd,Nd)
        COMPLEX (Kind=8)            :: S_minus ! output
        !----
        INTEGER (Kind=4) :: m
        REAL (Kind=8) :: S, Ms
        COMPLEX (Kind=8) :: sum, arg
        REAL (Kind=8) 	            :: coeff

        S =spinModulus(Nd)
        sum = (0.d0, 0.d0)
        do m = 1,Nd-1     ! ACHTUNG |M_s - 1> = | m+1>
            arg = S_minus_coeff(S, m)*dconjg(eigenVect(m+1,p))*eigenVect(m,q)
            sum = sum + arg
        end do
        S_minus = sum
    end function S_minus

    !---------------------------------------------------------

    pure function S_plus_sq(Nd,eigenVect,p,q) ! <p| ... |q>
        implicit none
        INTEGER (Kind=4),intent(in) :: Nd, p,q
        COMPLEX (Kind=8),intent(in) :: eigenVect(Nd,Nd)
        COMPLEX (Kind=8)            :: S_plus_sq ! output
        !----
        INTEGER (Kind=4) :: m
        REAL (Kind=8) :: S, Ms
        COMPLEX (Kind=8) :: sum, arg
        REAL (Kind=8) 	            :: coeff

        S =spinModulus(Nd)
        sum = (0.d0, 0.d0)
        do m = 3,Nd     ! ACHTUNG |M_s + 1> = | m-1>
            arg = S_plus_coeff(S, m)*S_plus_coeff(S, m-1)*dconjg(eigenVect(m-2,p))*eigenVect(m,q)
            sum = sum + arg
        end do
        S_plus_sq = sum
    end function S_plus_sq

    !---------------------------------------------------------

    pure function S_plus_S_z(Nd,eigenVect,p,q) ! <p| ... |q>  anticommutator: {S_+, S_z} = S_+ S_z + S_z S_+
        implicit none
        INTEGER (Kind=4),intent(in) :: Nd, p,q
        COMPLEX (Kind=8),intent(in) :: eigenVect(Nd,Nd)
        COMPLEX (Kind=8)            :: S_plus_S_z ! output
        !----
        INTEGER (Kind=4) :: m
        REAL (Kind=8) :: S, Ms
        COMPLEX (Kind=8) :: sum, arg
        REAL (Kind=8) 	            :: coeff

        S =spinModulus(Nd)
        sum = (0.d0, 0.d0)
        do m = 2,Nd     ! ACHTUNG |M_s + 1> = | m-1>
            arg = (2.*S_z_coeff(S, m) + 1.d0)*S_plus_coeff(S, m)*dconjg(eigenVect(m-1,p))*eigenVect(m,q)
            sum = sum + arg
        end do
        S_plus_S_z = sum
    end function S_plus_S_z

    !---------------------------------------------------------
    pure function S_minus_sq(Nd,eigenVect,p,q) ! <p| ... |q>
        implicit none
        INTEGER (Kind=4),intent(in) :: Nd, p,q
        COMPLEX (Kind=8),intent(in) :: eigenVect(Nd,Nd)
        COMPLEX (Kind=8)            :: S_minus_sq ! output
        !----
        INTEGER (Kind=4) :: m
        REAL (Kind=8) :: S, Ms
        COMPLEX (Kind=8) :: sum, arg
        REAL (Kind=8) 	            :: coeff

        S =spinModulus(Nd)
        sum = (0.d0, 0.d0)
        do m = 1,Nd-2     ! ACHTUNG |M_s - 1> = | m+1>
            arg = S_minus_coeff(S, m)*S_minus_coeff(S, m+1)*dconjg(eigenVect(m+2,p))*eigenVect(m,q)
            sum = sum + arg
        end do
        S_minus_sq = sum
    end function S_minus_sq

    !---------------------------------------------------------
    pure function S_minus_S_z(Nd,eigenVect,p,q) ! <p| ... |q> anticommutator: {S_-, S_z} = S_- S_z + S_z S_-
        implicit none
        INTEGER (Kind=4),intent(in) :: Nd, p,q
        COMPLEX (Kind=8),intent(in) :: eigenVect(Nd,Nd)
        COMPLEX (Kind=8)            :: S_minus_S_z ! output
        !----
        INTEGER (Kind=4) :: m
        REAL (Kind=8) :: S, Ms
        COMPLEX (Kind=8) :: sum, arg
        REAL (Kind=8) 	            :: coeff

        S =spinModulus(Nd)
        sum = (0.d0, 0.d0)
        do m = 1,Nd-1    ! ACHTUNG |M_s - 1> = | m+1>
            arg = (2.*S_z_coeff(S, m) - 1.d0)*S_minus_coeff(S, m)*dconjg(eigenVect(m+1,p))*eigenVect(m,q)
            sum = sum + arg
        end do
        S_minus_S_z = sum
    end function S_minus_S_z


    !---------------------------------------------------------
    function S_prog_real(S_proj_complex)
        REAL (Kind=8) :: S_prog_real,S_prog_imag
        COMPLEX (Kind=8) :: S_proj_complex
        S_prog_real = REAL(S_proj_complex, kind=8)
        S_proj_imag = AIMAG(S_proj_complex)
        if(abs(S_prog_imag).gt.1.d-13)write(*,*) "ACHUNG imaginary spin projection!!!"
    end function S_prog_real


End Module SpinAlgebra