Module TransitionProbabilities
    USE SpinAlgebra

Contains

    SUBROUTINE transition_rates(Nd, TIN, H_x, H_y, H_z, W_Ale, eig, S_proj_out)
        implicit none

        INTEGER (Kind = 4) :: Nd, I, p, q, mm, jm(Nd)
        COMPLEX (Kind = 8) :: tmp_C, tmp_x, tmp_y
        REAL (Kind = 8) :: lambda(Nd, Nd), dlambda(Nd, Nd)
        REAL (Kind = 8) :: DASH(Nd), IDASH(Nd), delpop(Nd), vh(2000), Vmagn(2000)
        REAL (Kind = 8) :: vmeq(2000)

        REAL (Kind = 8) :: proj_max(Nd), eigenvalue(Nd), tmp, H
        REAL (Kind = 8) :: T, TIN, H_x, H_y, H_z
        REAL (Kind = 8) :: Energy_Ale(Nd), W_Ale(Nd, Nd)
        REAL (Kind = 8) :: S_plus_coeff, S_minus_coeff
        REAL (Kind = 8) :: En_levels(Nd)
        COMPLEX (Kind = 8) :: component_plus, component_minus, test_norm
        REAL (Kind = 8) :: S_z_proj(Nd)
        COMPLEX (Kind = 8) :: S_x_proj(Nd), S_y_proj(Nd)
        COMPLEX (Kind = 8) :: S_plus_mat(Nd, Nd), S_minus_mat(Nd, Nd)
        COMPLEX (Kind = 8) :: S_x_mat(Nd, Nd), S_y_mat(Nd, Nd)
        COMPLEX (Kind = 8) :: S_proj_out(3, Nd)
        logical :: empty(Nd)
        ! - diag routine ZHEEV
        REAL (Kind = 8) :: eig(Nd)
        COMPLEX (Kind = 8) :: HamMat(Nd, Nd)
        COMPLEX (Kind = 8), allocatable :: WORK(:), RWORK(:)
        INTEGER (Kind = 4) :: LWORK, INFO
        ! -- commons
        REAL (Kind = 8) :: S, GI, D, E, B, B42, C, B43, B66, gamma_0, gamma_tunnel
        REAL (Kind = 8) :: g1, g2, gam1a, gam2a, gam1b, gam2b
        COMMON S,GI,D,E,B,B42,C,B43,B66,gamma_0,gamma_tunnel,g1,g2

        T = TIN

        call Hamiltonian(Nd, H_x, H_y, H_z, HamMat)

        LWORK = 363 ! 2*Nd
        allocate(WORK(LWORK), RWORK(3 * Nd - 2))

        call zheev('V', 'U', Nd, HamMat, Nd, eig, WORK, LWORK, RWORK, INFO)
        !        write(*,*)'info = ', INFO
        !        write(*,*)'optimal LWORK = ', WORK(1)

        deallocate(WORK, RWORK)

        call PhysicalLabeling(Nd, HamMat, eig)

        do q = 1,Nd
            do p = 1, Nd
                W_Ale(q, p) = transition_rate(T, gamma_tunnel, gamma_0, g1, g2, Nd, eig, HamMat, p, q)
            end do
         ! spin components
            S_proj_out(1,q) = (0.5d0,0.d0)*(S_plus(Nd, HamMat, q, q) + S_minus(Nd, HamMat, q, q))
            S_proj_out(2,q) = (0.d0,-0.5d0)*(S_plus(Nd, HamMat, q, q) - S_minus(Nd, HamMat, q, q))
            S_proj_out(3,q) = S_z(Nd, HamMat, q, q)
        end do

        call flush(6)

        RETURN
    END

    !---------------------------------------------------------

    pure function transition_rate(T, gamma_tunnel, gg0, gg1, gg2, Nd, eigenVal, eigenVect, p, q) ! <p| ... |q>
        implicit none
        REAL (Kind = 8), intent(in) :: T, gamma_tunnel, gg0, gg1, gg2
        INTEGER (Kind = 4), intent(in) :: Nd, p, q
        COMPLEX (Kind = 8), intent(in) :: eigenVect(Nd, Nd)
        REAL (Kind = 8), intent(in) :: eigenVal(Nd)
        REAL (Kind = 8) :: S_plus_sq_sq, S_minus_sq_sq, S_plus_S_z_sq, S_minus_S_z_sq
        REAL (Kind = 8) :: spin_phonon, thermal_weight, deltaE
        REAL (Kind = 8) :: transition_rate ! output
        !----
        INTEGER (Kind = 4) :: m
        REAL (Kind = 8) :: S, Ms
        COMPLEX (Kind = 8) :: sum, arg
        REAL (Kind = 8) :: coeff

        S_plus_sq_sq = modSquare(S_plus_sq(Nd, eigenVect, p, q))
        S_minus_sq_sq = modSquare(S_minus_sq(Nd, eigenVect, p, q))
        S_plus_S_z_sq = S_plus_S_z(Nd, eigenVect, p, q)**2.
        S_minus_S_z_sq = S_minus_S_z(Nd, eigenVect, p, q)**2.

        spin_phonon = gg2*gg2*(S_plus_sq_sq + S_minus_sq_sq) + gg1*gg1*(S_plus_S_z_sq + S_minus_S_z_sq)
        spin_phonon = gg0 * spin_phonon

        deltaE = eigenVal(p) - eigenVal(q)
        thermal_weight = (deltaE**3.) / (dexp(deltaE / T) - 1.d0)

        if (abs(deltaE).lt.1.d-1) then ! neighborhood of level crossing
            transition_rate = spin_phonon * (deltaE**2) * T
            !     -----    INTRODUCE PURE TUNNELING CHANNEL  ----
             if (abs(p - q).ge.9) transition_rate = transition_rate + gamma_tunnel
        else
            if (abs(deltaE) / T.lt.40.d0) then
                thermal_weight = (deltaE**3.) / (dexp(deltaE / T) - 1.d0)
            else
                if (deltaE.gt.0.) then
                    ! limit T -> 0 for absorption transitions
                    thermal_weight = 0.d0
                else
                    ! limit T -> 0 for emission transitions
                    thermal_weight = -deltaE**3.
                endif
            endif

            transition_rate = spin_phonon * thermal_weight

        endif

    end function transition_rate

    !---------------------------------------------------------

    pure function modSquare(cc)
        implicit none
        COMPLEX (Kind = 8), intent(in) :: cc
        REAL (Kind = 8) :: modSquare
        modSquare = abs(cc) * abs(cc)
    end function modSquare

    !---------------------------------------------------------

    subroutine PhysicalLabeling(Nd, eigenVect, eig)

        implicit none

        logical :: empty(Nd)
        REAL (Kind = 8) :: component_sq, sum
        REAL (Kind = 8) :: eig(Nd), eigRoutine(Nd), proj_max(Nd)
        COMPLEX (Kind = 8) :: eigenVect(Nd, Nd), eigenVectRoutine(Nd, Nd)
        INTEGER (Kind = 4) :: Nd, I, J, JM(Nd), II

        empty(:) = .true.  ! initialization
        ! we assume that |m = 0> has always the highest energy: it can be checked with J = MAXLOC(eig, Nd)
        JM(11) = 6
        empty(6) = .false.

        do J = 1, Nd - 1 ! i exclude the last one already assigned
            proj_max(J) = 0.
            do I = 1, Nd
                if(empty(I))then
                    component_sq = dconjg(eigenVect(I, J)) * eigenVect(I, J)
                    if(component_sq>proj_max(J))then
                        proj_max(J) = component_sq
                        JM(J) = I
                    endif
                endif
            enddo
            II = JM(J)
            empty(II) = .false.
        enddo

        eigRoutine(:) = eig(:)
        eigenVectRoutine(:, :) = eigenVect(:, :)

        do J = 1, Nd
            I = JM(J)
            eig(I) = eigRoutine(J)    ! reorder the eigenvalues according to max proj on |m>
            do II = 1, Nd
                eigenVect(II, I) = eigenVectRoutine(II, J)    ! reorder the egenvectors
            enddo
        enddo

    end subroutine PhysicalLabeling


    !-----------------------------

    SUBROUTINE ZEEM(Nd, H_x, H_y, H_z, DAT)
        implicit none
        COMPLEX (Kind = 8) :: DAT(Nd, Nd)
        REAL (Kind = 8), allocatable :: Z(:)
        COMPLEX (Kind = 8) :: C_x, C_y, C_z
        COMPLEX (Kind = 8), parameter :: PR = (1.d0, 0.d0)
        COMPLEX (Kind = 8), parameter :: PI = (0.d0, 1.d0)
        REAL (Kind = 8) :: H_x, H_y, H_z, SP, BETA
        INTEGER (Kind = 4) :: Nd, I, J, IS1, IDF
        ! -- commons
        REAL (Kind = 8) :: S, GI, D, E, B, B42, C, B43, B66, gamma_0, gamma_tunnel
        REAL (Kind = 8) :: g1, g2, gam1a, gam2a, gam1b, gam2b
        COMMON S,GI,D,E,B,B42,C,B43,B66,gamma_0,gamma_tunnel,g1,g2

        allocate(Z(Nd))
        I = 1
        do IS1 = int(S), -int(S), -1
            Z(I) = real(IS1)
            I = I + 1
        enddo

        BETA = (9.274078d-5) / 1.380662d0        ! Kelvin/Gauss
        C_x = H_x * PR
        C_y = H_y * PI
        C_z = H_z * PR
        do I = 1, Nd
            do J = I, Nd
                DAT(I, J) = (0.00, 0.00)
                IDF = J - I
                IF(Idf.eq.0) DAT(I, I) = BETA * GI * C_z * Z(I)

                IF(IDF.eq.1) then
                    SP = SQRT(S * (S + 1.0) - Z(J) * (Z(J) + 1.0))
                    DAT(I, J) = GI * BETA * (C_x - C_y) * .5 * SP
                endif
                DAT(J, I) = CONJG(DAT(I, J))
            enddo
        enddo
        RETURN
    END subroutine ZEEM
    !---------------------------------------------------------

    subroutine Hamiltonian(Nd, H_x, H_y, H_z, Ham)
        implicit none
        REAL (Kind = 8) :: MAT(Nd, Nd)
        COMPLEX (Kind = 8) :: ZEE(Nd, Nd)
        REAL (Kind = 8), allocatable :: Z(:)
        COMPLEX (Kind = 8), parameter :: PR = (1.d0, 0.d0)
        COMPLEX (Kind = 8) :: Ham
        dimension :: Ham(Nd, Nd)
        REAL (Kind = 8) :: SS
        REAL (Kind = 8) :: H_x, H_y, H_z
        INTEGER (Kind = 4) :: Nd, I, J, IS1, IDF
        ! -- commons
        REAL (Kind = 8) :: S, GI, D, E, B, B42, C, B43, B66, gamma_0, gamma_tunnel
        REAL (Kind = 8) :: g1, g2, gam1a, gam2a, gam1b, gam2b
        COMMON S,GI,D,E,B,B42,C,B43,B66,gamma_0,gamma_tunnel,g1,g2

        allocate(Z(Nd))
        I = 1
        do IS1 = int(S), -int(S), -1
            Z(I) = real(IS1, Kind=8)
            I = I + 1
        enddo

        do I = 1, Nd
            do J = I, Nd
                MAT(I, J) = 0.d0
                IDF = J - I
                if (IDF.eq.0) then
                    MAT(I, I) = D * Z(I)**2 + 35. * B * Z(I)**4 - B * (30 * S * (S + 1) - 25.) * z(I)**2
                endif
                if(IDF.eq.2) then
                    MAT(I, J) = E * 0.5 * DSQRT((S - Z(J)) * (S + Z(J) + 1.) * (S - Z(J) - 1.) * (S + Z(J) + 2.))
                    MAT(I, J) = MAT(I, J) + 0.25 * B42 * DSQRT(S * (S + 1) - z(J) * (z(j) + 1.)) * &
                    DSQRT(S * (S + 1.) - (z(J) + 1.) * (Z(J) + 2.)) * &
                            (7. * z(j) * z(J) + 14. * z(J) - S * (S + 1.) + 9.) * 2
                endif
                if (IDF.eq.3) then
                    SS = S * (S + 1.)
                    MAT(I, J) = MAT(I, J) + 0.25 * B43 * Z(I + 3) * DSQRT((SS - Z(I) * Z(I + 1)) * &
                            (SS - Z(I + 1) * Z(I + 2)) * (SS - Z(I + 2) * Z(I + 3))) + 0.25 * B43 * Z(I) * &
                            DSQRT((SS - Z(I) * Z(I + 1)) * (SS - Z(I + 1) * Z(I + 2)) * (SS - Z(I + 2) * Z(I + 3)))
                endif
                if (IDF.eq.4) then
                    SS = S * (S + 1.)
                    MAT(I, J) = MAT(I, J) + 0.5 * C * DSQRT((SS - Z(I) * Z(I + 1)) * (SS - Z(I + 1) * z(I + 2)) * &
                            (SS - Z(I + 2) * z(I + 3)) * (SS - Z(I + 3) * Z(I + 4)))
                endif
                if (IDF.eq.6) then
                    MAT(I, J) = MAT(I, J) + 0.5 * B66 * DSQRT((SS - Z(I) * Z(I + 1)) * (SS - Z(I + 1) * &
                            Z(I + 2)) * (SS - Z(I + 2) * Z(I + 3)) * (SS - Z(I + 3) * Z(I + 4)) * (SS - Z(I + 4) * &
                            Z(I + 5)) * (SS - Z(I + 5) * Z(I + 6)))
                end if
            end do
        end do

        deallocate(Z)

        call ZEEM(Nd, H_x, H_y, H_z, ZEE)
        do I = 1, Nd
            do J = I, Nd
                Ham(I, J) = (0.d0, 0.d0)
                Ham(I, J) = MAT(I, J) * PR + ZEE(I, J)
                Ham(J, I) = DCONJG(Ham(I, J))
            end do
        end do

        return

    end subroutine Hamiltonian

End Module TransitionProbabilities