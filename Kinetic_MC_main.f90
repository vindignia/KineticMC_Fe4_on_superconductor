! ordering of include is important aslo within
include 'SpinAlgebra.f90'
include 'TransitionProbModule.f90'
include 'FileNamesModule.f90'
include 'TestModule.f90'
include 'RandomModule.f90'
include 'Parameters.f90'
include 'VariablesMain.f90'

! I am changing al the types into REAL (Kind=8) consistently in every subrouitne and commons
program Kinetic_MC_Fe4

    USE SpinAlgebra
    USE TransitionProbabilities
    USE FileNames
    USE Tests
    USE Randomize
    USE ComputationParameters
    USE GlobalVariables

    S = spinModulus(matrix_size)

    call cpu_time(start)

    print '(/"kinetic MonteCarlo Fe4 on superconductor"/)'

    ! RND number generator
    idum=seed
    mmm = mzranset(521288629,362436069,16163801,idum)

    t_inc=time_max/dble(N_time_slot)
    h_step=dabs(h_f-h_i)/dble(N_time_slot)		! with h=0 in the range

    if(.NOT.read_fields) then
        print '("t_inc = " (F8.4) x "[s]")', t_inc
        print '("h_step = " (F8.4) x "[T]")', (1d-4)*h_step
        print '("v_sweep = " (F8.4) x "[T/s]")', (1d-4)*h_step/t_inc
    endif

    i=1

    WW(:,:,:)=0.d0

    theta=rad*theta_0

    if(read_fields)open(unit=2, file=filename_fields)


    field_critical = 0.d0
    do icritical = 0,0

        if (icritical == 1) field_critical =  h_crit_Pb

        ! centered in the interval
        do itime=1,N_time_slot

            if( read_fields ) then
                read(2,*)ih, field_nominal, field_true
                field_values(itime)=field_true
            elseif ( average_critical_field )then

                tmp = dble(itime-1)*h_step + h_i + h_step*0.5
                field_values(itime) = tmp
                if ( dabs(tmp) < field_critical ) field_values(itime) = 0.d0
            else
                field_values(itime)=dble(itime-1)*h_step + h_i + h_step*0.5
            endif

        enddo ! itime loop

        if(read_fields)close(2)

        tmp_h =  h_f

        ih = (tmp_h-h_i)/h_step
        ih = ih + 1

        print '(" time max = " (F8.4) x "[s]"/)', time_max

        ! write the eigenvalues and transition rates in a separate files
        if (save_eig_and_W) then
            open(unit=1, file=filename_eigenvalues)
            open(unit=2, file=filename_transition_rates)
            phi = rad*phi_0
            call labelsEnergyLevels(matrix_size, T,theta,phi,N_time_slot,field_values)
            close(1)
            close(2)
        end if


        !==================  MAIN BLOCK ====================================

        call random_num_generator	! intialize the seed for rnd number generator

        Prob(:,:)=0.d0
        Norm(:)=0.d0

        do iter=1,iter_max

            ! --- define transition rates differently for phi_rnd or not
            do itime=1,N_time_slot
                if(random_phi) then
                    rnd_phi = random()
                    phi = pi*rnd_phi
                else
                    phi = rad*phi_0
                endif

                H_x = dsin(theta)*dcos(phi)*field_values(itime)
                H_y = dsin(theta)*dsin(phi)*field_values(itime)
                H_z = dcos(theta)*field_values(itime)	! ACHTUNG field_values(itime)

                ! transition probabilities Gatteschi-Sessoli-Villain
                call transition_rates(matrix_size,T,H_x,H_y,H_z,W_single,En_levels,S_proj)

                do q=1,matrix_size
                    do p=1,matrix_size
                        if(p.ne.q)WW(p,q,itime)=W_single(p,q)
                    enddo
                enddo

            enddo   		!   itime transition rates

            call MonteCarlo_sweep(Prob,Norm)

            !----  write file
            if(mod(iter,iwrite)==0) then
                write(*,*)'write temportary output iter = ',iter
                call write_current_averages
            end if

        enddo !iter

    enddo ! h critical loop

    write(*,*)
    call cpu_time(finish)
    write(*,*)'CPU time [seconds]= ',finish-start


contains

    !---------------------------------------------------------------------

    subroutine Init_ferro
        INTEGER (Kind=4)	      :: ii

        do ii=1,number_of_spins
            conf(ii) = p_initial
        enddo

    end Subroutine Init_ferro

    !---------------------------------------------------------------------

    pure function Boltzmann_weight(Nd,TT,eigenVal)
        implicit none
        INTEGER (Kind=4),intent(in) :: Nd
        REAL (Kind=8),intent(in)    :: eigenVal(Nd),TT
        REAL (Kind=8)               :: Boltzmann_weight(Nd) ! output
        INTEGER (Kind=4)            :: pp
        REAL (Kind=8)               :: zeta_part

        zeta_part=0.
        do pp=1,Nd
            zeta_part = zeta_part + dexp(-eigenVal(pp)/TT)
        enddo
        do pp=1,Nd
            Boltzmann_weight(pp) = dexp(-eigenVal(pp)/TT)/zeta_part
        enddo
    end function

    !---------------------------------------------------------------------

    subroutine mag_equilibrium_z
        REAL (Kind=8)		:: Proj_eq,zeta_sum
        REAL (Kind=8)       :: P_eq(matrix_size)

        zeta_sum=0.
        do p=1,matrix_size
            zeta_sum = zeta_sum + dexp(-En_levels(p)/T)
        enddo

        do p=1,matrix_size
            P_eq(p)=dexp(-En_levels(p)/T)/zeta_sum
        enddo

        Proj_eq = 0.d0
        do p=1,matrix_size
            Proj_eq = Proj_eq - S_proj(3,p)*P_eq(p)
        enddo
        Mag_eq =  Proj_eq/S

    end subroutine mag_equilibrium_z

    !---------------------------------------------------------------------

    Subroutine Advance_time_squares(ii,t_i,q_final,Delta_t_min)

        INTEGER (Kind=4)	:: q_final,ii,q_max
        REAL (Kind=8)		:: tmp_int
        REAL (Kind=8)		:: delta_tt,control,rr,h_mod,tmp_min,control_min
        REAL (Kind=8)		:: int_tt(matrix_size),log_rr(matrix_size),t_step(matrix_size)
        REAL (Kind=8)		:: t_i,tt,Delta_t_min,W_max

        ! Implementation of the First-Reaction-Algorithm with discretized, "squarelike" transition rates.
        ! write(*,*)'in routine Advance time squares'

        p = conf(ii)
        W_max = 0.d0
        do q=1,matrix_size

            !R  rr = 0
            !R  call RANDOM_NUMBER(rr)
            rr = random()
            log_rr(q)=-dlog(rr)

            W_single(p,q)=WW(p,q,itime)

            if(W_single(p,q)>W_max)then
                W_max=W_single(p,q)
                q_max=q
            endif

        enddo !q loop

        int_tt(:)=0.
        t_step(:)=time_max

        do q=1,matrix_size

            tmp_int=WW(p,q,itime)*(dble(itime)*t_inc -t_i)

            if(tmp_int .ge. log_rr(q))then
                t_step(q)=log_rr(q)/WW(p,q,itime)
            else
                int_tt(q) = tmp_int
                ih = itime
                do while(int_tt(q) < log_rr(q))

                    ih = ih + 1
                    if(ih > N_time_slot) then
                        t_step(q)=time_max
                        exit
                    endif

                    tmp_int = int_tt(q) + WW(p,q,ih)*t_inc

                    if(tmp_int > log_rr(q))then
                        t_step(q) = (log_rr(q)-int_tt(q))/WW(p,q,ih)
                        t_step(q) = t_step(q) + (ih-1)*t_inc - t_i
                        exit
                    else
                        int_tt(q) = tmp_int
                    endif

                enddo ! while

            endif ! integral condition

        enddo ! q loop

        rr = random()

        q_final = 1 + matrix_size*rr

        Delta_t_min = time_max ! - t_i

        do q=1,matrix_size
            if(t_step(q) < Delta_t_min)then
                Delta_t_min = t_step(q)
                q_final=q
            endif
        enddo

        if(Delta_t_min<1.d-13)then
            write(*,*)'routine Achtung Delta_t_min<0!!!',Delta_t_min,q_max,log_rr(q_max)/W_max
            q_final=q_max
            Delta_t_min=log_rr(q_max)/W_max

        endif

    End Subroutine Advance_time_squares


    !------------------------------------------------------------

    subroutine MonteCarlo_sweep(PP,NN)

        REAL (Kind=8)   :: PP(number_of_spins*matrix_size,N_time_slot),NN(N_time_slot)
        REAL (Kind=8)	:: h_ext_x,h_ext_y,h_ext_z

        time = 0.d0
        Delta_t = time_max -time

        h_ext_x = dsin(theta)*dcos(phi)*field_values(1)
        h_ext_y = dsin(theta)*dsin(phi)*field_values(1)
        h_ext_z = dcos(theta)*field_values(1)

        !--- initialization
        call Init_ferro
        !-----------------------------------
        i=1
        time = 0.d0
        time_old = 0.d0
        itime = 1
        itime_old = 1
        icount = 0


        do while (i < number_of_spins*10 )	! A condition which is always fulfilled!

            ! the spin which would make the transition within the shortest tiem is the one which is actually flipped => i = ir

            iadvance=0
            time_left = time_max -time
            Delta_t = time_left

            do ir=1,number_of_spins

                call Advance_time_squares(ir,time,qr,Delta_t_local)

                if(Delta_t_local < Delta_t)then
                    iadvance = iadvance+1
                    Delta_t=Delta_t_local
                    i=ir
                    q_i=qr
                endif

            enddo

            !  exiting at this point does not seem to have an influence
            if(iadvance==0)Delta_t = time_max + 1.	!no event has happened

            if(Delta_t<1.d-13)write(*,*)time,'achtung Delta_t=0!',Delta_t,'iadvance=',iadvance
            ! Advance time and define probability

            time_old = time
            itime_old = itime
            time = time + Delta_t

            tmp_h = h_i + time*h_step/t_inc
            itime=int(time/t_inc)  + 1   			! N.B. Delta_t is allowed to be larger than t_inc with this trick

            if(itime == itime_old) then
                icount = icount + 1
                NN(itime) = NN(itime) + Delta_t 			! it has to be outside the loop on the lattice

                do ir=1,number_of_spins
                    p = conf(ir)	! 1 + (sigma(ir) + 1)/2			! vector of the states in which the each spin is found
                    ip = (ir-1)*matrix_size + p 					! label site + level before flipping
                    PP(ip,itime) = PP(ip,itime) + Delta_t		    ! probabilities are stored with global time
                enddo

            else

                if(mod(iter,iwrite)==0.and.mod(itime,10)==0)write(*,*)'iter=',iter,'itime=', itime, icount

                icount = 0

                dt_f = itime_old*t_inc - time_old			! the slice of delta t in the final time slot

                NN(itime_old) = NN(itime_old) + dt_f	    ! it has to be outside the loop on the lattice

                it_max=itime						        ! used for updating intermediate time slots
                if(itime > N_time_slot)then
                    it_max = N_time_slot
                endif

                do ir=1,number_of_spins

                    p = conf(ir)
                    ip = (ir-1)*matrix_size + p

                    PP(ip,itime_old) = PP(ip,itime_old) + dt_f

                    !	update intermediate time slots
                    do it = itime_old+1,it_max
                        dt = time - (it-1)*t_inc
                        if (dt>t_inc)dt = t_inc
                        PP(ip,it) = PP(ip,it) + dt
                        if(ir==1)NN(it) = NN(it) + dt
                    enddo

                enddo ! ir loop to define

            endif ! close (itime == itime_old)


            ! flip the spin if not time > time_max

            !-------------------------------------------
            if (time < time_max)then
                conf(i) = q_i
            else
                exit ! do while
            endif

            !	END UPDATE CONFIGURATION
            !-------------------------------------------

            if (itime > N_time_slot)exit	! do while =>  this happens also when Delta_t = 0, the following one not

        enddo   ! do while

    end subroutine MonteCarlo_sweep

    !------------------------------------------------------------

    subroutine write_current_averages

        open(unit=1, file=filename_1)
        if(read_fields) open(unit=2, file=filename_fields)

        !----------- output header ------------------
        if ( average_critical_field )then
            !   if ( icritical == 0)then
            !      write(1,*)
            !      write(1,*)
            !     write(1,*)'# H_c,rit H nominal ', 'H true ' ,'M stoch ','M eq ','time'
            !    else
            do i_skip = 1, (N_time_slot + 2)*icritical
                read(1,*)
            enddo
            !    endif  ! icritical loop

        else
            write(1,*)'# H nominal ', 'H true ' ,'M stoch ','M eq ','time', 'itime'
        endif
        !--------------------------------------------


        time_check = 0.d0
        do itime=1,N_time_slot

            if( read_fields ) then
                read(2,*)ih, field_nominal, field_true
                field_values(itime)=field_true
            elseif ( average_critical_field )then
                field_nominal = dble(itime-1)*h_step + h_i + h_step*0.5
                field_true =  field_values(itime)
            else
                field_true = field_values(itime)
                field_nominal = field_true
            endif

            Population(:,:)=0.d0

            do ir=1,number_of_spins
                if(Norm(itime)>1.d-14)then
                    do p=1,matrix_size
                        ip = (ir-1)*matrix_size + p !label site + level
                        Population(ir,p)=Prob(ip,itime)/Norm(itime)     ! needed in dipolar_fields
                    enddo
                endif
            enddo ! lattice sweep

            h_ext_x = dsin(theta)*dcos(phi)*field_values(itime)
            h_ext_y = dsin(theta)*dsin(phi)*field_values(itime)
            h_ext_z = dcos(theta)*field_values(itime)

            Mag_eq_average = 0.d0
            Mag(:) = 0.d0
            Proj(:) = 0.d0
            icount=0

            Population_states(:)=0.d0

            do ir=1,number_of_spins

                H_x = h_ext_x
                H_y = h_ext_y
                H_z = h_ext_z

                call transition_rates(matrix_size,T,H_x,H_y,H_z,W_single,En_levels,S_proj)
                !                        call transition_rates_H(matrix_size,T,H_x,H_y,H_z,W_single,En_levels,S_proj)

                if(Norm(itime)>1.d-14)then
                    icount=icount+1
                    do p=1,matrix_size
                        Proj(1) = Proj(1) - S_prog_real(S_proj(1,p))*Population(ir,p)
                        Proj(2) = Proj(2) - S_prog_real(S_proj(2,p))*Population(ir,p)
                        Proj(3) = Proj(3) - S_prog_real(S_proj(3,p))*Population(ir,p)

                        Population_states(p) = Population_states(p) + Population(ir,p)
                    enddo
                endif 		!Norm condition

                time_check = time_check + Norm(itime)

                call mag_equilibrium_z		! actung it is computed locally

                Mag_eq_average = Mag_eq_average +  Mag_eq

            enddo ! ir loop

            Mag_eq_average = Mag_eq_average/dble(number_of_spins)

            Proj(1)=Proj(1)/(S*real(icount, kind=8))
            Proj(2)=Proj(2)/(S*real(icount, kind=8))
            Proj(3)=Proj(3)/(S*real(icount, kind=8))

            Mag_stoch = dsin(theta)*dcos(phi)*Proj(1)
            Mag_stoch = Mag_stoch + dsin(theta)*dsin(phi)*Proj(2)
            Mag_stoch = Mag_stoch + dcos(theta)*Proj(3)

            if (average_critical_field) then
                write(1,*) field_critical, field_nominal, field_true, Mag_stoch, Mag_eq_average, t_inc*itime,itime
            else
                write(1,*)  field_nominal, field_true, Mag_stoch, Mag_eq_average, t_inc*itime,itime
            endif

            if(itime==1)write(*,*)icritical ,t_inc*itime,field_values(itime),Proj(3),Mag_stoch,Mag_eq_average


        enddo ! itime

        if(read_fields)close(2)

        write(*,*)'time_check =',time_check

        if (average_critical_field) then
            write(1,*)
            write(1,*)
        else

            !----------- output ending ------------------
            write(1,*)'# Source code: Kinetic_MC_Fe4.f90'
            write(1,*)'# computation parameters'
            write(1,*)'#'
            write(1,*)'# Spin Hamiltonian parameters'
            write(1,*)'# D,E,B,B42,C,B43,B66'
            write(1,*)'#',D,E,B,B42,C,B43,B66
            write(1,*)'#'
            write(1,*)'# number_of_spins=',number_of_spins
            write(1,*)'# iter_max=',iter_max,' iter=',iter
            write(1,*)'# iwrite=',iwrite
            write(1,*)'# N_time_slot=',N_time_slot
            write(1,*)'# time_max=',time_max
            write(1,*)'# h_i=',h_i
            write(1,*)'# h_f=',h_f
            ! write(1,*)'# p_initial=',p_initial
            write(1,*)'# T=',T
            write(1,*)'# theta=',theta/rad
            write(1,*)'# phi=',phi/rad
            write(1,*)'# Gamma_0 (rkk) =',Gamma_0
            write(1,*)'# t_inc=',t_inc,'# h_step=',h_step
            write(1,*)'# v sweep=',h_step/t_inc
            !--------------------------------------------

        endif ! write ending

        flush(1)
        close(1)


    end subroutine write_current_averages


end program Kinetic_MC_Fe4




