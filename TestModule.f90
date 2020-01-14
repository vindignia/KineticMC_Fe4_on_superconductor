Module Tests
    USE SpinAlgebra
    USE TransitionProbabilities
    USE FileNames

Contains

    subroutine LabelsEnergyLevels(matrix_size, T,theta,phi,N_time_slot,field_values)

        implicit none
        INTEGER (Kind=4) :: itime,N_time_slot, matrix_size,p,q
        REAL (Kind=8) :: H_x,H_y,H_z,field_values(N_time_slot)
        REAL (Kind=8) :: eig(matrix_size),W_single(matrix_size,matrix_size)
        REAL (Kind=8) :: theta, phi, T
        double precision :: En_levels(matrix_size)
        COMPLEX (Kind=8) :: S_proj(3,11) = (0.d0,0.d0)

        do itime=1,N_time_slot

            H_x = dsin(theta)*dcos(phi)*field_values(itime)
            H_y = dsin(theta)*dsin(phi)*field_values(itime)
            H_z = dcos(theta)*field_values(itime)	!ACHTUNG field_values(itime)

            call transition_rates_Ale(matrix_size,T,H_x,H_y,H_z,W_single,eig,S_proj)
            write(1,*)field_values(itime),eig(1),eig(11),eig(2),eig(10),eig(3),eig(9),eig(4),eig(8), &
                    eig(5),eig(7),eig(6)
!           old ones
!            call transition_rates_H(matrix_size,T,H_x,H_y,H_z,W_single,En_levels,S_proj)
!            write(1,*)field_values(itime),En_levels(1),En_levels(11),En_levels(2),En_levels(10),En_levels(3), &
!            En_levels(9),En_levels(4),En_levels(8),En_levels(5),En_levels(7),En_levels(6)

!           transition (+5 -> -5), (+5 -> +4), (+5 -> -4)
            write(2,*)field_values(itime),W_single(1,11),W_single(11,1)
            !,W_single(1,4),W_single(1,5), &
            !        W_single(1,6),W_single(1,7),W_single(1,8),W_single(1,9),W_single(1,10),W_single(1,11)

        end do
    end subroutine LabelsEnergyLevels

End Module Tests