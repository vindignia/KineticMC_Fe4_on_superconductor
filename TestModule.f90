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
            ! ACHTUNG field_values(itime)
            H_x = dsin(theta)*dcos(phi)*field_values(itime)
            H_y = dsin(theta)*dsin(phi)*field_values(itime)
            H_z = dcos(theta)*field_values(itime)

            call transition_rates(matrix_size,T,H_x,H_y,H_z,W_single,eig,S_proj)
            write(1,*)field_values(itime),eig(1),eig(11),eig(2),eig(10),eig(3),eig(9),eig(4),eig(8), &
                    eig(5),eig(7),eig(6)
!           transition (+5 -> -5), (+5 -> +4), (+5 -> -4)
            write(2,*)field_values(itime),W_single(1,11),W_single(11,1)

        end do
    end subroutine LabelsEnergyLevels

End Module Tests