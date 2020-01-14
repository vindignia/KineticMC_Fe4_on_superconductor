module GlobalVariables
    USE ComputationParameters

    Implicit none
    REAL (Kind=8)				:: theta, phi
    REAL (Kind=8)				:: rndReal,rnd_time,rnd_phi,x,y,time_check
    REAL (Kind=8)				:: Mag_eq,Mag_eq_average,Mag_stoch
    REAL (Kind=8)				:: tmp,tmp_rnd,tmp_h,Delta_t,Delta_t_local,time,time_old,dt,dt_f,t_inc,time_left,rho
    REAL (Kind=8)				:: time_local(number_of_spins),Mag(number_of_spins),Proj(3)
    REAL (Kind=8)				:: field_values(N_time_slot) ! h_ext_x,h_ext_y,h_ext_z,
    REAL (Kind=8)				:: field_nominal, field_true, field_critical
    INTEGER (Kind=4)			:: i,q_i	!i=molecule that makes a transition q_i final state: do not use for other loops!!!!
    INTEGER (Kind=4)			:: iter,N_points,icount,j,p,q,ip,pp,itime,it,it_max,itime_old,i_rnd,iadvance,ir,qr
    INTEGER (Kind=4)			:: ih, icritical, i_skip
    REAL (Kind=8)				:: W(number_of_spins*matrix_size,matrix_size)
    REAL (Kind=8)				:: Prob(number_of_spins*matrix_size,N_time_slot),Norm(N_time_slot),dipolar_field_ref(N_time_slot)
    REAL (Kind=8)				:: W_single(matrix_size,matrix_size),Population(number_of_spins,matrix_size)
    REAL (Kind=8)				:: Population_states(matrix_size)
    REAL (Kind=8)				:: WW(matrix_size,matrix_size,N_time_slot)
    REAL (Kind=8)				:: En_levels(matrix_size)
    COMPLEX (Kind=8)			:: S_proj(3,matrix_size)	! it used to be REAL (Kind=8)
    REAL (Kind=8)				:: X_coord(number_of_spins),Y_coord(number_of_spins),Z_coord(number_of_spins)
    INTEGER (Kind=4),dimension(1:number_of_spins)		:: conf
    INTEGER (Kind=4)            :: idum,mmm
    REAL (Kind=8) 				:: H_x,H_y,H_z
    logical                     :: save_eig_and_W=.true.
    logical 					:: random_phi=.false.
    logical 					:: heat_bath=.false.
    logical 					:: read_fields=.false.
    logical 					:: average_critical_field =.false.
    REAL (Kind=8)				:: start, finish
    REAL (Kind=8)				:: time_exp,field_exp,mag_exp
    REAL (Kind=8)				:: h_dip_x(number_of_spins),h_dip_y(number_of_spins),h_dip_z(number_of_spins)
    REAL (Kind=8)				:: h_dip_old_x(number_of_spins),h_dip_old_y(number_of_spins),h_dip_old_z(number_of_spins)

end module GlobalVariables