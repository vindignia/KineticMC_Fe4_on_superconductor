module ComputationParameters

    Implicit none
    COMMON S,GI,D,E,B,B42,C,B43,B66,gamma_0,gamma_tunnel,g1,g2
    ! ----------------- spin-Hamiltonian parameters
    REAL (Kind=8)       		    :: GI = 2.d0
    REAL (Kind=8)                   :: D = -0.6d0
    REAL (Kind=8)                   :: E = 2.d-2
    REAL (Kind=8)                   :: B = 0.d0
    REAL (Kind=8)                   :: B42 = 0.d0
    REAL (Kind=8)                   :: C = 0.d0
    REAL (Kind=8)                   :: B43 = 5.d-3
    REAL (Kind=8)                   :: B66 = 3.d-6
    ! ----------------- Transition-rate parameters -----------
    REAL (Kind=8)                   :: gamma_0 = 4.7987d3
    REAL (Kind=8)                   :: gamma_tunnel = 0.5d0
    REAL (Kind=8)                   :: g1 = 4.d-2
    REAL (Kind=8)                   :: g2 = 4.d-2
    ! ----------------- System parameters --------------------
    INTEGER (Kind=4),parameter      :: number_of_spins = 10		! spin can potentially have a different tilting
    INTEGER (Kind=4),parameter      :: matrix_size = 11		    ! 2*S+1 we define this instead of the spin S
    REAL (Kind=8)                   :: S
    ! --------------------------------------------------------
    INTEGER (Kind=4),parameter      :: iter_max = 200           ! Number of MC sweeps, iterations
    INTEGER (Kind=4) 	            :: iwrite = 20              ! partial averages are written in the output file every "iwrite" iterations
    INTEGER (Kind=4),parameter	    :: N_time_slot = 400        ! number of slots in which the elapsed time is discretized
    REAL (Kind=8),parameter		    :: time_max = 0.3*dble(N_time_slot)     ! time elapsed during the simulation in seconds
    REAL (Kind=8)	 			    :: h_step
    REAL (Kind=8),parameter	        :: h_i = -0.6d4			    ! initial field [Gauss]
    REAL (Kind=8),parameter		    :: h_f = 0.6d4			    ! final field [Gauss]
    INTEGER (Kind=4),parameter      :: seed = 980759715 		! random numbers seed
    INTEGER (Kind=4),parameter      :: p_initial = 1			! Initial state for each spin
    REAL (Kind=8),parameter         :: T = 4.d-1
    REAL (Kind=8),parameter         :: pi = 3.141592653589793238462643d0
    REAL (Kind=8),parameter         :: mu_B = (9.274078d-5)/1.380662d0  	! Bohr magneton [Kelvin/Gauss]
    REAL (Kind=8),parameter         :: eps = 1.d-13
    REAL (Kind=8),parameter         :: rad = pi/180.d0
    REAL (Kind=8)			 	    :: theta_0 = 0.d0			! tilting easy axis w.r.t. normal to the plane in deg
    REAL (Kind=8)	                :: phi_0 = 0.d0
    REAL (Kind=8)                   :: h_crit_Pb = 7.d2         ! critical field of the superconductor
    ! ----------------- logical parameters --------------------
    logical                         :: save_eig_and_W=.false.
    logical 					    :: random_phi=.false.
    logical 					    :: heat_bath=.false.
    logical 					    :: read_fields=.false.
    logical 					    :: average_critical_field =.false.
contains

    ! --------------------------------------------------------

    subroutine show_consts()
        print*, "Pi = ", pi
        print*,  "e = ", e
    end subroutine show_consts

end module ComputationParameters