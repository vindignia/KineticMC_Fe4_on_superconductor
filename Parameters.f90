module ComputationParameters

    Implicit none
    COMMON S,GI,D,E,B,B42,C,B43,B66,gamma_0,gamma_tunnel,g1,g2
    ! ----------------- spin-Hamiltonian parameters
    REAL (Kind=8)       		    :: GI = 2.d0  ! before g was a parameter
    REAL (Kind=8)                   :: D = -0.6d0
    REAL (Kind=8)                   :: E = 2.d-2
    REAL (Kind=8)                   :: B = 0.d0
    REAL (Kind=8)                   :: B42 = 0.d0
    REAL (Kind=8)                   :: C = 0.d0
    REAL (Kind=8)                   :: B43 = 5.d-3
    REAL (Kind=8)                   :: B66 = 3.d-6
    ! ----------------- Transition-rate parameters
    REAL (Kind=8)                   :: gamma_0 = 1.1563d3       ! 1156.3d-5		! Graphene
    REAL (Kind=8)                   :: gamma_tunnel = 0.5d0
    REAL (Kind=8)                   :: g1 = 4.d-2
    REAL (Kind=8)                   :: g2 = 4.d-2
    ! ----------------- system parameters
    INTEGER (Kind=4),parameter      :: number_of_spins = 10		! spin can potentially have a different tilting
    INTEGER (Kind=4),parameter      :: matrix_size = 11		    ! 2*S+1 we define this instead of the spin S
    REAL (Kind=8) :: S
    ! --------------------------------------------------------
    INTEGER (Kind=4),parameter      :: iter_max =200 !  50!0 !135 !999
    INTEGER (Kind=4) 	            :: iwrite=20
    INTEGER (Kind=4),parameter	    :: N_time_slot =400 ! 4000 !8000! 4800!2400!386!0
    REAL (Kind=8),parameter		    :: time_max=(1.d1)*dble(N_time_slot)/2.
    !REAL (Kind=8),parameter	    :: time_max=(7.5d0)*dble(N_time_slot)!real time slot
    REAL (Kind=8)	 			    :: h_step
    REAL (Kind=8),parameter	        :: h_i=-0.6d4			    ! Gauss
    REAL (Kind=8),parameter		    :: h_f= 0.6d4			    ! Gauss
    INTEGER (Kind=4),parameter      :: seed = 980759715 		! random numbers seed
    INTEGER (Kind=4),parameter      :: p_initial=1			    ! Initial state for each spin
    REAL (Kind=8),parameter         :: T=4.d-1
    REAL (Kind=8),parameter         :: pi=3.141592653589793238462643D0
    REAL (Kind=8),parameter         :: mu_B=(9.274078d-5)/1.380662d0  		! Kelvin/Gauss
    REAL (Kind=8),parameter         :: eps =1.d-13
    REAL (Kind=8),parameter         :: rad=pi/180.d0
    REAL (Kind=8)			 	    :: theta_0=22.d0			! transverse fields in deg
    REAL (Kind=8)	                :: phi_0=pi/2.
    REAL (Kind=8)                   :: h_crit_Pb = 7.d2

contains

    ! --------------------------------------

    subroutine show_consts()
        print*, "Pi = ", pi
        print*,  "e = ", e
    end subroutine show_consts

end module ComputationParameters