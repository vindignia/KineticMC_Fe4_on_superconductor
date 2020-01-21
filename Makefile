makefile_Fe4: Kinetic_MC_main.f90 ParametersModule.f90 FileNamesModule.f90
	gfortran Kinetic_MC_main.f90 rnd_subroutine.f90 -llapack  -lblas -o Fe4
