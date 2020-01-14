subroutine random_num_generator
!---- random number generator see file: random_number_check.f90 in Amade's folder 
integer                      	:: seedSize
integer                      	:: my_seed=9458
integer,dimension(:),allocatable:: seed,my_seed_array
integer, dimension(8)		:: dtVals
! 
 call DATE_AND_TIME(VALUES=dtVals)
 call RANDOM_SEED(SIZE=seedSize)
 allocate(my_seed_array(seedSize))
 allocate(seed(seedSize)) 
 call RANDOM_SEED(GET=seed)
 my_seed_array(1:8)=dtVals(1:8)
 my_seed_array(9:seedSize)=my_seed
 call RANDOM_SEED(PUT=my_seed_array)
 call RANDOM_SEED(GET=seed)
! 
end subroutine random_num_generator
