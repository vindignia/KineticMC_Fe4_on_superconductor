# Kinetic MonteCarlo for Fe4 molecules on a superconductor
### © 2020 ETH Zurich, [PD Dr. Alessandro Vindigni]
This project computes hysteresis curves for Fe4 molecules deposited on a superconducting substrate. 
For details, please refer to the manuscript "Quantum dynamics of a single molecule magnet on superconducting Pb(111)", G. Serrano et al., Nat. Mater. (2020). This branch is a copy of the master at the time the article was published. 
 
The present project can be compiled through the Makefile, with the command line
> make   

To execute the program type
> ./Fe4

The project is compatible with the GNU Fortran (GCC) 8.1.0 (or later versions) and uses the lapack libraries. 

The names of output files can be defined in the "FileNamesModule.f90". 
Computation parameters can be set in the module "ParametersModule.f90". 
Referring to the Supplementary Note 8 associated with the manuscript mentioned above, the meaning of computation parameters should be clear. 
Differently from the manuscript the spin-phonon couplings are here named g1 and g2; while gamma_tunneling represents the transition rate due to pure tunnelling between almost degenerate levels. The meaning of other computation parameters is explained in the module "ParametersModule.f90" itself.  

Changes of file names or parameters become effective only after having compiled the project again. 
