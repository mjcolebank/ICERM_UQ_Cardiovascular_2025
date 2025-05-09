
This directory contains a combination of C++, FORTRAN and MATLAB files that simulate blood pressure and flow in the 
systemic blood circulation. 

COMPILING:
Instructions on how to compile the C++ / FORTRAN code can be found in the file Makefile_Instructions.pdf 
accompanied by the video Makefile_Instructions.mp4.

Data Files: 

newDORV2q.dat: Inflow data that drives the system of PDEs.

Code Files:

Makefile:     Compiles C++/FORTRAN files together
sor06:	      The executable file. (sor06.exe on Windows machines)
sor06.h:      Header file that initiates the basic, global parameter values.
sor06.c:      A main solver that calls the artery initialization and solves the system of PDEs.
arteries.h:   Header file which declares the class "Tube" and all of its objects and functions.
arteries.c:   Main computational code that computes all variables associated with every vessel in the
	      system. Includes calls to tools.c, tools.h, and arteries.h.
tools.h:      Header file that initiates the solver tools.
tools.c:      Contains solvers for the nonlinear systems, including Newton-Raphson method.
junction.c: Contains the Jacobian for the type of junction i.e. single vessels, bifurcation, trifurcation, and converging vessels
junction.h: Header file which declares all of the arguments needed to define a junction point for all cases
impedance_init_sub.f90: Defines the initial impedance at the end of the structured tree
impedance_sub.f90: Calls the function to determine the impedance at the root of the structured tree using recursive relations
root_imp.f90: Contains subroutines and functions needed to compute the impedance at the root of a structured tree
f90_tools.f90: Contains subroutines and functions needed to compute the fast Fourier transform and to flip a vector up down
run_1D.m:	      MATLAB file that calls the C++ code via the "unix" command. Passes parameter values from
              MATLAB into sor06.C.
gnuplot.m:    MATLAB file that reads the output files and converts them to the appropriate format for plotting.


Output Files:
output_ID.2d:    The pressure, flow and area for file with index ID.


License:
This project is licensed under the MIT License -  see the LICENSE.txt file for details.
