#*************************************************************************#
#                                                                         #
#  Program: Makefile - ONE SIDED MODEL                                    #
#  Version: 3.0                                                           #
#  By: Mette Olufsen  									  #
#  Edited by: Alyssa Taylor, Mitchel Colebank                             #
#  Date: March 1st, 202                                                   #
#                                                                         # 
#  A makefile that ensures that all modules are linked together in the    #
#  right order.                                                           #
#*************************************************************************#

# The Makefile requires using the GNU available software "make", as 
# well as "gcc" and "gfortran".
#
# Users can download the make software from https://www.gnu.org/software/make/
# For gcc: https://gcc.gnu.org/install/download.html
# For gfortran: https://gcc.gnu.org/wiki/GFortranBinaries

# When you download these programs, you will have to identify which directory they are stored to.

# For the C++ compiler, edit the following line
#CXX=/PATH/TO/YOUR/g++
CXX =C:\cygwin64\bin\g++ # For Colebank's computer

# For the Fortran compiler, edit the following line
#FC=/PATH/TO/YOUR/gfortran
FC=/usr/bin/gfortran # For Colebank's computer
####################################################
# Do not change anything from here on
####################################################
CXXFLAGS=-O2 -g -Wall -D_REENTRANT -fPIC
FFLAGS=-O2 -g -Wall
FLIBS=-lgfortran -lquadmath


LIBS=$(FLIBS) -lm

LDFLAGS=-O2

OBJS1=tools.o junction.o sor06.o arteries.o
OBJS2=impedance_sub.o impedance_init_sub.o root_imp.o f90_tools.o

MAIN=sor06

all: $(MAIN)

$(MAIN): $(OBJS1) $(OBJS2) 
	$(CXX) -o $(MAIN) $(LDFLAGS) $(OBJS1) $(OBJS2) $(LIBS)
	
sor06.o: sor06.c sor06.h
	$(CXX) -c $(CXXFLAGS) sor06.c

junction.o: junction.c junction.h arteries.h tools.h
	$(CXX) -c $(CXXFLAGS) junction.c
	
arteries.o: arteries.c arteries.h junction.h tools.h sor06.h
	$(CXX) -c $(CXXFLAGS) arteries.c
	
tools.o: tools.c tools.h
	$(CXX) -c $(CXXFLAGS) tools.c
		
root_imp.o: root_imp.f90 f90_tools.o
	$(FC) -c $(FFLAGS) root_imp.f90
	
f90_tools.o: f90_tools.f90
	$(FC) -c $(FFLAGS) f90_tools.f90
	
impedance_sub.o: impedance_sub.f90 f90_tools.o root_imp.o
	$(FC) -c $(FFLAGS) impedance_sub.f90
	
impedance_init_sub.o: impedance_init_sub.f90 f90_tools.o root_imp.o
	$(FC) -c $(FFLAGS) impedance_init_sub.f90
	
clean:
	-rm -f *.o *.mod *.2d Zhat*
	
veryclean: clean
	-rm $(MAIN) a.out *~ STnetwork*