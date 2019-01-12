########################################################################
# Simple Makefile designed to compile Test_Psolver.                    # 
# User only needs to indicate the path to lapack and BigDFT Psolver    #
# libraries.                                                           #
# To run in parallel: mpirun -np nproc Test_Psolver.x                  #
########################################################################
# PLT@CFM(2018)                                                        #
########################################################################
#
# Name of executable 
#
EXE = Test_Psolver_BigDFT.x
#
FC = mpif90 
#
FCFLAGS = -O2 -L/userdefined/lapack-3.4.2-gfortran/lib -llapack -lblas \
              -I/userdefined/bigdft-1.8.1/build_solver/install/include \
              -L/userdefined/bigdft-1.8.1/build_solver/install/lib \
              -lPSolver-1 -lfutile-1 -framework Accelerate -lyaml -ldl
#
#
# Objects, order matters:
#
OBJECTS = timestamp.o m_utils.o m_psolver.o main.o
#
# Makefile
#
$(EXE): $(OBJECTS)
	$(FC) -o $(EXE) $(FCFLAGS) $(OBJECTS)
#
# Simple subroutines: 
# (Order does not matter)
timestamp.o: timestamp.o 
	$(FC) -c $(FCFLAGS) timestamp.f90 
# Modules:
# (Order does not matter)
main.o: main.f90
	$(FC) -c $(FCFLAGS) main.f90
m_utils.mod: m_utils.o m_utils.f90
	$(FC) -c $(FCFLAGS) m_utils.f90 
m_utils.o: m_utils.f90 
	$(FC) -c $(FCFLAGS) m_utils.f90 
m_psolver.mod: m_psolver.o m_psolver.f90 
	$(FC) -c $(FCFLAGS) m_psolver.f90 
m_psolver.o: m_psolver.f90 
	$(FC) -c $(FCFLAGS) m_psolver.f90 
#
%.o: %.f90
	$(FC) -c $(FCFLAGS) $<
# Cleaning everything
clean:
	rm -f *.o *.x *.mod *.L *~
# End of the makefile
