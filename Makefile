FC=gfortran
FFLAGS=-O3 -Wall -Wextra
#-- optional stuff for specifying dependencies automatically
QMLDIR=/home/elprof/QuantumModelLib
DEPFLAGS=-M -cpp
SRC= NumParameters_m.f90 \
 UtilLib_m.f90 \
 Molec_m.f90 NDindex_m.f90 Basis_m.f90 psi_m.f90 Op_m.f90 GWP1D_m.f90 GWPnD_m.f90 propa_m.f90 TD_SCHROD.f90
OBJ=${SRC:.f90=.o}
%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $<
	
TD_SCHROD: $(OBJ)
	$(FC) $(FFLAGS) -o $@ $(OBJ) $QMLDIR/libpot.a -lblas -llapack
clean:
	 rm *.o *.mod maths *.deps
	 
deps:$(SRC)
	 @$(FC) $(DEPFLAGS) $(SRC)>make.deps
make.deps: deps
#include make.deps
	 
