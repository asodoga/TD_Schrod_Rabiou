FC=gfortran
#FFLAGS=-O3 -Wall -Wextra -fopenmp -J$(MOD_DIR)
FFLAGS= -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan -fopenmp -J$(MOD_DIR)
OBJ_DIR=obj
MOD_DIR=obj
SRC_DIR=src
MAIN_DIR=app


QMLDIR=/home/elprof/QuantumModelLib
AD_dnSVMDIR=/home/elprof/QuantumModelLib
#QMLDIR=/Users/lauvergn/git/QuantumModelLib
#AD_dnSVMDIR=/Users/lauvergn/git/QuantumModelLib

MAIN=TD_SCHROD

LIBSRC= NumParameters_m.f90 UtilLib_m.f90 diago_m.f90 \
 Molec_m.f90 NDindex_m.f90 Basis_m.f90 psi_m.f90 Ana_psi_m.f90 lanczos_m.f90 Op_m.f90 param_WP0_m.f90 Propa_m.f90
OBJ0=${LIBSRC:.f90=.o}

OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
$(info ************ $(OBJ_DIR)/files: $(OBJ))

#LIB = -lblas -llapack $(AD_dnSVMDIR)/libAD_dnSVM.a $(QMLDIR)/libpot.a
LIB =$(QMLDIR)/libpot.a $(AD_dnSVMDIR)/libAD_dnSVM.a -llapack -lblas
$(MAIN).x: $(OBJ_DIR)/$(MAIN).o $(OBJ) $(QMLDIR)/libpot.a $(AD_dnSVMDIR)/libAD_dnSVM.a
	$(FC) $(FFLAGS) -o $(MAIN).x  $(OBJ_DIR)/$(MAIN).o $(OBJ) $(LIB)

$(OBJ_DIR)/$(MAIN).o: $(MAIN_DIR)/$(MAIN).f90
	$(FC) $(FFLAGS) -o $@ -c $<

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -o $@ -c $<


clean:
	rm -f $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod

#QML and AD_dnSVM
$(QMLDIR)/libpot.a:
	test -e $(QMLDIR)/libpot.a || exit 1
$(AD_dnSVMDIR)/libAD_dnSVM.a:
	test -e $(AD_dnSVMDIR)/libAD_dnSVM.a || exit 1


#dependencies
$(OBJ_DIR)/$(MAIN).o: $(OBJ)

$(OBJ_DIR)/Basis_m.o: $(OBJ_DIR)/NDindex_m.o $(OBJ_DIR)/UtilLib_m.o $(OBJ_DIR)/diago_m.o $(OBJ_DIR)/NumParameters_m.o
$(OBJ_DIR)/Molec_m.o: $(OBJ_DIR)/UtilLib_m.o $(OBJ_DIR)/diago_m.o $(OBJ_DIR)/NumParameters_m.o
$(OBJ_DIR)/param_WP0_m.o: $(OBJ_DIR)/UtilLib_m.o

$(OBJ_DIR)/psi_m.o: $(OBJ_DIR)/Basis_m.o $(OBJ_DIR)/param_WP0_m.o

$(OBJ_DIR)/Ana_psi_m.o: $(OBJ_DIR)/psi_m.o

$(OBJ_DIR)/Op_m.o: $(OBJ_DIR)/Molec_m.o $(OBJ_DIR)/Basis_m.o $(OBJ_DIR)/psi_m.o

$(OBJ_DIR)/lanczos_m.o: $(OBJ_DIR)/Op_m.o $(OBJ_DIR)/psi_m.o $(OBJ_DIR)/Basis_m.o $(OBJ_DIR)/diago_m.o
$(OBJ_DIR)/Propa_m.o: $(OBJ_DIR)/lanczos_m.o $(OBJ_DIR)/Op_m.o $(OBJ_DIR)/Ana_psi_m.o $(OBJ_DIR)/psi_m.o $(OBJ_DIR)/Basis_m.o

$(OBJ_DIR)/NDindex_m.o: $(OBJ_DIR)/UtilLib_m.o $(OBJ_DIR)/NumParameters_m.o
$(OBJ_DIR)/diago_m.o:   $(OBJ_DIR)/NumParameters_m.o
$(OBJ_DIR)/UtilLib_m.o: $(OBJ_DIR)/NumParameters_m.o