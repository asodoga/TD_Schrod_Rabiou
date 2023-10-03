#FC=gfortran
#FC=ifort
#FFLAGS=-O3 -Wall -Wextra -fopenmp -J$(MOD_DIR)
FFLAGS= -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan -fopenmp -J$(MOD_DIR)
OBJ_DIR=obj
MOD_DIR=obj
SRC_DIR=src
MAIN_DIR=app


QMLDIR= /home/rabiou_issa/rabiou.issa/TD_Schrod_Rabiou/EXT_Lib/QuantumModelLib
QDUDIR= /home/rabiou_issa/rabiou.issa/TD_Schrod_Rabiou/EXT_Lib/QDUtilLib
MODEXT =/home/rabiou_issa/rabiou.issa/TD_Schrod_Rabiou/EXT_Lib/QDUtilLib/OBJ/obj_gfortran_opt1_omp1_lapack1

FFLAGS += -I$(MODEXT)
#QMLDIR=/Users/lauvergn/git/QuantumModelLib

MAIN=TD_SCHROD

# source files ============================================================================================================================

LIBSRC=   Molec_m.f90 NDindex_m.f90 poly0rtho_m.f90 Basis_m.f90 psi_m.f90 Ana_psi_m.f90 lanczos_m.f90 Op_m.f90 param_WP0_m.f90 Auto_corr_m.f90 Propa_m.f90

# construct file.o frome file.f90 ==========================================================================================================

OBJ0=${LIBSRC:.f90=.o}

OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
$(info ************ $(OBJ_DIR)/files: $(OBJ))


# libriries====================================================================================================================
LIB =$(QMLDIR)/libQMLibFull_gfortran_opt1_omp1_lapack1.a  $(QDUDIR)/libQD_gfortran_opt1_omp1_lapack1.a 
# target and dependancy===========================================================================================================


$(MAIN).x: $(OBJ_DIR)/$(MAIN).o $(OBJ) $(LIB)
	$(FC) $(FFLAGS) -o $(MAIN).x  $(OBJ_DIR)/$(MAIN).o $(OBJ) $(LIB) -llapack -lblas

$(OBJ_DIR)/$(MAIN).o: $(MAIN_DIR)/$(MAIN).f90 
	$(FC) $(FFLAGS) -o $@ -c $< 

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 
	$(FC) $(FFLAGS) -o $@ -c $< 


clean:
	rm -f $(OBJ_DIR)/*.o $(MOD_DIR)/*.mod

#QML and QDU
$(QMLDIR)/libpot.a:
	test -e $(QMLDIR)/libpot.a || exit 1
$(QDUDIR)/libpot.a:
	test -e $(QDUDIR)/libutil.a || exit 1

#=============================================================================================== 
 #============= module dependencies ============================================================
#===============================================================================================

$(OBJ_DIR)/$(MAIN).o: $(OBJ)

$(OBJ_DIR)/Basis_m.o: $(OBJ_DIR)/NDindex_m.o  $(OBJ_DIR)/poly0rtho_m.o

$(OBJ_DIR)/psi_m.o: $(OBJ_DIR)/Basis_m.o $(OBJ_DIR)/param_WP0_m.o

$(OBJ_DIR)/Ana_psi_m.o: $(OBJ_DIR)/psi_m.o

$(OBJ_DIR)/Auto_corr_m.o: $(OBJ_DIR)/psi_m.o $(OBJ_DIR)/poly0rtho_m.o

$(OBJ_DIR)/Op_m.o: $(OBJ_DIR)/Molec_m.o $(OBJ_DIR)/Basis_m.o $(OBJ_DIR)/psi_m.o

$(OBJ_DIR)/lanczos_m.o: $(OBJ_DIR)/Op_m.o $(OBJ_DIR)/psi_m.o  

$(OBJ_DIR)/Propa_m.o: $(OBJ_DIR)/Auto_corr_m.o $(OBJ_DIR)/lanczos_m.o $(OBJ_DIR)/Op_m.o $(OBJ_DIR)/Ana_psi_m.o $(OBJ_DIR)/psi_m.o $(OBJ_DIR)/Basis_m.o

