 begining: max threads?           4
 pot_name,option: read_model           1
 ndim,nsurf     :            0           0
 =================================================
 =================================================
 == QML: Quantum Model Lib (E-CAM) ===============
 == QML version:       8.5                 
 == QML path:          /home/elprof/QuantumModelLib
 -------------------------------------------------
 == Compiled on       "elprof-Latitude-E6440" the jeu. 30 juin 2022 - 10:35:30
 == Compiler:         gfortran
 == Compiler version: GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0
 == Compiler options: -O5 -g -fbacktrace -fopenmp -funroll-loops -ftree-vectorize -falign-loops=16
 == Compiler libs:     -llapack -lblas
 -------------------------------------------------
 QML is under GNU LGPL3 license and 
   is written by David Lauvergnat [1]
   with contributions of
      Félix MOUHAT [2]
      Liang LIANG [3]
      Emanuele MARSILI [1,4]

 [1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
 [2]: Laboratoire PASTEUR, ENS-PSL-Sorbonne Université-CNRS, France
 [3]: Maison de la Simulation, CEA-CNRS-Université Paris-Saclay,France
 [4]: Durham University, Durham, UK
 =================================================
 =================================================
 == Initialization of the Model ==================
 Reading input file . . .
 Non-adiabatic potential . . .
 pot_name_loc: morse
 Non-adiabatic potential . . .
 =================================================
  Quantum Model
 -----------------------------------------------
 Output file for potential library
 Morse current parameters:

     V(R) = D.( 1 - exp(-a.(r-req)) )^2
   D:     0.22500000000000001     
   a:      1.1740999999999999     
   req:    1.7329000000000001     

 end Morse current parameters

 d0GGdef: 1        0.000573196

  Q0:       1.7329000000

 -----------------------------------------------
 Extra action(s):
 opt F
 irc F
 -----------------------------------------------
 =================================================
 =================================================
 ndim,nsurf           1           1
 pot_nameread_model      
 pot_nameread_model      
&BASIS_ND
 NAME="dp                  ",
 NB_BASIS=2          ,
 NB=0          ,
 NQ=0          ,
 A=  0.0000000000000000     ,
 B=  0.0000000000000000     ,
 SCALEQ=  1.0000000000000000     ,
 Q0=  0.0000000000000000     ,
 /
&BASIS_ND
 NAME="HO                  ",
 NB_BASIS=0          ,
 NB=5          ,
 NQ=5          ,
 A=  0.0000000000000000     ,
 B=  0.0000000000000000     ,
 SCALEQ=  1.0000000000000000     ,
 Q0=  0.0000000000000000     ,
 /
 Sii-1,Sij   3.1121771826292388E-012   1.0765240210083283E-016

&BASIS_ND
 NAME="el                  ",
 NB_BASIS=0          ,
 NB=1          ,
 NQ=0          ,
 A=  0.0000000000000000     ,
 B=  0.0000000000000000     ,
 SCALEQ=  1.0000000000000000     ,
 Q0=  0.0000000000000000     ,
 /
 Electronic basis. Electronic state number:           1
 This routine is .not. possible Basis el
 Initialization of a complex psi
 *************************************************
 *************************************************
 Q0 =   0.0000000000000000     
 K=    0.0000000000000000     
 phase=  0.28599999999999998     
 DQ=   1.4142135623700001     
 ndim =           1
 *************************************************
 *************************************************
  The basis is linked to psi.
 BEGINNING GridTOBasis_Basis
 END GridTOBasis_Basis
 BEGINNING GridTOBasis_Basis
 END GridTOBasis_Basis
  E = <Psi | H | Psi> =   47.734172888515353      <Psi | Psi> =   1.7724538508906056     
