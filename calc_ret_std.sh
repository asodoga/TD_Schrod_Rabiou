#!/bin/bash
#
#make clean
make

nb1=$1
nb2=$2

nq1=$nb1
nq2=$nb2
echo " Parameter for the first Basis" $nb1 $nq1
echo "Parameter for the second Basis" $nb2 $nq2

mkdir -p results_std_nb$nb1$nb2
cd results_std_nb$nb1$nb2

../TD_SCHROD.x << ** > resultat_ret_std.lua
&potential
    pot_name  = 'Retinal_JPCB2000'
    option=1
    adiabatic = f
    nsurf = 2
    ndim = 2
   /
   &basis_nd name='dp' nb_basis=3/
   &basis_nd name='fourier' nb=$nb1 nq=$nq1 SCALEQ=1.0  Imp_k=0.0  alpha =( 1.0,0.0)/
    &basis_nd name='HO' nb=$nb2 nq=$nq2  SCALEQ= 0.9592722 Imp_k=0.0  alpha =( 0.9202031537,0.00)  /
      &basis_nd name= 'el' nb=2/

 &defGWP ndim=2  Elecindex=2 Coef=(1.,0.)/
 &defWP0  sigma= 0.181095798 Beta=0.0 Qeq=0.0   imp_k=0.0 gamma=0.0 /
 &defWP0  sigma=1.474256632  Beta=0.0 Qeq=0.0   imp_k=0.0 gamma=0.0 /

&prop t0  = 0  tf  = 10000 delta_t =10.0  eps=1.e-20
      max_iter=500   propa_name='non_hagedorn'  propa_name2='taylor' Beta=T  P=T renorm=T  /
**
