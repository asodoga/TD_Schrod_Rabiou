#!/bin/bash
#
#make clean
make

nb1=$1
nb2=$2
nq1=$((nb1+10))
nq2=$((nb2+10))
echo $nb1 $nb2


mkdir -p results_H_nb$nb1$nb2
cd results_H_nb$nb1$nb2

../TD_SCHROD.x << ** > resultat_H.lua

&potential
  pot_name='TwoD_RJDI2014' ! potential surface name
  adiabatic=f ! diabatic basis
  nsurf = 2
  ndim = 2
  read_nml=t
/

&TwoD_RJDI2014 wX=1. wY=1. a=1. c=0.1 Delta=1. /

&basis_nd name='dp' nb_basis=3/
  &basis_nd name='HO' nb =$nb1 nq=$nq1 Q0 =1.0 SCALEQ = 1.0 Imp_k=0.0  alpha =( 1.,0.0) /
  &basis_nd name='HO' nb =$nb2 nq=$nq2 Q0 =0.0 SCALEQ = 1.0 Imp_k=0.0  alpha =( 1.,0.0) /
  &basis_nd name= 'el' nb=2 /

&defGWP ndim=2  Elecindex = 1 Coef=(1.,0.)/
  &defWP0  sigma=  1.4142135623730951  Beta=0.0  Qeq= 1.0   imp_k=0.0 gamma=0.0 /
  &defWP0  sigma=  1.4142135623730951  Beta=0.0  Qeq= 0.0   imp_k=0.0 gamma=0.0 /

&prop t0=0.0  tf=250. delta_t=0.1 eps=1.e-20
      max_iter=500   propa_name='hagedorn'  propa_name2='taylor' Beta=T  P=T renorm=T /
*