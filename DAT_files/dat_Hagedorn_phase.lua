
&potential
    pot_name  = 'Retinal_JPCB2000'
    option=1
   /
&basis_nd name='dp' nb_basis=3/
   &basis_nd name='HO' nb =3 nq=5 Q0 =0.0 Imp_k=0.0  alpha =( 1.,0.0) /
   &basis_nd name='HO' nb =3 nq=5 Q0 =2.0 Imp_k=0.0  alpha =( 1.,0.0) /
  
      &basis_nd name= 'el' nb=2/
 &defGWP ndim=2  Elecindex = 2 Coef=(1.,0.)/
 &defWP0  sigma=   1.41421356237  Qeq= 0.0   imp_k=0.0 phase=0.0 /
 &defWP0  sigma=   1.41421356237  Qeq= 2.0   imp_k=0.0 phase=0.0 /
 

   &prop t0  = 0.0  tf  = 6.0 delta_t =0.1 eps= 0.0000000001
        max_iter =  500   propa_name ='hagedorn'  propa_name2 = 'taylor'
        Beta = F  P = F /
