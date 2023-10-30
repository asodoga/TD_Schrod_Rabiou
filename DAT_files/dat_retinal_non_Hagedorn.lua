
&potential
    pot_name='Retinal_JPCB2000' ! potential surface name
    option=1
   /
&basis_nd name='dp' nb_basis=3/
   &basis_nd name='fourier' nb =156  nq=165  Imp_k=0.0  alpha =( 1.0,0.0) /
    &basis_nd name='HO' nb =50  nq=55  Imp_k=0.0 alpha =( 0.920203149,0.0) /
      &basis_nd name= 'el' nb=2/
 &defGWP ndim=2  Elecindex = 2 Coef=(1.,0.)/
  &defWP0  sigma= 0.181095798  Qeq= 0.0   imp_k=0.0 phase=0.0 /
 &defWP0  sigma= 1.474256632  Qeq= 0.0   imp_k=0.0 phase=0.0 /

   &prop t0  = 0  tf  = 10000 delta_t =10    eps= 0.000000001
        max_iter =  500 propa_name ='hagedorn'  propa_name2 = 'taylor' /
