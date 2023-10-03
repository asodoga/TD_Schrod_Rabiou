
&potential
    pot_name  = 'HenonHeiles'
    option=1
   /
&basis_nd name='dp' nb_basis=3/
   &basis_nd name='HO' nb =35 nq=35 Q0 =1.0 Imp_k=1.0  alpha =( 1.,0.0) /
     &basis_nd name='HO' nb =35 nq=35 Q0 =1.0 Imp_k=1.0  alpha =(1.,0.0) /
      &basis_nd name= 'el' nb=2/
 &defGWP ndim=2  Elecindex = 1 Coef=(1.,0.)/
 &defWP0  sigma=   1.41421356237  Qeq= 1.0  imp_k=1.0 phase=0.0 /
  &defWP0  sigma=  1.41421356237   Qeq=1.0   imp_k=1.0 phase=0.0 /
  

   &prop t0  = 0.0  tf  = 10.0 delta_t =0.01 eps= 0.000000000000001
        max_iter =  500  kmax=10 propa_name ='hagedorn'  propa_name2 = 'taylor' /
