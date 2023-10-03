
&potential
    pot_name  = 'HenonHeiles'
    option=1
   /
&basis_nd name='dp' nb_basis=2/
   &basis_nd name='HO' nb =15 nq=15 Q0 =1.0 Imp_k=1.0  alpha =( 1.,0.0) /
      &basis_nd name= 'el' nb=1/
 &defGWP ndim=1  Elecindex = 1 Coef=(1.,0.)/
 &defWP0  sigma=  1.41421356237  Qeq= 1.0   imp_k=1.0 phase=0.0 /

   &prop t0  = 0.0  tf  = 20.0 delta_t =0.1 eps= 0.000000000000001
        max_iter =  500  kmax=10 propa_name ='non_hagedorn'  propa_name2 = 'taylor' /
