
&potential
    pot_name  = 'HenonHeiles'
    option=1
   /
&basis_nd name='dp' nb_basis=3/
   &basis_nd name='fourier' nb = 256 nq= 256 SCALEQ=  1/
   &basis_nd name='HO' nb =30        nq=30 SCALEQ=   0.9592722/
      &basis_nd name= 'el' nb=1/
 &defGWP ndim=2  Elecindex = 1 Coef=(1.,0.)/
  &defWP0  sigma= 1.474256632  Qeq= 0.0   imp_k=0.0 phase=0.0 /
  &defWP0  sigma= 0.181095798  Qeq= 0.0   imp_k=0.0 phase=0.0 /

   &prop t0  = 0  tf  = 10 delta_t =0.1  eps= 0.00000001
        max_iter =  500 propa_name ='hagedorn'  propa_name2 = 'taylor' /
