
&potential
    pot_name  = 'H2'
    option=1
   /
&basis_nd name='dp' nb_basis=3/
   &basis_nd name='HO' nb =10  nq=10 Q0 = -1.0 SCALEQ=  1.0  /
    &basis_nd name='HO' nb =10  nq=10 Q0 = 0.0 SCALEQ=  1.0  /
      &basis_nd name= 'el' nb=1/
 &defGWP ndim=2  Elecindex = 1 Coef=(1.,0.)/
 &defWP0  sigma= 1.4142135623730951  Qeq= -1.0   imp_k=0.0 phase=0.0 /
 &defWP0  sigma= 1.4142135623730951  Qeq= 0.0   imp_k=0.0 phase=0.0 /

   &prop t0  = 0  tf  = 50 delta_t =0.1    eps= 0.000000001
        max_iter =  500 propa_name ='hagedorn'  propa_name2 = 'taylor' /
