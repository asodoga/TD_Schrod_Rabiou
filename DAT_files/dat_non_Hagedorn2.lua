
&potential
    pot_name  = 'HenonHeiles'
    option=1
   /
&basis_nd name='dp' nb_basis=4/
   &basis_nd name='HO' nb =10  nq= 10 Q0 = 0.0 SCALEQ=  1/
   &basis_nd name='HO' nb =10  nq= 10 Q0 = 0.0 SCALEQ=  1/
   &basis_nd name='HO' nb =10  nq= 10 Q0 = 0.0 SCALEQ=  1/
      &basis_nd name= 'el' nb=1/
 &defGWP ndim=3  Elecindex = 1 Coef=(1.,0.)/
 &defWP0  sigma= 1.41421356237  Qeq= -1.0   imp_k=0.0 phase=0.0 /
 &defWP0  sigma= 1.41421356237  Qeq= -1.0   imp_k=0.0 phase=0.0 /
 &defWP0  sigma= 1.41421356237  Qeq= -1.0   imp_k=0.0 phase=0.0 /
 
 
 
  

   &prop t0  = 0  tf  = 20 delta_t =0.1  eps= 0.00000001
        max_iter =  500 propa_name ='non_hagedorn'  propa_name2 = 'taylor' /
