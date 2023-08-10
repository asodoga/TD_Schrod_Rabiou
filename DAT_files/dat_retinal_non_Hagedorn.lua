
&potential
    pot_name='Retinal_JPCB2000' ! potential surface name
    option=1
   /
&basis_nd name='dp' nb_basis=3/
   &basis_nd name='fourier' nb =256  nq=256 /
    &basis_nd name='HO' nb =30  nq=30  SCALEQ= 0.9592722  /
      &basis_nd name= 'el' nb=2/
 &defGWP ndim=2  Elecindex = 2 Coef=(1.,0.)/
 &defWP0  sigma= 1.474256632  Qeq= 0.0   imp_k=0.0 phase=0.0 /
 &defWP0  sigma= 0.181095798  Qeq= 0.0   imp_k=0.0 phase=0.0 /

   &prop t0  = 0  tf  = 10000 delta_t =10    eps= 0.000000001
        max_iter =  500 propa_name ='non_hagedorn'  propa_name2 = 'taylor' /
