
&potential
  pot_name='Retinal_JPCB2000' ! potential surface name
  option=1
   /
&basis_nd name='dp' nb_basis=3/  
&basis_nd name='HO' nb = 20 nq=22   /
   &basis_nd name='HO' nb = 30 nq=32    /
      &basis_nd name= 'el' nb=2/


   &defGWP I_ElecChannel=2 Coef=(1.,0.) /
       &defWP0 DQ= 1.474256632  Q0=0.0  k=0.0 phase=0.0 /
        &defWP0 DQ=  0.181095798  Q0=0.0  k=0.0 phase=0.0 /


      &prop t0  = 0  tf  =1000 delta_t =10     eps= 0.00000000001
        max_iter =  5000 propa_name = 'taylor' kmax =  200/
