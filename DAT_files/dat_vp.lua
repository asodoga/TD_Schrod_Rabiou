
&potential
    pot_name  = 'Tully'
    option=1
   /
&basis_nd name='dp' nb_basis=2 /
 &basis_nd name='boxab' nb=5   nq=15 A=-25 B=25/
&basis_nd name= 'el' nb=2/
 &defGWP ndim=1  Elecindex = 1 Coef=(1.,0.)/
 &defWP0  sigma=   1.41421356237  Qeq= 0.0   imp_k=1.0 phase=0.0 / 

   &prop t0  = 0.0  tf  = 0.5 delta_t =0.001 eps= 0.0000000001
        max_iter =  500   propa_name ='non_hagedorn'  propa_name2 = 'VP'
        Beta = F  P = F /
