
&potential
    pot_name  = 'HenonHeiles'
    option=1
    adiabatic = f
    nsurf = 1
    ndim = 1
   /
&basis_nd name='dp' nb_basis=2/
&basis_nd name='HO' nb =50 nq=70 Q0 =0.0 SCALEQ = 1.0 Imp_k=0.0  alpha =( 1.,0.0) /
&basis_nd name= 'el' nb=1/

&defGWP ndim=1  Elecindex = 1 Coef=(1.,0.)/
&defWP0  sigma=  1.4142135623730951  Beta=0.0  Qeq= 0.0   imp_k=0.0 gamma=0.0 /

&prop t0  = 0.0  tf  = 60 delta_t =0.1 eps= 0.0000000000000001
max_iter =  20   propa_name ='non_hagedorn'  propa_name2 = 'taylor' Beta = T  P = T /
      