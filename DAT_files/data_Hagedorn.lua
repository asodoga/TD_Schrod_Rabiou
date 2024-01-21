
&potential
    pot_name  = 'HenonHeiles'
    option=1
    adiabatic = f
    nsurf = 1
    ndim = 2
   /
&basis_nd name='dp' nb_basis=3/
&basis_nd name='HO' nb =35 nq=45 Q0 =2.0 SCALEQ = 1.0 Imp_k=0.0  alpha =( 1.,0.0) /
&basis_nd name='HO' nb =35 nq=45 Q0 =0.0 SCALEQ = 1.0 Imp_k=0.0  alpha =( 1.,0.0) /
&basis_nd name= 'el' nb=1/

&defGWP ndim=2  Elecindex = 1 Coef=(1.,0.)/
&defWP0  sigma=  1.4142135623730951  Beta=0.0  Qeq= 2.0   imp_k=0.0 gamma=0.0 /
&defWP0  sigma=  1.4142135623730951  Beta=0.0  Qeq= 0.0   imp_k=0.0 gamma=0.0 /

&prop t0  = 0.0  tf  = 60 delta_t =0.1 eps= 0.0000000001
max_iter =  500   propa_name ='hagedorn'  propa_name2 = 'tayloR' Beta = T  P = T /
      