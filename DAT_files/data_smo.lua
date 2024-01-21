
&potential
    pot_name  = 'HenonHeiles'
    option=1
    adiabatic = f
    nsurf = 1
    ndim = 2
   /

&basis_nD name='smolyak' nb_basis=3 LB=2  /
&basis_nD name='HO' A_smol=1  B_smol=1  SCALEQ =1. Q0=1  Imp_k =0. alpha = (1.,0.)/
&basis_nD name='HO' A_smol=1  B_smol=1  SCALEQ =1. Q0=1  Imp_k =0. alpha = (1.,0.)/
&basis_nd name= 'el' nb=1/

