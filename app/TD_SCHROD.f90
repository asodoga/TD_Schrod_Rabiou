PROGRAM TD_SCHROD
   USE UtilLib_m
   USE Basis_m
   USE psi_m
   USE op_m
   USE Param_WP0_m
   USE Propa_m
   USE Ana_psi_m
   USE lanczos_m
   IMPLICIT NONE
   TYPE(Basis_t), target          :: Basis, Basis0
   !TYPE(Op_t)                    :: H
   TYPE(psi_t)                    :: psi0, psif, psi
   TYPE(propa_t)                  :: propa
   TYPE(GWP_t), allocatable       :: tab_GWP(:)
   real(Kind=Rk)                  :: E, Norm, K(1)
   real(kind=Rk)                  :: x(2), sx(2), Pt(2)
!====================================================================
! for QML
   integer :: ndim, nsurf, option
   logical :: adiabatic
   character(len=16)                  :: pot_name
   ndim = 2
   nsurf = 1
   pot_name = 'read_model'
   adiabatic = .false.
   option = 1
   call sub_Init_Qmodel(ndim, nsurf, pot_name, adiabatic, option)
   write (out_unitp, *) 'ndim,nsurf', ndim, nsurf
   write (out_unitp, *) 'pot_name'
   !====================================================================
   ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
   ! the basis/grid informations have to be put in a module
   call Read_Basis(Basis, nio=in_unitp)
   call init_Basis1_TO_Basis2(Basis0, Basis)
   call construct_primitive_basis(Basis)
   call construct_primitive_basis(Basis0)
   !Call Write_Basis(Basis)
!====================================================================
!print*,"Basis is allocated",Basis_IS_allocated(Basis)
   write (out_unitp, *) 'Initialization of  psi0'

   call init_psi(psi0, Basis, cplx=.TRUE., grid=.false.)
   call init_psi(psif, Basis, cplx=.TRUE., grid=.false.)
   call init_psi(psi, Basis0, cplx=.TRUE., grid=.false.)

   call Read_tab_GWP(tab_GWP=tab_GWP, nb_GWP=1, nio=in_unitp)
   !call test_basitogridgridtobasis(Basis)
   call psi_init_GWP0(psi=psi0, Tab_GWP=tab_GWP)
   psi%CVec= ZERO
   !call Calc_psi_tild_coeff(psi0,psi0)
   !call Calc_average_energy(psi0, E)
   !call TEST_D(psi2=psi,psi1=psi0)
   !call psi0_init(psi0)
   !call Calc_average_energy(psi0, E)
    !call Test_calc_S(S=S, nb=3, nq=15, x1=ZERO, x2=ONE, s1=ONE, s2=ONE, p1=ZERO, p2=ONE)
   !call Calc_Av_imp_k_nD(psi0, Pt)
   !call Calc_AVQ_nD0(psi0=psi0, SQ=X)
   !call Calc_AVQ_nD0(psi0=psi0, AVQ=x, SQ=sx)
   !call Set_Op(H,Basis)
   ! call Make_Mat_OP(H)
   !call  write_Op(H)
   STOP 'calcul de H|psi> est fait'

   call read_propa(propa)
   call propagation(psif, psi0, propa)
   ! CALL Write_psi(psif)
   write (out_unitp, *) 'deallocation'
   call dealloc_psi(psi0)
   call dealloc_psi(psif)

END PROGRAM TD_SCHROD
