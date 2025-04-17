PROGRAM TD_SCHROD
   USE QDUtil_m
   USE Basis_m
   USE psi_m
   USE op_m
   USE Param_WP0_m
   USE Propa_m
   USE Ana_psi_m
   USE lanczos_m
   USE Sub_Vp_m
   IMPLICIT NONE
   TYPE(Basis_t), target             :: Basis
   TYPE(Op_t)                        :: H
   TYPE(psi_t)                       :: psi0, psif
   TYPE(propa_t)                     :: propa
   TYPE(GWP_t), allocatable          :: tab_GWP(:)
   integer, allocatable              :: Tab_iq(:, :)
   real(Kind=Rkind), allocatable     :: V(:, :, :)
   real(Kind=Rkind)                  :: E, Norm ,Pt(2),Qt(2),SQt(2)
   complex(Kind=Rkind)               :: At(2)
   real(Kind=Rkind)                  :: t_deltat=0
 
!-------------------------------------------------------------------------------
! for QML
   integer :: ndim, nsurf, option
   logical :: adiabatic
   character(len=16)                  :: pot_name

   ndim = -1
   nsurf = -1
   pot_name = 'read_model'
   adiabatic = .false.
   option = 1
   call sub_Init_Qmodel(ndim, nsurf, pot_name, adiabatic, option)
   write (out_unit, *) 'ndim,nsurf', ndim, nsurf
   write (out_unit, *) 'pot_name'
   
   
   !---------------------------------------------------------------------------------
   ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
   ! the basis/grid informations have to be put in a module
   call Read_Basis(Basis, nio=in_unit)
   call construct_primitive_basis(Basis)
   !call Write_Basis(Basis) 
   !stop 'cc'
!-------------------------------------------------------------------------------------------
!print*,"Basis is allocated",Basis_IS_allocated(Basis)
   write (out_unit, *) "------------------- Initialization of  psi0 ---------------------"

   call init_psi(psi0, Basis, cplx=.true., grid=.false.)
   call init_psi(psif, Basis, cplx=.true., grid=.false.)
   call Read_tab_GWP(tab_GWP=tab_GWP, nb_GWP=1, nio=in_unit)
   !call test_basitogridgridtobasis(Basis)
   call Calc_tab_Iq0(Tab_Iq,psi0%Basis)
   call Set_Op(H, psi0%Basis,Tab_Iq)
   call psi_init_GWP(psi=psi0, Tab_GWP=tab_GWP)
   !call Ecrire_psi(psi0,nio=100,t=ZERO)

   !print*,"psi0",psi0%CVec(:)
   !psi0%CVec=CZERO
   !psi0%CVec(1) =CONE
   !call Calc_Basis_parameters(psi0,Qt,SQt,At,Pt)
   !call Get_Basis_Parameters(Basis,Qt,SQt,At,Pt)
   !call Calc_Avg_A_nD(psi0, At)
   !call Calc_Av_E(E,psi0,H)
   !call Calc_Norm_OF_Psi(psi0,Norm)
   !write (out_unit, *)'-------------Energy And Norme initial WP0-----------------------------'
   !write (out_unit, *) ' <psi|H|psi> ',E,'<psi|psi>',Norm
   !write (out_unit, *) '-------------End  Initialization of  psi0 ---------------------'
   !call test_op(Basis,psi0) ! diagonalisation
   ! call Make_Mat_H(H)
   !call test_openmp_op(Basis)
   !call calc_nac(Basis)
   !call calc_VV(Basis)
   call read_propa(propa)
   !call Calc_Scalar_Pot(V, Basis)
   !call test_taylor(psi0, psif, propa, H)
   !STOP 'calcul de H|psi> est fait'
   call propagation(psif, psi0, propa)
   print*, 'Fin de la propagationt'
   ! CALL Write_psi(psif)
   write (out_unit, *) 'deallocation'
   call dealloc_psi(psi0)
   call dealloc_psi(psif)

END PROGRAM TD_SCHROD
