PROGRAM TD_SCHROD
   USE QDUtil_m
   USE Basis_m
   USE psi_m
   USE op_m
   USE Param_WP0_m
   USE Propa_m
   USE Ana_psi_m
   USE lanczos_m
   Use Vp_m
   IMPLICIT NONE
   TYPE(Basis_t), target             :: Basis
   !TYPE(Op_t)                       :: H
   TYPE(psi_t)                       :: psi0, psif
   TYPE(propa_t)                     :: propa
   TYPE(GWP_t), allocatable          :: tab_GWP(:)
   real(Kind=Rkind)                  :: E, Norm 
   real(Kind=Rkind) ,allocatable     :: Mat(:,:)
 
!-------------------------------------------------------------------------------
! for QML
   integer :: ndim, nsurf, option
   logical :: adiabatic
   character(len=16)                  :: pot_name

   ndim = 2
   nsurf = 2
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
!-------------------------------------------------------------------------------------------
!print*,"Basis is allocated",Basis_IS_allocated(Basis)
   write (out_unit, *) "------------------- Initialization of  psi0 ---------------------"

   call init_psi(psi0, Basis, cplx=.TRUE., grid=.false.)
   call init_psi(psif, Basis, cplx=.TRUE., grid=.false.)

   call Read_tab_GWP(tab_GWP=tab_GWP, nb_GWP=1, nio=in_unit)
   !call test_basitogridgridtobasis(Basis)
   call psi_init_GWP(psi=psi0, Tab_GWP=tab_GWP)
   call Calc_average_energy(psi0, E)
   call Calc_Norm_OF_Psi(psi0,Norm)
   write (out_unit, *)'-------------Energy And Norme initial WP0-----------------------------'
   write (out_unit, *) ' <psi|H|psi> ',E,'<psi|psi>',Norm
   !call H_test(psi0)
  ! call Vp_test(psi0)
  !call Calc_GlobalOverlap_S(Mat,Basis)
   !call Set_Op(H,Basis)
   ! call Make_Mat_OP(H)
   !call  write_Op(H)
   !call  Calc_Avg_A_nDtemp(psi0,At)
   STOP 'calcul de H|psi> est fait'

   call read_propa(propa)
   call propagation(psif, psi0, propa)
   print*, 'Fin de la propagationt'
   ! CALL Write_psi(psif)
   write (out_unit, *) 'deallocation'
   call dealloc_psi(psi0)
   call dealloc_psi(psif)

END PROGRAM TD_SCHROD
