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
   TYPE(Basis_t), target          :: Basis,Basis0
   !TYPE(Op_t)                    :: H
   TYPE(psi_t)                    :: psi0, psif, psi
   TYPE(propa_t)                  :: propa
   TYPE(GWP_t), allocatable       :: tab_GWP(:)
   real(Kind=Rk)                  :: E, Norm !, x(2), y1(2), y2(2) ,p(1), S(1),x0(1)
   complex(kind=Rk)               :: Alpha(1)
   real(kind=Rk),allocatable  :: K(:,:)
!====================================================================
! for QML
   integer :: ndim, nsurf, option
   logical :: adiabatic
   character(len=16)                  :: pot_name

   ndim = 1
   nsurf = 1
   pot_name = 'read_model'
   adiabatic = .false.
   option = 1
   call sub_Init_Qmodel(ndim, nsurf, pot_name, adiabatic, option)
   write (out_unitp, *) 'ndim,nsurf', ndim, nsurf
   write (out_unitp, *) 'pot_name'

   !p(1)=ONETENTH; S(1)=ONE;x0(1)=ONETENTH
   !====================================================================
   ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
   ! the basis/grid informations have to be put in a module
   call Read_Basis(Basis, nio=in_unitp)
   call construct_primitive_basis(Basis)
   !call init_Basis1_TO_Basis2(Basis0, Basis)
   !write (out_unitp, *) 'p',p,'*********************************'
   !call construct_primitive_basis(Basis0, x=x0, sx=s,p=p)
   
   !Call Write_Basis(Basis)
!====================================================================
!print*,"Basis is allocated",Basis_IS_allocated(Basis)
   write (out_unitp, *) 'Initialization of  psi0'

   call init_psi(psi0, Basis, cplx=.TRUE., grid=.false.)
   call init_psi(psif, Basis, cplx=.TRUE., grid=.false.)
   call init_psi(psi, Basis, cplx=.TRUE., grid=.false.)

   call Read_tab_GWP(tab_GWP=tab_GWP, nb_GWP=1, nio=in_unitp)
   !call test_basitogridgridtobasis(Basis)
   call psi_init_GWP0(psi=psi0, Tab_GWP=tab_GWP)
   !call psi0_init(psi0)
   call Calc_average_energy(psi0, E)
   !call Calc_Integral_cplx(psi0, Alpha, 2)
   !call Construct_triband_matrix(Mat=K,psi=psi0,num_max=10)
   !call TEST_Lonaczos_cplx(psi0,10)
   !psi%CVec = CZERO
   !psi0%CVec = CZERO
   !psi%CVec(1) = CONE
   !call Projection(psi0, psi)
   !call TEST_S_cplx(Basis)
   !call Calc_AVQ_nD0(psi0=psi0,AVQ=y1, SQ=y2)
   !call Calc_AVQ_nD(psi0=psi0, AVQ=y1, SQ=y2)
   !call Set_Op(H,Basis)
   ! call Make_Mat_OP(H)
   !call  write_Op(H)
   !STOP 'calcu de H|psi> est fait'

   call read_propa(propa)
   call propagation(psif, psi0, propa)
   ! CALL Write_psi(psif)
   write (out_unitp, *) 'deallocation'
   call dealloc_psi(psi0)
   call dealloc_psi(psif)

END PROGRAM TD_SCHROD
