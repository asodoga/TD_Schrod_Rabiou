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
  TYPE (Basis_t), target         :: Basis,Basis0
    TYPE (Op_t)                  :: H
  TYPE (psi_t)                   :: psi0,psif,psi
  TYPE(propa_t)                  :: propa
  TYPE (GWP_t),allocatable       :: tab_GWP(:)
  real(Kind = Rk)                :: E,Norm,AVQ,DQ

!====================================================================
! for QML
integer :: ndim,nsurf,option,iq,ib
logical :: adiabatic
character (len=16)                  :: pot_name

ndim      = 1
nsurf     = 1
pot_name  = 'read_model'
adiabatic = .false.
option    = 1
CALL sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
 write(out_unitp,*)'ndim,nsurf',ndim,nsurf
write(out_unitp,*) 'pot_name'
  !====================================================================
  ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
  ! the basis/grid informations have to be put in a module
  CALL Read_Basis(Basis,nio=in_unitp)
  call  init_Basis1_TO_Basis2 (Basis0,Basis)
  CALL construct_primitive_basis(Basis,-ONE,0.70710678118_Rk)
  CALL construct_primitive_basis(Basis0,ZERO,ONE)
  !Call Write_Basis(Basis)
!====================================================================
!print*,"Basis is allocated",Basis_IS_allocated(Basis)
  write(out_unitp,*) 'Initialization of a complex psi'
  CALL init_psi(psi0,   Basis,    cplx=.TRUE.   ,grid =.false.)
  CALL init_psi(psif,   Basis,    cplx=.TRUE.   ,grid =.false.)
  CALL init_psi(psi,    Basis,    cplx=.TRUE.   ,grid =.false.)
   call Read_tab_GWP(tab_GWP=tab_GWP,nb_GWP=1,nio=in_unitp)
   !call psi_init_GWP0(psi=psi0,Tab_GWP=tab_GWP)
  ! call write_psi(psi=psi0,psi_cplx=.false.,print_psi_grid=.true.,print_basis=.false.,t=ZERO,int_print=100,real_part=.true.)
    !call Calc_average_energy(psi0,E)
  !call Calc_std_dev_AVQ_1D(psi0,1,AVQ,DQ)
  !call Set_Op(H,Basis)
  ! CALL Make_Mat_OP(H)
  !call  write_Op(H)
  STOP 'calcul de H|psi> est fait'

  CALL read_propa(propa)
  CALL propagation(Psif,Psi0,propa)
 ! CALL Write_psi(psif)
  write(out_unitp,*) 'deallocation'
  CALL dealloc_psi(psi0)
  CALL dealloc_psi(psif)


END PROGRAM TD_SCHROD
