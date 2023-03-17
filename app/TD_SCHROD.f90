PROGRAM TD_SCHROD
  USE UtilLib_m
  USE param_WP0_m
  USE Basis_m
  USE psi_m
  USE op_m
  USE Propa_m
  USE Ana_psi_m
  USE lanczos_m
  IMPLICIT NONE
  TYPE (Basis_t), target       :: Basis_i,Basis_f
    TYPE (Op_t)                :: H
  TYPE (psi_t)                 :: psi0,psif,psi
  TYPE(propa_t)                :: propa
  TYPE (GWP_t),allocatable     :: tab_GWP(:)
  real(Kind = Rk)              :: E,Norm,AVQ,SQ

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
call sub_Init_Qmodel(ndim,nsurf,pot_name,adiabatic,option)
write(out_unitp,*)'ndim,nsurf',ndim,nsurf
write(out_unitp,*) 'pot_name'
  !====================================================================
  ! read some informations (basis set/grid) : numbers of basis functions, grid points ...
  ! the basis/grid informations have to be put in a module
  call Read_Basis(Basis_i,nio=in_unitp)
  call construct_primitive_basis(Basis_i)
  !STOP 'calcul de H|psi> est fait'
  call init_Basis1_TO_Basis2 (Basis_f,Basis_i)
  call construct_primitive_basis(Basis_f)

  !====================================================================
!print*,"Basis is allocated",Basis_IS_allocated(Basis)
  write(out_unitp,*) 'Initialization of a complex psi'
  call init_psi(psi0,   Basis_i,    cplx=.TRUE.   ,grid =.false.)
  call init_psi(psif,   Basis_f,    cplx=.TRUE.   ,grid =.false.)
  call init_psi(psi,   Basis_i,    cplx=.TRUE.   ,grid =.true.)
  call Read_tab_GWP(tab_GWP=tab_GWP,nb_GWP=1,nio=in_unitp)
  call psi_init_GWP0(psi=psi0,Tab_GWP=tab_GWP)
 ! call write_psi(psi=psi0,psi_cplx=.false.,print_psi_grid=.true.,print_basis=.false.,t=ZERO,int_print=100,real_part=.true.)
  !call Write_Basis(Basis_0)


  !call  Buld_S(S=Basis_i%tab_basis(1)%S,d0gb1=Basis_i%tab_basis(1)%d0gb &
     ! ,d0gb2=Basis_i%tab_basis(1)%d0gb,nb=Basis_i%tab_basis(1)%nb,w1=Basis_i%tab_basis(1)%w)

  !CALL  Calc_basis(psi0%Basis, Basis_f,ONE,ONE)
  !call Hagedorn(psi,psi0,Basis_f)  
  
   !CALL  Calc_std_dev_AVQ_1D(psi0,1,AVQ,SQ)
   !call Calc_average_energy(psi0,E)
 !call diff()
  !call Set_Op(H,Basis)
  ! CALL Make_Mat_OP(H)
  !call  write_Op(H)
  !STOP 'calcul de H|psi> est fait'

  CALL read_propa(propa)
  CALL propagation(psif,psi0,propa)
 ! CALL Write_psi(psif)
  write(out_unitp,*) 'deallocation'
  CALL dealloc_psi(psi0)
  CALL dealloc_psi(psif)


END PROGRAM TD_SCHROD
